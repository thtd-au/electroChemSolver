#include "ButlerVolmer.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "PstreamReduceOps.H"
#include "token.H"
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <limits>

namespace Foam
{
namespace
{
// Numerical helpers are private to this translation unit.
namespace electroChemBV
{

template<class Real> struct Evaluation
{
    Real residual;
    Real derivative;
    Real logScale;
};

// Integrated BV current is exp(logA+a*x)-exp(logC-c*x).
// Scale all terms together before exponentiation to avoid overflow while
// bracketing. residual/derivative is still the unscaled Newton step.
template<class Real>
Evaluation<Real> evaluate
(
    Real logA, Real logC, Real a, Real c, Real x, Real target
)
{
    const Real positive = logA + a*x;
    const Real negative = logC - c*x;
    const Real logTarget = target == Real(0)
        ? -std::numeric_limits<Real>::infinity() : std::log(std::abs(target));
    const Real scale = std::max(std::max(positive, negative), logTarget);
    const Real ep = std::exp(positive - scale);
    const Real en = std::exp(negative - scale);
    const Real et = target == Real(0) ? Real(0)
        : std::copysign(std::exp(logTarget - scale), target);
    return {ep-en-et, a*ep+c*en, scale};
}

template<class Real> struct Root
{
    Real x;
    int iterations;
};

template<class Real>
Root<Real> solve
(
    Real logA, Real logC, Real a, Real c, Real target,
    Real currentTolerance, Real initial, Real initialSpan, int maxIterations
)
{
    Real span = initialSpan;
    Real lo = initial-span, hi = initial+span;
    int expanded = 0;
    while
    (
        evaluate(logA,logC,a,c,lo,target).residual > Real(0)
     || evaluate(logA,logC,a,c,hi,target).residual < Real(0)
    )
    {
        if (++expanded >= maxIterations)
            throw std::runtime_error("Could not bracket the electrode potential");
        span *= Real(2);
        lo = initial-span;
        hi = initial+span;
        if (!std::isfinite(lo) || !std::isfinite(hi))
            throw std::runtime_error("Non-finite electrode-potential bracket");
    }

    Real x = initial;
    Real lastStep = hi-lo;
    Real previousStep = lastStep;
    for (int iter=0; iter<maxIterations; ++iter)
    {
        const auto e = evaluate(logA,logC,a,c,x,target);
        const Real scaledTolerance = std::exp(std::log(currentTolerance)-e.logScale);
        if (std::abs(e.residual) <= scaledTolerance) return {x,iter+1};
        if (e.residual > Real(0)) hi=x;
        else lo=x;
        const Real newtonStep = e.derivative > Real(0)
            ? e.residual/e.derivative : std::numeric_limits<Real>::infinity();
        Real next = x-newtonStep;
        // Accept Newton while it stays in the bracket and makes sufficient
        // step-size progress. Otherwise bisect. Unlike periodic bisection,
        // this preserves rapid local Newton convergence (rtsafe criterion).
        if
        (
            !std::isfinite(next) || next<=lo || next>=hi
         || std::abs(newtonStep) > Real(0.5)*std::abs(previousStep)
        ) next = lo+(hi-lo)/Real(2);
        previousStep = lastStep;
        lastStep = next-x;
        x=next;
    }
    throw std::runtime_error("Electrode-potential root did not converge");
}

template<class Real>
Real localCurrent(Real i0, Real a, Real c, Real eta)
{
    const Real largest = std::max(a*eta, -c*eta);
    const Real limit = std::log(std::numeric_limits<Real>::max())-Real(4);
    if (!std::isfinite(eta) || largest+std::log(i0)>limit)
        throw std::runtime_error("Local BV exponential exceeds the scalar range");
    // expm1 retains precision close to equilibrium; both branches vanish
    // together at eta=0. No clipping of the physical current is performed.
    if (largest < Real(50))
        return i0*(std::expm1(a*eta)-std::expm1(-c*eta));
    return std::exp(std::log(i0)+a*eta)-std::exp(std::log(i0)-c*eta);
}

// Statistics use the area-weighted inner product over both electrodes.
template<class Real>
Real aitkenFactor
(
    Real omega, Real initial, Real lower, Real upper, bool history,
    Real numerator, Real denominator, Real oldNorm, Real newNorm
)
{
    if (!history) return initial;
    if (newNorm > Real(4)*oldNorm) return std::max(lower, Real(0.5)*omega);
    if
    (
        denominator > Real(100)*std::numeric_limits<Real>::epsilon()
            *std::max(oldNorm, newNorm)
    )
    {
        const Real proposed = -omega*numerator/denominator;
        return std::isfinite(proposed) && proposed > Real(0)
            ? std::max(lower, std::min(upper, proposed))
            : std::max(lower, Real(0.5)*omega);
    }
    return omega;
}

} // namespace electroChemBV

scalar electronCoefficient(const fvPatchScalarField& patchField)
{
    OStringStream output;
    output << token::BEGIN_BLOCK << nl;
    patchField.write(output);
    output << token::END_BLOCK << nl;
    IStringStream input(output.str());
    const dictionary dict(input);
    return dict.get<scalar>("electronStoichCoeff");
}

bool finitePositive(const scalar value)
{
    return std::isfinite(value) && value > 0;
}
}

defineTypeNameAndDebug(ButlerVolmer, 0);
addToRunTimeSelectionTable(currentModel, ButlerVolmer, dictionary);

ButlerVolmer::ButlerVolmer
(
    fvMesh& mesh,
    const dictionary& currentDict,
    const volScalarField& phiE
)
:
    currentModel(mesh, currentDict, phiE),
    appliedCurrentDensity_(currentDict.get<scalar>("appliedCurrentDensity")),
    referencePatchName_(currentDict.get<word>("currentDensityReferencePatch")),
    electrodeNames_(currentDict.subDict("electrodes").sortedToc()),
    patchIDs_(electrodeNames_.size(), -1),
    nuE_(electrodeNames_.size(), 0),
    i0_(electrodeNames_.size(), 0),
    alphaA_(electrodeNames_.size(), 0),
    alphaC_(electrodeNames_.size(), 0),
    Eeq_(electrodeNames_.size(), 0),
    phiS_(electrodeNames_.size(), 0),
    targetCurrent_(electrodeNames_.size(), 0),
    area_(electrodeNames_.size(), 0),
    actualCurrent_(electrodeNames_.size(), 0),
    bvResidual_(electrodeNames_.size(), 0),
    boundaryResidual_(electrodeNames_.size(), 0),
    rootIterations_(electrodeNames_.size(), 0),
    T_(0), F_(0), R_(8.31446261815324),
    absTol_(0), relTol_(0), rootAbsTol_(0), rootRelTol_(0),
    rootMaxIter_(100), currentRelaxation_(1),
    relaxationMethod_("fixed"),
    relaxationMin_(0.005), relaxationMax_(0.1), activeRelaxation_(1),
    previousResidual_(electrodeNames_.size()),
    relaxationTimeIndex_(-1), haveResidual_(false),
    initialized_(false), evaluated_(false), converged_(false)
{
    const IOdictionary& ec = mesh.lookupObject<IOdictionary>("electroChemistry");
    if (ec.subDict("potential").get<word>("model") == "none")
    {
        FatalIOErrorInFunction(currentDict_)
            << "ButlerVolmer requires an active electrolyte potential model. "
            << "Use constantCurrent for potential-isolation tests."
            << exit(FatalIOError);
    }
    T_ = ec.get<scalar>("temperature");
    F_ = ec.get<scalar>("FaradayConstant");
    const dictionary& controls = mesh.solutionDict().subDict("Picard").subDict("current");
    absTol_ = controls.get<scalar>("absTol");
    relTol_ = controls.get<scalar>("relTol");
    rootAbsTol_ = controls.getOrDefault<scalar>("rootAbsTol", 1e-10);
    rootRelTol_ = controls.getOrDefault<scalar>("rootRelTol", 1e-10);
    rootMaxIter_ = controls.getOrDefault<label>("rootMaxIter", 100);
    currentRelaxation_ = controls.getOrDefault<scalar>("relaxation", 1.0);
    relaxationMethod_ = controls.getOrDefault<word>("relaxationMethod", "fixed");
    relaxationMin_ = controls.getOrDefault<scalar>("relaxationMin", 0.005);
    relaxationMax_ = controls.getOrDefault<scalar>("relaxationMax", 0.1);
    activeRelaxation_ = currentRelaxation_;
    if (relaxationMethod_ != "fixed" && relaxationMethod_ != "Aitken")
        FatalIOErrorInFunction(controls)
            << "current/relaxationMethod must be fixed or Aitken"
            << exit(FatalIOError);
    if
    (
        relaxationMethod_ == "Aitken"
     && (!finitePositive(relaxationMin_) || !finitePositive(relaxationMax_)
      || relaxationMin_ > currentRelaxation_
      || currentRelaxation_ > relaxationMax_ || relaxationMax_ > 1)
    )
        FatalIOErrorInFunction(controls)
            << "Require 0 < relaxationMin <= relaxation <= relaxationMax <= 1"
            << exit(FatalIOError);
    if (!finitePositive(currentRelaxation_) || currentRelaxation_ > 1)
        FatalIOErrorInFunction(controls)
            << "current/relaxation must be finite and in (0,1]"
            << exit(FatalIOError);

    if
    (
        !std::isfinite(appliedCurrentDensity_) || appliedCurrentDensity_ < 0
     || !finitePositive(T_) || !finitePositive(F_)
     || !finitePositive(absTol_) || !std::isfinite(relTol_) || relTol_ < 0
     || !finitePositive(rootAbsTol_) || !std::isfinite(rootRelTol_) || rootRelTol_ < 0
     || rootMaxIter_ < 1
    )
    {
        FatalIOErrorInFunction(currentDict_)
            << "Invalid BV settings: require finite non-negative appliedCurrentDensity, "
            << "positive T/F/absolute tolerances, non-negative relative tolerances "
            << "and rootMaxIter>=1." << exit(FatalIOError);
    }
    if (electrodeNames_.size() != 2)
    {
        FatalIOErrorInFunction(currentDict_)
            << "This BV implementation requires exactly two electrode patches, "
            << "one anodic and one cathodic." << exit(FatalIOError);
    }

    forAll(electrodeNames_, e)
    {
        const label patchi = mesh.boundaryMesh().findPatchID(electrodeNames_[e]);
        if (patchi < 0)
            FatalIOErrorInFunction(currentDict_) << "Unknown electrode patch "
                << electrodeNames_[e] << exit(FatalIOError);
        patchIDs_[e] = patchi;
        const fvPatchScalarField& p = phiE.boundaryField()[patchi];
        if (p.type() != "galvanostaticPotentialFlux")
            FatalIOErrorInFunction(currentDict_) << "Electrode " << electrodeNames_[e]
                << " requires galvanostaticPotentialFlux on phiE" << exit(FatalIOError);
        nuE_[e] = electronCoefficient(p);
        const dictionary& electrode = currentDict.subDict("electrodes").subDict(electrodeNames_[e]);
        i0_[e] = electrode.get<scalar>("i0");
        alphaA_[e] = electrode.get<scalar>("alpha_a");
        alphaC_[e] = electrode.get<scalar>("alpha_c");
        const dictionary& eq = electrode.subDict("equilibriumPotential");
        if (eq.get<word>("model") != "constant")
            FatalIOErrorInFunction(eq)
                << "Only equilibriumPotential/model constant is implemented"
                << exit(FatalIOError);
        Eeq_[e] = eq.get<scalar>("Eeq");
        if
        (
            !finitePositive(i0_[e]) || !finitePositive(alphaA_[e])
         || !finitePositive(alphaC_[e]) || !std::isfinite(Eeq_[e])
         || !std::isfinite(nuE_[e]) || mag(nuE_[e]) <= SMALL
        )
            FatalIOErrorInFunction(electrode)
                << "Require finite positive i0/alpha_a/alpha_c, finite Eeq "
                << "and non-zero finite electronStoichCoeff" << exit(FatalIOError);
    }
    if ((nuE_[0] > 0) == (nuE_[1] > 0))
        FatalIOErrorInFunction(currentDict_) << "Electrode electron coefficients "
            << "must have opposite signs" << exit(FatalIOError);
    if (electrodeNames_.find(referencePatchName_) < 0)
        FatalIOErrorInFunction(currentDict_) << "currentDensityReferencePatch must "
            << "name one of the configured electrodes" << exit(FatalIOError);
    forAll(phiE.boundaryField(), patchi)
    {
        if
        (
            phiE.boundaryField()[patchi].type() == "galvanostaticPotentialFlux"
         && electrodeNames_.find(mesh.boundary()[patchi].name()) < 0
        )
            FatalIOErrorInFunction(currentDict_) << "Electrode patch "
                << mesh.boundary()[patchi].name() << " lacks a current/electrodes entry"
                << exit(FatalIOError);
    }
    // No kappa or diffusion-current field is needed for this initial root.
    correct();
}

scalar ButlerVolmer::currentTolerance(const scalar value) const
{
    return max(absTol_, relTol_*max(mag(value), appliedCurrentDensity_));
}

void ButlerVolmer::correct()
{
    const label refi = mesh_.boundaryMesh().findPatchID(referencePatchName_);
    if (refi < 0) FatalErrorInFunction << "Missing reference electrode patch"
        << exit(FatalError);
    const scalar refArea = gSum(mesh_.magSf().boundaryField()[refi]);
    if (!finitePositive(refArea)) FatalErrorInFunction << "Invalid reference electrode area"
        << exit(FatalError);
    const scalar systemCurrent = appliedCurrentDensity_*refArea;
    if (!std::isfinite(systemCurrent)) FatalErrorInFunction
        << "System current exceeds the scalar range" << exit(FatalError);
    electrodeCurrentDensity_.primitiveFieldRef() = 0.0;
    forAll(electrodeCurrentDensity_.boundaryField(), patchi)
    {
        // Retain the previous electrode distribution for Picard relaxation.
        if (electrodeNames_.find(mesh_.boundary()[patchi].name()) < 0)
            electrodeCurrentDensity_.boundaryFieldRef()[patchi] = 0.0;
    }

    // Use one area-weighted factor over both electrodes: their distributions
    // are coupled through the electrolyte potential. Restart history per step.
    if (relaxationTimeIndex_ != mesh_.time().timeIndex())
    {
        relaxationTimeIndex_ = mesh_.time().timeIndex();
        activeRelaxation_ = currentRelaxation_;
        haveResidual_ = false;
    }
    List<scalarField> baseCurrent(electrodeNames_.size());
    List<scalarField> residual(electrodeNames_.size());
    scalar numerator = 0, denominator = 0, oldNorm = 0, newNorm = 0;
    bool historyValid = haveResidual_;

    forAll(electrodeNames_, e)
    {
        const label patchi = mesh_.boundaryMesh().findPatchID(electrodeNames_[e]);
        if (patchi < 0) FatalErrorInFunction << "Missing electrode "
            << electrodeNames_[e] << exit(FatalError);
        patchIDs_[e] = patchi;
        const scalarField& areas = mesh_.magSf().boundaryField()[patchi];
        const fvPatchScalarField& p = phiE_.boundaryField()[patchi];
        area_[e] = gSum(areas);
        if (!finitePositive(area_[e])) FatalErrorInFunction << "Invalid area on "
            << electrodeNames_[e] << exit(FatalError);

        // Center the root at the area-mean electrolyte potential plus Eeq.
        const scalar phiRef = gSum(areas*scalarField(p))/area_[e] + Eeq_[e];
        const scalar a = alphaA_[e]*F_/(R_*T_);
        const scalar c = alphaC_[e]*F_/(R_*T_);
        scalarField logTermsA(p.size(), 0), logTermsC(p.size(), 0);
        scalar maxA = -GREAT, maxC = -GREAT;
        forAll(p, f)
        {
            if (!finitePositive(areas[f]) || !std::isfinite(p[f]))
                FatalErrorInFunction << "Invalid area/potential on "
                    << electrodeNames_[e] << " face " << f << exit(FatalError);
            const scalar delta = phiRef-p[f]-Eeq_[e];
            const scalar logWeight = std::log(i0_[e]) + std::log(areas[f]);
            logTermsA[f] = logWeight+a*delta;
            logTermsC[f] = logWeight-c*delta;
            maxA = max(maxA, logTermsA[f]);
            maxC = max(maxC, logTermsC[f]);
        }
        reduce(maxA, maxOp<scalar>());
        reduce(maxC, maxOp<scalar>());
        scalar sumA = 0, sumC = 0;
        forAll(p, f)
        {
            sumA += std::exp(logTermsA[f]-maxA);
            sumC += std::exp(logTermsC[f]-maxC);
        }
        reduce(sumA, sumOp<scalar>());
        reduce(sumC, sumOp<scalar>());
        const scalar logA = maxA+std::log(sumA);
        const scalar logC = maxC+std::log(sumC);
        targetCurrent_[e] = (nuE_[e] > 0 ? 1 : -1)*systemCurrent;
        // Always solve the controller more tightly than the outer test.
        const scalar rootTolerance = min
        (
            max(rootAbsTol_*area_[e], rootRelTol_*systemCurrent),
            scalar(0.01)*max(absTol_*area_[e], relTol_*systemCurrent)
        );

        try
        {
            const auto root = electroChemBV::solve
            (
                logA, logC, a, c, targetCurrent_[e], rootTolerance,
                phiS_[e]-phiRef, scalar(2)*R_*T_/F_, int(rootMaxIter_)
            );
            phiS_[e] = phiRef+root.x;
            rootIterations_[e] = root.iterations;
            fvsPatchScalarField& jn = electrodeCurrentDensity_.boundaryFieldRef()[patchi];
            const scalarField previousCurrent(jn);
            scalarField candidate(jn.size(), 0);
            forAll(jn, f)
                candidate[f] = -electroChemBV::localCurrent
                    (i0_[e], a, c, phiS_[e]-p[f]-Eeq_[e]);
            // Remove only the root/roundoff integral error. The final BV
            // acceptance check still compares against the unmodified BV law.
            candidate += (-targetCurrent_[e]-gSum(candidate*areas))/area_[e];
            baseCurrent[e] = previousCurrent;
            baseCurrent[e] +=
                (-targetCurrent_[e]-gSum(previousCurrent*areas))/area_[e];
            residual[e] = candidate-baseCurrent[e];
            if (previousResidual_[e].size() != residual[e].size())
                historyValid = false;
            if (historyValid)
            {
                forAll(jn, f)
                {
                    const scalar old = previousResidual_[e][f];
                    const scalar delta = residual[e][f]-old;
                    numerator += areas[f]*old*delta;
                    denominator += areas[f]*sqr(delta);
                    oldNorm += areas[f]*sqr(old);
                    newNorm += areas[f]*sqr(residual[e][f]);
                }
            }

        }
        catch (const std::exception& error)
        {
            FatalErrorInFunction << "BV electrode " << electrodeNames_[e]
                << ": " << error.what() << exit(FatalError);
        }
    }
    reduce(historyValid, andOp<bool>());
    reduce(numerator, sumOp<scalar>());
    reduce(denominator, sumOp<scalar>());
    reduce(oldNorm, sumOp<scalar>());
    reduce(newNorm, sumOp<scalar>());
    activeRelaxation_ = relaxationMethod_ == "Aitken"
        ? electroChemBV::aitkenFactor
          (
              activeRelaxation_, currentRelaxation_, relaxationMin_, relaxationMax_,
              historyValid, numerator, denominator, oldNorm, newNorm
          )
        : currentRelaxation_;
    forAll(electrodeNames_, e)
    {
        const label patchi = patchIDs_[e];
        const scalarField& areas = mesh_.magSf().boundaryField()[patchi];
        fvsPatchScalarField& jn = electrodeCurrentDensity_.boundaryFieldRef()[patchi];
        // The first initialization has no previous physical current to mix.
        const scalar omega = initialized_ ? activeRelaxation_ : scalar(1);
        forAll(jn, f) jn[f] = baseCurrent[e][f]+omega*residual[e][f];
        const scalar offset = (-targetCurrent_[e]-gSum(scalarField(jn)*areas))/area_[e];
        forAll(jn, f) jn[f] += offset;
        previousResidual_[e] = residual[e];
    }
    haveResidual_ = initialized_;
    initialized_ = true;
    evaluated_ = false;
}

void ButlerVolmer::evaluateConstraint()
{
    const surfaceScalarField& diffusion = mesh_.lookupObject<surfaceScalarField>("jDiffFlux");
    const surfaceScalarField& conductivity = mesh_.lookupObject<surfaceScalarField>("kappaFace");
    converged_ = true;
    forAll(electrodeNames_, e)
    {
        const label patchi = patchIDs_[e];
        const fvPatchScalarField& p = phiE_.boundaryField()[patchi];
        const scalarField gradient(p.snGrad());
        const fvsPatchScalarField& jn = electrodeCurrentDensity_.boundaryField()[patchi];
        const scalarField& areas = mesh_.magSf().boundaryField()[patchi];
        const scalar a = alphaA_[e]*F_/(R_*T_);
        const scalar c = alphaC_[e]*F_/(R_*T_);
        bvResidual_[e] = 0;
        boundaryResidual_[e] = 0;
        scalar total = 0;
        try
        {
            forAll(jn, f)
            {
                const scalar bv = electroChemBV::localCurrent
                    (i0_[e],a,c,phiS_[e]-p[f]-Eeq_[e]);
                const scalar kappa = conductivity.boundaryField()[patchi][f];
                if (!finitePositive(kappa) || !finitePositive(areas[f]))
                    FatalErrorInFunction << "Invalid conductivity/area on electrode "
                        << electrodeNames_[e] << exit(FatalError);
                const scalar ionic = diffusion.boundaryField()[patchi][f]/areas[f]
                    -kappa*gradient[f];
                const scalar bvError = mag(jn[f]+bv);
                const scalar boundaryError = mag(jn[f]-ionic);
                if
                (
                    !std::isfinite(bvError) || !std::isfinite(boundaryError)
                 || bvError > currentTolerance(jn[f])
                 || boundaryError > currentTolerance(jn[f])
                ) converged_ = false;
                bvResidual_[e] = max(bvResidual_[e],bvError);
                boundaryResidual_[e] = max(boundaryResidual_[e],boundaryError);
                total -= jn[f]*areas[f]; // anodic-positive integrated BV current
            }
        }
        catch (const std::exception& error)
        {
            FatalErrorInFunction << "BV residual on " << electrodeNames_[e]
                << ": " << error.what() << exit(FatalError);
        }
        reduce(total,sumOp<scalar>());
        reduce(bvResidual_[e],maxOp<scalar>());
        reduce(boundaryResidual_[e],maxOp<scalar>());
        actualCurrent_[e] = total;
        const scalar integratedTolerance = max(absTol_*area_[e],relTol_*mag(targetCurrent_[e]));
        if (!std::isfinite(total) || mag(total-targetCurrent_[e]) > integratedTolerance)
            converged_ = false;
    }
    reduce(converged_,andOp<bool>());
    evaluated_ = true;
}

bool ButlerVolmer::converged() const
{
    if (!evaluated_) FatalErrorInFunction << "BV convergence requested before evaluation"
        << exit(FatalError);
    return converged_;
}

void ButlerVolmer::report() const
{
    Info<< "    ButlerVolmer: converged=" << (evaluated_ && converged_ ? "yes" : "no") << nl;
    Info<< "    current relaxation=" << activeRelaxation_
        << " (" << relaxationMethod_ << ")" << nl;
    forAll(electrodeNames_, e)
    {
        Info<< "    " << electrodeNames_[e] << ": phiS=" << phiS_[e]
            << " V, I_applied/target=" << actualCurrent_[e] << " / " << targetCurrent_[e]
            << " A, root iterations=" << rootIterations_[e] << nl
            << "    max BV residual=" << bvResidual_[e]
            << ", max ionic boundary residual=" << boundaryResidual_[e]
            << " A/m2" << nl;
    }
}

} // namespace Foam

