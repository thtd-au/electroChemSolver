#include "concentrationButlerVolmer.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"
#include "dimensionSets.H"
#include "fixedGradientFvPatchFields.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "token.H"
#include "PstreamReduceOps.H"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Foam
{
namespace
{
namespace electroChemConcentration
{
// Returns log((c/cref)^order). Exact zero factors extinguish a branch;
// no positive concentration floor or artificial reaction is introduced.
template<class Real>
Real logFactor(Real concentration, Real reference, Real order)
{
    if (!std::isfinite(concentration) || concentration < Real(0))
        throw std::runtime_error("BV requires finite non-negative surface concentrations");
    if (!std::isfinite(reference) || reference <= Real(0) || !std::isfinite(order))
        throw std::runtime_error("Invalid concentration reference/order");
    if (order == Real(0)) return Real(0);
    if (concentration == Real(0))
    {
        if (order < Real(0))
            throw std::runtime_error("Negative kinetic order is singular at zero concentration");
        return -std::numeric_limits<Real>::infinity();
    }
    const Real value=order*(std::log(concentration)-std::log(reference));
    if (!std::isfinite(value))
        throw std::runtime_error("BV logarithmic concentration factor exceeds scalar range");
    return value;
}
// Direct partial at fixed cell value, gradient, geometry and material data.
template<class Real>
Real boundaryCurrentResponse
(
    Real nu, Real nue, Real faraday, Real diffusion, Real delta, Real beta
)
{
    const long double value=-static_cast<long double>(nu)
        /(static_cast<long double>(nue)*faraday*diffusion*(delta+beta));
    if (!std::isfinite(value) || std::abs(value)>std::numeric_limits<Real>::max())
        return std::numeric_limits<Real>::quiet_NaN();
    return static_cast<Real>(value);
}


struct DiagnosticBV
{
    long double ia, ic, current, slope;
};

template<class Real>
bool diagnosticBV
(
    Real i0, Real a, Real c, Real eta, Real logCa, Real logCc, DiagnosticBV& result
)
{
    const Real negInf=-std::numeric_limits<Real>::infinity();
    if (!std::isfinite(eta) || !std::isfinite(i0) || i0<=0
     || (!std::isfinite(logCa) && logCa!=negInf)
     || (!std::isfinite(logCc) && logCc!=negInf)) return false;
    const long double lp=std::log(static_cast<long double>(i0))+logCa
        +static_cast<long double>(a)*eta;
    const long double lc=std::log(static_cast<long double>(i0))+logCc
        -static_cast<long double>(c)*eta;
    const long double limit=std::log(std::numeric_limits<long double>::max())-8;
    if (lp>limit || lc>limit) return false;
    result.ia=std::exp(lp);
    result.ic=std::exp(lc);
    // expm1 retains accuracy in counterfactual states near branch cancellation.
    if (lp== -std::numeric_limits<long double>::infinity()
     && lc== -std::numeric_limits<long double>::infinity()) result.current=0;
    else if (lp>=lc) result.current=result.ia*(-std::expm1(lc-lp));
    else result.current=result.ic*std::expm1(lp-lc);
    result.slope=static_cast<long double>(a)*result.ia
        +static_cast<long double>(c)*result.ic;
    return std::isfinite(result.current) && std::isfinite(result.slope);
}

bool scalarDiagnostic(const long double value)
{
    return std::isfinite(value) && std::abs(value)<=std::numeric_limits<scalar>::max();
}

} // namespace electroChemConcentration
}

defineTypeNameAndDebug(concentrationButlerVolmer, 0);
addToRunTimeSelectionTable(currentModel, concentrationButlerVolmer, dictionary);

concentrationButlerVolmer::concentrationButlerVolmer
(
    fvMesh& mesh, const dictionary& currentDict, const volScalarField& phiE
)
:
    ButlerVolmer(mesh, currentDict, phiE, false, true),
    dependencyNames_(electrodeNames_.size()),
    referenceConcentrations_(electrodeNames_.size()),
    anodicOrders_(electrodeNames_.size()),
    cathodicOrders_(electrodeNames_.size()),
    speciesNu_(electrodeNames_.size()),
    diffusivity_(electrodeNames_.size()),
    chargeNumber_(electrodeNames_.size()),
    responseMethod_(mesh.solutionDict().subDict("Picard").subDict("current")
        .getOrDefault<word>("responseMethod", "boundaryTaylor")),
    responseDiagnostics_(mesh.solutionDict().subDict("Picard").subDict("current")
        .getOrDefault<Switch>("responseDiagnostics", false)),
    responseStatus_(electrodeNames_.size(), 0),
    denominatorMin_(electrodeNames_.size(), 1),
    denominatorMax_(electrodeNames_.size(), 1),
    responseStepScale_(electrodeNames_.size(), 1),
    potentialIncrement_(electrodeNames_.size(), 0),
    predictedFaces_(electrodeNames_.size()),
    predictionAbsError_(electrodeNames_.size(), 0),
    predictionRelError_(electrodeNames_.size(), 0),
    havePrediction_(false), predictionTested_(false),
    linearizationPhiS_(electrodeNames_.size(),0),
    beforeFaces_(electrodeNames_.size()),
    predictionPhiE_(electrodeNames_.size()),
    previousPhiE_(electrodeNames_.size()), previousCurrent_(electrodeNames_.size()),
    previousLogCa_(electrodeNames_.size()), previousLogCc_(electrodeNames_.size()),
    previousPhiS_(electrodeNames_.size(),0), havePreviousState_(false),
    diagnosticTimeIndex_(-1),
    predictedBVResidual_(electrodeNames_.size(),0),
    kineticPredictionError_(electrodeNames_.size(),0), tangentRemainder_(electrodeNames_.size(),0),
    potentialSensitivity_(electrodeNames_.size(),0), potentialTolerance_(electrodeNames_.size(),0),
    sweepPotentialEffect_(electrodeNames_.size(),0),
    electrolyteEffect_(electrodeNames_.size(),0), netEtaEffect_(electrodeNames_.size(),0),
    compositionEffect_(electrodeNames_.size(),0), currentStep_(electrodeNames_.size(),0),
    decompositionError_(electrodeNames_.size(),0),
    predictionCurrentValid_(electrodeNames_.size(),false),
    tangentValid_(electrodeNames_.size(),false), couplingValid_(electrodeNames_.size(),false)
{
    if (responseMethod_ != "none" && responseMethod_ != "boundaryTaylor")
        FatalIOErrorInFunction(currentDict_)
            << "current/responseMethod must be none or boundaryTaylor"
            << exit(FatalIOError);
    const IOdictionary& reactions=mesh.lookupObject<IOdictionary>("reactions");
    const wordList species(reactions.lookup("species"));
    const word inert(reactions.get<word>("inertSpecie"));
    forAll(electrodeNames_, e)
    {
        const dictionary& electrode=currentDict.subDict("electrodes").subDict(electrodeNames_[e]);
        // An empty dictionary is an explicit concentration-independent test.
        const dictionary& dependencies=electrode.subDict("concentrationDependence");
        dependencyNames_[e]=dependencies.sortedToc();
        referenceConcentrations_[e].setSize(dependencyNames_[e].size());
        anodicOrders_[e].setSize(dependencyNames_[e].size());
        cathodicOrders_[e].setSize(dependencyNames_[e].size());
        speciesNu_[e].setSize(dependencyNames_[e].size(), 0);
        diffusivity_[e].setSize(dependencyNames_[e].size(), 0);
        chargeNumber_[e].setSize(dependencyNames_[e].size(), 0);
        forAll(dependencyNames_[e], i)
        {
            const word& specie=dependencyNames_[e][i];
            const dictionary& entry=dependencies.subDict(specie);
            if (species.find(specie)<0)
                FatalIOErrorInFunction(entry) << "Unknown concentration species " << specie
                    << "; available: " << species << exit(FatalIOError);
            referenceConcentrations_[e][i]=entry.get<scalar>("c_ref");
            anodicOrders_[e][i]=entry.get<scalar>("anodicOrder");
            cathodicOrders_[e][i]=entry.get<scalar>("cathodicOrder");
            if (!std::isfinite(referenceConcentrations_[e][i])
             || referenceConcentrations_[e][i]<=0
             || !std::isfinite(anodicOrders_[e][i]) || !std::isfinite(cathodicOrders_[e][i]))
                FatalIOErrorInFunction(entry) << "Require positive finite c_ref and finite branch orders"
                    << exit(FatalIOError);
            if (specie==inert && (anodicOrders_[e][i]!=0 || cathodicOrders_[e][i]!=0))
                FatalIOErrorInFunction(entry)
                    << "Inert closure species " << specie << " is updated after Picard convergence. "
                    << "Use unit solvent activity (omit its dependency or set both orders zero)."
                    << exit(FatalIOError);
            const volScalarField& field=mesh.lookupObject<volScalarField>("c_"+specie);
            if (field.dimensions()!=dimMoles/dimVolume)
                FatalIOErrorInFunction(entry) << field.name() << " must have concentration dimensions"
                    << exit(FatalIOError);
            diffusivity_[e][i]=reactions.subDict("diffusivity").get<scalar>(specie);
            if (!std::isfinite(diffusivity_[e][i]) || diffusivity_[e][i]<=0)
                FatalIOErrorInFunction(entry) << "Require positive finite diffusivity for "
                    << specie << exit(FatalIOError);
            chargeNumber_[e][i]=reactions.found("chargeNumber")
                ? reactions.subDict("chargeNumber").getOrDefault<scalar>(specie,0) : 0;
            if (!std::isfinite(chargeNumber_[e][i]))
                FatalIOErrorInFunction(entry) << "Require finite chargeNumber for "
                    << specie << exit(FatalIOError);
            const fvPatchScalarField& patch=field.boundaryField()[patchIDs_[e]];
            if (patch.type()=="galvanostaticMolarFlux")
            {
                // Read existing metadata without changing or evaluating the BC.
                OStringStream output;
                output << token::BEGIN_BLOCK << nl;
                patch.write(output);
                output << token::END_BLOCK << nl;
                IStringStream input(output.str());
                const dictionary patchDict(input);
                speciesNu_[e][i]=patchDict.get<scalar>("speciesStoichCoeff");
                if (!std::isfinite(speciesNu_[e][i]))
                    FatalIOErrorInFunction(patchDict) << "Require finite speciesStoichCoeff"
                        << exit(FatalIOError);
            }

        }
    }
    // All derived state now exists; virtual factors participate in this root.
    correct();
}

void concentrationButlerVolmer::logConcentrationFactors
(
    const label e, scalarField& logCa, scalarField& logCc
) const
{
    ButlerVolmer::logConcentrationFactors(e,logCa,logCc);
    bool valid=true;
    word failedSpecies;
    forAll(dependencyNames_[e], i)
    {
        const volScalarField& field=mesh_.lookupObject<volScalarField>("c_"+dependencyNames_[e][i]);
        const fvPatchScalarField& patch=field.boundaryField()[patchIDs_[e]];
        // Use stored Picard face values: evaluating a BC here would change
        // the frozen state and could create a current/concentration recursion.
        try
        {
            forAll(patch, f)
            {
                logCa[f]+=electroChemConcentration::logFactor
                    (patch[f],referenceConcentrations_[e][i],anodicOrders_[e][i]);
                logCc[f]+=electroChemConcentration::logFactor
                    (patch[f],referenceConcentrations_[e][i],cathodicOrders_[e][i]);
            }
        }
        catch (const std::exception&)
        {
            valid=false;
            failedSpecies=dependencyNames_[e][i];
        }
    }
    // Reject invalid data collectively before the controller's global sums.
    reduce(valid,andOp<bool>());
    if (!valid)
        FatalErrorInFunction << "Invalid surface concentration or singular order on electrode "
            << electrodeNames_[e] << "; local species=" << failedSpecies << nl
            << "Require non-negative finite concentrations; negative orders need positive concentrations."
            << exit(FatalError);
}


scalar concentrationButlerVolmer::adjustCurrentCandidate
(
    const label e, const scalarField& baseCurrent, scalarField& candidate,
    const scalarField& logCa, const scalarField& logCc, bool& historyValid
)
{
    if (responseDiagnostics_) linearizationPhiS_[e]=phiS_[e];
    denominatorMin_[e]=1;
    denominatorMax_[e]=1;
    responseStepScale_[e]=1;
    potentialIncrement_[e]=0;
    const auto setStatus=[&](const label status)
    {
        if (responseStatus_[e]!=status) historyValid=false;
        responseStatus_[e]=status;
    };
    // Constructor and zero-factor startup use the existing robust controller.
    if (responseMethod_=="none" || !initialized_ || dependencyNames_[e].empty())
    {
        setStatus(0);
        return 0;
    }

    const label patchi=patchIDs_[e];
    const fvPatchScalarField& phiPatch=phiE_.boundaryField()[patchi];
    const fixedGradientFvPatchScalarField& phiBC=
        refCast<const fixedGradientFvPatchScalarField>(phiPatch);
    // Match the stored gradient consumed by galvanostaticMolarFlux exactly.
    const scalarField& gradient=phiBC.gradient();
    const scalarField& delta=mesh_.boundary()[patchi].deltaCoeffs();
    const scalarField& areas=mesh_.magSf().boundaryField()[patchi];
    const scalar a=alphaA_[e]*F_/(R_*T_);
    const scalar c=alphaC_[e]*F_/(R_*T_);
    scalarField gEta(candidate.size(),0), diagonal(candidate.size(),1);
    scalarField compositionOffset(candidate.size(),0);
    scalarField branchA(candidate.size(),0), branchC(candidate.size(),0);
    List<scalarField> anchorConcentration(dependencyNames_[e].size());
    List<scalarField> dcDj(dependencyNames_[e].size());
    bool valid=true;

    forAll(candidate,f)
    {
        const scalar eta=phiS_[e]-phiPatch[f]-Ebase_[e];
        const scalar ia=std::exp(std::log(i0_[e])+logCa[f]+a*eta);
        const scalar ic=std::exp(std::log(i0_[e])+logCc[f]-c*eta);
        branchA[f]=ia;
        branchC[f]=ic;
        const long double slope=static_cast<long double>(a)*ia
            +static_cast<long double>(c)*ic;
        if (!std::isfinite(slope) || slope<=0
         || slope>std::numeric_limits<scalar>::max()) valid=false;
        else gEta[f]=static_cast<scalar>(slope);
    }
    forAll(dependencyNames_[e],m)
    {
        dcDj[m].setSize(candidate.size(),0);
        anchorConcentration[m].setSize(candidate.size(),0);
        const scalar pa=anodicOrders_[e][m], pc=cathodicOrders_[e][m];
        if (pa==0 && pc==0) continue;
        const volScalarField& field=mesh_.lookupObject<volScalarField>("c_"+dependencyNames_[e][m]);
        const fvPatchScalarField& patch=field.boundaryField()[patchi];
        const scalarField cellValues(patch.patchInternalField());
        // Fixed-value and zeroGradient dependencies have zero direct current
        // partial. Their bulk/potential response is still resolved by Picard.
        if (patch.type()!="galvanostaticMolarFlux"
         && patch.type()!="zeroGradient" && !patch.fixesValue()) valid=false;
        forAll(candidate,f)
        {
            if (!std::isfinite(patch[f]) || patch[f]<=0)
            {
                valid=false; // No quotient derivative at zero; no artificial floor.
                continue;
            }
            anchorConcentration[m][f]=patch.type()=="zeroGradient"
                ? cellValues[f] : patch[f];
            if (patch.type()=="galvanostaticMolarFlux")
            {
                const scalar beta=chargeNumber_[e][m]*F_/(R_*T_)*gradient[f];
                const scalar denom=delta[f]+beta;
                if (!std::isfinite(denom) || denom<=SMALL*max(mag(delta[f]),scalar(1)))
                {
                    valid=false;
                    continue;
                }
                dcDj[m][f]=electroChemConcentration::boundaryCurrentResponse
                (
                    speciesNu_[e][m],nuE_[e],F_,diffusivity_[e][m],delta[f],beta
                );
                // Account for the known change in gradient since the last
                // species sweep: stored face values need not satisfy the new BC.
                anchorConcentration[m][f]=delta[f]*cellValues[f]/denom
                    +dcDj[m][f]*baseCurrent[f];
            }
            if (!std::isfinite(anchorConcentration[m][f]) || anchorConcentration[m][f]<0)
            {
                valid=false;
                continue;
            }
            const scalar ia=branchA[f], ic=branchC[f];
            const long double gC=(static_cast<long double>(pa)*ia
                -static_cast<long double>(pc)*ic)/patch[f];
            const long double nextDiagonal=diagonal[f]+gC*dcDj[m][f];
            const long double nextOffset=compositionOffset[f]
                +gC*(anchorConcentration[m][f]-patch[f]);
            if (!std::isfinite(nextDiagonal) || !std::isfinite(nextOffset)
             || std::abs(nextDiagonal)>std::numeric_limits<scalar>::max()
             || std::abs(nextOffset)>std::numeric_limits<scalar>::max()) valid=false;
            else
            {
                diagonal[f]=static_cast<scalar>(nextDiagonal);
                compositionOffset[f]=static_cast<scalar>(nextOffset);
            }
            if (!std::isfinite(dcDj[m][f])) valid=false;
        }
    }
    forAll(diagonal,f)
        if (!std::isfinite(diagonal[f]) || diagonal[f]<=SMALL) valid=false;
    reduce(valid,andOp<bool>());
    if (!valid)
    {
        setStatus(2);
        return 0; // Whole-electrode fallback preserves the base current constraint.
    }

    scalar numerator=0, denominator=0;
    scalar dmin=GREAT,dmax=-GREAT;
    forAll(candidate,f)
    {
        numerator+=areas[f]*(candidate[f]-baseCurrent[f]-compositionOffset[f])/diagonal[f];
        denominator+=areas[f]*gEta[f]/diagonal[f];
        dmin=min(dmin,diagonal[f]);
        dmax=max(dmax,diagonal[f]);
    }
    reduce(numerator,sumOp<scalar>());
    reduce(denominator,sumOp<scalar>());
    reduce(dmin,minOp<scalar>());
    reduce(dmax,maxOp<scalar>());
    if (!std::isfinite(numerator) || !std::isfinite(denominator) || denominator<=VSMALL)
    {
        setStatus(2);
        return 0;
    }
    const scalar dPhi=numerator/denominator;
    scalarField step(candidate.size(),0);
    scalar scale=1;
    if (!std::isfinite(dPhi)) valid=false;
    // Limit the additional tangent-potential shift to half an exponential unit.
    // This is a numerical safeguard, not a small-overpotential physical model.
    if (std::isfinite(dPhi) && max(a,c)*mag(dPhi)>0.5)
        scale=min(scale,scalar(0.5)/(max(a,c)*mag(dPhi)));
    forAll(step,f)
    {
        step[f]=(candidate[f]-baseCurrent[f]-compositionOffset[f]-gEta[f]*dPhi)/diagonal[f];
        if (!std::isfinite(step[f])) valid=false;
    }
    // One common electrode factor preserves the integral constraint.
    // It controls only the predicted direct response, not the full transport.
    forAll(dependencyNames_[e],m)
    {
        if (anodicOrders_[e][m]==0 && cathodicOrders_[e][m]==0) continue;
        forAll(step,f)
        {
            const scalar dc=dcDj[m][f]*step[f];
            if (!std::isfinite(dc)) valid=false;
            else if (dc<0) scale=min(scale,scalar(0.95)*anchorConcentration[m][f]/(-dc));
        }
    }
    reduce(valid,andOp<bool>());
    reduce(scale,minOp<scalar>());
    if (!valid || !std::isfinite(scale) || scale<=0)
    {
        setStatus(2);
        return 0;
    }
    forAll(candidate,f) candidate[f]=baseCurrent[f]+scale*step[f];
    denominatorMin_[e]=dmin;
    denominatorMax_[e]=dmax;
    responseStepScale_[e]=scale;
    potentialIncrement_[e]=scale*dPhi;
    setStatus(1);
    if (scale<1) historyValid=false;
    return potentialIncrement_[e];
}


void concentrationButlerVolmer::correct()
{
    if (responseDiagnostics_ && diagnosticTimeIndex_!=mesh_.time().timeIndex())
    {
        havePreviousState_=false;
        diagnosticTimeIndex_=mesh_.time().timeIndex();
    }
    ButlerVolmer::correct();
    predictionTested_=false;
    havePrediction_=responseMethod_=="boundaryTaylor";
    if (!havePrediction_) return;
    // Predict the NEXT face values at frozen current cell concentrations and
    // stored potential gradient, using the ACTUALLY relaxed new face current.
    // No BC evaluation and no physical field mutation occurs here.
    forAll(electrodeNames_,e)
    {
        const label patchi=patchIDs_[e];
        const fixedGradientFvPatchScalarField& phiBC=
            refCast<const fixedGradientFvPatchScalarField>(phiE_.boundaryField()[patchi]);
        const scalarField& delta=mesh_.boundary()[patchi].deltaCoeffs();
        const fvsPatchScalarField& jn=electrodeCurrentDensity_.boundaryField()[patchi];
        if (responseDiagnostics_) predictionPhiE_[e]=scalarField(phiE_.boundaryField()[patchi]);
        predictedFaces_[e].setSize(dependencyNames_[e].size());
        if (responseDiagnostics_) beforeFaces_[e].setSize(dependencyNames_[e].size());
        forAll(dependencyNames_[e],m)
        {
            const volScalarField& field=mesh_.lookupObject<volScalarField>("c_"+dependencyNames_[e][m]);
            const fvPatchScalarField& patch=field.boundaryField()[patchi];
            const scalarField cells(patch.patchInternalField());
            predictedFaces_[e][m]=scalarField(patch);
            if (responseDiagnostics_) beforeFaces_[e][m]=scalarField(patch);
            if (patch.type()=="zeroGradient") predictedFaces_[e][m]=cells;
            if (patch.type()!="galvanostaticMolarFlux") continue;
            forAll(jn,f)
            {
                const scalar beta=chargeNumber_[e][m]*F_/(R_*T_)*phiBC.gradient()[f];
                const scalar denom=delta[f]+beta;
                // The physical BC owns singularity errors. Do not create an
                // invalid diagnostic division before the ordinary transport.
                if (!std::isfinite(denom) || mag(denom)<=SMALL*max(mag(delta[f]),scalar(1)))
                {
                    havePrediction_=false;
                    continue;
                }
                predictedFaces_[e][m][f]=delta[f]*cells[f]/denom
                    +electroChemConcentration::boundaryCurrentResponse
                    (speciesNu_[e][m],nuE_[e],F_,diffusivity_[e][m],delta[f],beta)*jn[f];
                if (!std::isfinite(predictedFaces_[e][m][f])) havePrediction_=false;
            }
        }
    }
    reduce(havePrediction_,andOp<bool>());
}

void concentrationButlerVolmer::evaluateConstraint()
{
    ButlerVolmer::evaluateConstraint();
    if (!havePrediction_)
    {
        if (responseDiagnostics_) evaluateResponseDiagnostics();
        return;
    }
    forAll(electrodeNames_,e)
    {
        predictionAbsError_[e]=0;
        predictionRelError_[e]=0;
        forAll(dependencyNames_[e],m)
        {
            if (anodicOrders_[e][m]==0 && cathodicOrders_[e][m]==0) continue;
            const volScalarField& field=mesh_.lookupObject<volScalarField>("c_"+dependencyNames_[e][m]);
            const fvPatchScalarField& actual=field.boundaryField()[patchIDs_[e]];
            forAll(actual,f)
            {
                const scalar predicted=predictedFaces_[e][m][f];
                const scalar error=mag(actual[f]-predicted);
                const scalar scale=max(referenceConcentrations_[e][m],max(mag(actual[f]),mag(predicted)));
                predictionAbsError_[e]=max(predictionAbsError_[e],error);
                predictionRelError_[e]=max(predictionRelError_[e],error/scale);
            }
        }
        reduce(predictionAbsError_[e],maxOp<scalar>());
        reduce(predictionRelError_[e],maxOp<scalar>());
    }
    predictionTested_=true;
    if (responseDiagnostics_) evaluateResponseDiagnostics();
}

void concentrationButlerVolmer::evaluateResponseDiagnostics()
{
    using electroChemConcentration::DiagnosticBV;
    using electroChemConcentration::diagnosticBV;
    using electroChemConcentration::scalarDiagnostic;
    const auto maxMagnitude=[](const long double value, scalar& maximum, bool& valid)
    {
        if (!scalarDiagnostic(value)) valid=false;
        else maximum=max(maximum,scalar(std::abs(value)));
    };
    forAll(electrodeNames_,e)
    {
        const label patchi=patchIDs_[e];
        const fvPatchScalarField& phi=phiE_.boundaryField()[patchi];
        const fvsPatchScalarField& jn=electrodeCurrentDensity_.boundaryField()[patchi];
        const scalar a=alphaA_[e]*F_/(R_*T_), c=alphaC_[e]*F_/(R_*T_);
        scalarField actualA,actualC;
        logConcentrationFactors(e,actualA,actualC);
        scalarField predictedA(phi.size(),0),predictedC(phi.size(),0);
        scalarField beforeA(phi.size(),0),beforeC(phi.size(),0);
        bool predictionOK=havePrediction_ && predictionPhiE_[e].size()==phi.size();
        bool tangentOK=predictionOK && responseStatus_[e]==1;
        bool couplingOK=havePreviousState_ && previousPhiE_[e].size()==phi.size()
            && previousCurrent_[e].size()==phi.size()
            && previousLogCa_[e].size()==phi.size() && previousLogCc_[e].size()==phi.size();
        forAll(dependencyNames_[e],m)
        {
            if (anodicOrders_[e][m]==0 && cathodicOrders_[e][m]==0) continue;
            forAll(phi,f)
            {
                if (predictionOK)
                {
                    try
                    {
                        predictedA[f]+=electroChemConcentration::logFactor
                            (predictedFaces_[e][m][f],referenceConcentrations_[e][m],anodicOrders_[e][m]);
                        predictedC[f]+=electroChemConcentration::logFactor
                            (predictedFaces_[e][m][f],referenceConcentrations_[e][m],cathodicOrders_[e][m]);
                    }
                    catch (const std::exception&) { predictionOK=false; tangentOK=false; }
                }
                if (tangentOK)
                {
                    const scalar old=beforeFaces_[e][m][f];
                    if (!std::isfinite(old) || old<=0) { tangentOK=false; continue; }
                    beforeA[f]+=electroChemConcentration::logFactor
                        (old,referenceConcentrations_[e][m],anodicOrders_[e][m]);
                    beforeC[f]+=electroChemConcentration::logFactor
                        (old,referenceConcentrations_[e][m],cathodicOrders_[e][m]);
                }
            }
        }
        predictedBVResidual_[e]=0; kineticPredictionError_[e]=0;
        tangentRemainder_[e]=0; potentialSensitivity_[e]=0;
        potentialTolerance_[e]=GREAT; sweepPotentialEffect_[e]=0;
        electrolyteEffect_[e]=0; netEtaEffect_[e]=0;
        compositionEffect_[e]=0; currentStep_[e]=0; decompositionError_[e]=0;
        bool sensitivityOK=true;
        forAll(phi,f)
        {
            const scalar eta=phiS_[e]-phi[f]-Ebase_[e];
            DiagnosticBV actual;
            if (!diagnosticBV(i0_[e],a,c,eta,actualA[f],actualC[f],actual))
            {
                predictionOK=false; tangentOK=false; couplingOK=false; sensitivityOK=false;
                continue;
            }
            maxMagnitude(actual.slope,potentialSensitivity_[e],sensitivityOK);
            if (actual.slope>0)
            {
                const long double voltage=currentTolerance(jn[f])/actual.slope;
                if (scalarDiagnostic(voltage)) potentialTolerance_[e]=min(potentialTolerance_[e],scalar(voltage));
                else sensitivityOK=false;
            }
            if (predictionOK)
            {
                const scalar etaPred=phiS_[e]-predictionPhiE_[e][f]-Ebase_[e];
                DiagnosticBV predicted,actualAtPredictedEta;
                if (!diagnosticBV(i0_[e],a,c,etaPred,predictedA[f],predictedC[f],predicted)
                 || !diagnosticBV(i0_[e],a,c,etaPred,actualA[f],actualC[f],actualAtPredictedEta))
                { predictionOK=false; tangentOK=false; }
                else
                {
                    maxMagnitude(jn[f]+predicted.current,predictedBVResidual_[e],predictionOK);
                    maxMagnitude(actualAtPredictedEta.current-predicted.current,kineticPredictionError_[e],predictionOK);
                    maxMagnitude(actual.current-actualAtPredictedEta.current,sweepPotentialEffect_[e],predictionOK);
                    if (tangentOK)
                    {
                        DiagnosticBV root;
                        const scalar etaRoot=linearizationPhiS_[e]-predictionPhiE_[e][f]-Ebase_[e];
                        if (!diagnosticBV(i0_[e],a,c,etaRoot,beforeA[f],beforeC[f],root)) tangentOK=false;
                        else
                        {
                            long double tangent=root.current+root.slope*(phiS_[e]-linearizationPhiS_[e]);
                            forAll(dependencyNames_[e],m)
                            {
                                const scalar pa=anodicOrders_[e][m],pc=cathodicOrders_[e][m];
                                if (pa==0 && pc==0) continue;
                                const long double old=beforeFaces_[e][m][f];
                                tangent+=(pa*root.ia-pc*root.ic)/old*(predictedFaces_[e][m][f]-old);
                            }
                            maxMagnitude(predicted.current-tangent,tangentRemainder_[e],tangentOK);
                        }
                    }
                }
            }
            if (couplingOK)
            {
                DiagnosticBV old,crossEta,electrolyte;
                const scalar oldEta=previousPhiS_[e]-previousPhiE_[e][f]-Ebase_[e];
                const scalar electrolyteEta=previousPhiS_[e]-phi[f]-Ebase_[e];
                if (!diagnosticBV(i0_[e],a,c,oldEta,previousLogCa_[e][f],previousLogCc_[e][f],old)
                 || !diagnosticBV(i0_[e],a,c,eta,previousLogCa_[e][f],previousLogCc_[e][f],crossEta)
                 || !diagnosticBV(i0_[e],a,c,electrolyteEta,previousLogCa_[e][f],previousLogCc_[e][f],electrolyte))
                    couplingOK=false;
                else
                {
                    const long double dj=static_cast<long double>(jn[f])-previousCurrent_[e][f];
                    const long double potential=crossEta.current-old.current;
                    const long double composition=actual.current-crossEta.current;
                    const long double residualChange=(jn[f]+actual.current)-(previousCurrent_[e][f]+old.current);
                    maxMagnitude(dj,currentStep_[e],couplingOK);
                    maxMagnitude(electrolyte.current-old.current,electrolyteEffect_[e],couplingOK);
                    maxMagnitude(potential,netEtaEffect_[e],couplingOK);
                    maxMagnitude(composition,compositionEffect_[e],couplingOK);
                    maxMagnitude(residualChange-dj-potential-composition,decompositionError_[e],couplingOK);
                }
            }
        }
        reduce(predictionOK,andOp<bool>()); reduce(tangentOK,andOp<bool>());
        reduce(couplingOK,andOp<bool>()); reduce(sensitivityOK,andOp<bool>());
        reduce(predictedBVResidual_[e],maxOp<scalar>());
        reduce(kineticPredictionError_[e],maxOp<scalar>());
        reduce(tangentRemainder_[e],maxOp<scalar>());
        reduce(sweepPotentialEffect_[e],maxOp<scalar>());
        reduce(potentialSensitivity_[e],maxOp<scalar>());
        reduce(potentialTolerance_[e],minOp<scalar>());
        reduce(electrolyteEffect_[e],maxOp<scalar>());
        reduce(netEtaEffect_[e],maxOp<scalar>());
        reduce(compositionEffect_[e],maxOp<scalar>());
        reduce(currentStep_[e],maxOp<scalar>()); reduce(decompositionError_[e],maxOp<scalar>());
        predictionCurrentValid_[e]=predictionOK;
        tangentValid_[e]=tangentOK && predictionOK;
        couplingValid_[e]=couplingOK;
        if (!sensitivityOK) { potentialSensitivity_[e]=-1; potentialTolerance_[e]=-1; }
        previousPhiE_[e]=scalarField(phi); previousCurrent_[e]=scalarField(jn);
        previousLogCa_[e]=actualA; previousLogCc_[e]=actualC; previousPhiS_[e]=phiS_[e];
    }
    havePreviousState_=true;
}

void concentrationButlerVolmer::report() const
{
    ButlerVolmer::report();
    if (responseMethod_=="none") return;
    forAll(electrodeNames_,e)
    {
        const char* status=responseStatus_[e]==1 ? "active"
            : responseStatus_[e]==2 ? "fallback" : "startup/identity";
        Info<< "    boundaryTaylor " << electrodeNames_[e] << ": " << status;
        if (responseStatus_[e]==1)
            Info<< ", diagonal min/max=" << denominatorMin_[e] << '/' << denominatorMax_[e]
                << ", stepScale=" << responseStepScale_[e]
                << ", proposed dPhiS=" << potentialIncrement_[e] << " V";
        if (predictionTested_)
            Info<< ", max face prediction error=" << predictionAbsError_[e]
                << " mol/m3, scaled=" << predictionRelError_[e];
        Info<< nl;
        if (!responseDiagnostics_ || !evaluated_) continue;
        Info<< "    BV prediction " << electrodeNames_[e] << ": ";
        if (predictionCurrentValid_[e])
        {
            Info<< "predicted max residual=" << predictedBVResidual_[e]
                << " A/m2, concentration error=" << kineticPredictionError_[e]
                << " A/m2, sweep potential effect=" << sweepPotentialEffect_[e] << " A/m2";
            if (tangentValid_[e]) Info<< ", Taylor remainder=" << tangentRemainder_[e] << " A/m2";
            else Info<< ", Taylor remainder=n/a";
        }
        else Info<< "prediction=n/a";
        if (potentialSensitivity_[e]>=0)
        {
            Info<< ", max gEta=" << potentialSensitivity_[e] << " A/m2/V";
            if (potentialSensitivity_[e]>0)
                Info<< ", min voltage allowance=" << potentialTolerance_[e] << " V";
        }
        else Info<< ", sensitivity=n/a";
        Info<< nl;
        Info<< "    BV coupling " << electrodeNames_[e] << ": ";
        if (couplingValid_[e])
            Info<< "electrolyte effect=" << electrolyteEffect_[e]
                << " A/m2, net eta effect=" << netEtaEffect_[e]
                << " A/m2, composition effect=" << compositionEffect_[e]
                << " A/m2, current step=" << currentStep_[e]
                << " A/m2, closure=" << decompositionError_[e] << " A/m2";
        else Info<< "previous tested state unavailable";
        Info<< nl;
    }
}

} // namespace Foam

