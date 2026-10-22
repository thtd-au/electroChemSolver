#include "concentrationButlerVolmer.H"
#include "potentialManager.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"
#include "dimensionSets.H"
#include "fixedGradientFvPatchFields.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "token.H"
#include "mixedFvPatchFields.H"
#include "fvmLaplacian.H"
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


// Validate extended-precision intermediate values before scalar conversion.
bool scalarRepresentable(const long double value)
{
    return std::isfinite(value) && std::abs(value)<=std::numeric_limits<scalar>::max();
}

// Mixed coefficients for k*L*snGrad(psi)+g*psi=-h (zero metal increment).
template<class Real>
bool potentialRobinCoefficients
(
    Real k, Real L, Real delta, Real g, Real h, Real minimumFraction,
    Real& fraction, Real& value
)
{
    if (!std::isfinite(k) || !std::isfinite(L) || !std::isfinite(delta)
     || !std::isfinite(g) || !std::isfinite(h) || k<=0 || L<=0 || delta<=0 || g<=0)
        return false;
    const long double denom=static_cast<long double>(k)*L*delta+g;
    const long double w=static_cast<long double>(g)/denom;
    const long double v=-static_cast<long double>(h)/g;
    if (!std::isfinite(w) || !std::isfinite(v)
     || w<=minimumFraction || w>1 || std::abs(v)>std::numeric_limits<Real>::max())
        return false;
    fraction=static_cast<Real>(w);
    value=static_cast<Real>(v);
    return true;
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
    responseStatus_(electrodeNames_.size(), 0),
    responseStepScale_(electrodeNames_.size(), 1),
    potentialDiagonal_(electrodeNames_.size()), potentialSlope_(electrodeNames_.size()),
    potentialRhs_(electrodeNames_.size()), potentialAnchor_(electrodeNames_.size()),
    potentialDcDj_(electrodeNames_.size()), potentialResponseStatus_(0),
    potentialSolveIterations_(0), potentialCurrentDefect_(0), potentialRowResidual_(0), potentialStepScale_(1),
    predictedPotentialMax_(electrodeNames_.size(),0)
{
    if (responseMethod_ != "none" && responseMethod_ != "boundaryTaylor"
     && responseMethod_ != "boundaryPotential")
        FatalIOErrorInFunction(currentDict_)
            << "current/responseMethod must be none, boundaryTaylor or boundaryPotential"
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
    responseStepScale_[e]=1;
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
    forAll(candidate,f)
    {
        numerator+=areas[f]*(candidate[f]-baseCurrent[f]-compositionOffset[f])/diagonal[f];
        denominator+=areas[f]*gEta[f]/diagonal[f];
    }
    reduce(numerator,sumOp<scalar>());
    reduce(denominator,sumOp<scalar>());
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
    if (responseMethod_=="boundaryPotential")
    {
        potentialDiagonal_[e]=diagonal;
        potentialSlope_[e]=gEta;
        potentialRhs_[e]=candidate-baseCurrent-compositionOffset;
        potentialAnchor_[e]=anchorConcentration;
        potentialDcDj_[e]=dcDj;
    }
    forAll(candidate,f) candidate[f]=baseCurrent[f]+scale*step[f];
    responseStepScale_[e]=scale;
    setStatus(1);
    if (scale<1) historyValid=false;
    return scale*dPhi;
}


void concentrationButlerVolmer::adjustCurrentCandidates
(
    const List<scalarField>& baseCurrent, List<scalarField>& candidates,
    scalarList& increments, bool& historyValid
)
{
    const label oldStatus=potentialResponseStatus_;
    potentialResponseStatus_=0;
    potentialSolveIterations_=0;
    potentialCurrentDefect_=0;
    potentialRowResidual_=0;
    potentialStepScale_=1;
    predictedPotentialMax_=scalar(0);
    const auto finish=[&]()
    {
        if (oldStatus!=potentialResponseStatus_) historyValid=false;
    };
    if (responseMethod_!="boundaryPotential" || !initialized_)
    { finish(); return; }
    // Local boundaryTaylor candidates are already available for safe fallback.
    potentialResponseStatus_=2;
    if (Pstream::parRun() || electrodeNames_.size()!=2
     || !mesh_.foundObject<surfaceScalarField>("kappaFace"))
    { finish(); return; }
    forAll(electrodeNames_,e)
        if (responseStatus_[e]!=1) { finish(); return; }

    // This two-solve Schur construction uses the all-flux boundary topology
    // of the supplied case. A fixed-potential reservoir needs a larger Schur
    // system; retain boundaryTaylor rather than impose an incorrect closure.
    forAll(phiE_.boundaryField(),patchi)
    {
        const fvPatchScalarField& p=phiE_.boundaryField()[patchi];
        if (!isA<fixedGradientFvPatchScalarField>(p) && p.type()!="zeroGradient")
        { finish(); return; }
    }
    const surfaceScalarField& conductivity=
        mesh_.lookupObject<surfaceScalarField>("kappaFace");
    bool valid=true;
    forAll(conductivity.primitiveField(),f)
        if (!std::isfinite(conductivity.primitiveField()[f])
         || conductivity.primitiveField()[f]<=0) valid=false;
    forAll(conductivity.boundaryField(),p)
        forAll(conductivity.boundaryField()[p],f)
            if (!std::isfinite(conductivity.boundaryField()[p][f])
             || conductivity.boundaryField()[p][f]<=0) valid=false;
    if (!valid) { finish(); return; }

    const IOdictionary& ec=mesh_.lookupObject<IOdictionary>("electroChemistry");
    const word reference=ec.subDict("potential").get<word>("referencePatch");
    const label refPatch=mesh_.boundaryMesh().findPatchID(reference);
    if (refPatch<0 || mesh_.boundary()[refPatch].faceCells().empty())
    { finish(); return; }
    const label refCell=mesh_.boundary()[refPatch].faceCells()[0];
    wordList types(mesh_.boundary().size(),"fixedGradient");
    forAll(patchIDs_,e) types[patchIDs_[e]]="mixed";
    volScalarField response0
    (
        IOobject("bvPotentialResponse0",mesh_.time().timeName(),mesh_,
            IOobject::NO_READ,IOobject::NO_WRITE,false),
        mesh_,dimensionedScalar("zero",phiE_.dimensions(),0),types
    );
    volScalarField response1
    (
        IOobject("bvPotentialResponse1",mesh_.time().timeName(),mesh_,
            IOobject::NO_READ,IOobject::NO_WRITE,false),
        mesh_,dimensionedScalar("zero",phiE_.dimensions(),0),types
    );
    forAll(types,patchi)
    {
        if (types[patchi]=="fixedGradient")
        {
            refCast<fixedGradientFvPatchScalarField>(response0.boundaryFieldRef()[patchi]).gradient()=0;
            refCast<fixedGradientFvPatchScalarField>(response1.boundaryFieldRef()[patchi]).gradient()=0;
        }
    }
    forAll(electrodeNames_,e)
    {
        const label patchi=patchIDs_[e];
        const scalarField& delta=mesh_.boundary()[patchi].deltaCoeffs();
        mixedFvPatchScalarField& p0=
            refCast<mixedFvPatchScalarField>(response0.boundaryFieldRef()[patchi]);
        mixedFvPatchScalarField& p1=
            refCast<mixedFvPatchScalarField>(response1.boundaryFieldRef()[patchi]);
        p0.refGrad()=0; p1.refGrad()=0;
        forAll(delta,f)
        {
            const scalar g=potentialSlope_[e][f],L=potentialDiagonal_[e][f];
            const scalar k=conductivity.boundaryField()[patchi][f];
            scalar fraction=0,value=0;
            if (!electroChemConcentration::potentialRobinCoefficients
                (k,L,delta[f],g,potentialRhs_[e][f],scalar(SMALL),fraction,value))
            { valid=false; continue; }
            p0.valueFraction()[f]=fraction;
            p1.valueFraction()[f]=fraction;
            p0.refValue()[f]=value;
            // Basis for a +1 V increment of electrode 1's metal potential.
            p1.refValue()[f]=e==1 ? scalar(1) : scalar(0);
        }
    }
    if (!valid) { finish(); return; }
    response0.correctBoundaryConditions();
    response1.correctBoundaryConditions();
    fvScalarMatrix eqn0(fvm::laplacian(conductivity,response0));
    dictionary controls(eqn0.solverDict(phiE_.name()));
    controls.set("log",label(0));
    controls.set("relTol",scalar(0));
    controls.set("tolerance",min(controls.get<scalar>("tolerance"),scalar(1e-12)));
    controls.set("minIter",label(1));
    const auto perf0=eqn0.solve(controls);
    potentialSolveIterations_+=perf0.nIterations();
    if (!perf0.converged()) { finish(); return; }
    response0.correctBoundaryConditions();
    fvScalarMatrix eqn1(fvm::laplacian(conductivity,response1));
    const auto perf1=eqn1.solve(controls);
    potentialSolveIterations_+=perf1.nIterations();
    if (!perf1.converged()) { finish(); return; }
    response1.correctBoundaryConditions();

    const label constrainedPatch=patchIDs_[1];
    const scalarField& constraintAreas=mesh_.magSf().boundaryField()[constrainedPatch];
    const scalarField& constraintK=conductivity.boundaryField()[constrainedPatch];
    const scalarField grad0(response0.boundaryField()[constrainedPatch].snGrad());
    const scalarField grad1(response1.boundaryField()[constrainedPatch].snGrad());
    const scalar I0=-gSum(constraintAreas*constraintK*grad0);
    const scalar I1=-gSum(constraintAreas*constraintK*grad1);
    if (!std::isfinite(I0) || !std::isfinite(I1) || mag(I1)<=VSMALL)
    { finish(); return; }
    const long double ratio=-static_cast<long double>(I0)/I1;
    if (!electroChemConcentration::scalarRepresentable(ratio)) { finish(); return; }
    const scalar metal1=scalar(ratio);
    const scalar gauge=response0[refCell]+metal1*response1[refCell];
    if (!std::isfinite(gauge)) { finish(); return; }
    scalarList metal(electrodeNames_.size(),-gauge);
    metal[1]+=metal1;
    List<scalarField> step(electrodeNames_.size());
    scalar scale=1;
    forAll(electrodeNames_,e)
    {
        const label patchi=patchIDs_[e];
        const scalarField& areas=mesh_.magSf().boundaryField()[patchi];
        const scalarField& k=conductivity.boundaryField()[patchi];
        const scalarField G0(response0.boundaryField()[patchi].snGrad());
        const scalarField G1(response1.boundaryField()[patchi].snGrad());
        const fvPatchScalarField& p0=response0.boundaryField()[patchi];
        const fvPatchScalarField& p1=response1.boundaryField()[patchi];
        step[e].setSize(areas.size(),0);
        const scalar exponent=max(alphaA_[e],alphaC_[e])*F_/(R_*T_);
        if (!std::isfinite(metal[e])) valid=false;
        else if (exponent*mag(metal[e])>0.5)
            scale=min(scale,scalar(0.5)/(exponent*mag(metal[e])));
        forAll(areas,f)
        {
            const long double dj=-static_cast<long double>(k[f])*(G0[f]
                +static_cast<long double>(metal1)*G1[f]);
            const long double psi=p0[f]+static_cast<long double>(metal1)*p1[f]-gauge;
            if (!electroChemConcentration::scalarRepresentable(dj)
             || !electroChemConcentration::scalarRepresentable(psi))
            { valid=false; continue; }
            step[e][f]=scalar(dj);
            const long double row=static_cast<long double>(potentialDiagonal_[e][f])*dj
                +static_cast<long double>(potentialSlope_[e][f])*(metal[e]-psi)
                -potentialRhs_[e][f];
            if (!electroChemConcentration::scalarRepresentable(row)) valid=false;
            else
            {
                potentialRowResidual_=max(potentialRowResidual_,scalar(std::abs(row)));
                if (std::abs(row)>scalar(0.01)*currentTolerance(baseCurrent[e][f])) valid=false;
            }
            predictedPotentialMax_[e]=max(predictedPotentialMax_[e],scalar(std::abs(psi)));
            const long double etaChange=metal[e]-psi;
            if (exponent*std::abs(etaChange)>0.5)
                scale=min(scale,scalar(0.5L/(exponent*std::abs(etaChange))));
        }
        const scalar defect=gSum(areas*step[e]);
        potentialCurrentDefect_=max(potentialCurrentDefect_,mag(defect));
        const scalar currentScale=max(appliedCurrentDensity_*area_[e],gSum(areas*mag(step[e])));
        const scalar allowed=max(rootAbsTol_*area_[e],rootRelTol_*currentScale);
        if (!std::isfinite(defect) || mag(defect)>allowed) valid=false;
        // Roundoff projection only; a larger conservation defect triggers fallback.
        if (std::isfinite(defect)) step[e]-=defect/area_[e];
        forAll(dependencyNames_[e],m)
        {
            if (anodicOrders_[e][m]==0 && cathodicOrders_[e][m]==0) continue;
            forAll(step[e],f)
            {
                const long double dc=static_cast<long double>(potentialDcDj_[e][m][f])*step[e][f];
                if (!electroChemConcentration::scalarRepresentable(dc)) valid=false;
                else if (dc<0)
                    scale=min(scale,scalar(0.95L*potentialAnchor_[e][m][f]/(-dc)));
            }
        }
    }
    if (!valid || !std::isfinite(scale) || scale<=0) { finish(); return; }
    forAll(electrodeNames_,e)
    {
        candidates[e]=baseCurrent[e]+scale*step[e];
        increments[e]=scale*metal[e];
        responseStepScale_[e]=scale;
        predictedPotentialMax_[e]*=scale;
    }
    potentialResponseStatus_=1;
    potentialStepScale_=scale;
    if (scale<1) historyValid=false;
    finish();
}

void concentrationButlerVolmer::report() const
{
    ButlerVolmer::report();
    if (picardOutputLevel(mesh_.time().controlDict()) == 0) return;
    if (!evaluated_ || converged_ || responseMethod_=="none") return;
    if (responseMethod_=="boundaryPotential")
    {
        Info<< "    Potential response: "
            << (potentialResponseStatus_==1 ? "active" : potentialResponseStatus_==2 ? "fallback to boundaryTaylor" : "startup")
            << ", linear iterations=" << potentialSolveIterations_;
        if (potentialResponseStatus_==1)
        {
            scalar predicted=0;
            forAll(predictedPotentialMax_,e) predicted=max(predicted,predictedPotentialMax_[e]);
            Info<< ", stepScale=" << potentialStepScale_
                << ", max|deltaPhiE_predicted|=" << predicted << " V"
                << ", max Robin residual=" << potentialRowResidual_ << " A/m2"
                << ", max current defect=" << potentialCurrentDefect_ << " A";
        }
        Info<< nl;
        if (potentialResponseStatus_==1) return;
    }
    Info<< "    Concentration response:";
    forAll(electrodeNames_,e)
    {
        Info<< " " << electrodeNames_[e] << '=';
        if (responseStatus_[e]==1)
            Info<< "boundaryTaylor (stepScale=" << responseStepScale_[e] << ')';
        else Info<< (responseStatus_[e]==2 ? "frozen-state fallback" : "startup/identity");
    }
    Info<< nl;
}

} // namespace Foam

