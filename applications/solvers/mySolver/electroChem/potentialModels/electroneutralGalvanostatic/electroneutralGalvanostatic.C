#include "electroneutralGalvanostatic.H"
#include "addToRunTimeSelectionTable.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "fvc.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvmLaplacian.H"
#include "PstreamReduceOps.H"

namespace Foam
{

defineTypeNameAndDebug(electroneutralGalvanostatic, 0);
addToRunTimeSelectionTable
(
    potentialModel,
    electroneutralGalvanostatic,
    dictionary
);

electroneutralGalvanostatic::electroneutralGalvanostatic
(
    fvMesh& mesh,
    const dictionary& electroChemistry,
    const dictionary& reactions,
    const wordList& species,
    PtrList<volScalarField>& c,
    volScalarField& phiE,
    const dimensionedScalar& T,
    const dimensionedScalar& Rgas,
    const dimensionedScalar& F
)
:
    potentialModel(mesh, electroChemistry, reactions, species, c, phiE, T, Rgas, F),
    potentialDict_(electroChemistry.subDict("potential")),
    appliedCurrentDensity_(potentialDict_.get<scalar>("appliedCurrentDensity")),
    referencePotential_(potentialDict_.get<scalar>("referencePotential")),
    referencePatchName_(potentialDict_.get<word>("referencePatch")),
    referencePatchID_(mesh.boundaryMesh().findPatchID(referencePatchName_)),
    referenceCell_(-1),
    minConductivity_(potentialDict_.getOrDefault<scalar>("minConductivity", SMALL)),
    absTol_(0.0),
    relTol_(0.0),
    relaxation_(1.0),
    z_(species.size(), 0.0),
    D_(species.size(), 0.0),
    chargeResidual_(mesh.nCells(), 0.0),
    chargeScale_(mesh.nCells(), 0.0),
    evaluated_(false),
    converged_(false),
    kappa_
    (
        IOobject("kappaModelInternal", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1,-3,3,0,0,2,0), 0.0)
    ),
    jDiff_
    (
        IOobject("jDiffModelInternal", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh,
        dimensionedVector("zero", dimensionSet(0,-2,0,0,0,1,0), vector::zero)
    ),
    kappaFace_(nullptr),
    jDiffFlux_(nullptr)
{
    if (appliedCurrentDensity_ < 0.0)
    {
        FatalIOErrorInFunction(potentialDict_)
            << "appliedCurrentDensity must be non-negative"
            << exit(FatalIOError);
    }
    if (minConductivity_ <= 0.0)
    {
        FatalIOErrorInFunction(potentialDict_)
            << "minConductivity must be positive"
            << exit(FatalIOError);
    }
    if (referencePatchID_ < 0)
    {
        FatalIOErrorInFunction(potentialDict_)
            << "Unknown referencePatch " << referencePatchName_
            << exit(FatalIOError);
    }

    const dictionary& picard = mesh.solutionDict().subDict("Picard");
    const dictionary& en = picard.subDict("electroneutrality");
    absTol_ = en.get<scalar>("absTol");
    relTol_ = en.get<scalar>("relTol");
    relaxation_ = en.get<scalar>("relaxation");

    if (absTol_ <= 0 || relTol_ < 0 || relaxation_ <= 0 || relaxation_ > 1)
    {
        FatalIOErrorInFunction(en)
            << "Require absTol>0, relTol>=0 and 0<relaxation<=1"
            << exit(FatalIOError);
    }

    const dictionary& dDict = reactions.subDict("diffusivity");
    const bool hasZ = reactions.found("chargeNumber");
    const dictionary* zDict = hasZ ? &reactions.subDict("chargeNumber") : nullptr;

    forAll(species_, i)
    {
        D_[i] = dDict.get<scalar>(species_[i]);
        if (D_[i] <= 0)
        {
            FatalIOErrorInFunction(dDict)
                << "diffusivity for " << species_[i] << " must be positive"
                << exit(FatalIOError);
        }
        z_[i] = zDict ? zDict->getOrDefault<scalar>(species_[i], 0.0) : 0.0;
    }

    selectReferenceCell();

    // Register the face fields required by galvanostaticPotentialFlux, but do
    // not evaluate the coupled electrode BC during construction.  The first
    // species transport uses the gradient already stored in phiE (phi^(1)).
    // After transport, correct() refreshes kappa/jDiff and deliberately updates
    // the absolute potential BC for the following Picard state.
    kappaFace_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "kappaFace",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(kappa_)
        )
    );

    jDiffFlux_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "jDiffFlux",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (fvc::interpolate(jDiff_) & mesh.Sf())
        )
    );

    // Preserve the supplied phiE field as the first Picard state.  Only apply
    // a constant gauge shift, which leaves all electric-field gradients
    // unchanged, so the selected reference cell matches referencePotential.
    bool hasFixedPotential = false;
    forAll(phiE_.boundaryField(), patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(phiE_.boundaryField()[patchi]))
        {
            hasFixedPotential = true;
        }
    }
    if (!hasFixedPotential)
    {
        const scalar gaugeShift = referencePotential_ - phiE_[referenceCell_];
        phiE_.primitiveFieldRef() += gaugeShift;
        forAll(phiE_.boundaryField(), patchi)
        {
            phiE_.boundaryFieldRef()[patchi] += gaugeShift;
        }
    }
    else
    {
        Info<< "Fixed-potential boundary detected; initial gauge shift disabled" << nl;
    }

}

void electroneutralGalvanostatic::selectReferenceCell()
{
    const labelUList& faceCells = mesh_.boundary()[referencePatchID_].faceCells();
    label localCandidate = faceCells.empty() ? -1 : faceCells[0];

    // v1 source implementation: deterministic local selection. Parallel-global
    // reference selection must be compile/runtime verified in the v2512 build.
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Parallel-safe reference-cell selection for referencePatch "
            << referencePatchName_ << " still requires v2512 integration verification."
            << exit(FatalError);
    }
    if (localCandidate < 0)
    {
        FatalErrorInFunction
            << "referencePatch " << referencePatchName_ << " has no faces"
            << exit(FatalError);
    }
    referenceCell_ = localCandidate;

}

void electroneutralGalvanostatic::updateCurrentFields()
{
    kappa_ = dimensionedScalar("zero", kappa_.dimensions(), 0.0);
    jDiff_ = dimensionedVector("zero", jDiff_.dimensions(), vector::zero);

    forAll(c_, i)
    {
        if (mag(z_[i]) <= SMALL)
        {
            continue;
        }

        const dimensionedScalar Di
        (
            "D_" + species_[i],
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            D_[i]
        );

        // Full-field algebra is intentional here.  This routine is called only
        // after the species transport step, whose electrode BCs have already
        // been evaluated using the frozen phiE gradient from Picard state k.
        // Hence boundary concentrations/gradients now belong to c^(k+1).
        kappa_ +=
            sqr(F_)*sqr(z_[i])*Di*c_[i]
           /(Rgas_*T_);

        jDiff_ -= F_*z_[i]*Di*fvc::grad(c_[i]);
    }

    kappaFace_() = fvc::interpolate(kappa_);
    jDiffFlux_() = fvc::interpolate(jDiff_) & mesh_.Sf();


    // Use species-face diffusion discretization consistent with transport.
    surfaceScalarField jDiffSpeciesFaces
    (
        IOobject("jDiffSpeciesFaces", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh_, dimensionedScalar("zero", jDiffFlux_().dimensions(), 0.0)
    );
    forAll(c_, i)
    {
        if (mag(z_[i]) <= SMALL) continue;
        const dimensionedScalar Di
        (
            "D_" + species_[i],
            dimensionSet(0, 2, -1, 0, 0, 0, 0), D_[i]
        );
        jDiffSpeciesFaces -=
            F_*z_[i]*Di*fvc::snGrad(c_[i])*mesh_.magSf();
    }
    jDiffFlux_() = jDiffSpeciesFaces;

}

void electroneutralGalvanostatic::evaluateConstraint()
{
    evaluated_ = true;
    converged_ = true;

    forAll(chargeResidual_, celli)
    {
        scalar E = 0.0;
        scalar ionic = 0.0;
        forAll(c_, i)
        {
            const scalar ci = c_[i][celli];
            E += z_[i]*ci;
            ionic += mag(z_[i])*ci;
        }

        chargeResidual_[celli] = E;
        chargeScale_[celli] = max(ionic, absTol_);

        if (!(mag(E) <= absTol_ || mag(E)/chargeScale_[celli] <= relTol_))
        {
            converged_ = false;
        }
    }

    reduce(converged_, andOp<bool>());
}

bool electroneutralGalvanostatic::converged() const
{
    if (!evaluated_)
    {
        FatalErrorInFunction
            << "electroneutrality requested before evaluateConstraint()"
            << exit(FatalError);
    }
    return converged_;
}

void electroneutralGalvanostatic::setCorrectionBoundaryConditions
(
    volScalarField& deltaPhi
) const
{
    forAll(deltaPhi.boundaryField(), patchi)
    {
        fvPatchScalarField& dp = deltaPhi.boundaryFieldRef()[patchi];
        const fvPatchScalarField& p = phiE_.boundaryField()[patchi];

        if (isA<fixedValueFvPatchScalarField>(p))
        {
            fixedValueFvPatchScalarField& dpv =
                refCast<fixedValueFvPatchScalarField>(dp);
            // The absolute potential boundary is already prescribed.
            // Its correction must vanish at the reservoir.
            dpv == scalar(0);
        }
        else
        {
            fixedGradientFvPatchScalarField& dpg =
                refCast<fixedGradientFvPatchScalarField>(dp);
            if (isA<fixedGradientFvPatchScalarField>(p))
            {
                const fixedGradientFvPatchScalarField& pg =
                    refCast<const fixedGradientFvPatchScalarField>(p);
                dpg.gradient() = pg.gradient() - p.snGrad();
            }
            else
            {
                dpg.gradient() = 0.0;
            }
        }
    }
}

void electroneutralGalvanostatic::correct()
{
    if (!evaluated_)
    {
        FatalErrorInFunction
            << "correct() requires a transported state evaluated first"
            << exit(FatalError);
    }
    if (converged_) return;

    updateCurrentFields();

    // Update the absolute phiE galvanostatic gradients now that jDiffFlux and
    // kappaFace correspond to the current transported concentrations.
    phiE_.boundaryFieldRef().evaluate();

    wordList correctionPatchTypes(mesh_.boundary().size(), "fixedGradient");
    bool hasFixedPotential = false;
    forAll(phiE_.boundaryField(), patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(phiE_.boundaryField()[patchi]))
        {
            correctionPatchTypes[patchi] = "fixedValue";
            hasFixedPotential = true;
        }
    }

    volScalarField deltaPhi
    (
        IOobject("deltaPhiE", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar("zero", phiE_.dimensions(), 0.0),
        correctionPatchTypes
    );

    setCorrectionBoundaryConditions(deltaPhi);
    deltaPhi.correctBoundaryConditions();

    // Residual correction of the frozen-concentration absolute equation:
    // div(kappa grad(deltaPhi)) = div(jDiffSpeciesFaces)
    //                                - div(kappa grad(phiE)).
    fvScalarMatrix deltaEqn
    (
        fvm::laplacian(kappa_, deltaPhi)
      ==
        fvc::div(jDiffFlux_())
      - fvc::laplacian(kappa_, phiE_)
    );


    if (!hasFixedPotential)
    {
        deltaEqn.setReference(referenceCell_, 0.0);
    }
    deltaEqn.solve(phiE_.name());

    const scalar maxAbsDeltaPhi = gMax(mag(deltaPhi.primitiveField()));
    Info<< "Potential correction: max|deltaPhiE|="
        << maxAbsDeltaPhi << " V"
        << ", relaxation=" << relaxation_
        << ", applied max|alpha*deltaPhiE|="
        << relaxation_*maxAbsDeltaPhi << " V"
        << nl;

    phiE_.primitiveFieldRef() += relaxation_*deltaPhi.primitiveField();
    phiE_.correctBoundaryConditions();

    evaluated_ = false; // correction must be followed by transport before retest
}

void electroneutralGalvanostatic::report() const
{
    scalar maxAbs = 0.0;
    scalar maxRel = 0.0;
    forAll(chargeResidual_, celli)
    {
        maxAbs = max(maxAbs, mag(chargeResidual_[celli]));
        maxRel = max(maxRel,
            mag(chargeResidual_[celli])/chargeScale_[celli]);
    }
    reduce(maxAbs, maxOp<scalar>());
    reduce(maxRel, maxOp<scalar>());
    Info<< "Electroneutrality: max|E|=" << maxAbs
        << " mol/m3, max relative=" << maxRel
        << ", converged=" << (evaluated_ && converged_ ? "yes" : "no")
        << nl;
}

} // End namespace Foam

