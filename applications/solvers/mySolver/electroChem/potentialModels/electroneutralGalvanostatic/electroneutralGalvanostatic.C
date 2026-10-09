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
        IOobject("kappa", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1,-3,3,0,0,2,0), 0.0)
    ),
    jDiff_
    (
        IOobject("jDiff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
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
        Info<< "Test1F: fixed-potential reservoir detected; initial gauge shift disabled" << nl;
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
    // TEST 1C: optional gauge-location-only A/B override.
    const label overrideCell = potentialDict_.getOrDefault<label>
    (
        "test1CReferenceCell", -1
    );
    if (overrideCell >= mesh_.nCells())
    {
        FatalIOErrorInFunction(potentialDict_)
            << "test1CReferenceCell out of range: " << overrideCell
            << exit(FatalIOError);
    }
    referenceCell_ = overrideCell >= 0 ? overrideCell : localCandidate;
    Info<< "Test1C selected reference cell=" << referenceCell_
        << " (patch candidate=" << localCandidate << ")" << nl;
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


    // DIAGNOSTIC ONLY: compare two diffusion-current face discretizations.
    // The first uses the snGrad operator of the species diffusion Laplacian;
    // the second is the existing reconstructed-gradient potential flux.
    surfaceScalarField jDiffSpeciesFaces
    (
        IOobject
        (
            "jDiffSpeciesFacesAudit", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", jDiffFlux_().dimensions(), 0.0)
    );

    forAll(c_, i)
    {
        if (mag(z_[i]) <= SMALL) continue;
        const dimensionedScalar Di
        (
            "D_audit_" + species_[i],
            dimensionSet(0, 2, -1, 0, 0, 0, 0), D_[i]
        );
        jDiffSpeciesFaces -=
            F_*z_[i]*Di*fvc::snGrad(c_[i])*mesh_.magSf();
    }

    surfaceScalarField jDiffMismatch
    (
        IOobject
        (
            "jDiffMismatchAudit", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        jDiffSpeciesFaces - jDiffFlux_()
    );

    tmp<volScalarField> tDivMismatch = fvc::div(jDiffMismatch);
    const volScalarField& divMismatch = tDivMismatch();

    Info<< "Diffusion-flux audit: max internal |face mismatch|="
        << gMax(mag(jDiffMismatch.primitiveField())) << " A" << nl
        << "  max cell |div(mismatch)|="
        << gMax(mag(divMismatch.primitiveField())) << " A/m3" << nl;

    for (label celli = 2597; celli <= 2599; ++celli)
    {
        if (celli < mesh_.nCells())
        {
            Info<< "  diffusion mismatch divergence cell " << celli
                << " = " << divMismatch[celli] << " A/m3" << nl;
        }
    }

    forAll(mesh_.boundary(), patchi)
    {
        const scalarField& mismatchPatch = jDiffMismatch.boundaryField()[patchi];
        const scalarField& speciesPatch = jDiffSpeciesFaces.boundaryField()[patchi];
        const scalarField& potentialPatch = jDiffFlux_().boundaryField()[patchi];
        Info<< "  patch " << mesh_.boundary()[patchi].name()
            << ": sum species diffusion current=" << sum(speciesPatch)
            << " A, sum potential diffusion current=" << sum(potentialPatch)
            << " A, sum mismatch=" << sum(mismatchPatch) << " A" << nl;
    }

    // A/B experiment: use the species-Laplacian face-normal diffusion current
    // for BOTH the potential equation RHS and galvanostatic potential BC.
    // The audit above intentionally retains the original reconstruction
    // for comparison before this assignment.
    jDiffFlux_() = jDiffSpeciesFaces;
    Info<< "Diffusion A/B: potential current now uses species snGrad face flux" << nl;
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
    Info<< "Test1F correction BC mode="
        << (hasFixedPotential ? "mixed Dirichlet-Neumann" : "all-Neumann") << nl;

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


    // TEST 1E: diagnostic only. Evaluate species diffusive and migration
    // currents on every boundary, using the frozen transported state.
    // No boundary conditions or equation coefficients are modified.
    Info<< "Test1E boundary ionic-current audit (outward A):" << nl;
    forAll(mesh_.boundary(), patchi)
    {
        const word& patchName = mesh_.boundary()[patchi].name();
        if
        (
            patchName != "inlet" && patchName != "outlet"
         && patchName != "cathode" && patchName != "anode"
        ) continue;

        scalar diffSum = 0.0;
        scalar migSum = 0.0;
        const scalarField& faceAreas = mesh_.magSf().boundaryField()[patchi];
        const scalarField potentialNormalGrad =
            phiE_.boundaryField()[patchi].snGrad();
        const scalarField kappaFaces(kappaFace_().boundaryField()[patchi]);

        forAll(c_, i)
        {
            if (mag(z_[i]) <= SMALL) continue;
            const scalarField speciesNormalGrad
            (
                c_[i].boundaryField()[patchi].snGrad()
            );
            scalar speciesDiff = 0.0;
            forAll(faceAreas, facei)
            {
                speciesDiff -=
                    F_.value()*z_[i]*D_[i]
                   *speciesNormalGrad[facei]*faceAreas[facei];
            }
            diffSum += speciesDiff;
            if (patchName == "inlet" || patchName == "outlet")
            {
                Info<< "  Test1E " << patchName << " species "
                    << species_[i]
                    << " bc=" << c_[i].boundaryField()[patchi].type()
                    << " diffCurrent=" << speciesDiff << " A" << nl;
            }
        }

        forAll(faceAreas, facei)
        {
            migSum -= kappaFaces[facei]
                *potentialNormalGrad[facei]*faceAreas[facei];
        }

        const scalar potentialDiff =
            sum(jDiffFlux_().boundaryField()[patchi]);
        Info<< "  Test1E patch " << patchName
            << " phiBC=" << phiE_.boundaryField()[patchi].type()
            << " speciesDiff=" << diffSum
            << " potentialDiff=" << potentialDiff
            << " migration=" << migSum
            << " totalIonic=" << diffSum + migSum
            << " A" << nl;
    }

    // TEST 1C: capture the unreferenced equation BEFORE gauge modification.
    // fvMatrix::residual() returns a per-cell algebraic residual for the
    // CURRENT (pre-solve) deltaPhi, including its explicit source.
    const scalarField preReferenceResidual(deltaEqn.residual());
    // TEST 1D: for a pure-Neumann Laplacian, sum of the unreferenced
    // algebraic residual at deltaPhi=0 is the discrete compatibility defect.
    // At this point deltaPhi is identically zero and BCs are evaluated.
    const scalar preReferenceGlobalSum = gSum(preReferenceResidual);
    Info<< "Test1D pre-reference global sum residual="
        << preReferenceGlobalSum << " A (integrated)" << nl;

    if (!hasFixedPotential)
    {
        const scalar preReferenceDiag = deltaEqn.diag()[referenceCell_];
        const scalar preReferenceSource = deltaEqn.source()[referenceCell_];
        Info<< "Test1C before reference: cell=" << referenceCell_
            << " diag=" << preReferenceDiag
            << " source=" << preReferenceSource
            << " residual=" << preReferenceResidual[referenceCell_] << nl;
        deltaEqn.setReference(referenceCell_, 0.0);
        Info<< "Test1C reference modification: deltaDiag="
            << deltaEqn.diag()[referenceCell_] - preReferenceDiag
            << " deltaSource="
            << deltaEqn.source()[referenceCell_] - preReferenceSource << nl;
    }
    else
    {
        Info<< "Test1F: skipping setReference; Dirichlet reservoir fixes gauge" << nl;
    }
    deltaEqn.solve(phiE_.name());

    // Evaluate ORIGINAL unreferenced discrete balance at the solved deltaPhi.
    // Keep the original equation as a separate matrix object so the
    // reference-row modification cannot hide the physical residual.
    fvScalarMatrix unreferencedEqn
    (
        fvm::laplacian(kappa_, deltaPhi)
      ==
        fvc::div(jDiffFlux_())
      - fvc::laplacian(kappa_, phiE_)
    );
    const scalarField physicalResidual(unreferencedEqn.residual());
    const scalarField referencedResidual(deltaEqn.residual());
    // TEST 1D: the sum of unreferenced residuals should equal the
    // global Neumann compatibility defect, regardless of gauge location.
    const scalar globalPhysicalResidual = gSum(physicalResidual);
    const scalar globalReferencedResidual = gSum(referencedResidual);
    const scalar refResidual = physicalResidual[referenceCell_]; // diagnostic only for Dirichlet
    Info<< "Test1D post-solve global unreferenced residual="
        << globalPhysicalResidual
        << " A, referenceCellResidual=" << refResidual
        << " A, global referenced residual=" << globalReferencedResidual
        << " A, (global-refCell)="
        << globalPhysicalResidual-refResidual << " A" << nl;


    scalar maxPhysical = 0.0;
    label maxPhysicalCell = -1;
    forAll(physicalResidual, celli)
    {
        if (mag(physicalResidual[celli]) > maxPhysical)
        {
            maxPhysical = mag(physicalResidual[celli]);
            maxPhysicalCell = celli;
        }
    }
    Info<< "Test1C original unreferenced equation: max|algebraic residual|="
        << maxPhysical << " (integrated equation units), cell="
        << maxPhysicalCell << nl;
    for (label celli = 0; celli < mesh_.nCells(); ++celli)
    {
        if (celli != referenceCell_ && celli > 3 && celli < 2597) continue;
        if (celli > 2599 && celli != referenceCell_) continue;
        Info<< "  Test1C cell " << celli
            << ": originalResidual=" << physicalResidual[celli]
            << ", referencedResidual=" << referencedResidual[celli]
            << ", originalResidual/V="
            << physicalResidual[celli]/mesh_.V()[celli]
            << " (volume-normalized)" << nl;
    }

    // TEST 1B diagnostic: frozen-state species migration versus the
    // conductivity-weighted potential migration, before applying deltaPhi.
    // This is a face-flux compatibility check, not an exact matrix residual.
    surfaceScalarField speciesMigrationFrozenAudit
    (
        IOobject("speciesMigrationFrozenAudit", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh_, dimensionedScalar("zero", jDiffFlux_().dimensions(), 0.0)
    );
    const surfaceScalarField gradPhiFrozenAudit
    (
        IOobject("gradPhiFrozenAudit", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        fvc::snGrad(phiE_)*mesh_.magSf()
    );
    forAll(c_, i)
    {
        if (mag(z_[i]) <= SMALL) continue;
        const dimensionedScalar Di
        (
            "D_audit_" + species_[i],
            dimensionSet(0, 2, -1, 0, 0, 0, 0), D_[i]
        );
        speciesMigrationFrozenAudit -=
            F_*z_[i]
           *fvc::interpolate(Di*z_[i]*F_/(Rgas_*T_)*c_[i])
           *gradPhiFrozenAudit;
    }
    const surfaceScalarField potentialMigrationFrozenAudit
    (
        IOobject("potentialMigrationFrozenAudit", mesh_.time().timeName(),
            mesh_, IOobject::NO_READ, IOobject::NO_WRITE, false),
        -fvc::interpolate(kappa_)*gradPhiFrozenAudit
    );
    const surfaceScalarField migrationMismatchFrozenAudit
    (
        IOobject("migrationMismatchFrozenAudit", mesh_.time().timeName(),
            mesh_, IOobject::NO_READ, IOobject::NO_WRITE, false),
        speciesMigrationFrozenAudit - potentialMigrationFrozenAudit
    );
    const volScalarField migrationMismatchDivAudit
    (
        fvc::div(migrationMismatchFrozenAudit)
    );
    Info<< "Test1B frozen migration algebra: max internal |difference|="
        << gMax(mag(migrationMismatchFrozenAudit.primitiveField()))
        << " A, max cell |div difference|="
        << gMax(mag(migrationMismatchDivAudit.primitiveField()))
        << " A/m3" << nl;

    // The correction matrix is the actual potential operator solved here.
    // Its flux contains ONLY the implicit deltaPhi Laplacian; the absolute
    // phiE gradient and jDiff are separate explicit RHS terms.
    const surfaceScalarField deltaMatrixFluxAudit
    (
        IOobject("deltaMatrixFluxAudit", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        deltaEqn.flux()
    );
    const surfaceScalarField deltaSnGradFluxAudit
    (
        IOobject("deltaSnGradFluxAudit", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        fvc::interpolate(kappa_)*fvc::snGrad(deltaPhi)*mesh_.magSf()
    );
    const surfaceScalarField deltaFluxMismatchAudit
    (
        IOobject("deltaFluxMismatchAudit", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        deltaMatrixFluxAudit - deltaSnGradFluxAudit
    );
    Info<< "Test1B potential correction matrix flux: max internal |matrix-snGrad|="
        << gMax(mag(deltaFluxMismatchAudit.primitiveField()))
        << " A" << nl;

    const surfaceScalarField absoluteCurrentFrozenAudit
    (
        IOobject("absoluteCurrentFrozenAudit", mesh_.time().timeName(),
            mesh_, IOobject::NO_READ, IOobject::NO_WRITE, false),
        jDiffFlux_() + potentialMigrationFrozenAudit
    );
    const surfaceScalarField correctedCurrentUnrelaxedAudit
    (
        IOobject("correctedCurrentUnrelaxedAudit", mesh_.time().timeName(),
            mesh_, IOobject::NO_READ, IOobject::NO_WRITE, false),
        absoluteCurrentFrozenAudit - deltaMatrixFluxAudit
    );
    // TEST 1D: boundary currents of the UNRELAXED correction state.
    // Outward patch sums must agree with the volume-integrated current
    // divergence, including the contribution at the reference cell.
    scalar boundaryNetCurrent = 0.0;
    forAll(mesh_.boundary(), patchi)
    {
        const scalar patchCurrent =
            sum(correctedCurrentUnrelaxedAudit.boundaryField()[patchi]);
        boundaryNetCurrent += patchCurrent;
        Info<< "Test1D boundary current "
            << mesh_.boundary()[patchi].name()
            << "=" << patchCurrent << " A" << nl;
    }
    Info<< "Test1D net outward boundary current="
        << boundaryNetCurrent << " A" << nl;

    const volScalarField correctedCurrentDivAudit
    (
        fvc::div(correctedCurrentUnrelaxedAudit)
    );
    scalar volumeIntegratedDiv = 0.0;
    forAll(mesh_.V(), celli)
    {
        volumeIntegratedDiv +=
            correctedCurrentDivAudit[celli]*mesh_.V()[celli];
    }
    Info<< "Test1D volume integrated corrected-current divergence="
        << volumeIntegratedDiv << " A, boundary-minus-volume="
        << boundaryNetCurrent-volumeIntegratedDiv << " A" << nl;

    Info<< "Test1B unrelaxed correction current: max cell |div|="
        << gMax(mag(correctedCurrentDivAudit.primitiveField()))
        << " A/m3 (matrix-flux reconstruction, not assembled residual)" << nl;

    forAll(mesh_.boundary(), patchi)
    {
        Info<< "  Test1B patch " << mesh_.boundary()[patchi].name()
            << ": speciesMig="
            << sum(speciesMigrationFrozenAudit.boundaryField()[patchi])
            << " A, potentialMig="
            << sum(potentialMigrationFrozenAudit.boundaryField()[patchi])
            << " A, deltaMatrix="
            << sum(deltaMatrixFluxAudit.boundaryField()[patchi])
            << " A, correctedCurrent="
            << sum(correctedCurrentUnrelaxedAudit.boundaryField()[patchi])
            << " A" << nl;
    }
    for (label celli = 0; celli < mesh_.nCells(); ++celli)
    {
        if (celli > 3 && celli < 2597) continue;
        if (celli > 2599) break;
        Info<< "  Test1B cell " << celli
            << ": migrationMismatchDiv=" << migrationMismatchDivAudit[celli]
            << ", correctedCurrentDiv=" << correctedCurrentDivAudit[celli]
            << " A/m3" << nl;
    }

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
    label worstCell = -1;

    forAll(chargeResidual_, celli)
    {
        const scalar absE = mag(chargeResidual_[celli]);
        const scalar relE = absE/chargeScale_[celli];

        if (absE > maxAbs)
        {
            maxAbs = absE;
            worstCell = celli;
        }

        maxRel = max(maxRel, relE);
    }

    reduce(maxAbs, maxOp<scalar>());
    reduce(maxRel, maxOp<scalar>());

    Info<< "Electroneutrality: max|E|=" << maxAbs
        << " mol/m3, max relative=" << maxRel
        << ", converged="
        << (evaluated_ && converged_ ? "yes" : "no")
        << nl;

    if (worstCell >= 0)
    {
        Info<< "  worst cell=" << worstCell
            << ", centre=" << mesh_.C()[worstCell]
            << ", E=" << chargeResidual_[worstCell]
            << " mol/m3"
            << ", relative="
            << mag(chargeResidual_[worstCell])
              /chargeScale_[worstCell]
            << nl;

        Info<< "  charged-species contributions:" << nl;

        forAll(c_, i)
        {
            if (mag(z_[i]) > SMALL)
            {
                Info<< "    " << species_[i]
                    << ": z=" << z_[i]
                    << ", c=" << c_[i][worstCell]
                    << ", z*c=" << z_[i]*c_[i][worstCell]
                    << " mol/m3"
                    << nl;
            }
        }

        // Global charge diagnostic for the currently evaluated transported state.
        scalar integratedE = 0.0;
        scalar integratedAbsE = 0.0;
        scalar totalVolume = 0.0;

        const scalarField& cellVolumes = mesh_.V();
        forAll(chargeResidual_, celli)
        {
            const scalar dQ = cellVolumes[celli]*chargeResidual_[celli];
            integratedE += dQ;
            integratedAbsE += mag(dQ);
            totalVolume += cellVolumes[celli];
        }
        reduce(integratedE, sumOp<scalar>());
        reduce(integratedAbsE, sumOp<scalar>());
        reduce(totalVolume, sumOp<scalar>());

        Info<< "  Integrated charge: sum(V*E)=" << integratedE
            << " mol charge-equivalents"
            << ", F*sum(V*E)=" << F_.value()*integratedE << " C" << nl
            << "  Absolute charge integral: sum(V*|E|)="
            << integratedAbsE << " mol charge-equivalents" << nl
            << "  Volume-mean signed E=" << integratedE/totalVolume
            << " mol/m3"
            << ", cancellation ratio="
            << (integratedAbsE > VSMALL
                ? mag(integratedE)/integratedAbsE : 0.0)
            << nl;

        // Diagnostic: electroneutrality profile through the worst-cell column.
        // Select cells whose centres have the same x,z location as worstCell,
        // then report the cells closest to the worst cell in y.

        const vector& worstCentre = mesh_.C()[worstCell];

        DynamicList<label> columnCells;

        const scalar geomTol = 1e-9;

        forAll(mesh_.C(), celli)
        {
            const vector& cc = mesh_.C()[celli];

            if
            (
                mag(cc.x() - worstCentre.x()) < geomTol
            && mag(cc.z() - worstCentre.z()) < geomTol
            )
            {
                columnCells.append(celli);
            }
        }

        // Sort manually by distance in y from the worst cell.
        // Diagnostic only, so keep this independent of mesh numbering.
        for (label i = 0; i < columnCells.size(); ++i)
        {
            for (label j = i + 1; j < columnCells.size(); ++j)
            {
                const scalar di =
                    mag(mesh_.C()[columnCells[i]].y() - worstCentre.y());

                const scalar dj =
                    mag(mesh_.C()[columnCells[j]].y() - worstCentre.y());

                if (dj < di)
                {
                    const label tmp = columnCells[i];
                    columnCells[i] = columnCells[j];
                    columnCells[j] = tmp;
                }
            }
        }

        Info<< "  electroneutrality profile from worst cell:" << nl;
        Info<< "    cell        y [m]        E [mol/m3]        relative" << nl;

        const label nReport = min(label(10), columnCells.size());

        for (label n = 0; n < nReport; ++n)
        {
            const label celli = columnCells[n];

            Info<< "    "
                << celli << "    "
                << mesh_.C()[celli].y() << "    "
                << chargeResidual_[celli] << "    "
                << mag(chargeResidual_[celli])/chargeScale_[celli]
                << nl;
        }
    }
}

} // End namespace Foam

