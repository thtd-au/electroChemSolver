// cEqn.H -- PotentialReaction reactive Nernst-Planck transport

const dictionary& diffusivityDict = reactions.subDict("diffusivity");
const bool hasChargeNumber = reactions.found("chargeNumber");

surfaceScalarField gradPhiEFlux
(
    IOobject
    (
        "gradPhiEFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::snGrad(phiE)*mesh.magSf()
);

// Kinetic rates are frozen for the complete species sweep and are updated
// exactly once by electroChem.C before this include is entered.

const scalar domainVolume = gSum(mesh.V());
scalarList cMinAfterTransport(c.size(), 0.0);
scalarList cMeanAfterTransport(c.size(), 0.0);
scalarList cMaxAfterTransport(c.size(), 0.0);

forAll(c, i)
{
    cMinAfterTransport[i] = gMin(c[i].primitiveField());
    cMeanAfterTransport[i] =
        gSum(c[i].primitiveField()*mesh.V())/domainVolume;
    cMaxAfterTransport[i] = gMax(c[i].primitiveField());
}

// TEST 1 / stage A: charge-weighted implicit matrix fluxes.
// fvMatrix::flux() contains the face flux from implicit spatial terms;
// explicit migration is recorded separately, with its original sign.
surfaceScalarField chargeImplicitMatrixFluxAudit
(
    IOobject("chargeImplicitMatrixFluxAudit", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE, false),
    mesh, dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 1, 0), 0.0)
);
surfaceScalarField chargeExplicitMigrationAudit
(
    IOobject("chargeExplicitMigrationAudit", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE, false),
    mesh, dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 1, 0), 0.0)
);
surfaceScalarField chargeReconstructedImplicitAudit
(
    IOobject("chargeReconstructedImplicitAudit", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE, false),
    mesh, dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 1, 0), 0.0)
);

// TEST 1G: charge-weighted, frozen-sweep bookkeeping (diagnostic only).
// These are not guaranteed to be the same-time-state fluxes because the
// species equations are solved sequentially.  Matrix residuals are recorded
// immediately after each solve, before fvOptions.correct and clipping.
volScalarField test1GReactionCharge
(
    IOobject("test1GReactionCharge", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE, false),
    mesh, dimensionedScalar("zero", dimCurrent/dimVolume, 0.0)
);
scalar test1GMaxMatrixResidual = 0.0;
scalar test1GSumMatrixResidual = 0.0;
scalar test1GMaxPostSolveChange = 0.0;
scalar test1GMaxNegativeBeforeClip = 0.0;

forAll(c, i)
{
    if (i == inertIndex)
    {
        continue;
    }

    volScalarField& ci = c[i];

    dimensionedScalar Di
    (
        "D_" + species[i],
        dimArea/dimTime,
        diffusivityDict.get<scalar>(species[i])
    );

    scalar zi = 0.0;
    if (hasChargeNumber)
    {
        zi = reactions.subDict("chargeNumber")
            .getOrDefault<scalar>(species[i], 0.0);
    }

    surfaceScalarField migrationFlux
    (
        IOobject
        (
            "migrationFlux_" + species[i],
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(Di*zi*Fdim/(Rgas*T)*ci)*gradPhiEFlux
    );

    tmp<volScalarField> tReactionSource = reactionMgr.source(i);

    // Store the total reaction source for post-convergence diagnostics.
    // No kinetic/equilibrium split is inferred from the transported state.
    Rc[i] = tReactionSource();

    fvScalarMatrix ciEqn
    (
        fvm::ddt(ci)
      + fvm::div(phi, ci)
     ==
        fvm::laplacian(Di, ci)
      + fvc::div(migrationFlux)
      + tReactionSource()
      + fvOptions(ci)
    );

    ciEqn.relax();
    fvOptions.constrain(ciEqn);

    // S_i = V/aP_i must come from the final constrained transport matrix.
    // Exact v2512 fvScalarMatrix::diag() semantics remain compile-verification item.
    reactionMgr.setTransportResponse(i, ciEqn);

    ciEqn.solve(ci.name());

    // Sample before fvOptions.correct or concentration clipping modifies ci.
    // Flux extraction depends on fvSchemes/fluxRequired; the case has default yes.
    if (mag(zi) > SMALL)
    {
        chargeImplicitMatrixFluxAudit += Fdim*zi*ciEqn.flux();
        chargeExplicitMigrationAudit -= Fdim*zi*migrationFlux;
        chargeReconstructedImplicitAudit += Fdim*zi
        *(phi*fvc::interpolate(ci) - Di*fvc::snGrad(ci)*mesh.magSf());
    }

    // TEST 1G: exact algebraic residual of this species' solved matrix.
    // Residual units are integrated mol/s; multiply by F*z for amperes.
    if (mag(zi) > SMALL)
    {
        const scalarField speciesMatrixResidual(ciEqn.residual());
        forAll(speciesMatrixResidual, celli)
        {
            const scalar weighted = Fdim.value()*zi*speciesMatrixResidual[celli];
            test1GMaxMatrixResidual = max(test1GMaxMatrixResidual, mag(weighted));
            test1GSumMatrixResidual += weighted;
        }
        test1GReactionCharge += Fdim*zi*Rc[i];
    }
    const scalarField beforePostSolve(ci.primitiveField());
    test1GMaxNegativeBeforeClip = min
    (
        test1GMaxNegativeBeforeClip, gMin(beforePostSolve)
    );

    fvOptions.correct(ci);
    ci.correctBoundaryConditions();
    if (mag(zi) > SMALL)
    {
        forAll(beforePostSolve, celli)
        {
            test1GMaxPostSolveChange = max
            (
                test1GMaxPostSolveChange,
                mag(zi*(ci[celli] - beforePostSolve[celli]))
            );
        }
    }

    // cMinAfterTransport[i] = gMin(ci.primitiveField());
    // cMeanAfterTransport[i] =
    //     gSum(ci.primitiveField()*mesh.V())/domainVolume;
    // cMaxAfterTransport[i] = gMax(ci.primitiveField());

    // if (cMinAfterTransport[i] < 0.0)
    // {
    //     FatalErrorInFunction
    //         << "Negative concentration after transport for " << species[i]
    //         << ": min(c)=" << cMinAfterTransport[i] << nl
    //         << "PotentialReaction does not clip concentrations because that "
    //         << "would alter the transported physical state."
    //         << exit(FatalError);
    // }

    const scalar minCi = gMin(c[i]);

    if (minCi < -negativeConcentrationTolerance)
    {
        FatalErrorInFunction
            << "Negative concentration after transport for "
            << species[i] << ": min(c)=" << minCi << " mol/m3" << nl
            << "negativeConcentrationTolerance="
            << negativeConcentrationTolerance << " mol/m3" << nl
            << "The undershoot exceeds the permitted numerical tolerance."
            << exit(FatalError);
    }

    if (minCi < 0.0)
    {
        Info<< "Clipping numerical concentration undershoot for "
            << species[i]
            << ": min(c)=" << minCi << " mol/m3"
            << " -> 0"
            << nl;

        c[i].max(0.0);
    }
}

// No fast-equilibrium concentration projection is permitted here.
// lambda corrections are calculated by reactionManager after this complete
// transported state has been evaluated, and affect concentrations only through
// the next cEqn solve.

// TEST 1G: frozen-sweep charge budget.  The face fluxes below are taken
// from the exact implicit species matrices and explicit migration fluxes.
// The transient term is a backward-difference PROXY: for non-Euler ddt
// schemes it is not the assembled transient operator.  It is also evaluated
// after any fvOptions corrections/clipping, unlike the matrix residual.
{
    volScalarField test1GChargeNow
    (
        IOobject("test1GChargeNow", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh, dimensionedScalar("zero", dimMoles/dimVolume, 0.0)
    );
    volScalarField test1GChargeOld
    (
        IOobject("test1GChargeOld", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh, dimensionedScalar("zero", dimMoles/dimVolume, 0.0)
    );
    forAll(c, i)
    {
        const scalar zi = hasChargeNumber
            ? reactions.subDict("chargeNumber").getOrDefault<scalar>
                (species[i], 0.0) : 0.0;
        if (mag(zi) <= SMALL) continue;
        test1GChargeNow += zi*c[i];
        test1GChargeOld += zi*c[i].oldTime();
    }

    const surfaceScalarField test1GTotalChargeFlux
    (
        IOobject("test1GTotalChargeFlux", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        chargeImplicitMatrixFluxAudit + chargeExplicitMigrationAudit
    );
    const volScalarField test1GFluxDivergence
    (
        fvc::div(test1GTotalChargeFlux)
    );
    const scalar dt1G = runTime.deltaTValue();
    scalar maxE = 0.0;
    label worstCell1G = -1;
    forAll(mesh.V(), celli)
    {
        if (mag(test1GChargeNow[celli]) > maxE)
        {
            maxE = mag(test1GChargeNow[celli]);
            worstCell1G = celli;
        }
    }
    Info<< "Test1G: max matrix residual (charge weighted, per cell)="
        << test1GMaxMatrixResidual << " A; summed="
        << test1GSumMatrixResidual << " A" << nl
        << "  max post-solve concentration change (charge weighted)="
        << test1GMaxPostSolveChange << " mol/m3; min before clipping="
        << test1GMaxNegativeBeforeClip << " mol/m3" << nl
        << "  worst electroneutrality cell=" << worstCell1G
        << ", max|sum(z*c)|=" << maxE << " mol/m3" << nl;

    forAll(mesh.boundary(), patchi)
    {
        Info<< "  Test1G boundary " << mesh.boundary()[patchi].name()
            << " matrix+explicit-migration current="
            << sum(test1GTotalChargeFlux.boundaryField()[patchi])
            << " A" << nl;
    }
    const label focusCells[] = {0, 1, 2, 49, 50, 51, 2598, 2599};
    for (label k = 0; k < 8; ++k)
    {
        const label celli = focusCells[k];
        if (celli >= mesh.nCells()) continue;
        const scalar accumulation = Fdim.value()
            *(test1GChargeNow[celli]-test1GChargeOld[celli])/dt1G;
        const scalar fluxDiv = test1GFluxDivergence[celli];
        const scalar reaction = test1GReactionCharge[celli];
        Info<< "  Test1G cell=" << celli
            << " E=" << test1GChargeNow[celli]
            << " oldE=" << test1GChargeOld[celli]
            << " ddtChargeProxy=" << accumulation
            << " divChargeMatrixFlux=" << fluxDiv
            << " reactionCharge=" << reaction
            << " proxyBalance=" << accumulation + fluxDiv - reaction
            << " [A/m3]" << nl;
    }
    if (worstCell1G >= 0)
    {
        const label celli = worstCell1G;
        const scalar accumulation = Fdim.value()
            *(test1GChargeNow[celli]-test1GChargeOld[celli])/dt1G;
        Info<< "  Test1G WORST cell=" << celli
            << " E=" << test1GChargeNow[celli]
            << " ddtChargeProxy=" << accumulation
            << " divChargeMatrixFlux=" << test1GFluxDivergence[celli]
            << " reactionCharge=" << test1GReactionCharge[celli]
            << " proxyBalance=" << accumulation
               + test1GFluxDivergence[celli] - test1GReactionCharge[celli]
            << " [A/m3]" << nl;
    }
}


// TEST 1H: all fields below are reconstructed AFTER the complete species
// sweep. Unlike the stored matrix fluxes, they use one common transported
// concentration state and the same frozen phiE used in this sweep.
// This is a face-operator comparison, not a second potential solve.
{
    surfaceScalarField test1HSpeciesDiff
    (
        IOobject("test1HSpeciesDiff", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh, dimensionedScalar("zero", dimCurrent, 0.0)
    );
    surfaceScalarField test1HSpeciesMig
    (
        IOobject("test1HSpeciesMig", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh, dimensionedScalar("zero", dimCurrent, 0.0)
    );
    surfaceScalarField test1HSpeciesAdv
    (
        IOobject("test1HSpeciesAdv", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh, dimensionedScalar("zero", dimCurrent, 0.0)
    );
    volScalarField test1HKappa
    (
        IOobject("test1HKappa", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh, dimensionedScalar("zero", dimensionSet(-1,-3,3,0,0,2,0), 0.0)
    );
    const surfaceScalarField test1HGradPhi
    (
        IOobject("test1HGradPhi", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        fvc::snGrad(phiE)*mesh.magSf()
    );

    forAll(c, i)
    {
        const scalar zi = hasChargeNumber
            ? reactions.subDict("chargeNumber").getOrDefault<scalar>
                (species[i], 0.0) : 0.0;
        if (mag(zi) <= SMALL) continue;
        const dimensionedScalar Di
        (
            "D_test1H_" + species[i], dimArea/dimTime,
            diffusivityDict.get<scalar>(species[i])
        );
        test1HSpeciesDiff -= Fdim*zi*Di*fvc::snGrad(c[i])*mesh.magSf();
        test1HSpeciesMig -= Fdim*zi
            *fvc::interpolate(Di*zi*Fdim/(Rgas*T)*c[i])
            *test1HGradPhi;
        test1HSpeciesAdv += Fdim*zi*phi*fvc::interpolate(c[i]);
        test1HKappa += sqr(Fdim)*sqr(zi)*Di*c[i]/(Rgas*T);
    }

    const surfaceScalarField test1HPotentialMig
    (
        IOobject("test1HPotentialMig", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        -fvc::interpolate(test1HKappa)*test1HGradPhi
    );
    const surfaceScalarField test1HMigrationDifference
    (
        IOobject("test1HMigrationDifference", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        test1HSpeciesMig - test1HPotentialMig
    );
    const surfaceScalarField test1HPotentialCurrent
    (
        IOobject("test1HPotentialCurrent", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        test1HSpeciesDiff + test1HPotentialMig
    );
    const surfaceScalarField test1HSpeciesIonicCurrent
    (
        IOobject("test1HSpeciesIonicCurrent", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        test1HSpeciesDiff + test1HSpeciesMig
    );
    const surfaceScalarField test1HSpeciesTotalCurrent
    (
        IOobject("test1HSpeciesTotalCurrent", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        test1HSpeciesIonicCurrent + test1HSpeciesAdv
    );
    const surfaceScalarField test1HMatrixVsFrozen
    (
        IOobject("test1HMatrixVsFrozen", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        chargeImplicitMatrixFluxAudit + chargeExplicitMigrationAudit
        - test1HSpeciesTotalCurrent
    );
    const volScalarField test1HDivPotential(fvc::div(test1HPotentialCurrent));
    const volScalarField test1HDivIonic(fvc::div(test1HSpeciesIonicCurrent));
    const volScalarField test1HDivAdv(fvc::div(test1HSpeciesAdv));
    const volScalarField test1HDivTotal(fvc::div(test1HSpeciesTotalCurrent));
    const volScalarField test1HDivMatrixDifference(fvc::div(test1HMatrixVsFrozen));

    Info<< "Test1H FROZEN-STATE: max internal migration face difference="
        << gMax(mag(test1HMigrationDifference.primitiveField()))
        << " A; max |div(potential ionic current)|="
        << gMax(mag(test1HDivPotential.primitiveField()))
        << " A/m3; max |div(matrix-minus-frozen total)|="
        << gMax(mag(test1HDivMatrixDifference.primitiveField()))
        << " A/m3" << nl;

    scalar netPotential = 0.0;
    forAll(mesh.boundary(), patchi)
    {
        const scalar patchPotential =
            sum(test1HPotentialCurrent.boundaryField()[patchi]);
        netPotential += patchPotential;
        Info<< "  Test1H patch " << mesh.boundary()[patchi].name()
            << " diff=" << sum(test1HSpeciesDiff.boundaryField()[patchi])
            << " mig=" << sum(test1HPotentialMig.boundaryField()[patchi])
            << " adv=" << sum(test1HSpeciesAdv.boundaryField()[patchi])
            << " potentialI=" << patchPotential
            << " matrixMinusFrozen="
            << sum(test1HMatrixVsFrozen.boundaryField()[patchi])
            << " A" << nl;
    }
    Info<< "  Test1H net outward potential ionic current="
        << netPotential << " A (unsolved frozen-state residual)" << nl;

    const label focusCells1H[] = {0, 1, 49, 50, 51, 2598, 2599};
    for (label k=0; k<7; ++k)
    {
        const label celli = focusCells1H[k];
        if (celli >= mesh.nCells()) continue;
        Info<< "  Test1H cell=" << celli
            << " divPotentialI=" << test1HDivPotential[celli]
            << " divSpeciesIonic=" << test1HDivIonic[celli]
            << " divAdvection=" << test1HDivAdv[celli]
            << " divSpeciesTotal=" << test1HDivTotal[celli]
            << " divMatrixMinusFrozen=" << test1HDivMatrixDifference[celli]
            << " [A/m3]" << nl;
    }
}

// Test 1 stage A output.  Matrix and reconstructed fluxes have identical
// face orientation and electrical-current units; they need not be equal
// under nonorthogonal corrections, convection scheme differences, etc.
{
    surfaceScalarField implicitFluxDifferenceAudit
    (
        IOobject("implicitFluxDifferenceAudit", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false),
        chargeImplicitMatrixFluxAudit - chargeReconstructedImplicitAudit
    );
    const volScalarField divergenceDifferenceAudit
    (
        fvc::div(implicitFluxDifferenceAudit)
    );
    Info<< "TEST1 MATRIX-FLUX audit: max internal |matrix-reconstructed|="
        << gMax(mag(implicitFluxDifferenceAudit.primitiveField()))
        << " A; max |div difference|="
        << gMax(mag(divergenceDifferenceAudit.primitiveField()))
        << " A/m3" << nl;
    for (label celli=0; celli<=3; ++celli)
    {
        if (celli < mesh.nCells())
            Info<< "  matrix-flux divergence mismatch cell " << celli
                << " = " << divergenceDifferenceAudit[celli]
                << " A/m3" << nl;
    }
    for (label celli=2597; celli<=2599; ++celli)
    {
        if (celli < mesh.nCells())
            Info<< "  matrix-flux divergence mismatch cell " << celli
                << " = " << divergenceDifferenceAudit[celli]
                << " A/m3" << nl;
    }
    forAll(mesh.boundary(), patchi)
    {
        Info<< "  matrix-flux patch " << mesh.boundary()[patchi].name()
            << ": implicit="
            << sum(chargeImplicitMatrixFluxAudit.boundaryField()[patchi])
            << " A; reconstructed="
            << sum(chargeReconstructedImplicitAudit.boundaryField()[patchi])
            << " A; migration-used="
            << sum(chargeExplicitMigrationAudit.boundaryField()[patchi])
            << " A; difference="
            << sum(implicitFluxDifferenceAudit.boundaryField()[patchi])
            << " A" << nl;
    }
}

