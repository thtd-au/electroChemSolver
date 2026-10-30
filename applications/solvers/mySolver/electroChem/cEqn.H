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
    mesh,
    dimensionedScalar("zero", phiE.dimensions()*dimLength, 0.0)
);

if (potentialMgr.active())
{
    gradPhiEFlux = fvc::snGrad(phiE)*mesh.magSf();
}

// Kinetic rates are frozen for the complete species sweep and are updated
// exactly once by electroChem.C before this include is entered.

reactionMgr.beginTransportDiagnostics();

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

    const bool captureReactionSplit = diagTransportActive || diagConservationActive;
    tmp<volScalarField> tReactionSource = reactionMgr.source
    (
        i,
        captureReactionSplit ? &RcKinetic[i] : nullptr,
        captureReactionSplit ? &RcEquilibrium[i] : nullptr
    );

    // Store the total reaction source for post-convergence diagnostics.
    // Split snapshots use the same frozen kinetic rates and lambda as this solve.
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
    reactionMgr.setTransportResponse(i, ciEqn);

    dictionary transportControls(ciEqn.solverDict(ci.name()));
    transportControls.set("log", reactionMgr.transportLogLevel());
    const auto transportPerformance = ciEqn.solve(transportControls);
    if (reactionMgr.equilibriumDiagnostics())
    {
        Pout<< "EQDIAG linearSolve species=" << species[i]
            << " initialResidual=" << transportPerformance.initialResidual()
            << " finalResidual=" << transportPerformance.finalResidual()
            << " iterations=" << transportPerformance.nIterations()
            << " converged=" << transportPerformance.converged() << nl;
    }

    fvOptions.correct(ci);
    ci.correctBoundaryConditions();

    // Inspect before clipping or the existing fatal negativity check.
    reactionMgr.reportTransportDiagnostics(i);

    if (printCFieldsBeforeClipping)
    {
        cPreMin[i] = gMin(ci.primitiveField());
        cPreMean[i] = gSum(ci.primitiveField()*mesh.V())/gSum(mesh.V());
        cPreMax[i] = gMax(ci.primitiveField());
    }

    const scalar minCi = gMin(c[i]);

    if (minCi < -negativeConcentrationTolerance)
    {
        if (permitInitialGuessRetry)
        {
            Info<< "    Rejecting initial equilibrium guess: "
                << species[i] << " min(c)=" << minCi << " mol/m3" << nl;
            transportGuessRejected = true;
            break; // Caller restores all species and reduces the lambda guess.
        }
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

