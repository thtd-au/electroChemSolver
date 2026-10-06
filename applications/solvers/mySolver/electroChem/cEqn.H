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
    fvOptions.correct(ci);
    ci.correctBoundaryConditions();

    cMinAfterTransport[i] = gMin(ci.primitiveField());
    cMeanAfterTransport[i] =
        gSum(ci.primitiveField()*mesh.V())/domainVolume;
    cMaxAfterTransport[i] = gMax(ci.primitiveField());

    if (cMinAfterTransport[i] < 0.0)
    {
        FatalErrorInFunction
            << "Negative concentration after transport for " << species[i]
            << ": min(c)=" << cMinAfterTransport[i] << nl
            << "PotentialReaction does not clip concentrations because that "
            << "would alter the transported physical state."
            << exit(FatalError);
    }
}

// No fast-equilibrium concentration projection is permitted here.
// lambda corrections are calculated by reactionManager after this complete
// transported state has been evaluated, and affect concentrations only through
// the next cEqn solve.
