// cEqn.H
label inertIndex = species.find(inertSpecie);

const dictionary& diffusivityDict = reactions.subDict("diffusivity");
const bool hasChargeNumber = reactions.found("chargeNumber");

dimensionedScalar T
(
    "T",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    reactions.getOrDefault<scalar>("temperature", 298.15)
);

dimensionedScalar Rgas
(
    "R",
    dimensionSet(1, 2, -2, -1, -1, 0, 0),
    8.31446261815324
);

dimensionedScalar Fdim
(
    "FaradayConstant",
    dimensionSet(0, 0, 1, 0, -1, 1, 0),
    reactions.getOrDefault<scalar>("FaradayConstant", 96485.3329)
);

forAll(Rc, i)
{
    tmp<volScalarField> tSource
    (
        reactionMgr.kineticSource(species[i])
    );

    Rc[i] = tSource();
}

// Store concentration statistics before positivity clipping.
// These are printed later in diagnostics.H so that solver output
// is not interrupted by diagnostic messages.
scalarList cMinBeforeClipping(c.size(), 0.0);
scalarList cMeanBeforeClipping(c.size(), 0.0);
scalarList cMaxBeforeClipping(c.size(), 0.0);

const scalar domainVolume = gSum(mesh.V());

// Initialise all species, including the inert species.
// Active species values are overwritten immediately after their solves.
forAll(c, i)
{
    cMinBeforeClipping[i] =
        gMin(c[i].primitiveField());

    cMeanBeforeClipping[i] =
        gSum(c[i].primitiveField()*mesh.V())
       /domainVolume;

    cMaxBeforeClipping[i] =
        gMax(c[i].primitiveField());
}

// Solve concentration transport
forAll(c, i)
{
    if (i != inertIndex)
    {
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
            const dictionary& zDict = reactions.subDict("chargeNumber");
            zi = zDict.getOrDefault<scalar>(species[i], 0.0);
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
            fvc::interpolate(Di*zi*Fdim/(Rgas*T)*ci)
          * fvc::snGrad(phiE)
          * mesh.magSf()
        );

        fvScalarMatrix ciEqn
        (
            fvm::ddt(ci)
          + fvm::div(phi, ci)
         ==
            fvm::laplacian(Di, ci)
          + fvc::div(migrationFlux)
          + Rc[i]
          + fvOptions(ci)
        );

        ciEqn.relax();
        fvOptions.constrain(ciEqn);
        ciEqn.solve(ci.name());
        fvOptions.correct(ci);

        // Save the raw transport solution before positivity clipping.
        if (printCFieldsBeforeClipping)
        {
            cMinBeforeClipping[i] =
                gMin(ci.primitiveField());

            cMeanBeforeClipping[i] =
                gSum(ci.primitiveField()*mesh.V())
            /domainVolume;

            cMaxBeforeClipping[i] =
                gMax(ci.primitiveField());
        }
        
        // Do not clip here. Equilibrium reactions must see the raw transported
        // concentrations so deficits (negative values) are retained in their
        // conserved quantities, e.g. A = [H+] - [OH-].
        // ci.max(0.0);
        ci.correctBoundaryConditions();
    }
}
// Apply all fast equilibrium projections only after all species transport solves.
reactionMgr.applyEquilibriumReactions();

// Clip only after the equilibrium projection has had a chance to account for
// any local species deficit. Equilibrium species should already be positive;
// this also retains the existing non-negative safeguard for all other species.
forAll(c, i)
{
    if (i != inertIndex)
    {
        c[i].max(0.0);
        c[i].correctBoundaryConditions();
    }
}

if (pHActive)
{
    volScalarField& pH =
        pHPtr();

    const volScalarField& cHplus =
        c[HplusIndex];

    const dimensionedScalar standardConcentration
    (
        "standardConcentration",
        cHplus.dimensions(),
        1000.0
    );

    const dimensionedScalar concentrationFloor
    (
        "concentrationFloor",
        cHplus.dimensions(),
        SMALL
    );

    pH =
        -log10
        (
            max(cHplus, concentrationFloor)
           /standardConcentration
        );

    pH.correctBoundaryConditions();

    if (printPH)
    {
        Info<< "pH min/mean/max:"
            << gMin(pH)
            << "/" << gAverage(pH)
            << "/" << gMax(pH)
            << nl;
    }
}

