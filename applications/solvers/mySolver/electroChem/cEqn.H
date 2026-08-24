// cEqn.H
label inertIndex = species.find(inertSpecie);

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
    fvc::interpolate(gradPhiE) & mesh.Sf()
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
          * gradPhiEFlux
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

        // ci.max(0.0);
        ci.correctBoundaryConditions();
    }
}

// Store internal concentrations immediately before the fast equilibrium projection. 
List<scalarField> cBeforeEquilibrium(c.size());

forAll(c, i)
{
    cBeforeEquilibrium[i] = c[i].primitiveField();

    Req[i] =
        dimensionedScalar
        (
            "zero",
            dimMoles/dimVolume/dimTime,
            0.0
        );
}

// Store the water-equilibrium invariant before projection.
// Only allocate/populate this when the diagnostic is requested and a waterEquilibrium model actually exists.
List<scalarField> waterAPreclip;

if(diagWaterEquilibriumA  && reactionMgr.hasWaterEquilibrium())
{
    waterAPreclip.setSize
    (
        reactionMgr.nWaterEquilibriumReactions()
    );

    for
    (
        label reactioni = 0;
        reactioni < reactionMgr.nWaterEquilibriumReactions();
        ++reactioni
    )
    {
        const waterEquilibrium& waterEq =
            reactionMgr.waterEquilibriumReaction(reactioni);

        const volScalarField& cH =
            mesh.lookupObject<volScalarField>
            (
                "c_" + waterEq.HplusName()
            );

        const volScalarField& cOH =
            mesh.lookupObject<volScalarField>
            (
                "c_" + waterEq.OHminusName()
            );

        waterAPreclip[reactioni] = cH.primitiveField()- cOH.primitiveField();
    }
}


// Apply all fast equilibrium projections only after all species transport solves.
reactionMgr.applyEquilibriumReactions();

const scalar invDeltaT = 1.0/runTime.deltaTValue();

forAll(c, i)
{
    scalarField& ReqInternal = Req[i].primitiveFieldRef();

    ReqInternal =
        (c[i].primitiveField() - cBeforeEquilibrium[i])*invDeltaT;

    // Equilibrium source bookkeeping is intentionally internal-cell only.
    // Boundary changes caused by correctBoundaryConditions() belong to the
    // boundary transport balance, not the volume equilibrium source.
    forAll(Req[i].boundaryField(), patchi)
    {
        Req[i].boundaryFieldRef()[patchi] == 0.0;
    }
}

// Positivity enforcement only AFTER equilibrium projection.
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
    Info<< "pH fields:" << nl;

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
        Info<< "  pH min/mean/max:"
            << gMin(pH)
            << "/" << gAverage(pH)
            << "/" << gMax(pH)
            << nl;
    }
    Info<< endl;
}

