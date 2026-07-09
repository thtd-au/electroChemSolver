if (runWaterEquilibrium)
{
    const dictionary& molarMassDict =
        reactions.subDict("molarMass");

    dimensionedScalar MH
    (
        "MH",
        dimMass/dimMoles,
        molarMassDict.get<scalar>(species[HplusIndex])
    );

    dimensionedScalar MOH
    (
        "MOH",
        dimMass/dimMoles,
        molarMassDict.get<scalar>(species[OHminusIndex])
    );

    volScalarField YHBefore
    (
        IOobject
        (
            "YHBefore_waterEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Y[HplusIndex]
    );

    volScalarField YOHBefore
    (
        IOobject
        (
            "YOHBefore_waterEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Y[OHminusIndex]
    );

    volScalarField cH
    (
        IOobject
        (
            "cH_waterEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*Y[HplusIndex]/MH
    );

    volScalarField cOH
    (
        IOobject
        (
            "cOH_waterEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*Y[OHminusIndex]/MOH
    );

    volScalarField A
    (
        IOobject
        (
            "A_waterEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cH - cOH
    );

    volScalarField root
    (
        IOobject
        (
            "root_waterEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqrt(sqr(A) + 4.0*Kw)
    );

    cH = 0.5*(A + root);

    cOH = 0.5*(-A + root);

    Y[HplusIndex] =
        MH*cH/rho;

    Y[OHminusIndex] =
        MOH*cOH/rho;

    Y[HplusIndex].max(0.0);
    Y[OHminusIndex].max(0.0);

    Y[HplusIndex].correctBoundaryConditions();
    Y[OHminusIndex].correctBoundaryConditions();

    RYWaterEq[HplusIndex] =
        (Y[HplusIndex] - YHBefore)/runTime.deltaT();

    RYWaterEq[OHminusIndex] =
        (Y[OHminusIndex] - YOHBefore)/runTime.deltaT();

    Info<< "Water equilibrium correction applied." << endl;
}

