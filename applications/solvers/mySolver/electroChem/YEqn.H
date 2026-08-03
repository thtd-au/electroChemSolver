// YEqn.H
volScalarField Yt
(
    IOobject
    (
        "Yt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

forAll(Y, i)
{
    dimensionedScalar Mi
    (
        "M_" + species[i],
        dimMass/dimMoles,
        molarMassDict.get<scalar>(species[i])
    );

    if (i != inertIndex)
    {
        Y[i] = c[i]*Mi/rho;
        Y[i].max(0.0);
        Y[i].correctBoundaryConditions();

        Yt += Y[i];
    }
}

if (inertIndex >= 0)
{
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
    Y[inertIndex].correctBoundaryConditions();
}

// reconstruct inert concentration field
dimensionedScalar Minert
(
    "M_" + species[inertIndex],
    dimMass/dimMoles,
    molarMassDict.get<scalar>(species[inertIndex])
);

c[inertIndex] = rho*Y[inertIndex]/Minert;

c[inertIndex].correctBoundaryConditions();

