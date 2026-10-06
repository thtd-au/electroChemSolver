// YEqn.H -- post-convergence mass-fraction/inert closure
//
// Concentrations are the primary electrochemical state.  This file derives
// mass fractions only after joint Picard convergence and reconstructs the
// neutral inert concentration from the remaining mass fraction.  No clipping
// is permitted: an invalid closure is reported instead of silently changing
// the accepted electrochemical state.

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
        Y[i].correctBoundaryConditions();

        const scalar minYi = gMin(Y[i].primitiveField());
        if (minYi < 0.0)
        {
            FatalErrorInFunction
                << "Negative derived mass fraction for " << species[i]
                << ": min(Y)=" << minYi << nl
                << "PotentialReaction does not clip derived mass fractions."
                << exit(FatalError);
        }

        Yt += Y[i];
    }
}

if (inertIndex >= 0)
{
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].correctBoundaryConditions();

    const scalar minYInert = gMin(Y[inertIndex].primitiveField());
    if (minYInert < 0.0)
    {
        FatalErrorInFunction
            << "Inert mass-fraction closure is negative for "
            << species[inertIndex] << ": min(Y_inert)=" << minYInert << nl
            << "The non-inert mass fractions exceed unity; the state will not "
            << "be clipped or renormalized."
            << exit(FatalError);
    }
}

// Reconstruct the inert concentration from the accepted mass-fraction closure.
dimensionedScalar Minert
(
    "M_" + species[inertIndex],
    dimMass/dimMoles,
    molarMassDict.get<scalar>(species[inertIndex])
);

c[inertIndex] = rho*Y[inertIndex]/Minert;
c[inertIndex].correctBoundaryConditions();

