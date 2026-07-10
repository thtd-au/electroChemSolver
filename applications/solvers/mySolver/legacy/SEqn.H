Info<< "Updating Henry saturation diagnostics" << endl;

if (henrySpecies.size())
{
    dimensionedScalar pSmall
    (
        "pSmall",
        dimPressure,
        SMALL
    );

    volScalarField pAbs
    (
        IOobject
        (
            "pAbs_tmp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pAbsRef + rho*p
    );

    pEq_tot = dimensionedScalar("zero", dimPressure, 0.0);

    forAll(henrySpecies, hi)
    {
        const label speciei = henrySpeciesIndex[hi];
        const word& specieName = henrySpecies[hi];

        dimensionedScalar Hi
        (
            "H_" + specieName,
            dimMoles/dimVolume/dimPressure,
            henryConstants[hi]
        );

        pEqHenry[hi] = c[speciei]/Hi;
        pEqHenry[hi].correctBoundaryConditions();

        saturation[hi] = pEqHenry[hi]/max(pAbs, pSmall);
        saturation[hi].correctBoundaryConditions();

        pEq_tot += pEqHenry[hi];

        Info<< saturation[hi].name()
            << " min/max = "
            << min(saturation[hi]).value() << " / "
            << max(saturation[hi]).value()
            << endl;
    }
    S_tot = pEq_tot/max(pAbs, pSmall);

    pEq_tot.correctBoundaryConditions();    
    S_tot.correctBoundaryConditions();

    // Info<< pEq_tot.name()
    //     << " min/max = "
    //     << min(pEq_tot).value() << " / "
    //     << max(pEq_tot).value()
    //     << " Pa" << endl;

    Info<< S_tot.name()
        << " min/max = "
        << min(S_tot).value() << " / "
        << max(S_tot).value()
        << endl;
}

