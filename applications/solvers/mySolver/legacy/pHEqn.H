if (runPH)
{
    dimensionedScalar cRef
    (
        "cRef",
        dimMoles/dimVolume,
        1000.0
    );

    pH =
       -log10(
            max
            (
                c[HplusIndex],
                dimensionedScalar
                (
                    "small",
                    dimMoles/dimVolume,
                    SMALL
                )
            )
           /cRef
        );

    pH.correctBoundaryConditions();

    const scalar pHMean = gSum(pH.primitiveField()*mesh.V())/gSum(mesh.V());

    Info<< "pH min/mean/max = "
        << min(pH).value()
        << " / "
        << pHMean
        << " / "
        << max(pH).value()
        << endl;

    // Info<< "internal pH min/max = "
    //     << gMin(pH.primitiveField()) << " / "
    //     << gMax(pH.primitiveField()) << nl;

    // Info<< "field pH min/max = "
    //     << min(pH).value() << " / "
    //     << max(pH).value() << endl;
}

