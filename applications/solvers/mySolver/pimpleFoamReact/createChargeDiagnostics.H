word chargeDiagnosticsMode("none");
word chargeDiagnosticsOutput("molarCharge");

Switch chargeDiagnosticsLogToInfo(true);
Switch chargeDiagnosticsWriteToFile(true);
Switch chargeDiagnosticsWriteField(false);

scalar chargeDiagnosticsWarningTolerance(SMALL);

if (controlDict.found("chargeDiagnostics"))
{
    const dictionary& cdDict = controlDict.subDict("chargeDiagnostics");

    chargeDiagnosticsMode =
        cdDict.getOrDefault<word>("mode", word("none"));

    chargeDiagnosticsOutput =
        cdDict.getOrDefault<word>("output", word("molarCharge"));

    chargeDiagnosticsLogToInfo =
        cdDict.getOrDefault<Switch>("logToInfo", true);

    chargeDiagnosticsWriteToFile =
        cdDict.getOrDefault<Switch>("writeToFile", true);

    chargeDiagnosticsWriteField =
        cdDict.getOrDefault<Switch>("writeField", false);

    chargeDiagnosticsWarningTolerance =
        cdDict.getOrDefault<scalar>("warningTolerance", SMALL);
}

if
(
    chargeDiagnosticsMode != "none"
 && chargeDiagnosticsMode != "minimal"
 && chargeDiagnosticsMode != "total"
 && chargeDiagnosticsMode != "full"
)
{
    FatalIOErrorInFunction(controlDict)
        << "chargeDiagnostics/mode must be one of: "
        << "none, minimal, total, full"
        << exit(FatalIOError);
}

if
(
    chargeDiagnosticsOutput != "molarCharge"
 && chargeDiagnosticsOutput != "chargeDensity"
 && chargeDiagnosticsOutput != "both"
)
{
    FatalIOErrorInFunction(controlDict)
        << "chargeDiagnostics/output must be one of: "
        << "molarCharge, chargeDensity, both"
        << exit(FatalIOError);
}

const bool chargeDiagnosticsActive =
    chargeDiagnosticsMode != "none";

List<scalar> chargeNumber(Y.size(), 0.0);

if (chargeDiagnosticsActive && reactions.found("chargeNumber"))
{
    const dictionary& zDict = reactions.subDict("chargeNumber");

    forAll(species, i)
    {
        if (zDict.found(species[i]))
        {
            chargeNumber[i] = zDict.get<scalar>(species[i]);
        }
    }
}

volScalarField c_charge
(
    IOobject
    (
        "c_charge",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        (
            chargeDiagnosticsWriteField
         && (
                chargeDiagnosticsOutput == "molarCharge"
             || chargeDiagnosticsOutput == "both"
            )
        )
      ? IOobject::AUTO_WRITE
      : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimMoles/dimVolume, 0.0)
);

volScalarField rho_e
(
    IOobject
    (
        "rho_e",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        (
            chargeDiagnosticsWriteField
         && (
                chargeDiagnosticsOutput == "chargeDensity"
             || chargeDiagnosticsOutput == "both"
            )
        )
      ? IOobject::AUTO_WRITE
      : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimCurrent*dimTime/dimVolume,
        0.0
    )
);

autoPtr<OFstream> chargeDiagnosticsFilePtr;

if
(
    chargeDiagnosticsActive
 && chargeDiagnosticsWriteToFile
 && Pstream::master()
)
{
    fileName cdDir(runTime.globalPath()/"postProcessing"/"chargeDiagnostics");
    mkDir(cdDir);

    chargeDiagnosticsFilePtr.reset
    (
        new OFstream(cdDir/"chargeDiagnostics.dat")
    );

    OFstream& os = chargeDiagnosticsFilePtr();

    os  << "# time "
        << "Ccharge_min Ccharge_mean Ccharge_max Ccharge_integral "
        << "rhoE_min rhoE_mean rhoE_max rhoE_integral"
        << nl;
}

