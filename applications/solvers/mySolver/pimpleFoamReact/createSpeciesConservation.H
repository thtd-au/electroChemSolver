// Species conservation monitor setup

word speciesConservationMode("none");     // none | total | full
Switch speciesConservationLogToInfo(true);
Switch speciesConservationWriteToFile(true);
Switch speciesConservationPrintPatchBreakdown(false);
Switch speciesConservationAllPatchDiagnostics(true);
scalar speciesConservationWarningTolerance(SMALL);

// Patch roles used for the residual accounting
// 0 = ignored/wall/unclassified, 1 = flow, 2 = electrode
List<label> speciesConservationPatchRole(mesh.boundary().size(), 0);

if (controlDict.found("speciesConservation"))
{
    const dictionary& scDict = controlDict.subDict("speciesConservation");

    speciesConservationMode = scDict.getOrDefault<word>
    (
        "mode",
        word("none")
    );

    speciesConservationLogToInfo = scDict.getOrDefault<Switch>
    (
        "logToInfo",
        true
    );

    speciesConservationWriteToFile = scDict.getOrDefault<Switch>
    (
        "writeToFile",
        true
    );

    speciesConservationPrintPatchBreakdown = scDict.getOrDefault<Switch>
    (
        "printPatchBreakdown",
        false
    );

    speciesConservationAllPatchDiagnostics = scDict.getOrDefault<Switch>
    (
        "allPatchDiagnostics",
        true
    );

    speciesConservationWarningTolerance = scDict.getOrDefault<scalar>
    (
        "warningTolerance",
        SMALL
    );

    if
    (
        speciesConservationMode != "none"
     && speciesConservationMode != "total"
     && speciesConservationMode != "full"
    )
    {
        FatalIOErrorInFunction(scDict)
            << "speciesConservation/mode must be one of: none, total, full"
            << exit(FatalIOError);
    }

    if (scDict.found("flowPatches"))
    {
        wordList flowPatchNames(scDict.lookup("flowPatches"));

        forAll(flowPatchNames, i)
        {
            const label patchi = mesh.boundaryMesh().findPatchID(flowPatchNames[i]);

            if (patchi < 0)
            {
                FatalIOErrorInFunction(scDict)
                    << "Unknown speciesConservation flow patch "
                    << flowPatchNames[i] << exit(FatalIOError);
            }

            speciesConservationPatchRole[patchi] = 1;
        }
    }

    if (scDict.found("electrodePatches"))
    {
        wordList electrodePatchNames(scDict.lookup("electrodePatches"));

        forAll(electrodePatchNames, i)
        {
            const label patchi = mesh.boundaryMesh().findPatchID(electrodePatchNames[i]);

            if (patchi < 0)
            {
                FatalIOErrorInFunction(scDict)
                    << "Unknown speciesConservation electrode patch "
                    << electrodePatchNames[i] << exit(FatalIOError);
            }

            speciesConservationPatchRole[patchi] = 2;
        }
    }
}

const bool speciesConservationActive = speciesConservationMode != "none";

List<scalar> speciesConservationN0(Y.size(), 0.0);
List<scalar> speciesConservationAdvIn(Y.size(), 0.0);
List<scalar> speciesConservationAdvOut(Y.size(), 0.0);
List<scalar> speciesConservationDiffIn(Y.size(), 0.0);
List<scalar> speciesConservationDiffOut(Y.size(), 0.0);
List<scalar> speciesConservationFlowIn(Y.size(), 0.0);
List<scalar> speciesConservationFlowOut(Y.size(), 0.0);
List<scalar> speciesConservationGen(Y.size(), 0.0);
List<scalar> speciesConservationWaterEqGen(Y.size(), 0.0);

autoPtr<OFstream> speciesConservationFilePtr;

if (speciesConservationActive)
{
    const dimensionedScalar rho
    (
        "rho",
        dimDensity,
        laminarTransport
    );

    const dictionary& molarMassDict = reactions.subDict("molarMass");

    forAll(Y, i)
    {
        const scalar Mi = molarMassDict.get<scalar>(species[i]);

        speciesConservationN0[i] =
            rho.value()/Mi*gSum(Y[i].primitiveField()*mesh.V());
    }

    if (speciesConservationLogToInfo)
    {
        Info<< "Species conservation monitor: " << speciesConservationMode
            << nl << endl;
    }

    if (speciesConservationWriteToFile && Pstream::master())
    {
        fileName scDir(runTime.globalPath()/"postProcessing"/"speciesConservation");
        mkDir(scDir);

        speciesConservationFilePtr.reset
        (
            new OFstream(scDir/"speciesConservation.dat")
        );

        OFstream& os = speciesConservationFilePtr();

        if (speciesConservationMode == "full")
        {
            os  << "# time species N0 Nt dN "
                << "flowIn flowOut advIn advOut diffIn diffOut gen waterEqGen"
                << "residual relativeResidual" << nl;
        }
        else
        {
            os  << "# time species dN flowIn flowOut gen waterEqGen"
                << "residual relativeResidual" << nl;
        }
    }
}
else
{
    Info<< "Species conservation monitor: none" << nl << endl;
}

