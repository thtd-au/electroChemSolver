#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

//
// ############################################ create reaction dict ############################################
//

Info<< "Reading reactions dictionary\n" << endl;

IOdictionary reactions
(
    IOobject
    (
        "reactions",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dictionary& controlDict = runTime.controlDict();

Switch writeSources(false);

if (runTime.controlDict().found("reaction"))
{
    const dictionary& reactionDict =
        runTime.controlDict().subDict("reaction");

    writeSources =
        reactionDict.getOrDefault<Switch>
        (
            "writeSources",
            false
        );
}

Switch runConcentration(false);
Switch writeCFields(false);

if (controlDict.found("concentration"))
{
    const dictionary& cDict = controlDict.subDict("concentration");

    runConcentration = cDict.getOrDefault<Switch>("runConcentration", false);
    writeCFields = cDict.getOrDefault<Switch>("writeCFields", false);
}

Switch runSaturation(false);
Switch writeSFields(false);
Switch writeSTotField(false);
Switch writePEqFields(false);
Switch writePEqTotField(false);

if (controlDict.found("saturation"))
{
    const dictionary& sDict = controlDict.subDict("saturation");

    runSaturation = sDict.getOrDefault<Switch>("runSaturation", false);
    writeSFields = sDict.getOrDefault<Switch>("writeSFields", false);
    writeSTotField = sDict.getOrDefault<Switch>("writeSTotField", false);
    writePEqFields = sDict.getOrDefault<Switch>("writePEqFields", false);
    writePEqTotField = sDict.getOrDefault<Switch>("writePEqTotField", false);
}

if (runSaturation && !runConcentration)
{
    FatalIOErrorInFunction(controlDict)
        << "saturation/runSaturation is true, but "
        << "concentration/runConcentration is false."
        << exit(FatalIOError);
}

Switch printYFields(true);
Switch printDiagnostics(true);

if (controlDict.found("species"))
{
    const dictionary& speciesDict = controlDict.subDict("species");

    printYFields = speciesDict.getOrDefault<Switch>("printYFields", true);
    printDiagnostics = speciesDict.getOrDefault<Switch>("printDiagnostics", true);
}

//
// ############################################ create species mass-fraction field ############################################
//

wordList species(reactions.lookup("species"));
word inertSpecie(reactions.lookup("inertSpecie"));

PtrList<volScalarField> Y(species.size());
PtrList<volScalarField> RY(species.size());
PtrList<volScalarField> RYWaterEq(Y.size());

forAll(species, i)
{
    Info<< "Reading species field Y_" << species[i] << endl;
    
    Y.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Y_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    RY.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "RY_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                writeSources
                ? IOobject::AUTO_WRITE
                : IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless/dimTime, 0.0)
        )
    );

    RYWaterEq.set // Pseudo source from water equilibrium used for species conservation
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "RYWaterEq_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless/dimTime, 0.0)
        )
    );
}

//
// ############################################ convert to concentration ############################################
//

const dictionary& molarMassDict = reactions.subDict("molarMass");

const dimensionedScalar rho
(
    "rho",
    dimDensity,
    laminarTransport
);

PtrList<volScalarField> c(Y.size());

if (runConcentration)
{
    forAll(Y, i)
    {
        const volScalarField& Yi = Y[i];

        dimensionedScalar Mi
        (
            "M_" + species[i],
            dimMass/dimMoles,
            molarMassDict.get<scalar>(species[i])
        );

        c.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "c_" + species[i],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    writeCFields
                      ? IOobject::AUTO_WRITE
                      : IOobject::NO_WRITE
                ),
                rho*Yi/Mi
            )
        );
    }
}

//
// ############################################ water equilibrium ############################################
//

const bool runWaterEquilibrium = 
    reactions.found("reactions")
 && reactions.subDict("reactions").found("waterEquilibrium");

label HplusIndex(-1);
label OHminusIndex(-1);

dimensionedScalar Kw
(
    "Kw",
    sqr(dimMoles/dimVolume),
    0.0
);

if (runWaterEquilibrium)
{
    const dictionary& waterDict =
        reactions.subDict("reactions").subDict("waterEquilibrium");

    word HplusName
    (
        waterDict.lookup("Hplus")
    );

    word OHminusName
    (
        waterDict.lookup("OHminus")
    );

    HplusIndex = species.find(HplusName);
    OHminusIndex = species.find(OHminusName);

    if (HplusIndex < 0)
    {
        FatalIOErrorInFunction(waterDict)
            << "Unknown Hplus species "
            << HplusName
            << exit(FatalIOError);
    }

    if (OHminusIndex < 0)
    {
        FatalIOErrorInFunction(waterDict)
            << "Unknown OHminus species "
            << OHminusName
            << exit(FatalIOError);
    }

    Kw =
        dimensionedScalar
        (
            "Kw",
            sqr(dimMoles/dimVolume),
            waterDict.get<scalar>("Kw")
        );

    Info<< "Water equilibrium enabled" << nl
        << "    Hplus   : " << HplusName << nl
        << "    OHminus : " << OHminusName << nl
        << "    Kw      : " << Kw.value()
        << " (mol/m3)^2" << endl;
}
else
{
    Info<< "Water equilibrium disabled" << endl;
}

Switch runPH(false);
Switch writePHField(false);

if (runPH && !runConcentration)
{
    FatalIOErrorInFunction(controlDict)
        << "pH/runPH is true, but concentration/runConcentration is false."
        << exit(FatalIOError);
}

if (controlDict.found("pH"))
{
    const dictionary& pHDict =
        controlDict.subDict("pH");

    runPH =
        pHDict.getOrDefault<Switch>
        (
            "runPH",
            false
        );

    writePHField =
        pHDict.getOrDefault<Switch>
        (
            "writePHField",
            false
        );
}

volScalarField pH
(
    IOobject
    (
        "pH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writePHField
      ? IOobject::AUTO_WRITE
      : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("pH", dimless, 7.0)
);

//
// ############################################ saturation diagnostics  ############################################
//

// Henry-law saturation diagnostics setup
dimensionedScalar pAbsRef
(
    "absolutePressureReference",
    dimPressure,
    reactions.getOrDefault<scalar>("absolutePressureReference", 101325.0)
);

volScalarField pEq_tot
(
    IOobject
    (
        "pEq_tot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writePEqTotField
          ? IOobject::AUTO_WRITE
          : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimPressure, 0.0)
);

volScalarField S_tot
(
    IOobject
    (
        "S_tot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeSTotField
            ? IOobject::AUTO_WRITE
            : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

wordList henrySpecies(0);
labelList henrySpeciesIndex(0);
scalarList henryConstants(0);
if (runSaturation)
{
    if (reactions.found("henryConstants") || reactions.found("henryConstant"))
    {
        const dictionary& henryDict = reactions.found("henryConstants")
        ? reactions.subDict("henryConstants")
        : reactions.subDict("henryConstant");

        henrySpecies = henryDict.toc();
        henrySpeciesIndex.setSize(henrySpecies.size(), -1);
        henryConstants.setSize(henrySpecies.size(), 0.0);

        Info<< "Creating Henry saturation diagnostics for species:" << nl;

        forAll(henrySpecies, hi)
        {
            const word& specieName = henrySpecies[hi];
            const label speciei = species.find(specieName);

            if (speciei < 0)
            {
                FatalIOErrorInFunction(reactions)
                    << "Species " << specieName
                    << " is listed in henryConstants but not in species"
                    << exit(FatalIOError);
            }

            const scalar H = henryDict.get<scalar>(specieName);

            if (H <= SMALL)
            {
                FatalIOErrorInFunction(reactions)
                    << "Henry constant for species " << specieName
                    << " must be positive"
                    << exit(FatalIOError);
            }

            henrySpeciesIndex[hi] = speciei;
            henryConstants[hi] = H;

            Info<< "  " << specieName << ": H = " << H
                << " mol/(m3 Pa)" << nl;
        }

        Info<< "absolutePressureReference = " << pAbsRef.value()
            << " Pa" << nl << endl;
    }
    else
    {
        Info<< "Henry saturation diagnostics: none" << nl << endl;
    }
}

PtrList<volScalarField> pEqHenry(henrySpecies.size());
PtrList<volScalarField> saturation(henrySpecies.size());

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

    pEqHenry.set
    (
        hi,
        new volScalarField
        (
            IOobject
            (
                "pEq_" + specieName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                writePEqFields
                ? IOobject::AUTO_WRITE
                : IOobject::NO_WRITE
            ),
            c[speciei]/Hi
        )
    );

    saturation.set
    (
        hi,
        new volScalarField
        (
            IOobject
            (
                "S_" + specieName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                writeSFields
                ? IOobject::AUTO_WRITE
                : IOobject::NO_WRITE
            ),
            pEqHenry[hi]/max(pAbs, pSmall)
        )
    );

    pEq_tot += pEqHenry[hi];
}
S_tot = pEq_tot/max(pAbs, pSmall);

pEq_tot.correctBoundaryConditions();
S_tot.correctBoundaryConditions();


#include "createMRF.H"
#include "createFvOptions.H"

