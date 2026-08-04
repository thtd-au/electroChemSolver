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
// ############################################ read reactions dict ############################################
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

const dictionary& molarMassDict = reactions.subDict("molarMass");

const dimensionedScalar rho
(
    "rho",
    dimDensity,
    laminarTransport
);

//
// ############################################ read controlDict ############################################
//

const dictionary& controlDict = runTime.controlDict();

// Species control
Switch printCFields(false);
Switch writeSourceFields(false);
Switch printYFields(false);
Switch writeYFields(false);

if (controlDict.found("species"))
{
    const dictionary& speciesControlDict =
        controlDict.subDict("species");

    if (speciesControlDict.found("c"))
    {
        const dictionary& cControlDict =
            speciesControlDict.subDict("c");

        printCFields =
            cControlDict.getOrDefault<Switch>("print", false);

        writeSourceFields =
            cControlDict.getOrDefault<Switch>("writeSources", false);
    }

    if (speciesControlDict.found("Y"))
    {
        const dictionary& YControlDict =
            speciesControlDict.subDict("Y");

        printYFields =
            YControlDict.getOrDefault<Switch>("print", false);

        writeYFields =
            YControlDict.getOrDefault<Switch>("write", false);
    }
}

// pH controls
const bool pHConfigured = controlDict.found("pH");

Switch printPH(false);
Switch writePH(false);

if (pHConfigured)
{
    const dictionary& pHDict =
        controlDict.subDict("pH");

    printPH =
        pHDict.getOrDefault<Switch>
        (
            "print",
            false
        );

    writePH =
        pHDict.getOrDefault<Switch>
        (
            "write",
            false
        );
}

const bool pHActive =
    pHConfigured && (printPH || writePH);

// Diagnostics controls
Switch diagTransportPrint(false);
Switch diagTransportLog(false);
wordList diagTransportSpecies(0);

Switch diagConservationPrint(false);
Switch diagConservationLog(false);
word diagConservationMode("full");
wordList diagFlowPatches(0);
wordList diagElectrodePatches(0);

word diagChargeOutput("molarCharge");  // molarCharge | chargeDensity | both
Switch diagChargePrint(false);
Switch diagChargeLog(false);
Switch diagChargeWrite(false);

Switch diagSatPrintSi(false);
Switch diagSatPrintSTot(false);
Switch diagSatPrintPEqi(false);
Switch diagSatPrintPEqTot(false);

Switch diagSatWriteSi(false);
Switch diagSatWriteSTot(false);
Switch diagSatWritePEqi(false);
Switch diagSatWritePEqTot(false);

Switch diagInertClosurePrint(false);

Switch diagDebyePrint(false);
Switch diagDebyeLog(false);
label diagDebyeDirection(1);
scalar diagDebyeMinIonicSum(0.0);
Switch diagDebyeWriteIonicSum(false);
Switch diagDebyeWriteLambda(false);
Switch diagDebyeWriteE(false);

if (controlDict.found("diagnostics"))
{
    const dictionary& diagDict = controlDict.subDict("diagnostics");

    if (diagDict.found("transport"))
    {
        const dictionary& d = diagDict.subDict("transport");

        diagTransportPrint = d.getOrDefault<Switch>("print", false);
        diagTransportLog   = d.getOrDefault<Switch>("log", false);

        if (d.found("species"))
        {
            diagTransportSpecies = wordList(d.lookup("species"));
        }
    }

    if (diagDict.found("conservation"))
    {
        const dictionary& d =
            diagDict.subDict("conservation");

        diagConservationPrint =
            d.getOrDefault<Switch>("print", false);

        diagConservationLog =
            d.getOrDefault<Switch>("log", false);

        diagConservationMode =
            d.getOrDefault<word>("mode", word("full"));

        if
        (
            diagConservationMode != "full"
        && diagConservationMode != "residuals"
        )
        {
            FatalIOErrorInFunction(d)
                << "diagnostics/conservation/mode must be "
                << "full or residuals"
                << exit(FatalIOError);
        }

        if (d.found("flowPatches"))
        {
            diagFlowPatches =
                wordList(d.lookup("flowPatches"));
        }

        if (d.found("electrodePatches"))
        {
            diagElectrodePatches =
                wordList(d.lookup("electrodePatches"));
        }
    }

    if (diagDict.found("charge"))
    {
        const dictionary& d = diagDict.subDict("charge");

        diagChargeOutput = d.getOrDefault<word>("output", word("molarCharge"));
        diagChargePrint  = d.getOrDefault<Switch>("print", false);
        diagChargeLog    = d.getOrDefault<Switch>("log", false);
        diagChargeWrite  = d.getOrDefault<Switch>("write", false);

        if
        (
            diagChargeOutput != "molarCharge"
         && diagChargeOutput != "chargeDensity"
         && diagChargeOutput != "both"
        )
        {
            FatalIOErrorInFunction(d)
                << "diagnostics/charge/output must be "
                << "molarCharge, chargeDensity, or both"
                << exit(FatalIOError);
        }
    }

    if (diagDict.found("saturation"))
    {
        const dictionary& d = diagDict.subDict("saturation");

        if (d.found("print"))
        {
            const dictionary& p = d.subDict("print");

            diagSatPrintSi     = p.getOrDefault<Switch>("Si", false);
            diagSatPrintSTot   = p.getOrDefault<Switch>("STot", false);
            diagSatPrintPEqi   = p.getOrDefault<Switch>("pEqi", false);
            diagSatPrintPEqTot = p.getOrDefault<Switch>("pEqTot", false);
        }

        if (d.found("write"))
        {
            const dictionary& w = d.subDict("write");

            diagSatWriteSi     = w.getOrDefault<Switch>("Si", false);
            diagSatWriteSTot   = w.getOrDefault<Switch>("STot", false);
            diagSatWritePEqi   = w.getOrDefault<Switch>("pEqi", false);
            diagSatWritePEqTot = w.getOrDefault<Switch>("pEqTot", false);
        }
    }

    if (diagDict.found("inertClosure"))
    {
        const dictionary& d =
            diagDict.subDict("inertClosure");

        diagInertClosurePrint =
            d.getOrDefault<Switch>("print", false);
    }

    if (diagDict.found("DebyeNumber"))
    {
        const dictionary& d =
            diagDict.subDict("DebyeNumber");

        diagDebyePrint =
            d.getOrDefault<Switch>("print", false);

        diagDebyeLog =
            d.getOrDefault<Switch>("log", false);

        diagDebyeDirection =
            d.getOrDefault<label>("direction", 1);

        if
        (
            diagDebyeDirection < 0
        || diagDebyeDirection > 2
        )
        {
            FatalIOErrorInFunction(d)
                << "diagnostics/DebyeNumber/direction must be "
                << "0 for x, 1 for y, or 2 for z"
                << exit(FatalIOError);
        }

        if (d.found("write"))
        {
            const dictionary& wd =
                d.subDict("write");

            diagDebyeWriteIonicSum =
                wd.getOrDefault<Switch>
                (
                    "ionicSum",
                    false
                );

            diagDebyeWriteLambda =
                wd.getOrDefault<Switch>
                (
                    "lambda_D",
                    false
                );

            diagDebyeWriteE =
                wd.getOrDefault<Switch>
                (
                    "E_D",
                    false
                );
        }
        
        diagDebyeMinIonicSum =
            d.getOrDefault<scalar>
            (
                "minIonicSum",
                0.0
            );

        if (diagDebyeMinIonicSum < 0.0)
        {
            FatalIOErrorInFunction(d)
                << "diagnostics/DebyeNumber/minIonicSum must be non-negative"
                << exit(FatalIOError);
        }
    }
}

const bool diagTransportActive =
    diagTransportPrint || diagTransportLog;

const bool diagConservationActive =
    diagConservationPrint || diagConservationLog;

const bool diagChargeActive =
    diagChargePrint || diagChargeLog || diagChargeWrite;

const bool diagSaturationActive =
    diagSatPrintSi
 || diagSatPrintSTot
 || diagSatPrintPEqi
 || diagSatPrintPEqTot
 || diagSatWriteSi
 || diagSatWriteSTot
 || diagSatWritePEqi
 || diagSatWritePEqTot;

 const bool diagDebyeActive =
    diagDebyePrint
 || diagDebyeLog
 || diagDebyeWriteLambda
 || diagDebyeWriteE;

(void)diagTransportActive;
(void)diagConservationActive;
(void)diagChargeActive;
(void)diagSaturationActive;
(void)diagDebyeActive;
(void)diagDebyeDirection;


//
// ############################################ create species fields ############################################
//

wordList species(reactions.lookup("species"));
word inertSpecie(reactions.lookup("inertSpecie"));

PtrList<volScalarField> c(species.size());
PtrList<volScalarField> Y(species.size());
PtrList<volScalarField> Rc(species.size());

forAll(species, i)
{

    const word& specieName = species[i];

    dimensionedScalar Mi
    (
        "M_" + specieName,
        dimMass/dimMoles,
        molarMassDict.get<scalar>(specieName)
    );

    Info<< "Reading species field c_" << species[i] << endl;

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
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

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
                IOobject::NO_READ,
                writeYFields
                ? IOobject::AUTO_WRITE
                : IOobject::NO_WRITE
            ),
            c[i]*Mi/rho
        )
    );

    Rc.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Rc_" + specieName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                writeSourceFields 
                ? IOobject::AUTO_WRITE 
                : IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimMoles/dimVolume/dimTime, 0.0)
        )
    );

}

// Optional pH field
label HplusIndex = -1;

if (pHActive)
{
    HplusIndex = species.find("Hplus");

    if (HplusIndex < 0)
    {
        FatalIOErrorInFunction(reactions)
            << "pH calculation is enabled, but species Hplus "
            << "is not listed in species." << nl
            << "Available species are: " << species
            << exit(FatalIOError);
    }
}

autoPtr<volScalarField> pHPtr;

if (pHActive)
{
    pHPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "pH",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                writePH
              ? IOobject::AUTO_WRITE
              : IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "pH",
                dimless,
                0.0
            )
        )
    );
}

//
// ############################################ constant gradient potential ############################################
//

dimensionSet dimVoltage(1, 2, -3, 0, 0, -1, 0);

dimensionedScalar phiERef
(
    "phiERef",
    dimVoltage,
    reactions.getOrDefault<scalar>("phiERef", 0.0)
);

dimensionedVector gradPhiE
(
    "gradPhiE",
    dimVoltage/dimLength,
    reactions.lookupOrDefault<vector>("gradPhiE", vector::zero)
);

dimensionedVector phiEOrigin
(
    "phiEOrigin",
    dimLength,
    reactions.lookupOrDefault<vector>("phiEOrigin", vector::zero)
);

volScalarField phiE
(
    IOobject
    (
        "phiE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phiERef + ((mesh.C() - phiEOrigin) & gradPhiE)
);

volVectorField gradPhiEField
(
    IOobject
    (
        "gradPhiEField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    gradPhiE
);

surfaceScalarField snGradPhiE
(
    IOobject
    (
        "snGradPhiE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (mesh.Sf()/mesh.magSf()) & gradPhiE
);

//
// ############################################ diagnostics fields  ############################################
//

const bool writeMolarCharge =
    diagChargeWrite
 && (
        diagChargeOutput == "molarCharge"
     || diagChargeOutput == "both"
    );

const bool writeChargeDensity =
    diagChargeWrite
 && (
        diagChargeOutput == "chargeDensity"
     || diagChargeOutput == "both"
    );

dimensionSet dimElectricCharge(0, 0, 1, 0, 0, 1, 0);

volScalarField molarCharge
(
    IOobject
    (
        "molarCharge",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeMolarCharge
        ? IOobject::AUTO_WRITE
        : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimMoles/dimVolume, 0.0)
);

volScalarField chargeDensity
(
    IOobject
    (
        "chargeDensity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeChargeDensity
        ? IOobject::AUTO_WRITE
        : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimElectricCharge/dimVolume,0.0)
);

volScalarField ionicSumField
(
    IOobject
    (
        "ionicSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        diagDebyeWriteIonicSum
        ? IOobject::AUTO_WRITE
        : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimMoles/dimVolume,
        0.0
    )
);

volScalarField lambdaD
(
    IOobject
    (
        "lambda_D",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        diagDebyeWriteLambda
        ? IOobject::AUTO_WRITE
        : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "lambda_D",
        dimLength,
        0.0
    )
);

volScalarField DebyeNumber
(
    IOobject
    (
        "E_D",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        diagDebyeWriteE
        ? IOobject::AUTO_WRITE
        : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "E_D",
        dimless,
        0.0
    )
);


// Saturation diagnostic fields

dimensionedScalar pAbsRef
(
    "absolutePressureReference",
    dimPressure,
    reactions.getOrDefault<scalar>("absolutePressureReference", 101325.0)
);

dimensionedScalar pSmall
(
    "pSmall",
    dimPressure,
    SMALL
);

volScalarField pEq_tot
(
    IOobject
    (
        "pEq_tot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        diagSatWritePEqTot ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
        diagSatWriteSTot ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

wordList henrySpecies(0);
labelList henrySpeciesIndex(0);
scalarList henryConstants(0);

if (diagSaturationActive)
{
    if (reactions.found("henryConstants") || reactions.found("henryConstant"))
    {
        const dictionary& henryDict =
            reactions.found("henryConstants")
          ? reactions.subDict("henryConstants")
          : reactions.subDict("henryConstant");

        henrySpecies = henryDict.toc();
        henrySpeciesIndex.setSize(henrySpecies.size(), -1);
        henryConstants.setSize(henrySpecies.size(), 0.0);

        forAll(henrySpecies, hi)
        {
            const word& specieName = henrySpecies[hi];
            const label speciei = species.find(specieName);

            if (speciei < 0)
            {
                FatalIOErrorInFunction(reactions)
                    << "Species " << specieName
                    << " is listed in henryConstant but not in species"
                    << exit(FatalIOError);
            }

            const scalar H = henryDict.get<scalar>(specieName);

            if (H <= SMALL)
            {
                FatalIOErrorInFunction(reactions)
                    << "Henry constant for species "
                    << specieName << " must be positive"
                    << exit(FatalIOError);
            }

            henrySpeciesIndex[hi] = speciei;
            henryConstants[hi] = H;
        }
    }
}

PtrList<volScalarField> pEqHenry(henrySpecies.size());
PtrList<volScalarField> saturation(henrySpecies.size());

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
                diagSatWritePEqi
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
                diagSatWriteSi
              ? IOobject::AUTO_WRITE
              : IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
}

#include "createMRF.H"
#include "createFvOptions.H"

