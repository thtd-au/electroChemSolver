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
// ############################################ read constant dictionaries ############################################
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

Info<< "Reading electroChemistry dictionary\n" << endl;

IOdictionary electroChemistry
(
    IOobject
    (
        "electroChemistry",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dictionary& molarMassDict = reactions.subDict("molarMass");

dimensionSet dimVoltage(1, 2, -3, 0, 0, -1, 0);
dimensionSet dimCurrentDensity(0, -2, 0, 0, 0, 1, 0);
dimensionSet dimConductivity(-1, -3, 3, 0, 0, 2, 0);
dimensionSet dimPermittivity(-1, -3, 4, 0, 0, 2, 0);

dimensionedScalar T
(
    "temperature",
    dimTemperature,
    electroChemistry.get<scalar>("temperature")
);

dimensionedScalar Rgas
(
    "R",
    dimensionSet(1, 2, -2, -1, -1, 0, 0),
    8.31446261815324
);

dimensionedScalar Fdim
(
    "FaradayConstant",
    dimensionSet(0, 0, 1, 0, -1, 1, 0),
    electroChemistry.get<scalar>("FaradayConstant")
);

dimensionedScalar epsilon
(
    "permittivity",
    dimPermittivity,
    electroChemistry
);

if (T.value() <= SMALL)
{
    FatalIOErrorInFunction(electroChemistry)
        << "temperature must be positive"
        << exit(FatalIOError);
}

if (Fdim.value() <= SMALL)
{
    FatalIOErrorInFunction(electroChemistry)
        << "FaradayConstant must be positive"
        << exit(FatalIOError);
}

if (epsilon.value() <= SMALL)
{
    FatalIOErrorInFunction(electroChemistry)
        << "permittivity must be positive"
        << exit(FatalIOError);
}

Info<< "Reading transportProperties dictionary\n" << endl;

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
Switch printCFieldsBeforeClipping(false);
Switch writeSourceFields(false);
Switch printYFields(false);
Switch writeYFields(false);

if (controlDict.found("species"))
{
    const dictionary& speciesControlDict =
        controlDict.subDict("species");

    if (speciesControlDict.found("c"))
    {
        const dictionary& cControlDict = speciesControlDict.subDict("c");

        printCFields = cControlDict.getOrDefault<Switch>("print", false);

        printCFieldsBeforeClipping = cControlDict.getOrDefault<Switch>("printBeforeClipping",false);

        writeSourceFields = cControlDict.getOrDefault<Switch>("writeSources", false);
    }

    if (speciesControlDict.found("Y"))
    {
        const dictionary& YControlDict = speciesControlDict.subDict("Y");

        printYFields = YControlDict.getOrDefault<Switch>("print", false);

        writeYFields = YControlDict.getOrDefault<Switch>("write", false);
    }
}

// Potential controls
Switch printPhiE(false);
Switch printGradPhiE(false);
Switch printKappa(false);
Switch printJDiff(false);

Switch writePhiE(false);
Switch writeGradPhiE(false);
Switch writeKappa(false);
Switch writeJDiff(false);

if (controlDict.found("potential"))
{
    const dictionary& potentialControlDict =
        controlDict.subDict("potential");

    if (potentialControlDict.found("print"))
    {
        const dictionary& d = potentialControlDict.subDict("print");

        printPhiE     = d.getOrDefault<Switch>("phiE", false);
        printGradPhiE = d.getOrDefault<Switch>("gradPhiE", false);
        printKappa    = d.getOrDefault<Switch>("kappa", false);
        printJDiff    = d.getOrDefault<Switch>("jDiff", false);
    }

    if (potentialControlDict.found("write"))
    {
        const dictionary& d = potentialControlDict.subDict("write");

        writePhiE     = d.getOrDefault<Switch>("phiE", false);
        writeGradPhiE = d.getOrDefault<Switch>("gradPhiE", false);
        writeKappa    = d.getOrDefault<Switch>("kappa", false);
        writeJDiff    = d.getOrDefault<Switch>("jDiff", false);
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

    printPH = pHDict.getOrDefault<Switch>("print",false);
    
    writePH = pHDict.getOrDefault<Switch>("write",false);
}

const bool pHActive =
    pHConfigured && (printPH || writePH);


// Diagnostics controls

scalar diagLogInterval(1.0);

Switch diagTransportPrint(false);
Switch diagTransportLog(false);
wordList diagTransportSpecies(0);
word diagTransportReactionMode("total"); // total | split

Switch diagConservationPrint(false);
Switch diagConservationLog(false);
word diagConservationMode("residuals"); // residuals | breakdown
wordList diagFlowPatches(0);
wordList diagElectrodePatches(0);

Switch diagElectrodeFluxStabilityPrint(false);
Switch diagElectrodeFluxStabilityLog(false);

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

    diagLogInterval =
        diagDict.getOrDefault<scalar>("logInterval", diagLogInterval);

    if (diagLogInterval <= SMALL)
    {
        FatalIOErrorInFunction(diagDict)
            << "diagnostics/logInterval must be positive"
            << exit(FatalIOError);
    }

    if (diagDict.found("transport"))
    {
        const dictionary& d = diagDict.subDict("transport");

        diagTransportPrint = d.getOrDefault<Switch>("print", false);
        diagTransportLog   = d.getOrDefault<Switch>("log", false);
        diagTransportReactionMode =
            d.getOrDefault<word>("reactionMode", word("total"));

        if
        (
            diagTransportReactionMode != "total"
         && diagTransportReactionMode != "split"
        )
        {
            FatalIOErrorInFunction(d)
                << "diagnostics/transport/reactionMode must be "
                << "total or split"
                << exit(FatalIOError);
        }

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
            d.getOrDefault<word>("mode", word("residuals"));

        if
        (
            diagConservationMode != "residuals"
         && diagConservationMode != "breakdown"
        )
        {
            FatalIOErrorInFunction(d)
                << "diagnostics/conservation/mode must be "
                << "residuals or breakdown"
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

    if (diagDict.found("electrodeFluxStability"))
    {
        const dictionary& d =
            diagDict.subDict("electrodeFluxStability");

        diagElectrodeFluxStabilityPrint =
            d.getOrDefault<Switch>("print", false);

        diagElectrodeFluxStabilityLog =
            d.getOrDefault<Switch>("log", false);
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

const bool diagElectrodeFluxStabilityActive =
    diagElectrodeFluxStabilityPrint
 || diagElectrodeFluxStabilityLog;

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
(void)diagElectrodeFluxStabilityActive;
(void)diagChargeActive;
(void)diagSaturationActive;
(void)diagDebyeActive;
(void)diagDebyeDirection;


//
// ############################################ create species fields ############################################
//

wordList species(reactions.lookup("species"));
word inertSpecie(reactions.lookup("inertSpecie"));
label inertIndex = species.find(inertSpecie);

if (inertIndex < 0)
{
    FatalIOErrorInFunction(reactions)
        << "inertSpecie " << inertSpecie
        << " is not listed in species." << nl
        << "Available species are: " << species
        << exit(FatalIOError);
}

if (reactions.found("chargeNumber"))
{
    const dictionary& zDict = reactions.subDict("chargeNumber");
    const scalar zInert = zDict.getOrDefault<scalar>(inertSpecie, 0.0);

    if (mag(zInert) > SMALL)
    {
        FatalIOErrorInFunction(reactions)
            << "Charged inert species are not supported by the "
            << "current mass-fraction closure." << nl
            << "inertSpecie " << inertSpecie
            << " has charge number " << zInert << '.' << nl
            << "Choose an electrically neutral inert species."
            << exit(FatalIOError);
    }
}

PtrList<volScalarField> c(species.size());
PtrList<volScalarField> Y(species.size());
PtrList<volScalarField> Rc(species.size());
PtrList<volScalarField> Req(species.size());


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

    Req.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Req_" + specieName,
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
// Initial species inventories used by conservation diagnostics.
// Store these before the first time step so the conservation balance
// is anchored to the actual case initial condition.
scalarList initialSpeciesInventory(c.size(), 0.0);

forAll(c, i)
{
    initialSpeciesInventory[i] =
        gSum(c[i].primitiveField()*mesh.V());
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
// ############################################ simple-current potential ############################################
//

if (!electroChemistry.found("potential"))
{
    FatalIOErrorInFunction(electroChemistry)
        << "Missing required potential dictionary in constant/electroChemistry"
        << exit(FatalIOError);
}

const dictionary& potentialDict =
    electroChemistry.subDict("potential");

word potentialModel
(
    potentialDict.getOrDefault<word>("model", word("simpleCurrent"))
);

if (potentialModel != "simpleCurrent")
{
    FatalIOErrorInFunction(potentialDict)
        << "Unsupported potential model " << potentialModel << nl
        << "Currently supported model: simpleCurrent"
        << exit(FatalIOError);
}

dimensionedScalar jApplied
(
    "appliedCurrentDensity",
    dimCurrentDensity,
    potentialDict.get<scalar>("appliedCurrentDensity")
);

vector electrodeNormal
(
    potentialDict.get<vector>("electrodeNormal")
);

const scalar electrodeNormalMag = mag(electrodeNormal);

if (electrodeNormalMag <= SMALL)
{
    FatalIOErrorInFunction(potentialDict)
        << "electrodeNormal must be non-zero"
        << exit(FatalIOError);
}

electrodeNormal /= electrodeNormalMag;

label potentialDirection = 0;
scalar largestNormalComponent = mag(electrodeNormal[0]);

for (label dir = 1; dir < vector::nComponents; ++dir)
{
    if (mag(electrodeNormal[dir]) > largestNormalComponent)
    {
        potentialDirection = dir;
        largestNormalComponent = mag(electrodeNormal[dir]);
    }
}

vector alignedNormal(vector::zero);
alignedNormal[potentialDirection] =
    electrodeNormal[potentialDirection] >= 0.0 ? 1.0 : -1.0;

const scalar directionTolerance = 1e-8;

if (mag(electrodeNormal - alignedNormal) > directionTolerance)
{
    FatalIOErrorInFunction(potentialDict)
        << "simpleCurrent currently requires electrodeNormal to be "
        << "aligned with a Cartesian mesh direction." << nl
        << "Specified electrodeNormal = " << electrodeNormal
        << exit(FatalIOError);
}

electrodeNormal = alignedNormal;

word phiEReferencePatchName
(
    potentialDict.get<word>("referencePatch")
);

label phiEReferencePatch =
    mesh.boundaryMesh().findPatchID(phiEReferencePatchName);

if (phiEReferencePatch < 0)
{
    FatalIOErrorInFunction(potentialDict)
        << "Unknown potential referencePatch "
        << phiEReferencePatchName
        << exit(FatalIOError);
}

dimensionedScalar phiERef
(
    "referencePotential",
    dimVoltage,
    potentialDict.get<scalar>("referencePotential")
);

dimensionedScalar minConductivity
(
    "minConductivity",
    dimConductivity,
    potentialDict.getOrDefault<scalar>("minConductivity", SMALL)
);

if (minConductivity.value() < 0.0)
{
    FatalIOErrorInFunction(potentialDict)
        << "minConductivity must be non-negative"
        << exit(FatalIOError);
}

if (Pstream::parRun())
{
    FatalErrorInFunction
        << "The first simpleCurrent phiE reconstruction implementation "
        << "supports serial structured meshes only." << nl
        << "Parallel column reconstruction has not yet been implemented."
        << exit(FatalError);
}

const fvsPatchVectorField& refSf =
    mesh.Sf().boundaryField()[phiEReferencePatch];

vector refAreaVector(vector::zero);

forAll(refSf, facei)
{
    refAreaVector += refSf[facei];
}

if (mag(refAreaVector) <= SMALL)
{
    FatalIOErrorInFunction(potentialDict)
        << "referencePatch " << phiEReferencePatchName
        << " has zero resultant area vector"
        << exit(FatalIOError);
}

const vector referencePatchNormal =
    refAreaVector/mag(refAreaVector);

if
(
    mag(mag(referencePatchNormal & electrodeNormal) - 1.0)
  > directionTolerance
)
{
    FatalIOErrorInFunction(potentialDict)
        << "referencePatch " << phiEReferencePatchName
        << " is not normal to electrodeNormal." << nl
        << "reference patch normal = " << referencePatchNormal << nl
        << "electrodeNormal = " << electrodeNormal
        << exit(FatalIOError);
}

volScalarField phiE
(
    IOobject
    (
        "phiE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writePhiE 
        ? IOobject::AUTO_WRITE 
        : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimVoltage, 0.0)
);

volVectorField gradPhiE
(
    IOobject
    (
        "gradPhiE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeGradPhiE 
        ? IOobject::AUTO_WRITE 
        : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVoltage/dimLength, vector::zero)
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeKappa 
        ? IOobject::AUTO_WRITE 
        : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimConductivity, 0.0)
);

volVectorField jDiff
(
    IOobject
    (
        "jDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeJDiff
         ? IOobject::AUTO_WRITE
         : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimCurrentDensity, vector::zero)
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

