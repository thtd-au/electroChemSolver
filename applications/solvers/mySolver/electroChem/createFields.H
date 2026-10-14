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
Switch printKappa(false);
boolList printJDiffDirections(3, false);
bool printJDiff = false;
boolList printGradPhiEDirections(3, false);
bool printGradPhiE = false;

Switch writePhiE(false);
Switch writeKappa(false);
Switch writeJDiff(false);
Switch writeGradPhiE(false);
Switch writeE(false);
Switch writeRelE(false);

if (controlDict.found("potential"))
{
    const dictionary& potentialControlDict =
        controlDict.subDict("potential");

    if (potentialControlDict.found("print"))
    {
        const dictionary& d = potentialControlDict.subDict("print");

        printPhiE     = d.getOrDefault<Switch>("phiE", false);
        printKappa    = d.getOrDefault<Switch>("kappa", false);
        printJDiffDirections = potentialPrintDirections(d, "jDiff");
        printGradPhiEDirections = potentialPrintDirections(d, "gradPhiE");
        forAll(printJDiffDirections, cmpt)
        {
            printJDiff = printJDiff || printJDiffDirections[cmpt];
            printGradPhiE = printGradPhiE || printGradPhiEDirections[cmpt];
        }
    }

    if (potentialControlDict.found("write"))
    {
        const dictionary& d = potentialControlDict.subDict("write");

        writePhiE     = d.getOrDefault<Switch>("phiE", false);
        writeKappa    = d.getOrDefault<Switch>("kappa", false);
        writeJDiff    = d.getOrDefault<Switch>("jDiff", false);
        writeGradPhiE = d.getOrDefault<Switch>("gradPhiE", false);
        writeE = d.getOrDefault<Switch>("E", false);
        writeRelE = d.getOrDefault<Switch>("relE", false);
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

Switch diagConservationPrint(false);
Switch diagConservationLog(false);
word diagConservationMode("residuals"); // residuals | breakdown
wordList diagFlowPatches(0);
wordList diagElectrodePatches(0);

Switch diagSatLog(false);
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
        if (d.found("reactionMode"))
        {
            WarningInFunction
                << "diagnostics/transport/reactionMode is obsolete and ignored. "
                << "Transport diagnostics report the total transported reaction source. "
                << "Remove reactionMode from controlDict." << nl;
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

    if (diagDict.found("saturation"))
    {
        const dictionary& d = diagDict.subDict("saturation");

        diagSatLog = d.getOrDefault<Switch>("log", false);

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

// postProcess.H includes createFields without the normal diagnostic call sites.
(void)diagTransportActive;
(void)diagConservationActive;

const bool diagSaturationActive =
    diagSatLog
 || diagSatPrintSi
 || diagSatPrintSTot
 || diagSatPrintPEqi
 || diagSatPrintPEqTot
 || diagSatWriteSi
 || diagSatWriteSTot
 || diagSatWritePEqi
 || diagSatWritePEqTot;

bool diagDebyeActive =
    diagDebyePrint
 || diagDebyeWriteIonicSum
 || diagDebyeLog
 || diagDebyeWriteLambda
 || diagDebyeWriteE;



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
// Snapshots from the source used in the final accepted transport sweep.
List<scalarField> RcKinetic(species.size());
List<scalarField> RcEquilibrium(species.size());


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
// All dependencies now exist. Balance prescribed inlet patches before
// initial constraint evaluation; this does not project internal concentrations.
forAll(c, i)
{
    forAll(c[i].boundaryField(), patchi)
    {
        if (c[i].boundaryField()[patchi].type() == "chargeBalancedConcentration")
            c[i].boundaryFieldRef()[patchi].evaluate();
    }
}

// Reaction chemistry is constructed only after all concentration fields exist.
// reactionManager owns reaction definitions/stoichiometry/model instances, while
// the physical concentration fields remain owned by createFields/the solver.
reactionManager reactionMgr(mesh, reactions, species, c, inertIndex);

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
// ############################################ PotentialReaction potential ############################################
//

if (!electroChemistry.found("potential"))
{
    FatalIOErrorInFunction(electroChemistry)
        << "Missing required potential dictionary in constant/electroChemistry"
        << exit(FatalIOError);
}

const bool potentialActive =
    electroChemistry.subDict("potential").get<word>("model") != "none";

Info<< "Potential model: "
    << electroChemistry.subDict("potential").get<word>("model") << nl;

volScalarField phiE
(
    IOobject
    (
        "phiE",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        writePhiE ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh
);

if (!potentialActive)
{
    // Preserve patch types and electrode metadata while clearing the field.
    phiE = dimensionedScalar("zero", phiE.dimensions(), 0.0);
    phiE.correctBoundaryConditions();
}

// potentialManager performs runtime selection from electroChemistry/potential/model.
// The selected model owns electroneutrality evaluation and potential correction.
potentialManager potentialMgr
(
    mesh,
    electroChemistry,
    reactions,
    species,
    c,
    phiE,
    T,
    Rgas,
    Fdim
);


//
// ############################################ diagnostics fields  ############################################
//

// Debye resolution output is restricted to active-potential runs.
if (diagDebyeActive && !potentialMgr.active())
{
    WarningInFunction
        << "DebyeNumber diagnostics requested with potential/model none; "
        << "disabling Debye print, log and field output for this run." << nl;
    diagDebyeActive = false;
    diagDebyePrint = false;
    diagDebyeLog = false;
    diagDebyeWriteIonicSum = false;
    diagDebyeWriteLambda = false;
    diagDebyeWriteE = false;
}

// Canonical electroneutrality fields: same residual and scale as potentialMgr.
volScalarField E
(
    IOobject("E", runTime.timeName(), mesh, IOobject::NO_READ,
        writeE ? IOobject::AUTO_WRITE : IOobject::NO_WRITE),
    mesh, dimensionedScalar("zero", dimMoles/dimVolume, 0.0)
);
volScalarField relE
(
    IOobject("relE", runTime.timeName(), mesh, IOobject::NO_READ,
        writeRelE ? IOobject::AUTO_WRITE : IOobject::NO_WRITE),
    mesh, dimensionedScalar("zero", dimless, 0.0)
);
volVectorField gradPhiE
(
    IOobject("gradPhiE", runTime.timeName(), mesh, IOobject::NO_READ,
        writeGradPhiE ? IOobject::AUTO_WRITE : IOobject::NO_WRITE),
    mesh, dimensionedVector("zero", dimVoltage/dimLength, vector::zero)
);
volScalarField kappaOutput
(
    IOobject("kappa", runTime.timeName(), mesh, IOobject::NO_READ,
        writeKappa ? IOobject::AUTO_WRITE : IOobject::NO_WRITE),
    mesh, dimensionedScalar("zero", dimensionSet(-1,-3,3,0,0,2,0), 0.0)
);
volVectorField jDiffOutput
(
    IOobject("jDiff", runTime.timeName(), mesh, IOobject::NO_READ,
        writeJDiff ? IOobject::AUTO_WRITE : IOobject::NO_WRITE),
    mesh, dimensionedVector("zero", dimensionSet(0,-2,0,0,0,1,0), vector::zero)
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

