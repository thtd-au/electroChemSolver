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
Switch printE(false);
Switch printKappa(false);
Switch printJDiff(false);

Switch writePhiE(true);
Switch writeE(false);
Switch writeKappa(false);
Switch writeJDiff(false);

if (controlDict.found("potential"))
{
    const dictionary& potentialControlDict = controlDict.subDict("potential");

    if (potentialControlDict.found("print"))
    {
        const dictionary& d = potentialControlDict.subDict("print");

        printPhiE     = d.getOrDefault<Switch>("phiE", false);
        printE        = d.getOrDefault<Switch>("E", false);
        printKappa    = d.getOrDefault<Switch>("kappa", false);
        printJDiff    = d.getOrDefault<Switch>("jDiff", false);
    }

    if (potentialControlDict.found("write"))
    {
        const dictionary& d = potentialControlDict.subDict("write");

        writePhiE     = d.getOrDefault<Switch>("phiE", true);
        writeE        = d.getOrDefault<Switch>("E", false);
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
    const dictionary& pHDict = controlDict.subDict("pH");
    printPH = pHDict.getOrDefault<Switch>("print",false);
    writePH = pHDict.getOrDefault<Switch>("write",false);
}

const bool pHActive = pHConfigured && (printPH || writePH);


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

Switch diagWaterEquilibriumA(false);

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

    if (diagDict.found("waterEquilibrium"))
    {
        const dictionary& d =
            diagDict.subDict("waterEquilibrium");

        diagWaterEquilibriumA =
            d.getOrDefault<Switch>("diagA", false);
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

PtrList<volScalarField> cPicardLagged(species.size());

forAll(c, i)
{
    cPicardLagged.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "cPicardLagged_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "zero",
                dimMoles/dimVolume,
                0.0
            )
        )
    );

    // Copy internal values without triggering the physical BC.
    cPicardLagged[i].primitiveFieldRef() =
        c[i].primitiveField();

    // Copy the CURRENT boundary values explicitly.  These are the values
    // to be frozen for this Picard iteration.
    forAll(c[i].boundaryField(), patchi)
    {
        cPicardLagged[i].boundaryFieldRef()[patchi] ==
            c[i].boundaryField()[patchi];
    }
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

// Store concentration statistics before positivity clipping.
// These are printed later in diagnostics.H so that solver output
// is not interrupted by diagnostic messages.
scalarList cMinBeforeClipping(c.size(), 0.0);
scalarList cMeanBeforeClipping(c.size(), 0.0);
scalarList cMaxBeforeClipping(c.size(), 0.0);

// Store the water-equilibrium invariant before projection.
// Only allocate/populate this when the diagnostic is requested and a waterEquilibrium model actually exists.
List<scalarField> waterAPreclip;

//
// ############################################ electroneutralGalvanostatic galvanostatic potential ############################################
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
    potentialDict.getOrDefault<word>("model", word("electroneutralGalvanostatic"))
);

if (potentialModel != "electroneutralGalvanostatic")
{
    FatalIOErrorInFunction(potentialDict)
        << "Unsupported potential model " << potentialModel << nl
        << "Currently supported model: electroneutralGalvanostatic"
        << exit(FatalIOError);
}

dimensionedScalar jApplied
(
    "appliedCurrentDensity",
    dimCurrentDensity,
    potentialDict.get<scalar>("appliedCurrentDensity")
);

if (jApplied.value() < 0.0)
{
    FatalIOErrorInFunction(potentialDict)
        << "appliedCurrentDensity is defined as a magnitude and must be "
        << "non-negative"
        << exit(FatalIOError);
}

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

// Authoritative face fields for the ENPhi current discretisation.
// They are persistent/registered because galvanostaticPotentialFlux
// accesses them through the mesh objectRegistry.
surfaceScalarField jDiffFlux
(
    IOobject
    (
        "jDiffFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimCurrentDensity*dimArea,
        0.0
    )
);

surfaceScalarField kappaFace
(
    IOobject
    (
        "kappaFace",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimConductivity,
        0.0
    )
);

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

// The electroneutralGalvanostatic equation with flux boundary conditions determines phiE only up
// to an additive constant. Anchor this gauge at a physically named boundary
// patch. The patch itself keeps its normal-flux boundary condition; only one
// owner cell adjacent to the patch is used as the algebraic reference.
if (!potentialDict.found("referencePatch"))
{
    FatalIOErrorInFunction(potentialDict)
        << "Missing required referencePatch entry." << nl
        << "For example: referencePatch referenceElectrode;"
        << exit(FatalIOError);
}
const word phiEReferencePatchName
(
    potentialDict.get<word>("referencePatch")
);

const label phiEReferencePatchi = mesh.boundaryMesh().findPatchID(phiEReferencePatchName);

if (phiEReferencePatchi < 0)
{
    FatalIOErrorInFunction(potentialDict)
        << "referencePatch '" << phiEReferencePatchName
        << "' was not found in the mesh boundary." << nl
        << "Available patches are: " << mesh.boundaryMesh().names()
        << exit(FatalIOError);
}

const fvPatch& phiEReferencePatch = mesh.boundary()[phiEReferencePatchi];

// Determine the global area-weighted patch centroid. This gives a
// decomposition-independent physical location on the named reference patch.
const fvsPatchVectorField& phiEReferenceCf = mesh.Cf().boundaryField()[phiEReferencePatchi];

const fvsPatchScalarField& phiEReferenceMagSf = mesh.magSf().boundaryField()[phiEReferencePatchi];

scalar phiEReferenceArea = 0.0;
vector phiEReferenceAreaMoment(vector::zero);

forAll(phiEReferenceCf, facei)
{
    const scalar area = phiEReferenceMagSf[facei];
    phiEReferenceArea += area;
    phiEReferenceAreaMoment += area*phiEReferenceCf[facei];
}

reduce(phiEReferenceArea, sumOp<scalar>());
reduce(phiEReferenceAreaMoment, sumOp<vector>());

if (phiEReferenceArea <= SMALL)
{
    FatalErrorInFunction
        << "referencePatch '" << phiEReferencePatchName
        << "' has zero global area."
        << exit(FatalError);
}

const point phiEReferencePatchCentre
(
    phiEReferenceAreaMoment/phiEReferenceArea
);

// Select the patch face whose centre is closest to the global patch centre,
// then use its owner cell as the single algebraic gauge cell. In parallel,
// only the processor owning the selected face receives a non-negative cell.
scalar phiEReferenceLocalDistSqr = GREAT;
label phiEReferenceLocalFace = -1;

forAll(phiEReferenceCf, facei)
{
    const scalar distSqr =
        magSqr(phiEReferenceCf[facei] - phiEReferencePatchCentre);

    if (distSqr < phiEReferenceLocalDistSqr)
    {
        phiEReferenceLocalDistSqr = distSqr;
        phiEReferenceLocalFace = facei;
    }
}

scalar phiEReferenceGlobalDistSqr = phiEReferenceLocalDistSqr;
reduce(phiEReferenceGlobalDistSqr, minOp<scalar>());

label phiEReferenceProc = Pstream::nProcs();

if
(
    phiEReferenceLocalFace >= 0
 && mag(phiEReferenceLocalDistSqr - phiEReferenceGlobalDistSqr)
    <= SMALL*max(scalar(1), mag(phiEReferenceGlobalDistSqr))
)
{
    phiEReferenceProc = Pstream::myProcNo();
}

reduce(phiEReferenceProc, minOp<label>());

label phiERefCell = -1;

if (Pstream::myProcNo() == phiEReferenceProc)
{
    phiERefCell = phiEReferencePatch.faceCells()[phiEReferenceLocalFace];
}

const scalar phiERefValue = potentialDict.getOrDefault<scalar>("referencePotential", 0.0);

Info<< "electroneutralGalvanostatic potential reference: patch " << phiEReferencePatchName
    << ", referencePotential = " << phiERefValue << " V" << nl << endl;

volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeE ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    -fvc::grad(phiE)
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeKappa ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
        writeJDiff ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimCurrentDensity, vector::zero)
);

// surfaceScalarField jDiffFlux
// (
//     IOobject
//     (
//         "jDiffFlux",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar
//     (
//         "zero",
//         dimCurrentDensity*dimArea,
//         0.0
//     )
// );

// electroneutralGalvanostatic nonlinear corrector controls. These are deliberately separate from
// equilibriumProjection controls, which govern the local chemistry Newton
// iterations.
const dictionary& solutionDict = mesh.solutionDict();

if (!solutionDict.found("electroneutralGalvanostatic"))
{
    FatalIOErrorInFunction(solutionDict)
        << "Missing electroneutralGalvanostatic dictionary in system/fvSolution"
        << exit(FatalIOError);
}

const dictionary& enPhiDict = solutionDict.subDict("electroneutralGalvanostatic");

const label enPhiMaxIter = enPhiDict.get<label>("maxIter");
const label enPhiMinIter = enPhiDict.getOrDefault<label>("minIter", 1);

if (enPhiMaxIter <= 0 || enPhiMinIter <= 0 || enPhiMinIter > enPhiMaxIter)
{
    FatalIOErrorInFunction(enPhiDict)
        << "Require 0 < minIter <= maxIter for electroneutralGalvanostatic"
        << exit(FatalIOError);
}

if (!enPhiDict.found("concentration")
 || !enPhiDict.found("potential")
 || !enPhiDict.found("electroneutrality"))
{
    FatalIOErrorInFunction(enPhiDict)
        << "electroneutralGalvanostatic requires concentration, potential and electroneutrality "
        << "sub-dictionaries"
        << exit(FatalIOError);
}

const dictionary& enPhiCDict = enPhiDict.subDict("concentration");
const dictionary& enPhiPhiDict = enPhiDict.subDict("potential");
const dictionary& enPhiENDict = enPhiDict.subDict("electroneutrality");

const scalar enPhiCAbsTol = enPhiCDict.get<scalar>("absTol");
const scalar enPhiCRelTol = enPhiCDict.get<scalar>("relTol");
const scalar enPhiCScale = enPhiCDict.getOrDefault<scalar>("scale", SMALL);

const scalar enPhiPhiAbsTol = enPhiPhiDict.get<scalar>("absTol");
const scalar enPhiPhiRelTol = enPhiPhiDict.get<scalar>("relTol");
const scalar enPhiPhiScale = enPhiPhiDict.getOrDefault<scalar>("scale", 1.0);

const scalar enPhiENAbsTol = enPhiENDict.get<scalar>("absTol");
const scalar enPhiENRelTol = enPhiENDict.get<scalar>("relTol");
const scalar enPhiENScale = enPhiENDict.getOrDefault<scalar>("scale", SMALL);

if
(
    enPhiCAbsTol < 0.0 || enPhiCRelTol < 0.0 || enPhiCScale <= 0.0
 || enPhiPhiAbsTol < 0.0 || enPhiPhiRelTol < 0.0 || enPhiPhiScale <= 0.0
 || enPhiENAbsTol < 0.0 || enPhiENRelTol < 0.0 || enPhiENScale <= 0.0
)
{
    FatalIOErrorInFunction(enPhiDict)
        << "electroneutralGalvanostatic tolerances must be non-negative and scaling floors must "
        << "be positive"
        << exit(FatalIOError);
}

// Current-conservation monitoring values are updated by phiEqn.H.
scalar enPhiCurrentLocalError = GREAT;
scalar enPhiCurrentGlobalError = GREAT;

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

