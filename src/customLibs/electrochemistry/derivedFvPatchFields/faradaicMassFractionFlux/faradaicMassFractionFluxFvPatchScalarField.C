#include "faradaicMassFractionFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "Time.H"
#include "IOdictionary.H"

namespace Foam
{

defineTypeNameAndDebug(faradaicMassFractionFluxFvPatchScalarField, 0);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    faradaicMassFractionFluxFvPatchScalarField,
    dictionary
);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    faradaicMassFractionFluxFvPatchScalarField,
    patchMapper
);


faradaicMassFractionFluxFvPatchScalarField::
faradaicMassFractionFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    stoichCoeff_(0),
    nElectrons_(1),
    currentDensity_(0)
{}


faradaicMassFractionFluxFvPatchScalarField::
faradaicMassFractionFluxFvPatchScalarField
(
    const faradaicMassFractionFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    stoichCoeff_(ptf.stoichCoeff_),
    nElectrons_(ptf.nElectrons_),
    currentDensity_(ptf.currentDensity_)
{}


faradaicMassFractionFluxFvPatchScalarField::
faradaicMassFractionFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    stoichCoeff_(dict.get<scalar>("stoichCoeff")),
    nElectrons_(dict.get<scalar>("nElectrons")),
    currentDensity_(dict.get<scalar>("currentDensity"))
{
    // The base fixedGradient BC normally reads a user-supplied gradient.
    // Here the gradient is calculated later in updateCoeffs().  
    // So dummy gradient is created initialy and will be updated to actual gradient 
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );

    gradient() = 0.0;

    if (nElectrons_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "nElectrons must be positive" << exit(FatalIOError);
    }
}


faradaicMassFractionFluxFvPatchScalarField::
faradaicMassFractionFluxFvPatchScalarField
(
    const faradaicMassFractionFluxFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    stoichCoeff_(ptf.stoichCoeff_),
    nElectrons_(ptf.nElectrons_),
    currentDensity_(ptf.currentDensity_)
{}


faradaicMassFractionFluxFvPatchScalarField::
faradaicMassFractionFluxFvPatchScalarField
(
    const faradaicMassFractionFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    stoichCoeff_(ptf.stoichCoeff_),
    nElectrons_(ptf.nElectrons_),
    currentDensity_(ptf.currentDensity_)
{}


void faradaicMassFractionFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // create species name by removing "Y_" from sepceis field 
    const word fieldName = internalField().name();
    word specieName = fieldName;
    const word prefix("Y_");
    if (specieName.startsWith(prefix))
    {
        specieName = specieName.substr(prefix.size());
    }

    IOdictionary reactions
    (
        IOobject
        (
            "reactions",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    IOdictionary transportProperties    
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar F = reactions.getOrDefault<scalar>("FaradayConstant", 96485.3329);

    const dictionary& diffusivityDict = reactions.subDict("diffusivity");

    const dictionary& molarMassDict = reactions.subDict("molarMass");
    
    const scalar Di = diffusivityDict.get<scalar>(specieName);

    const scalar Mi = molarMassDict.get<scalar>(specieName);

    const dimensionedScalar rho
    (
        "rho",
        dimDensity,
        transportProperties
    );

    if (!diffusivityDict.found(specieName))
    {
        FatalIOErrorInFunction(reactions)
            << "No diffusivity entry for species " << specieName
            << " in constant/reactions/diffusivity" << exit(FatalIOError);
    }

    if (!molarMassDict.found(specieName))
    {
        FatalIOErrorInFunction(reactions)
            << "No molarMass entry for species " << specieName
            << " in constant/reactions/molarMass" << exit(FatalIOError);
    }

    if (Di <= SMALL)
    {
        FatalIOErrorInFunction(reactions)
            << "Diffusivity for species " << specieName
            << " must be positive" << exit(FatalIOError);
    }

    if (Mi <= SMALL)
    {
        FatalIOErrorInFunction(reactions)
            << "molarMass for species " << specieName
            << " must be positive" << exit(FatalIOError);
    }

    if (rho.value() <= SMALL)
    {
        FatalIOErrorInFunction(transportProperties)
            << "rho must be positive" << exit(FatalIOError);
    }

    dimensionSet dimVoltage(1, 2, -3, 0, 0, -1, 0);

    dimensionedVector gradPhiE
    (
        "gradPhiE",
        dimVoltage/dimLength,
        reactions.lookupOrDefault<vector>("gradPhiE", vector::zero)
    );

    tmp<vectorField> tnf = patch().nf();
    const vectorField& nf = tnf();

    scalarField snGradPhiE(patch().size(), 0.0);

    forAll(snGradPhiE, facei)
    {
        snGradPhiE[facei] = nf[facei] & gradPhiE.value();
    }

    scalar zi = 0.0;
    if (reactions.found("chargeNumber"))
    {
        const dictionary& zDict = reactions.subDict("chargeNumber");
        zi = zDict.getOrDefault<scalar>(specieName, 0.0);
    }

    const scalar T =
        reactions.getOrDefault<scalar>("temperature", 298.15);

    const scalar R = 8.31446261815324;

    const scalar faradaicGradient =
        stoichCoeff_*Mi*currentDensity_
    /(rho.value()*nElectrons_*F*Di);

    const scalarField YiInternal = patchInternalField();

    const scalarField beta =
        zi*F/(R*T)*snGradPhiE;

    gradient() =
        (faradaicGradient - beta*YiInternal)
    /(scalar(1.0) + beta/patch().deltaCoeffs());
    
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void faradaicMassFractionFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("stoichCoeff", stoichCoeff_);
    os.writeEntry("nElectrons", nElectrons_);
    os.writeEntry("currentDensity", currentDensity_);
    writeEntry("value", os);
}

}

