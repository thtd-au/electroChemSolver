#include "pHMassFractionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "Time.H"
#include "IOdictionary.H"

namespace Foam
{

defineTypeNameAndDebug(pHMassFractionFvPatchScalarField, 0);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    pHMassFractionFvPatchScalarField,
    dictionary
);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    pHMassFractionFvPatchScalarField,
    patchMapper
);


pHMassFractionFvPatchScalarField::
pHMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    pH_(7.0)
{}


pHMassFractionFvPatchScalarField::
pHMassFractionFvPatchScalarField
(
    const pHMassFractionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    pH_(ptf.pH_)
{}


pHMassFractionFvPatchScalarField::
pHMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    pH_(dict.get<scalar>("pH"))
{
    // No user-supplied value entry is required.
    // Initialise with zero; updateCoeffs() overwrites it.
    fvPatchField<scalar>::operator=(scalarField(p.size(), 0.0));

    updateCoeffs();
}


pHMassFractionFvPatchScalarField::
pHMassFractionFvPatchScalarField
(
    const pHMassFractionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    pH_(ptf.pH_)
{}


pHMassFractionFvPatchScalarField::
pHMassFractionFvPatchScalarField
(
    const pHMassFractionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    pH_(ptf.pH_)
{}


void pHMassFractionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

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

    if (!reactions.found("reactions"))
    {
        FatalIOErrorInFunction(reactions)
            << "Missing reactions sub-dictionary in constant/reactions"
            << exit(FatalIOError);
    }

    const dictionary& rxnDicts = reactions.subDict("reactions");

    if (!rxnDicts.found("waterEquilibrium"))
    {
        FatalIOErrorInFunction(reactions)
            << "pHMassFraction requires reactions/waterEquilibrium"
            << exit(FatalIOError);
    }

    const dictionary& waterDict = rxnDicts.subDict("waterEquilibrium");

    const word HplusName(waterDict.lookup("Hplus"));
    const word OHminusName(waterDict.lookup("OHminus"));

    const scalar Kw = waterDict.get<scalar>("Kw");

    if (Kw <= SMALL)
    {
        FatalIOErrorInFunction(waterDict)
            << "waterEquilibrium/Kw must be positive"
            << exit(FatalIOError);
    }

    const dictionary& molarMassDict = reactions.subDict("molarMass");

    if (!molarMassDict.found(specieName))
    {
        FatalIOErrorInFunction(reactions)
            << "No molarMass entry for species " << specieName
            << exit(FatalIOError);
    }

    const scalar Mi = molarMassDict.get<scalar>(specieName);

    const dimensionedScalar rho
    (
        "rho",
        dimDensity,
        transportProperties
    );

    if (rho.value() <= SMALL)
    {
        FatalIOErrorInFunction(transportProperties)
            << "rho must be positive"
            << exit(FatalIOError);
    }

    // pH uses mol/L, solver concentration uses mol/m3.
    const scalar cH = 1000.0*pow(10.0, -pH_);
    const scalar cOH = Kw/cH;

    scalar ci = 0.0;

    if (specieName == HplusName)
    {
        ci = cH;
    }
    else if (specieName == OHminusName)
    {
        ci = cOH;
    }
    else
    {
        FatalIOErrorInFunction(waterDict)
            << "pHMassFraction can only be used for "
            << HplusName << " or " << OHminusName
            << ". It was used for " << specieName
            << exit(FatalIOError);
    }

    const scalar Yi = Mi*ci/rho.value();

    fvPatchField<scalar>::operator=(scalarField(patch().size(), Yi));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void pHMassFractionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("pH", pH_);
    writeEntry("value", os);
}

}

