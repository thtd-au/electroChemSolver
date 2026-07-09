#include "pHConcentration.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "Time.H"
#include "IOdictionary.H"

namespace Foam
{

defineTypeNameAndDebug(pHConcentration, 0);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    pHConcentration,
    dictionary
);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    pHConcentration,
    patchMapper
);

pHConcentration::
pHConcentration
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    pH_(7.0)
{}

pHConcentration::
pHConcentration
(
    const pHConcentration& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    pH_(ptf.pH_)
{}

pHConcentration::
pHConcentration
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    pH_(dict.get<scalar>("pH"))
{
    fvPatchField<scalar>::operator=(scalarField(p.size(), 0.0));
    updateCoeffs();
}

pHConcentration::
pHConcentration
(
    const pHConcentration& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    pH_(ptf.pH_)
{}

pHConcentration::
pHConcentration
(
    const pHConcentration& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    pH_(ptf.pH_)
{}

void pHConcentration::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    word specieName = internalField().name();
    const word prefix("c_");

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
            << "pHConcentration requires reactions/waterEquilibrium"
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
            << "pHConcentration can only be used for "
            << HplusName << " or " << OHminusName
            << ". It was used for " << specieName
            << exit(FatalIOError);
    }

    fvPatchField<scalar>::operator=(scalarField(patch().size(), ci));

    fixedValueFvPatchScalarField::updateCoeffs();
}

void pHConcentration::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("pH", pH_);
}

}

