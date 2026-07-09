#include "faradaicMolarFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"
#include "IOdictionary.H"

namespace Foam
{

defineTypeNameAndDebug(faradaicMolarFlux, 0);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    faradaicMolarFlux,
    dictionary
);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    faradaicMolarFlux,
    patchMapper
);

faradaicMolarFlux::
faradaicMolarFlux
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

faradaicMolarFlux::
faradaicMolarFlux
(
    const faradaicMolarFlux& ptf,
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

faradaicMolarFlux::
faradaicMolarFlux
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
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );

    gradient() = 0.0;

    if (nElectrons_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "nElectrons must be positive"
            << exit(FatalIOError);
    }
}

faradaicMolarFlux::
faradaicMolarFlux
(
    const faradaicMolarFlux& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    stoichCoeff_(ptf.stoichCoeff_),
    nElectrons_(ptf.nElectrons_),
    currentDensity_(ptf.currentDensity_)
{}

faradaicMolarFlux::
faradaicMolarFlux
(
    const faradaicMolarFlux& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    stoichCoeff_(ptf.stoichCoeff_),
    nElectrons_(ptf.nElectrons_),
    currentDensity_(ptf.currentDensity_)
{}

void faradaicMolarFlux::updateCoeffs()
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

    const dictionary& diffusivityDict = reactions.subDict("diffusivity");

    if (!diffusivityDict.found(specieName))
    {
        FatalIOErrorInFunction(reactions)
            << "No diffusivity entry for species " << specieName
            << exit(FatalIOError);
    }

    const scalar Di = diffusivityDict.get<scalar>(specieName);

    if (Di <= SMALL)
    {
        FatalIOErrorInFunction(reactions)
            << "Diffusivity for species " << specieName
            << " must be positive"
            << exit(FatalIOError);
    }

    scalar zi = 0.0;

    if (reactions.found("chargeNumber"))
    {
        const dictionary& zDict = reactions.subDict("chargeNumber");
        zi = zDict.getOrDefault<scalar>(specieName, 0.0);
    }

    const scalar T = reactions.getOrDefault<scalar>("temperature", 298.15);
    const scalar F = reactions.getOrDefault<scalar>("FaradayConstant", 96485.3329);
    const scalar R = 8.31446261815324;

    if (T <= SMALL)
    {
        FatalIOErrorInFunction(reactions)
            << "temperature must be positive"
            << exit(FatalIOError);
    }

    const label patchi = patch().index();

    const volScalarField& phiE =
        mesh.lookupObject<volScalarField>("phiE");

    const scalarField cPatch(*this);

    tmp<scalarField> tsnGradPhiE =
        phiE.boundaryField()[patchi].snGrad();

    const scalarField& snGradPhiE = tsnGradPhiE();

    gradient() =
        stoichCoeff_*currentDensity_/(nElectrons_*F*Di)
      - zi*F/(R*T)*cPatch*snGradPhiE;

    fixedGradientFvPatchScalarField::updateCoeffs();
}

void faradaicMolarFlux::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("stoichCoeff", stoichCoeff_);
    os.writeEntry("nElectrons", nElectrons_);
    os.writeEntry("currentDensity", currentDensity_);
    writeEntry("value", os);
}

}

