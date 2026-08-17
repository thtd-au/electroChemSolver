#include "robinMolarFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"

namespace Foam
{

defineTypeNameAndDebug(robinMolarFlux, 0);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    robinMolarFlux,
    dictionary
);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    robinMolarFlux,
    patchMapper
);


robinMolarFlux::robinMolarFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    stoichCoeff_(0.0),
    nElectrons_(1.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


robinMolarFlux::robinMolarFlux
(
    const robinMolarFlux& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    stoichCoeff_(ptf.stoichCoeff_),
    nElectrons_(ptf.nElectrons_)
{}


robinMolarFlux::robinMolarFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    stoichCoeff_(dict.get<scalar>("stoichCoeff")),
    nElectrons_(dict.get<scalar>("nElectrons"))
{
    if (nElectrons_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "nElectrons must be positive"
            << exit(FatalIOError);
    }

    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=
        (
            patchInternalField()
        );
    }
}


robinMolarFlux::robinMolarFlux
(
    const robinMolarFlux& ptf
)
:
    mixedFvPatchScalarField(ptf),
    stoichCoeff_(ptf.stoichCoeff_),
    nElectrons_(ptf.nElectrons_)
{}


robinMolarFlux::robinMolarFlux
(
    const robinMolarFlux& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    stoichCoeff_(ptf.stoichCoeff_),
    nElectrons_(ptf.nElectrons_)
{}


void robinMolarFlux::updateCoeffs()
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

    const IOdictionary& reactions =
        mesh.lookupObject<IOdictionary>("reactions");

    const IOdictionary& electroChemistry =
        mesh.lookupObject<IOdictionary>("electroChemistry");

    if (!reactions.found("diffusivity"))
    {
        FatalIOErrorInFunction(reactions)
            << "Missing diffusivity dictionary"
            << exit(FatalIOError);
    }

    const dictionary& diffusivityDict =
        reactions.subDict("diffusivity");

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
        const dictionary& zDict =
            reactions.subDict("chargeNumber");

        zi = zDict.getOrDefault<scalar>(specieName, 0.0);
    }

    if (!electroChemistry.found("potential"))
    {
        FatalIOErrorInFunction(electroChemistry)
            << "Missing required potential dictionary"
            << exit(FatalIOError);
    }

    const dictionary& potentialDict =
        electroChemistry.subDict("potential");

    const scalar currentDensity =
        potentialDict.get<scalar>("appliedCurrentDensity");

    const scalar T =
        electroChemistry.get<scalar>("temperature");

    const scalar F =
        electroChemistry.get<scalar>("FaradayConstant");

    const scalar R = 8.31446261815324;

    if (T <= SMALL)
    {
        FatalIOErrorInFunction(electroChemistry)
            << "temperature must be positive"
            << exit(FatalIOError);
    }

    if (F <= SMALL)
    {
        FatalIOErrorInFunction(electroChemistry)
            << "FaradayConstant must be positive"
            << exit(FatalIOError);
    }

    if (!mesh.foundObject<volVectorField>("gradPhiE"))
    {
        FatalErrorInFunction
            << "Required field gradPhiE was not found in the mesh "
            << "object registry"
            << exit(FatalError);
    }

    const label patchi = patch().index();

    const volVectorField& gradPhiE =
        mesh.lookupObject<volVectorField>("gradPhiE");

    const fvPatchVectorField& gradPhiEPatch =
        gradPhiE.boundaryField()[patchi];

    const vectorField nf(patch().nf());
    const scalarField& delta = patch().deltaCoeffs();

    scalarField beta(patch().size(), 0.0);

    forAll(beta, facei)
    {
        const scalar normalGradPhiE =
            gradPhiEPatch[facei] & nf[facei];

        beta[facei] =
            zi*F/(R*T)*normalGradPhiE;
    }

    // Faradaic part of the concentration gradient
    //
    //     dc/dn = A - beta*c_b
    //
    // with
    //
    //     A = nu*j/(n*F*D)
    //
    // and OpenFOAM's boundary discretisation
    //
    //     dc/dn = delta*(c_b - c_P).
    //
    // Combining them gives the Robin relation
    //
    //     c_b = (delta*c_P + A)/(delta + beta).
    //
    // mixedFvPatchScalarField represents
    //
    //     c_b = w*refValue
    //         + (1-w)*(c_P + refGrad/delta),
    //
    // so choosing
    //
    //     refValue     = 0
    //     refGrad      = A
    //     valueFraction= beta/(delta + beta)
    //
    // reproduces the Robin relation exactly.
    const scalar A =
        stoichCoeff_*currentDensity
       /(nElectrons_*F*Di);

    scalarField denominator(delta + beta);

    forAll(denominator, facei)
    {
        const scalar scale = max(mag(delta[facei]), scalar(1.0));

        if (mag(denominator[facei]) <= SMALL*scale)
        {
            FatalErrorInFunction
                << "robinMolarFlux is singular or nearly singular." << nl
                << "species = " << specieName << nl
                << "patch = " << patch().name() << nl
                << "face = " << facei << nl
                << "deltaCoeffs = " << delta[facei] << " 1/m" << nl
                << "beta = " << beta[facei] << " 1/m" << nl
                << "deltaCoeffs + beta = "
                << denominator[facei] << " 1/m"
                << exit(FatalError);
        }
    }

    refValue() = 0.0;
    refGrad() = A;
    valueFraction() = beta/denominator;

    mixedFvPatchScalarField::updateCoeffs();
}


void robinMolarFlux::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("stoichCoeff", stoichCoeff_);
    os.writeEntry("nElectrons", nElectrons_);
    writeEntry("value", os);
}

} // End namespace Foam

