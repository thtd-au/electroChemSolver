#include "galvanostaticMolarFlux.H"
#include "galvanostaticPotentialFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"

namespace Foam
{

namespace
{
word concentrationSpeciesName(const word& fieldName)
{
    if
    (
        fieldName.size() <= 2
     || fieldName[0] != 'c'
     || fieldName[1] != '_'
    )
    {
        FatalErrorInFunction
            << "galvanostaticMolarFlux must be applied to a concentration "
            << "field named c_<species>. Field is " << fieldName
            << exit(FatalError);
    }

    return word(fieldName.substr(2));
}
}


defineTypeNameAndDebug(galvanostaticMolarFlux, 0);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    galvanostaticMolarFlux,
    dictionary
);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    galvanostaticMolarFlux,
    patchMapper
);


galvanostaticMolarFlux::galvanostaticMolarFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    speciesName_(concentrationSpeciesName(iF.name())),
    speciesStoichCoeff_(0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


galvanostaticMolarFlux::galvanostaticMolarFlux
(
    const galvanostaticMolarFlux& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    // Reconstruction may supply an unnamed/null destination internal field.
    // Mapping must preserve the source patch metadata without inspecting it.
    speciesName_(ptf.speciesName_),
    speciesStoichCoeff_(ptf.speciesStoichCoeff_)
{}


galvanostaticMolarFlux::galvanostaticMolarFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    speciesName_(concentrationSpeciesName(iF.name())),
    speciesStoichCoeff_(dict.get<scalar>("speciesStoichCoeff"))
{

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
        fvPatchField<scalar>::operator=(patchInternalField());
    }
}


galvanostaticMolarFlux::galvanostaticMolarFlux
(
    const galvanostaticMolarFlux& ptf
)
:
    mixedFvPatchScalarField(ptf),
    speciesName_(ptf.speciesName_),
    speciesStoichCoeff_(ptf.speciesStoichCoeff_)
{}


galvanostaticMolarFlux::galvanostaticMolarFlux
(
    const galvanostaticMolarFlux& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    speciesName_(ptf.speciesName_),
    speciesStoichCoeff_(ptf.speciesStoichCoeff_)
{}


void galvanostaticMolarFlux::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Enforce the concentration-field naming rule when the BC is actually
    // evaluated, after mapping/cloning has attached it to its real field.
    speciesName_ = concentrationSpeciesName(internalField().name());

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const IOdictionary& reactions =
        mesh.lookupObject<IOdictionary>("reactions");

    const IOdictionary& electroChemistry =
        mesh.lookupObject<IOdictionary>("electroChemistry");

    const wordList species(reactions.lookup("species"));

    if (species.find(speciesName_) < 0)
    {
        FatalIOErrorInFunction(reactions)
            << "Invalid concentration field species '" << speciesName_
            << "' for galvanostaticMolarFlux." << nl
            << "Available species are: " << species
            << exit(FatalIOError);
    }

    const dictionary& diffusivityDict = reactions.subDict("diffusivity");

    if (!diffusivityDict.found(speciesName_))
    {
        FatalIOErrorInFunction(reactions)
            << "No diffusivity entry for species " << speciesName_
            << exit(FatalIOError);
    }

    const scalar Di = diffusivityDict.get<scalar>(speciesName_);

    if (Di <= SMALL)
    {
        FatalIOErrorInFunction(reactions)
            << "Diffusivity for species " << speciesName_
            << " must be positive"
            << exit(FatalIOError);
    }

    scalar zi = 0.0;

    if (reactions.found("chargeNumber"))
    {
        zi = reactions.subDict("chargeNumber")
            .getOrDefault<scalar>(speciesName_, 0.0);
    }

    const dictionary& potentialDict =
        electroChemistry.subDict("potential");

    const bool potentialActive = potentialDict.get<word>("model") != "none";

    const scalar T = electroChemistry.get<scalar>("temperature");
    const scalar F = electroChemistry.get<scalar>("FaradayConstant");
    const scalar R = 8.31446261815324;

    if (T <= SMALL || F <= SMALL)
    {
        FatalIOErrorInFunction(electroChemistry)
            << "temperature and FaradayConstant must be positive"
            << exit(FatalIOError);
    }

    // Electrode stoichiometry always comes from the matching phiE patch.
    // In isolation mode only the metadata is consumed, never its gradient.
    const volScalarField& phiE = mesh.lookupObject<volScalarField>("phiE");
    const fvPatchScalarField& phiEPatch = phiE.boundaryField()[patch().index()];
    if (!isA<galvanostaticPotentialFlux>(phiEPatch))
    {
        FatalErrorInFunction
            << "galvanostaticMolarFlux requires galvanostaticPotentialFlux "
            << "on the matching phiE patch " << patch().name()
            << exit(FatalError);
    }
    const galvanostaticPotentialFlux& potentialFluxBC =
        refCast<const galvanostaticPotentialFlux>(phiEPatch);
    const scalar electronStoichCoeff = potentialFluxBC.electronStoichCoeff();
    scalarField snGradPhiE(patch().size(), 0.0);
    if (potentialActive) snGradPhiE = potentialFluxBC.gradient();
    if (mag(electronStoichCoeff) <= SMALL)
    {
        FatalErrorInFunction << "Zero electrode electronStoichCoeff on patch "
            << patch().name() << exit(FatalError);
    }
    const scalarField& delta = patch().deltaCoeffs();

    scalarField beta(patch().size(), 0.0);

    forAll(beta, facei)
    {
        beta[facei] = zi*F/(R*T)*snGradPhiE[facei];
    }

    if (!mesh.foundObject<surfaceScalarField>("electrodeCurrentDensity"))
    {
        FatalErrorInFunction
            << "galvanostaticMolarFlux requires the registered "
            << "electrodeCurrentDensity field supplied by currentManager."
            << exit(FatalError);
    }
    const surfaceScalarField& electrodeCurrentDensity =
        mesh.lookupObject<surfaceScalarField>("electrodeCurrentDensity");
    const fvsPatchScalarField& jnApplied =
        electrodeCurrentDensity.boundaryField()[patch().index()];

    // N_i dot n = nu_i*j_n/(nu_e*F).
    // dc_i/dn = A - beta*c_b, A = -nu_i*j_n/(nu_e*F*D_i).
    // This remains valid for signed, spatially varying model currents.
    scalarField A(patch().size(), 0.0);
    forAll(A, facei)
    {
        A[facei] = -speciesStoichCoeff_*jnApplied[facei]
            /(electronStoichCoeff*F*Di);
    }

    scalarField denominator(delta + beta);

    forAll(denominator, facei)
    {
        const scalar scale = max(mag(delta[facei]), scalar(1.0));

        if (mag(denominator[facei]) <= SMALL*scale)
        {
            FatalErrorInFunction
                << "galvanostaticMolarFlux is singular or nearly singular."
                << nl
                << "species = " << speciesName_ << nl
                << "patch = " << patch().name() << nl
                << "face = " << facei << nl
                << "deltaCoeffs = " << delta[facei] << " 1/m" << nl
                << "beta = " << beta[facei] << " 1/m" << nl
                << "deltaCoeffs + beta = " << denominator[facei]
                << " 1/m"
                << exit(FatalError);
        }
    }

    refValue() = 0.0;
    refGrad() = A;
    valueFraction() = beta/denominator;

    mixedFvPatchScalarField::updateCoeffs();
}


void galvanostaticMolarFlux::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("speciesStoichCoeff", speciesStoichCoeff_);
    writeEntry("value", os);
}

} // End namespace Foam

