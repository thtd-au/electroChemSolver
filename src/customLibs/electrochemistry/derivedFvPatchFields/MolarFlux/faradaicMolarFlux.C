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


// ************************************************************************* //
// Constructors
// ************************************************************************* //

faradaicMolarFlux::
faradaicMolarFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    stoichCoeff_(0.0),
    nElectrons_(1.0)
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
    nElectrons_(ptf.nElectrons_)
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
    nElectrons_(dict.get<scalar>("nElectrons"))
{
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );

    // Gradient is calculated in updateCoeffs().
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
    nElectrons_(ptf.nElectrons_)
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
    nElectrons_(ptf.nElectrons_)
{}


// ************************************************************************* //
// Update boundary coefficients
// ************************************************************************* //

void faradaicMolarFlux::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh =
        patch().boundaryMesh().mesh();


    // --------------------------------------------------------------------- //
    // Determine species name from concentration field name
    //
    // Expected field names:
    //
    //     c_Hplus
    //     c_OHminus
    //     c_Acetate
    //
    // giving species names:
    //
    //     Hplus
    //     OHminus
    //     Acetate
    // --------------------------------------------------------------------- //

    word specieName = internalField().name();

    const word prefix("c_");

    if (specieName.startsWith(prefix))
    {
        specieName =
            specieName.substr(prefix.size());
    }


    // --------------------------------------------------------------------- //
    // Access solver dictionaries
    // --------------------------------------------------------------------- //

    const IOdictionary& reactions =
        mesh.lookupObject<IOdictionary>
        (
            "reactions"
        );

    const IOdictionary& electroChemistry =
        mesh.lookupObject<IOdictionary>
        (
            "electroChemistry"
        );


    // --------------------------------------------------------------------- //
    // Species diffusivity
    // --------------------------------------------------------------------- //

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
            << "No diffusivity entry for species "
            << specieName
            << exit(FatalIOError);
    }

    const scalar Di =
        diffusivityDict.get<scalar>(specieName);

    if (Di <= SMALL)
    {
        FatalIOErrorInFunction(reactions)
            << "Diffusivity for species "
            << specieName
            << " must be positive"
            << exit(FatalIOError);
    }


    // --------------------------------------------------------------------- //
    // Species charge number
    //
    // Species absent from chargeNumber are treated as neutral.
    // --------------------------------------------------------------------- //

    scalar zi = 0.0;

    if (reactions.found("chargeNumber"))
    {
        const dictionary& zDict =
            reactions.subDict("chargeNumber");

        zi =
            zDict.getOrDefault<scalar>
            (
                specieName,
                0.0
            );
    }


    // --------------------------------------------------------------------- //
    // Electrochemical properties
    //
    // Use the same source as the main electroChem solver:
    //
    //     constant/electroChemistry
    //
    // --------------------------------------------------------------------- //

    if (!electroChemistry.found("potential"))
    {
        FatalIOErrorInFunction(electroChemistry)
            << "Missing required potential dictionary"
            << exit(FatalIOError);
    }

    const dictionary& potentialDict =
        electroChemistry.subDict("potential");

    if (!potentialDict.found("appliedCurrentDensity"))
    {
        FatalIOErrorInFunction(potentialDict)
            << "Missing appliedCurrentDensity"
            << exit(FatalIOError);
    }

    const scalar currentDensity =
        potentialDict.get<scalar>
        (
            "appliedCurrentDensity"
        );

    if (!electroChemistry.found("temperature"))
    {
        FatalIOErrorInFunction(electroChemistry)
            << "Missing temperature"
            << exit(FatalIOError);
    }

    const scalar T =
        electroChemistry.get<scalar>
        (
            "temperature"
        );

    if (!electroChemistry.found("FaradayConstant"))
    {
        FatalIOErrorInFunction(electroChemistry)
            << "Missing FaradayConstant"
            << exit(FatalIOError);
    }

    const scalar F =
        electroChemistry.get<scalar>
        (
            "FaradayConstant"
        );

    const scalar R =
        8.31446261815324;


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


    // --------------------------------------------------------------------- //
    // Resolved electric potential gradient
    // --------------------------------------------------------------------- //

    if (!mesh.foundObject<volVectorField>("gradPhiE"))
    {
        FatalErrorInFunction
            << "Required field gradPhiE was not found in the "
            << "mesh object registry."
            << exit(FatalError);
    }

    const label patchi =
        patch().index();

    const volVectorField& gradPhiE =
        mesh.lookupObject<volVectorField>
        (
            "gradPhiE"
        );

    const fvPatchVectorField& gradPhiEPatch =
        gradPhiE.boundaryField()[patchi];


    // --------------------------------------------------------------------- //
    // Patch geometry and concentrations
    // --------------------------------------------------------------------- //

    const vectorField nf =
        patch().nf();

    const scalarField cPatch(*this);

    scalarField normalGradPhiE
    (
        patch().size(),
        0.0
    );

    forAll(normalGradPhiE, facei)
    {
        normalGradPhiE[facei] =
            gradPhiEPatch[facei]
          & nf[facei];
    }


    // --------------------------------------------------------------------- //
    // Faradaic + migration boundary gradient
    //
    // Species molar flux normal to the electrode:
    //
    //     J_i,n = nu_i j/(n F)
    //
    // Nernst-Planck:
    //
    //     J_i,n =
    //       -D_i dc_i/dn
    //       -D_i z_i F/(R T) c_i d(phi)/dn
    //
    // Rearranging for the OpenFOAM fixed-gradient condition:
    //
    //     dc_i/dn =
    //         nu_i j/(n F D_i)
    //       - z_i F/(R T) c_i d(phi)/dn
    //
    // The existing solver sign convention is retained here.
    // --------------------------------------------------------------------- //

    gradient() =
        stoichCoeff_
       *currentDensity
       /(nElectrons_*F*Di)

      - zi
       *F
       /(R*T)
       *cPatch
       *normalGradPhiE;


    fixedGradientFvPatchScalarField::updateCoeffs();
}


// ************************************************************************* //
// Write
// ************************************************************************* //

void faradaicMolarFlux::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    os.writeEntry
    (
        "stoichCoeff",
        stoichCoeff_
    );

    os.writeEntry
    (
        "nElectrons",
        nElectrons_
    );

    writeEntry
    (
        "value",
        os
    );
}


} // End namespace Foam