#include "inertMolarConcentration.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "dimensionedScalar.H"

namespace Foam
{

defineTypeNameAndDebug(inertMolarConcentration, 0);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    inertMolarConcentration,
    dictionary
);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    inertMolarConcentration,
    patchMapper
);


inertMolarConcentration::inertMolarConcentration
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


inertMolarConcentration::inertMolarConcentration
(
    const inertMolarConcentration& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


inertMolarConcentration::inertMolarConcentration
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{
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


inertMolarConcentration::inertMolarConcentration
(
    const inertMolarConcentration& ptf
)
:
    fixedValueFvPatchScalarField(ptf)
{}


inertMolarConcentration::inertMolarConcentration
(
    const inertMolarConcentration& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF)
{}


void inertMolarConcentration::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (!mesh.foundObject<IOdictionary>("reactions"))
    {
        FatalErrorInFunction
            << "inertMolarConcentration requires the registered "
            << "constant/reactions dictionary." << nl
            << "Patch: " << patch().name() << nl
            << "Field: " << internalField().name()
            << exit(FatalError);
    }

    const IOdictionary& reactions =
        mesh.lookupObject<IOdictionary>("reactions");

    if (!reactions.found("species") || !reactions.found("inertSpecie"))
    {
        FatalIOErrorInFunction(reactions)
            << "inertMolarConcentration requires species and inertSpecie "
            << "entries in constant/reactions."
            << exit(FatalIOError);
    }

    if (!reactions.found("molarMass"))
    {
        FatalIOErrorInFunction(reactions)
            << "inertMolarConcentration requires the molarMass dictionary "
            << "in constant/reactions."
            << exit(FatalIOError);
    }

    const wordList species(reactions.lookup("species"));
    const word inertSpecie(reactions.lookup("inertSpecie"));

    if (species.find(inertSpecie) < 0)
    {
        FatalIOErrorInFunction(reactions)
            << "Configured inertSpecie '" << inertSpecie
            << "' is not listed in species." << nl
            << "Available species are: " << species
            << exit(FatalIOError);
    }

    const word fieldName = internalField().name();
    const word expectedField("c_" + inertSpecie);

    if (fieldName != inertSpecie && fieldName != expectedField)
    {
        FatalErrorInFunction
            << "inertMolarConcentration may only be applied to the "
            << "configured inert species." << nl
            << "Patch field: " << fieldName << nl
            << "Configured inert species: " << inertSpecie << nl
            << "Expected concentration field for the current solver: "
            << expectedField
            << exit(FatalError);
    }

    const dictionary& molarMassDict = reactions.subDict("molarMass");

    if (!molarMassDict.found(inertSpecie))
    {
        FatalIOErrorInFunction(molarMassDict)
            << "Missing molarMass entry for inert species "
            << inertSpecie
            << exit(FatalIOError);
    }

    const scalar MInert = molarMassDict.get<scalar>(inertSpecie);

    if (MInert <= SMALL)
    {
        FatalIOErrorInFunction(molarMassDict)
            << "Molar mass for inert species " << inertSpecie
            << " must be positive."
            << exit(FatalIOError);
    }

    // rho is currently supplied through constant/transportProperties in the
    // same way as the solver's singlePhaseTransportModel. Read an unregistered
    // dictionary here so the patch field remains independent of solver-local
    // variables.
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const dimensionedScalar rho
    (
        "rho",
        dimDensity,
        transportProperties
    );

    if (rho.value() <= SMALL)
    {
        FatalIOErrorInFunction(transportProperties)
            << "rho must be positive for inertMolarConcentration."
            << exit(FatalIOError);
    }

    const label patchi = patch().index();
    scalarField nonInertMassFraction(patch().size(), 0.0);

    forAll(species, i)
    {
        const word& specieName = species[i];

        if (specieName == inertSpecie)
        {
            continue;
        }

        if (!molarMassDict.found(specieName))
        {
            FatalIOErrorInFunction(molarMassDict)
                << "Missing molarMass entry for species " << specieName
                << exit(FatalIOError);
        }

        const scalar Mi = molarMassDict.get<scalar>(specieName);

        if (Mi <= SMALL)
        {
            FatalIOErrorInFunction(molarMassDict)
                << "Molar mass for species " << specieName
                << " must be positive."
                << exit(FatalIOError);
        }

        const word concentrationFieldName("c_" + specieName);

        if (!mesh.foundObject<volScalarField>(concentrationFieldName))
        {
            FatalErrorInFunction
                << "Could not find concentration field '"
                << concentrationFieldName << "' for species "
                << specieName << "." << nl
                << "inertMolarConcentration uses the same c_<species> "
                << "mapping as the current electroChem solver."
                << exit(FatalError);
        }

        volScalarField& concentration =
            mesh.lookupObjectRef<volScalarField>(concentrationFieldName);

        fvPatchScalarField& concentrationPatch =
            concentration.boundaryFieldRef()[patchi];

        // Make calculated/custom inlet values current before evaluating the
        // inert remainder. The inert field itself is skipped above, preventing
        // recursion through this boundary condition.
        if (!concentrationPatch.updated())
        {
            concentrationPatch.updateCoeffs();
        }

        const scalarField cPatch(concentrationPatch);

        nonInertMassFraction +=
            cPatch*Mi/rho.value();
    }

    scalarField inertMassFraction
    (
        scalar(1.0) - nonInertMassFraction
    );

    label worstFace = -1;
    scalar minInertMassFraction = GREAT;

    forAll(inertMassFraction, facei)
    {
        if (inertMassFraction[facei] < minInertMassFraction)
        {
            minInertMassFraction = inertMassFraction[facei];
            worstFace = facei;
        }
    }

    if (minInertMassFraction < 0.0)
    {
        FatalErrorInFunction
            << "inertMolarConcentration produced a negative inert mass "
            << "fraction." << nl
            << "Patch: " << patch().name() << nl
            << "Inert species: " << inertSpecie << nl
            << "Face: " << worstFace << nl
            << "Minimum inert mass fraction: "
            << minInertMassFraction << nl
            << "Non-inert mass-fraction sum at that face: "
            << nonInertMassFraction[worstFace] << nl
            << "The prescribed non-inert concentrations correspond to a "
            << "total mass fraction greater than one."
            << exit(FatalError);
    }

    const scalarField inertConcentration
    (
        rho.value()*inertMassFraction/MInert
    );

    operator==(inertConcentration);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void inertMolarConcentration::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}

} // End namespace Foam

