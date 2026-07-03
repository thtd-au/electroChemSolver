#include "inertMassFractionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "Time.H"
#include "IOdictionary.H"

namespace Foam
{

defineTypeNameAndDebug(inertMassFractionFvPatchScalarField, 0);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    inertMassFractionFvPatchScalarField,
    dictionary
);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    inertMassFractionFvPatchScalarField,
    patchMapper
);


inertMassFractionFvPatchScalarField::
inertMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


inertMassFractionFvPatchScalarField::
inertMassFractionFvPatchScalarField
(
    const inertMassFractionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


inertMassFractionFvPatchScalarField::
inertMassFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    // No user-supplied value is required.
    // Initialise to zero; updateCoeffs() computes 1 - sum(other Y).
    fvPatchField<scalar>::operator=(scalarField(p.size(), 0.0));
}


inertMassFractionFvPatchScalarField::
inertMassFractionFvPatchScalarField
(
    const inertMassFractionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf)
{}


inertMassFractionFvPatchScalarField::
inertMassFractionFvPatchScalarField
(
    const inertMassFractionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF)
{}


void inertMassFractionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const label patchi = patch().index();

    const word thisFieldName = internalField().name();

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

    wordList species(reactions.lookup("species"));

    scalarField Ypatch(patch().size(), 1.0);

    forAll(species, i)
    {
        const word fieldName = "Y_" + species[i];

        if (fieldName == thisFieldName)
        {
            continue;
        }

        if (!mesh.foundObject<volScalarField>(fieldName))
        {
            FatalIOErrorInFunction(reactions)
                << "Cannot find field " << fieldName
                << " needed by inertMassFraction on patch "
                << patch().name()
                << exit(FatalIOError);
        }

        const volScalarField& Yj =
            mesh.lookupObject<volScalarField>(fieldName);

        Ypatch -= Yj.boundaryField()[patchi];
    }

    if (min(Ypatch) < -SMALL)
    {
        FatalErrorInFunction
            << "Mass fractions exceed unity on patch "
            << patch().name()
            << " while applying inertMassFraction to "
            << thisFieldName
            << ". Minimum closure value = "
            << min(Ypatch)
            << exit(FatalError);
    }

    Ypatch = max(Ypatch, scalar(0.0));

    fvPatchField<scalar>::operator=(Ypatch);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void inertMassFractionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}

}

