#ifndef inertMassFractionFvPatchScalarField_H
#define inertMassFractionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class inertMassFractionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
public:

    TypeName("inertMassFraction");

    inertMassFractionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    inertMassFractionFvPatchScalarField
    (
        const inertMassFractionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    inertMassFractionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    inertMassFractionFvPatchScalarField
    (
        const inertMassFractionFvPatchScalarField&
    );

    inertMassFractionFvPatchScalarField
    (
        const inertMassFractionFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new inertMassFractionFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new inertMassFractionFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

}

#endif

