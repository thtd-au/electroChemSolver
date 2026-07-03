#ifndef pHMassFractionFvPatchScalarField_H
#define pHMassFractionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class pHMassFractionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    scalar pH_;

public:

    TypeName("pHMassFraction");

    pHMassFractionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    pHMassFractionFvPatchScalarField
    (
        const pHMassFractionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    pHMassFractionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    pHMassFractionFvPatchScalarField
    (
        const pHMassFractionFvPatchScalarField&
    );

    pHMassFractionFvPatchScalarField
    (
        const pHMassFractionFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new pHMassFractionFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new pHMassFractionFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

}

#endif

