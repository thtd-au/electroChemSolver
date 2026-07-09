#ifndef pHConcentration_H
#define pHConcentration_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class pHConcentration
:
    public fixedValueFvPatchScalarField
{
    scalar pH_;

public:

    TypeName("pHConcentration");

    pHConcentration
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    pHConcentration
    (
        const pHConcentration&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    pHConcentration
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    pHConcentration
    (
        const pHConcentration&
    );

    pHConcentration
    (
        const pHConcentration&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new pHConcentration(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new pHConcentration(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

}

#endif

