#ifndef inertMolarConcentration_H
#define inertMolarConcentration_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class inertMolarConcentration
:
    public fixedValueFvPatchScalarField
{
public:

    TypeName("inertMolarConcentration");

    inertMolarConcentration
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    inertMolarConcentration
    (
        const inertMolarConcentration&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    inertMolarConcentration
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    inertMolarConcentration
    (
        const inertMolarConcentration&
    );

    inertMolarConcentration
    (
        const inertMolarConcentration&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new inertMolarConcentration(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new inertMolarConcentration(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

