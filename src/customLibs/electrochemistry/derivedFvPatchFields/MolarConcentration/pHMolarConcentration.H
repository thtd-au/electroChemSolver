#ifndef pHMolarConcentration_H
#define pHMolarConcentration_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class pHMolarConcentration
:
    public fixedValueFvPatchScalarField
{
    // Name of the equilibrium reaction block in constant/reactions.
    word equilibriumReaction_;

    // Prescribed pH.
    scalar pH_;

public:

    TypeName("pHMolarConcentration");

    pHMolarConcentration
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    pHMolarConcentration
    (
        const pHMolarConcentration&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    pHMolarConcentration
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    pHMolarConcentration
    (
        const pHMolarConcentration&
    );

    pHMolarConcentration
    (
        const pHMolarConcentration&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new pHMolarConcentration(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new pHMolarConcentration(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

