#ifndef galvanostaticPotentialFlux_H
#define galvanostaticPotentialFlux_H

#include "fixedGradientFvPatchFields.H"

namespace Foam
{

class galvanostaticPotentialFlux
:
    public fixedGradientFvPatchScalarField
{
    // Signed electron stoichiometric coefficient of the electrode
    // half-reaction. Its sign determines the normal applied-current
    // direction. The magnitude is used by galvanostaticMolarFlux.
    scalar electronStoichCoeff_;

public:

    TypeName("galvanostaticPotentialFlux");

    galvanostaticPotentialFlux
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    galvanostaticPotentialFlux
    (
        const galvanostaticPotentialFlux&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    galvanostaticPotentialFlux
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    galvanostaticPotentialFlux
    (
        const galvanostaticPotentialFlux&
    );

    galvanostaticPotentialFlux
    (
        const galvanostaticPotentialFlux&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new galvanostaticPotentialFlux(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new galvanostaticPotentialFlux(*this, iF)
        );
    }

    // Used by galvanostaticMolarFlux on the same electrode patch.
    scalar electronStoichCoeff() const
    {
        return electronStoichCoeff_;
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

