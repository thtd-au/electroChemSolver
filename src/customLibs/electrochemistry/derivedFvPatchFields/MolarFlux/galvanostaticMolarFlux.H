#ifndef galvanostaticMolarFlux_H
#define galvanostaticMolarFlux_H

#include "mixedFvPatchFields.H"

namespace Foam
{

class galvanostaticMolarFlux
:
    public mixedFvPatchScalarField
{
    // Species name inferred from the concentration field name c_<species>.
    word speciesName_;

    // Stoichiometric coefficient of this aqueous species in the electrode
    // half-reaction. Positive means production, negative consumption.
    scalar speciesStoichCoeff_;

public:

    TypeName("galvanostaticMolarFlux");

    galvanostaticMolarFlux
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    galvanostaticMolarFlux
    (
        const galvanostaticMolarFlux&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    galvanostaticMolarFlux
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    galvanostaticMolarFlux
    (
        const galvanostaticMolarFlux&
    );

    galvanostaticMolarFlux
    (
        const galvanostaticMolarFlux&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new galvanostaticMolarFlux(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new galvanostaticMolarFlux(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

