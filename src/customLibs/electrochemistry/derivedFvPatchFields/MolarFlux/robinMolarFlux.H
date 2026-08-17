#ifndef robinMolarFlux_H
#define robinMolarFlux_H

#include "mixedFvPatchFields.H"

namespace Foam
{

class robinMolarFlux
:
    public mixedFvPatchScalarField
{
    // Stoichiometric coefficient in the electrode half reaction.
    scalar stoichCoeff_;

    // Number of electrons transferred in the electrode half reaction.
    scalar nElectrons_;

public:

    TypeName("robinMolarFlux");

    robinMolarFlux
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    robinMolarFlux
    (
        const robinMolarFlux&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    robinMolarFlux
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    robinMolarFlux
    (
        const robinMolarFlux&
    );

    robinMolarFlux
    (
        const robinMolarFlux&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new robinMolarFlux(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new robinMolarFlux(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

