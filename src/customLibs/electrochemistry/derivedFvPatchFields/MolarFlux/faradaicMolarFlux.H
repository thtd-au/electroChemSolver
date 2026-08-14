#ifndef faradaicMolarFlux_H
#define faradaicMolarFlux_H

#include "fixedGradientFvPatchFields.H"

namespace Foam
{

class faradaicMolarFlux
:
    public fixedGradientFvPatchScalarField
{
    scalar stoichCoeff_;
    scalar nElectrons_;

public:

    TypeName("faradaicMolarFlux");

    faradaicMolarFlux
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    faradaicMolarFlux
    (
        const faradaicMolarFlux&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    faradaicMolarFlux
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    faradaicMolarFlux
    (
        const faradaicMolarFlux&
    );

    faradaicMolarFlux
    (
        const faradaicMolarFlux&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new faradaicMolarFlux(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new faradaicMolarFlux(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

}

#endif

