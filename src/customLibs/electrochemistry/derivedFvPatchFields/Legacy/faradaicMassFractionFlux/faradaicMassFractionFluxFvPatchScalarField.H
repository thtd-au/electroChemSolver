#ifndef faradaicMassFractionFluxFvPatchScalarField_H
#define faradaicMassFractionFluxFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

namespace Foam
{

class faradaicMassFractionFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    scalar stoichCoeff_;
    scalar nElectrons_;
    scalar currentDensity_;

public:

    TypeName("faradaicMassFractionFlux");       // name used in type field for BC in 0/species 

    faradaicMassFractionFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    faradaicMassFractionFluxFvPatchScalarField
    (
        const faradaicMassFractionFluxFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    faradaicMassFractionFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    faradaicMassFractionFluxFvPatchScalarField
    (
        const faradaicMassFractionFluxFvPatchScalarField&
    );

    faradaicMassFractionFluxFvPatchScalarField
    (
        const faradaicMassFractionFluxFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new faradaicMassFractionFluxFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new faradaicMassFractionFluxFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

}

#endif

