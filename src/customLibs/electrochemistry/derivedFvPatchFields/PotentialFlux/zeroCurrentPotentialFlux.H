#ifndef zeroCurrentPotentialFlux_H
#define zeroCurrentPotentialFlux_H
#include "fixedGradientFvPatchFields.H"
namespace Foam
{
class zeroCurrentPotentialFlux : public fixedGradientFvPatchScalarField
{
public:
    TypeName("zeroCurrentPotentialFlux");
    zeroCurrentPotentialFlux(const fvPatch&, const DimensionedField<scalar, volMesh>&);
    zeroCurrentPotentialFlux(const fvPatch&, const DimensionedField<scalar, volMesh>&, const dictionary&);
    zeroCurrentPotentialFlux(const zeroCurrentPotentialFlux&, const fvPatch&, const DimensionedField<scalar, volMesh>&, const fvPatchFieldMapper&);
    zeroCurrentPotentialFlux(const zeroCurrentPotentialFlux&);
    zeroCurrentPotentialFlux(const zeroCurrentPotentialFlux&, const DimensionedField<scalar, volMesh>&);
    virtual tmp<fvPatchScalarField> clone() const
    { return tmp<fvPatchScalarField>(new zeroCurrentPotentialFlux(*this)); }
    virtual tmp<fvPatchScalarField> clone(const DimensionedField<scalar, volMesh>& f) const
    { return tmp<fvPatchScalarField>(new zeroCurrentPotentialFlux(*this, f)); }
    virtual void updateCoeffs();
    virtual void write(Ostream&) const;
};
}
#endif
