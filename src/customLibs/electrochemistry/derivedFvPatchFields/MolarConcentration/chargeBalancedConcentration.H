#ifndef chargeBalancedConcentration_H
#define chargeBalancedConcentration_H
#include "fixedValueFvPatchFields.H"
namespace Foam
{
class chargeBalancedConcentration : public fixedValueFvPatchScalarField
{
public:
    TypeName("chargeBalancedConcentration");
    chargeBalancedConcentration(const fvPatch&, const DimensionedField<scalar, volMesh>&);
    chargeBalancedConcentration(const fvPatch&, const DimensionedField<scalar, volMesh>&, const dictionary&);
    chargeBalancedConcentration(const chargeBalancedConcentration&, const fvPatch&, const DimensionedField<scalar, volMesh>&, const fvPatchFieldMapper&);
    chargeBalancedConcentration(const chargeBalancedConcentration&);
    chargeBalancedConcentration(const chargeBalancedConcentration&, const DimensionedField<scalar, volMesh>&);
    virtual tmp<fvPatchScalarField> clone() const
    { return tmp<fvPatchScalarField>(new chargeBalancedConcentration(*this)); }
    virtual tmp<fvPatchScalarField> clone(const DimensionedField<scalar, volMesh>& f) const
    { return tmp<fvPatchScalarField>(new chargeBalancedConcentration(*this, f)); }
    virtual void updateCoeffs();
    virtual void write(Ostream&) const;
};
}
#endif

