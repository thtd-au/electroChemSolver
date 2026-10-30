#include "chargeBalancedConcentration.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include <cmath>
namespace Foam
{
defineTypeNameAndDebug(chargeBalancedConcentration, 0);
addToRunTimeSelectionTable(fvPatchScalarField, chargeBalancedConcentration, dictionary);
addToRunTimeSelectionTable(fvPatchScalarField, chargeBalancedConcentration, patchMapper);
chargeBalancedConcentration::chargeBalancedConcentration(const fvPatch& p, const DimensionedField<scalar, volMesh>& f)
: fixedValueFvPatchScalarField(p, f) { fvPatchField<scalar>::operator=(patchInternalField()); }
chargeBalancedConcentration::chargeBalancedConcentration(const fvPatch& p, const DimensionedField<scalar, volMesh>& f, const dictionary& dict)
: fixedValueFvPatchScalarField(p, f)
{
    fvPatchField<scalar>::operator=(patchInternalField());
    
    if (dict.found("value")) fvPatchField<scalar>::operator=(scalarField("value", dict, p.size()));
    // Defer coupled lookup until all species/current fields are registered.
}
chargeBalancedConcentration::chargeBalancedConcentration(const chargeBalancedConcentration& b, const fvPatch& p, const DimensionedField<scalar, volMesh>& f, const fvPatchFieldMapper& m)
: fixedValueFvPatchScalarField(b, p, f, m) {}
chargeBalancedConcentration::chargeBalancedConcentration(const chargeBalancedConcentration& b) : fixedValueFvPatchScalarField(b) {}
chargeBalancedConcentration::chargeBalancedConcentration(const chargeBalancedConcentration& b, const DimensionedField<scalar, volMesh>& f)
: fixedValueFvPatchScalarField(b, f) {}
void chargeBalancedConcentration::updateCoeffs()
{
    if (updated()) return;
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const IOdictionary& reactions = mesh.lookupObject<IOdictionary>("reactions");
    const wordList species(reactions.lookup("species"));
    const dictionary& charges = reactions.subDict("chargeNumber");
    const word fieldName(internalField().name());
    if (fieldName.size() <= 2 || fieldName.substr(0, 2) != "c_")
        FatalErrorInFunction << "chargeBalancedConcentration requires c_<species>: "
            << fieldName << exit(FatalError);
    const word balancingSpecies(fieldName.substr(2));
    const scalar zb = charges.getOrDefault<scalar>(balancingSpecies, 0.0);
    if (species.find(balancingSpecies) < 0 || mag(zb) <= SMALL)
        FatalErrorInFunction << "Invalid or uncharged balancing species "
            << balancingSpecies << exit(FatalError);
    const label patchi = patch().index();
    // Validate all dependencies before evaluating any: exactly one balancing
    // ion, and other charged inlet concentrations prescribed independently.
    forAll(species, i)
    {
        const scalar zi = charges.getOrDefault<scalar>(species[i], 0.0);
        if (species[i] == balancingSpecies || mag(zi) <= SMALL) continue;
        const volScalarField& ci = mesh.lookupObject<volScalarField>("c_" + species[i]);
        const fvPatchScalarField& other = ci.boundaryField()[patchi];
        if (isA<chargeBalancedConcentration>(other)
         || !isA<fixedValueFvPatchScalarField>(other))
            FatalErrorInFunction << "Other charged species must have independent fixed-value inlet conditions. "
                << "Patch " << patch().name() << ", species " << species[i] << exit(FatalError);
    }
    scalarField charge(patch().size(), 0.0);
    forAll(species, i)
    {
        const scalar zi = charges.getOrDefault<scalar>(species[i], 0.0);
        if (species[i] == balancingSpecies || mag(zi) <= SMALL) continue;
        volScalarField& ci = const_cast<volScalarField&>
            (mesh.lookupObject<volScalarField>("c_" + species[i]));
        fvPatchScalarField& other = ci.boundaryFieldRef()[patchi];
        other.evaluate();
        forAll(other, facei)
            if (!std::isfinite(other[facei]) || other[facei] < 0)
                FatalErrorInFunction << "Invalid inlet concentration for " << species[i]
                    << " on " << patch().name() << exit(FatalError);
        charge += zi*other;
    }
    const scalarField balanced(-charge/zb);
    forAll(balanced, facei)
        if (!std::isfinite(balanced[facei]) || balanced[facei] < 0)
            FatalErrorInFunction << "Charge balance requires invalid concentration "
                << balanced[facei] << " for " << balancingSpecies << " on "
                << patch().name() << " face " << facei << exit(FatalError);
    operator==(balanced);
    fixedValueFvPatchScalarField::updateCoeffs();
}
void chargeBalancedConcentration::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}
}

