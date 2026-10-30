#include "zeroCurrentPotentialFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include <cmath>
namespace Foam
{
defineTypeNameAndDebug(zeroCurrentPotentialFlux, 0);
addToRunTimeSelectionTable(fvPatchScalarField, zeroCurrentPotentialFlux, dictionary);
addToRunTimeSelectionTable(fvPatchScalarField, zeroCurrentPotentialFlux, patchMapper);
zeroCurrentPotentialFlux::zeroCurrentPotentialFlux(const fvPatch& p, const DimensionedField<scalar, volMesh>& f)
: fixedGradientFvPatchScalarField(p, f) { gradient() = 0.0; }
zeroCurrentPotentialFlux::zeroCurrentPotentialFlux(const fvPatch& p, const DimensionedField<scalar, volMesh>& f, const dictionary& dict)
: fixedGradientFvPatchScalarField(p, f)
{
    gradient() = 0.0;
    if (dict.found("gradient")) gradient() = scalarField("gradient", dict, p.size());
    if (dict.found("value")) fvPatchField<scalar>::operator=(scalarField("value", dict, p.size()));
    // Defer coupled lookup until all species/current fields are registered.
}
zeroCurrentPotentialFlux::zeroCurrentPotentialFlux(const zeroCurrentPotentialFlux& b, const fvPatch& p, const DimensionedField<scalar, volMesh>& f, const fvPatchFieldMapper& m)
: fixedGradientFvPatchScalarField(b, p, f, m) {}
zeroCurrentPotentialFlux::zeroCurrentPotentialFlux(const zeroCurrentPotentialFlux& b) : fixedGradientFvPatchScalarField(b) {}
zeroCurrentPotentialFlux::zeroCurrentPotentialFlux(const zeroCurrentPotentialFlux& b, const DimensionedField<scalar, volMesh>& f)
: fixedGradientFvPatchScalarField(b, f) {}
void zeroCurrentPotentialFlux::updateCoeffs()
{
    if (updated()) return;
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const IOdictionary& ec = mesh.lookupObject<IOdictionary>("electroChemistry");
    const dictionary& controls = ec.subDict("potential");
    if (controls.get<word>("model") == "none")
    {
        gradient() = 0.0;
        fixedGradientFvPatchScalarField::updateCoeffs();
        return;
    }
    const surfaceScalarField& diffusion = mesh.lookupObject<surfaceScalarField>("jDiffFlux");
    const surfaceScalarField& conductivity = mesh.lookupObject<surfaceScalarField>("kappaFace");
    const label patchi = patch().index();
    const scalar minimum = controls.getOrDefault<scalar>("minConductivity", SMALL);
    forAll(gradient(), facei)
    {
        const scalar area = mesh.magSf().boundaryField()[patchi][facei];
        const scalar kappa = conductivity.boundaryField()[patchi][facei];
        if (!std::isfinite(kappa) || kappa < minimum || kappa <= 0 || area <= VSMALL)
            FatalErrorInFunction << "Invalid area/conductivity on " << patch().name()
                << " face " << facei << ": kappa=" << kappa << exit(FatalError);
        // Registered diffusion flux is integrated current [A]. No species
        // evaluation here: preserve the potential model's Picard lag.
        gradient()[facei] = diffusion.boundaryField()[patchi][facei]/(area*kappa);
    }
    fixedGradientFvPatchScalarField::updateCoeffs();
}
void zeroCurrentPotentialFlux::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
}
}
