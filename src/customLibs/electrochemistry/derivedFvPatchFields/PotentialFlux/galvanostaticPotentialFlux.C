#include "galvanostaticPotentialFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"

namespace Foam
{

defineTypeNameAndDebug(galvanostaticPotentialFlux, 0);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    galvanostaticPotentialFlux,
    dictionary
);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    galvanostaticPotentialFlux,
    patchMapper
);


galvanostaticPotentialFlux::galvanostaticPotentialFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    electronStoichCoeff_(1.0)
{
    gradient() = 0.0;
}


galvanostaticPotentialFlux::galvanostaticPotentialFlux
(
    const galvanostaticPotentialFlux& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    electronStoichCoeff_(ptf.electronStoichCoeff_)
{}


galvanostaticPotentialFlux::galvanostaticPotentialFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    electronStoichCoeff_
    (
        dict.get<scalar>("electronStoichCoeff")
    )
{
    if (mag(electronStoichCoeff_) <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "electronStoichCoeff must be non-zero"
            << exit(FatalIOError);
    }

    // The gradient is the Picard-lagged electrode state consumed by
    // galvanostaticMolarFlux.  A fresh case may omit it (zero is then the
    // initial guess); restart files written by this BC retain it.
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        gradient() = 0.0;
    }

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=
        (
            patchInternalField()
        );
    }
}


galvanostaticPotentialFlux::galvanostaticPotentialFlux
(
    const galvanostaticPotentialFlux& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    electronStoichCoeff_(ptf.electronStoichCoeff_)
{}


galvanostaticPotentialFlux::galvanostaticPotentialFlux
(
    const galvanostaticPotentialFlux& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    electronStoichCoeff_(ptf.electronStoichCoeff_)
{}


void galvanostaticPotentialFlux::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh =
        patch().boundaryMesh().mesh();

    const IOdictionary& electroChemistry =
        mesh.lookupObject<IOdictionary>
        (
            "electroChemistry"
        );

    const dictionary& potentialDict =
        electroChemistry.subDict("potential");

    if (potentialDict.get<word>("model") == "none")
    {
        // Keep electronStoichCoeff available to the species flux BC without
        // computing electric fields or requiring conductivity/current fields.
        gradient() = 0.0;
        fixedGradientFvPatchScalarField::updateCoeffs();
        return;
    }

    const scalar minKappa =
        potentialDict.getOrDefault<scalar>
        (
            "minConductivity",
            SMALL
        );

    if (!mesh.foundObject<surfaceScalarField>("electrodeCurrentDensity"))
    {
        FatalErrorInFunction
            << "galvanostaticPotentialFlux requires the registered "
            << "electrodeCurrentDensity field supplied by currentManager."
            << exit(FatalError);
    }
    const surfaceScalarField& electrodeCurrentDensity =
        mesh.lookupObject<surfaceScalarField>("electrodeCurrentDensity");

    // Diffusion-current and conductivity fields belong to the potential model.
    if
    (
        !mesh.foundObject<surfaceScalarField>("jDiffFlux")
     || !mesh.foundObject<surfaceScalarField>("kappaFace")
    )
    {
        FatalErrorInFunction
            << "galvanostaticPotentialFlux requires registered "
            << "surfaceScalarFields jDiffFlux and kappaFace." << nl
            << "They must be registered by the active potential model before boundary evaluation."
            << exit(FatalError);
    }

    const surfaceScalarField& jDiffFlux =
        mesh.lookupObject<surfaceScalarField>
        (
            "jDiffFlux"
        );

    const surfaceScalarField& kappaFace =
        mesh.lookupObject<surfaceScalarField>
        (
            "kappaFace"
        );

    const label patchi =
        patch().index();

    const fvsPatchScalarField& jDiffFluxPatch =
        jDiffFlux.boundaryField()[patchi];

    const fvsPatchScalarField& kappaFacePatch =
        kappaFace.boundaryField()[patchi];

    const fvsPatchScalarField& magSfPatch =
        mesh.magSf().boundaryField()[patchi];

    // Signed outward ionic current density [A/m2], shared with species BCs.
    const fvsPatchScalarField& jnApplied =
        electrodeCurrentDensity.boundaryField()[patchi];

    scalarField patchGradient
    (
        patch().size(),
        0.0
    );

    forAll(patchGradient, facei)
    {
        const scalar area =
            magSfPatch[facei];

        if (area <= VSMALL)
        {
            FatalErrorInFunction
                << "galvanostaticPotentialFlux encountered "
                << "zero-area face." << nl
                << "patch = " << patch().name() << nl
                << "face = " << facei
                << exit(FatalError);
        }

        const scalar kappaF =
            kappaFacePatch[facei];

        if (kappaF < minKappa)
        {
            FatalErrorInFunction
                << "galvanostaticPotentialFlux encountered "
                << "face conductivity below minConductivity."
                << nl
                << "patch = " << patch().name() << nl
                << "face = " << facei << nl
                << "kappaFace = " << kappaF
                << " S/m" << nl
                << "minConductivity = "
                << minKappa << " S/m"
                << exit(FatalError);
        }

        // jDiffFlux is integrated current [A].
        const scalar jDiffNormal =
            jDiffFluxPatch[facei]/area;

        // j_n = jDiff_n - kappa*d(phi)/dn
        //
        // => d(phi)/dn
        //      = (jDiff_n - j_n,app)/kappa
        patchGradient[facei] =
            (jDiffNormal - jnApplied[facei])
           /kappaF;
    }

    gradient() = patchGradient;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void galvanostaticPotentialFlux::write
(
    Ostream& os
) const
{
    // Use the fixedGradient writer so the current Picard-lagged gradient is
    // retained in restart files together with the patch value.
    fixedGradientFvPatchScalarField::write(os);
    os.writeEntry
    (
        "electronStoichCoeff",
        electronStoichCoeff_
    );
}

} // End namespace Foam

