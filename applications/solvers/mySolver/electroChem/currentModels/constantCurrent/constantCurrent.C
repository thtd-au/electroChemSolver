#include "constantCurrent.H"
#include "addToRunTimeSelectionTable.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "token.H"
#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(constantCurrent, 0);
addToRunTimeSelectionTable(currentModel, constantCurrent, dictionary);

constantCurrent::constantCurrent
(
    fvMesh& mesh,
    const dictionary& currentDict,
    const volScalarField& phiE
)
:
    currentModel(mesh, currentDict, phiE),
    appliedCurrentDensity_(currentDict.get<scalar>("appliedCurrentDensity"))
{
    if (!std::isfinite(appliedCurrentDensity_) || appliedCurrentDensity_ < 0.0)
    {
        FatalIOErrorInFunction(currentDict_)
            << "appliedCurrentDensity must be finite and non-negative"
            << exit(FatalIOError);
    }

    // Populate before potential initialization can evaluate electrode BCs.
    correct();
}

scalar constantCurrent::electronStoichCoeff
(
    const fvPatchScalarField& patchField
)
{
    const word& patchName = patchField.patch().name();
    if (!electronStoichCoeffs_.found(patchName))
    {
        // Virtual dispatch reaches the loaded BC's writer without including
        // or linking its concrete C++ class. Its existing write() retains
        // electronStoichCoeff. This does not evaluate/update the patch field.
        OStringStream patchOutput;
        patchOutput << token::BEGIN_BLOCK << nl;
        patchField.write(patchOutput);
        patchOutput << token::END_BLOCK << nl;
        IStringStream patchInput(patchOutput.str());
        const dictionary patchDict(patchInput);
        const scalar nuE = patchDict.get<scalar>("electronStoichCoeff");

        if (!std::isfinite(nuE) || mag(nuE) <= SMALL)
        {
            FatalErrorInFunction
                << "Non-finite or zero electronStoichCoeff on electrode patch "
                << patchName << exit(FatalError);
        }
        electronStoichCoeffs_.insert(patchName, nuE);
    }
    return electronStoichCoeffs_[patchName];
}

void constantCurrent::correct()
{
    electrodeCurrentDensity_.primitiveFieldRef() = 0.0;

    forAll(phiE_.boundaryField(), patchi)
    {
        const fvPatchScalarField& potentialPatch = phiE_.boundaryField()[patchi];
        fvsPatchScalarField& currentPatch =
            electrodeCurrentDensity_.boundaryFieldRef()[patchi];

        currentPatch = 0.0;
        if (potentialPatch.type() != "galvanostaticPotentialFlux") continue;

        const scalar nuE = electronStoichCoeff(potentialPatch);

        // Positive outward ionic current at the cathode; negative at anode.
        // Every electrode retains the old uniform density magnitude. There
        // is deliberately no area normalization or sign forcing in the BCs.
        currentPatch = -nuE/mag(nuE)*appliedCurrentDensity_;
    }
}

} // End namespace Foam

