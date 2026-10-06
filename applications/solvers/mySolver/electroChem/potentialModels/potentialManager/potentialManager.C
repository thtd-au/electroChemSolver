#include "potentialManager.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(potentialModel, 0);
defineRunTimeSelectionTable(potentialModel, dictionary);

potentialModel::potentialModel
(
    fvMesh& mesh,
    const dictionary& electroChemistry,
    const dictionary& reactions,
    const wordList& species,
    PtrList<volScalarField>& c,
    volScalarField& phiE,
    const dimensionedScalar& T,
    const dimensionedScalar& Rgas,
    const dimensionedScalar& F
)
:
    mesh_(mesh),
    electroChemistry_(electroChemistry),
    reactions_(reactions),
    species_(species),
    c_(c),
    phiE_(phiE),
    T_(T),
    Rgas_(Rgas),
    F_(F)
{}

autoPtr<potentialModel> potentialModel::New
(
    fvMesh& mesh,
    const dictionary& electroChemistry,
    const dictionary& reactions,
    const wordList& species,
    PtrList<volScalarField>& c,
    volScalarField& phiE,
    const dimensionedScalar& T,
    const dimensionedScalar& Rgas,
    const dimensionedScalar& F
)
{
    if (!electroChemistry.found("potential"))
    {
        FatalIOErrorInFunction(electroChemistry)
            << "Missing required potential dictionary"
            << exit(FatalIOError);
    }

    const dictionary& dict = electroChemistry.subDict("potential");
    const word modelType(dict.get<word>("model"));

    auto* ctorPtr = dictionaryConstructorTable(modelType);
    if (!ctorPtr)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown potential model " << modelType << nl
            << "Valid potential models are: "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<potentialModel>
    (
        ctorPtr(mesh, electroChemistry, reactions, species, c, phiE, T, Rgas, F)
    );
}

potentialManager::potentialManager
(
    fvMesh& mesh,
    const dictionary& electroChemistry,
    const dictionary& reactions,
    const wordList& species,
    PtrList<volScalarField>& c,
    volScalarField& phiE,
    const dimensionedScalar& T,
    const dimensionedScalar& Rgas,
    const dimensionedScalar& F
)
:
    model_
    (
        potentialModel::New
        (
            mesh, electroChemistry, reactions, species, c, phiE, T, Rgas, F
        )
    )
{}

} // End namespace Foam

