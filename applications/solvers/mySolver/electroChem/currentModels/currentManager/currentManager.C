#include "currentManager.H"

namespace Foam
{

defineTypeNameAndDebug(currentModel, 0);
defineRunTimeSelectionTable(currentModel, dictionary);

currentModel::currentModel
(
    fvMesh& mesh,
    const dictionary& currentDict,
    const volScalarField& phiE
)
:
    mesh_(mesh),
    currentDict_(currentDict),
    phiE_(phiE),
    electrodeCurrentDensity_
    (
        IOobject
        (
            "electrodeCurrentDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero", dimensionSet(0,-2,0,0,0,1,0), 0.0
        )
    )
{}

autoPtr<currentModel> currentModel::New
(
    fvMesh& mesh,
    const dictionary& electroChemistry,
    const volScalarField& phiE
)
{
    if (!electroChemistry.found("current"))
    {
        FatalIOErrorInFunction(electroChemistry)
            << "Missing required current dictionary in constant/electroChemistry."
            << nl << "Move appliedCurrentDensity from potential to current "
            << "and select model constantCurrent."
            << exit(FatalIOError);
    }

    const dictionary& dict = electroChemistry.subDict("current");
    const word modelType(dict.get<word>("model"));
    auto* ctorPtr = dictionaryConstructorTable(modelType);
    if (!ctorPtr)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown current model " << modelType << nl
            << "Valid current models are: "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    Info<< "Current model: " << modelType << nl;
    return autoPtr<currentModel>(ctorPtr(mesh, dict, phiE));
}

currentManager::currentManager
(
    fvMesh& mesh,
    const dictionary& electroChemistry,
    const volScalarField& phiE
)
:
    model_(currentModel::New(mesh, electroChemistry, phiE))
{}

} // End namespace Foam

