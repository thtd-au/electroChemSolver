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

// No potential constraints, currents or corrections in isolation mode.
class noPotentialModel : public potentialModel
{
    scalarField residual_, scale_;
    volScalarField kappa_;
    volVectorField jDiff_;
public:
    TypeName("none");
    noPotentialModel
    (
        fvMesh& mesh, const dictionary& ec, const dictionary& reactions,
        const wordList& species, PtrList<volScalarField>& c,
        volScalarField& phiE, const dimensionedScalar& T,
        const dimensionedScalar& Rgas, const dimensionedScalar& F
    )
    : potentialModel(mesh, ec, reactions, species, c, phiE, T, Rgas, F),
      residual_(mesh.nCells(), 0.0), scale_(mesh.nCells(), 1.0),
      kappa_(IOobject("kappaModelInternal", mesh.time().timeName(), mesh,
          IOobject::NO_READ, IOobject::NO_WRITE, false), mesh,
          dimensionedScalar("zero", dimensionSet(-1,-3,3,0,0,2,0), 0.0)),
      jDiff_(IOobject("jDiffModelInternal", mesh.time().timeName(), mesh,
          IOobject::NO_READ, IOobject::NO_WRITE, false), mesh,
          dimensionedVector("zero", dimensionSet(0,-2,0,0,0,1,0), vector::zero))
    {}
    void evaluateConstraint() override {}
    bool converged() const override { return true; }
    void correct() override {}
    void report() const override {}
    const scalarField& chargeResidual() const override { return residual_; }
    const scalarField& chargeScale() const override { return scale_; }
    void refreshDiagnosticCurrent() override {}
    const volScalarField& conductivityField() const override { return kappa_; }
    const volVectorField& diffusionCurrentField() const override { return jDiff_; }
};
defineTypeNameAndDebug(noPotentialModel, 0);

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

    if (modelType == "none")
    {
        return autoPtr<potentialModel>(new noPotentialModel
        (mesh, electroChemistry, reactions, species, c, phiE, T, Rgas, F));
    }

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

