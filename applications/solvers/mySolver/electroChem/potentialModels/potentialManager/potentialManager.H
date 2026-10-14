#ifndef potentialManager_H
#define potentialManager_H

#include "fvMesh.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "PtrList.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "error.H"
#include "boolList.H"
#include "token.H"

namespace Foam
{

// Shared output policy for the driver and both electrochemical managers.
inline label picardOutputLevel(const dictionary& controls)
{
    const word level = controls.lookupOrDefault<word>("picardOutput", "picard");
    if (level == "minimal") return 0;
    if (level == "picard") return 1;
    if (level == "full") return 2;
    if (level == "debug") return 3;
    FatalIOErrorInFunction(controls)
        << "picardOutput must be minimal, picard, full, or debug"
        << exit(FatalIOError);
    return 0;
}

// true/false or an unparenthesized sequence of component indices.
inline boolList potentialPrintDirections(const dictionary& controls, const word& key)
{
    boolList selected(3, false);
    if (!controls.found(key)) return selected;
    const ITstream& input = controls.lookup(key);
    if (input.size() == 1 && input[0].isWord())
    {
        const word value(input[0].wordToken());
        if (value == "true") return boolList(3, true);
        if (value == "false") return selected;
    }
    if (input.empty())
        FatalIOErrorInFunction(controls)
            << key << " requires true, false, or indices 0=x, 1=y, 2=z"
            << exit(FatalIOError);
    forAll(input, i)
    {
        if (!input[i].isLabel() || input[i].labelToken() < 0 || input[i].labelToken() > 2)
            FatalIOErrorInFunction(controls)
                << key << " requires true, false, or indices 0=x, 1=y, 2=z"
                << exit(FatalIOError);
        selected[input[i].labelToken()] = true;
    }
    return selected;
}

class potentialModel
{
protected:
    fvMesh& mesh_;
    const dictionary& electroChemistry_;
    const dictionary& reactions_;
    const wordList& species_;
    PtrList<volScalarField>& c_;
    volScalarField& phiE_;
    const dimensionedScalar& T_;
    const dimensionedScalar& Rgas_;
    const dimensionedScalar& F_;

public:
    TypeName("potentialModel");

    declareRunTimeSelectionTable
    (
        autoPtr, potentialModel, dictionary,
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
        ),
        (mesh, electroChemistry, reactions, species, c, phiE, T, Rgas, F)
    );

    potentialModel
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
    );

    virtual ~potentialModel() = default;

    static autoPtr<potentialModel> New
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
    );

    virtual void evaluateConstraint() = 0;
    virtual bool converged() const = 0;
    virtual void correct() = 0;
    virtual void report() const = 0;
    virtual const scalarField& chargeResidual() const = 0;
    virtual const scalarField& chargeScale() const = 0;
    virtual void refreshDiagnosticCurrent() = 0;
    virtual const volScalarField& conductivityField() const = 0;
    virtual const volVectorField& diffusionCurrentField() const = 0;
};


class potentialManager
{
    autoPtr<potentialModel> model_;

public:
    potentialManager
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
    );

    bool active() const { return model_->type() != "none"; }

    void evaluateConstraint() { model_->evaluateConstraint(); }
    bool converged() const { return model_->converged(); }
    void correct() { model_->correct(); }
    void report() const { model_->report(); }
    const scalarField& chargeResidual() const { return model_->chargeResidual(); }
    const scalarField& chargeScale() const { return model_->chargeScale(); }
    void refreshDiagnosticCurrent() { model_->refreshDiagnosticCurrent(); }
    const volScalarField& conductivityField() const { return model_->conductivityField(); }
    const volVectorField& diffusionCurrentField() const { return model_->diffusionCurrentField(); }
};

} // End namespace Foam

#endif

