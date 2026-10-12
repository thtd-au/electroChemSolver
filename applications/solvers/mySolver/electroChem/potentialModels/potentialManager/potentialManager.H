#ifndef potentialManager_H
#define potentialManager_H

#include "fvMesh.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "PtrList.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

namespace Foam
{

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

