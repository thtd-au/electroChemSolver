#ifndef reactionManager_H
#define reactionManager_H

#include "fvMesh.H"
#include "fvMatrices.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "PtrList.H"
#include "runTimeSelectionTables.H"
#include "scalarMatrices.H"
#include "volFields.H"

namespace Foam
{

class reactionModel
{
public:
    enum categoryType { kinetic, equilibrium };

protected:
    word name_;

public:
    TypeName("reactionModel");

    declareRunTimeSelectionTable
    (
        autoPtr, reactionModel, dictionary,
        (const word& name, const dictionary& dict, const wordList& species),
        (name, dict, species)
    );

    explicit reactionModel(const word& name) : name_(name) {}
    virtual ~reactionModel() = default;

    static autoPtr<reactionModel> New
    (
        const word& modelType,
        const word& name,
        const dictionary& dict,
        const wordList& species
    );

    const word& name() const { return name_; }
    virtual categoryType category() const = 0;

    virtual scalar kineticRate
    (
        const label celli,
        const scalarField& concentrations,
        const scalarField& nu
    ) const;

    virtual void equilibriumResidualAndJacobian
    (
        const label celli,
        const scalarField& concentrations,
        const scalarField& nu,
        scalar& residual,
        scalarField& dFdc,
        scalar& productTerm,
        scalar& reactantTerm
    ) const;
};


class reactionManager
{
    struct reactionRecord
    {
        word name;
        word type;
        labelList speciesIndices;
        scalarList stoichCoeffs;
        autoPtr<reactionModel> model;
        reactionRecord();
    };

    fvMesh& mesh_;
    const dictionary& reactionsDict_;
    const wordList& species_;
    PtrList<volScalarField>& c_;
    label inertIndex_;

    List<reactionRecord> reactions_;
    labelList kineticReactionIndices_;
    labelList equilibriumReactionIndices_;

    List<scalarField> kineticRates_;       // [kinetic reaction][cell]
    List<scalarField> lambda_;             // [equilibrium reaction][cell]
    List<scalarField> transportResponse_;  // [species][cell], V/aP
    List<scalarField> residual_;           // [equilibrium reaction][cell]
    List<scalarField> residualScale_;      // [equilibrium reaction][cell]

    scalarList equilibriumAbsTol_;
    scalarList equilibriumRelTol_;
    scalar equilibriumRelaxation_;
    bool equilibriumFBL_;
    scalar equilibriumFractionToBoundary_;

    // Dense equilibrium stoichiometry N(species,equilibrium reaction).
    scalarSquareMatrix equilibriumGram_; // startup rank validation workspace/result
    List<scalarField> equilibriumN_;      // rows by species

    bool equilibriumEvaluated_;
    bool equilibriumConverged_;

    label speciesIndex(const word& name, const dictionary& dict) const;
    void readReactions();
    void readEquilibriumControls();
    void buildAndValidateEquilibriumStoichiometry();
    void validateStoichiometry
    (
        const word& reactionName,
        const dictionary& reactionDict,
        reactionRecord& record
    ) const;

    scalar reactionStoich(const reactionRecord& record, const label speciesI) const;
    bool cellEquilibriumConverged(const label celli) const;

    // Small dense solve with scaled partial pivoting. Returns false if singular.
    bool solveLocalSystem
    (
        scalarSquareMatrix A,
        scalarField b,
        scalarField& x,
        scalar& minRelativePivot
    ) const;

public:
    reactionManager
    (
        fvMesh& mesh,
        const dictionary& reactionsDict,
        const wordList& species,
        PtrList<volScalarField>& c,
        const label inertIndex
    );

    void resetEquilibriumSources();
    void updateKineticRates();

    tmp<volScalarField> source
    (
        const label speciesI,
        scalarField* kineticSource = nullptr,
        scalarField* equilibriumSource = nullptr
    ) const;

    void setTransportResponse
    (
        const label speciesI,
        const fvScalarMatrix& cEqn
    );

    void evaluateEquilibrium();
    bool equilibriumConverged() const;
    void correctEquilibrium();
    void reportEquilibrium() const;

    label nKineticReactions() const { return kineticReactionIndices_.size(); }
    label nEquilibriumReactions() const { return equilibriumReactionIndices_.size(); }
};

} // End namespace Foam

#endif

