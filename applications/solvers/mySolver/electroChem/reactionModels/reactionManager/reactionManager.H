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
#include "boolList.H"

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
    List<scalarField> lastAcceptedLambda_; // Accepted physical-time state only
    bool haveAcceptedEquilibriumSources_;
    word equilibriumInitialGuess_;        // zero | lastAccepted
    bool usingAcceptedEquilibriumGuess_;
    List<scalarField> transportResponse_;  // [species][cell], V/aP
    List<scalarField> residual_;           // [equilibrium reaction][cell]
    List<scalarField> residualScale_;      // [equilibrium reaction][cell]

    scalarList equilibriumAbsTol_;
    scalarList equilibriumRelTol_;
    scalar equilibriumRelaxation_;
    bool equilibriumFBL_;
    scalar equilibriumFractionToBoundary_;
    dictionary responseSolver_;
    PtrList<fvScalarMatrix> transportMatrices_;
    label coupledRestart_;
    label coupledMaxIter_;
    scalar coupledRelTol_;
    scalar coupledAbsTol_;
    label coupledMaxBacktracks_;
    scalar coupledReduction_;
    scalar coupledSufficientDecrease_;
    // All response calculations are read-only predictions. Actual kinetic
    // rates are still refreshed only by updateKineticRates() in the outer loop.
    List<scalarField> kineticRefreshResponse() const;
    List<scalarField> transportMatrixResponse
    (
        const List<scalarField>& direction,
        const List<scalarField>* speciesSource = nullptr
    ) const;
    List<scalarField> coupledNewtonDirection
    (
        const List<scalarField>& baseline
    ) const;
    void applyTransportMatrixCorrection
    (
        const List<scalarField>& direction,
        const List<scalarField>& baseline,
        const List<scalarField>& kineticResponse,
        const scalar relaxationCap
    );
    label outputLevel_; // 0 minimal, 1 picard, 2 full, 3 debug
    bool equilibriumDiagnostics_;
    bool diagnosticStepValid_;
    List<scalarField> diagnosticBefore_;
    List<scalarField> diagnosticDeltaC_;
    List<scalarField> diagnosticResponse_;
    List<scalarField> diagnosticKineticBefore_;
    List<scalarField> diagnosticKineticDeltaC_;
    List<scalarField> diagnosticDeltaLambda_;
    scalarField diagnosticAlpha_;
    labelList diagnosticLimiter_;
    boolList diagnosticTransported_;

    // Dense equilibrium stoichiometry N(species,equilibrium reaction).
    scalarSquareMatrix equilibriumGram_; // startup rank validation workspace/result
    List<scalarField> equilibriumN_;      // rows by species

    bool equilibriumEvaluated_;
    bool equilibriumConverged_;

    label speciesIndex(const word& name, const dictionary& dict) const;
    void readReactions();
    void readEquilibriumControls();
    void readEquilibriumInitialGuess(const dictionary& controls);
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

    // Initialize a new timestep and invalidate all old acceptance flags.
    void resetEquilibriumSources();
    // Called only after every electrochemical constraint has passed.
    void acceptEquilibriumSources();
    bool usingAcceptedEquilibriumGuess() const
    { return usingAcceptedEquilibriumGuess_; }
    void scaleInitialEquilibriumGuess(const scalar scale);
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
    bool equilibriumDiagnostics() const { return equilibriumDiagnostics_; }
    void readOutputControls();
    label transportLogLevel() const { return outputLevel_ == 0 ? 0 : 1; }
    void reportNoEquilibriumCorrection(const char* reason) const;
    void beginTransportDiagnostics();
    void reportTransportDiagnostics(const label speciesI);

    label nKineticReactions() const { return kineticReactionIndices_.size(); }
    label nEquilibriumReactions() const { return equilibriumReactionIndices_.size(); }
};

} // End namespace Foam

#endif

