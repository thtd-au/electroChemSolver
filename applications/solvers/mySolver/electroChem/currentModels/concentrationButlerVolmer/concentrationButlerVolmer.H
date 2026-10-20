#ifndef concentrationButlerVolmer_H
#define concentrationButlerVolmer_H

#include "ButlerVolmer.H"

namespace Foam
{

// Reference-state kinetics with arbitrary configured species on each electrode.
// Concentration kinetics and a direct boundary-response Taylor correction.
// The controller, relaxation and exact acceptance tests remain shared with BV.
class concentrationButlerVolmer : public ButlerVolmer
{
    List<wordList> dependencyNames_;
    List<scalarList> referenceConcentrations_, anodicOrders_, cathodicOrders_;

    List<scalarList> speciesNu_, diffusivity_, chargeNumber_;
    word responseMethod_;
    bool responseDiagnostics_;
    labelList responseStatus_; // 0=none/startup, 1=Taylor, 2=undefined response
    scalarList denominatorMin_, denominatorMax_, responseStepScale_;
    scalarList potentialIncrement_;
    List<List<scalarField>> predictedFaces_;
    scalarList predictionAbsError_, predictionRelError_;
    bool havePrediction_, predictionTested_;
    scalarList linearizationPhiS_;
    List<List<scalarField>> beforeFaces_;
    List<scalarField> predictionPhiE_, previousPhiE_, previousCurrent_;
    List<scalarField> previousLogCa_, previousLogCc_;
    scalarList previousPhiS_;
    bool havePreviousState_;
    label diagnosticTimeIndex_;
    scalarList predictedBVResidual_, kineticPredictionError_, tangentRemainder_;
    scalarList potentialSensitivity_, potentialTolerance_, sweepPotentialEffect_;
    scalarList electrolyteEffect_, netEtaEffect_, compositionEffect_, currentStep_;
    scalarList decompositionError_;
    List<bool> predictionCurrentValid_, tangentValid_, couplingValid_;

    void evaluateResponseDiagnostics();

    virtual scalar adjustCurrentCandidate
    (
        const label electrode, const scalarField& baseCurrent,
        scalarField& candidate, const scalarField& logCa,
        const scalarField& logCc, bool& historyValid
    );

    virtual void logConcentrationFactors
    (
        const label electrode, scalarField& logCa, scalarField& logCc
    ) const;

public:
    TypeName("concentrationButlerVolmer");

    concentrationButlerVolmer
    (
        fvMesh& mesh, const dictionary& currentDict, const volScalarField& phiE
    );

    virtual void correct();
    virtual void evaluateConstraint();
    virtual void report() const;
};

} // namespace Foam
#endif

