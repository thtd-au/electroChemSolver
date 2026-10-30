#ifndef concentrationButlerVolmer_H
#define concentrationButlerVolmer_H

#include "ButlerVolmer.H"

namespace Foam
{

// Reference-state kinetics with arbitrary configured species on each electrode.
// Concentration kinetics with conservative concentration/potential response.
// The controller, relaxation and exact acceptance tests remain shared with BV.
class concentrationButlerVolmer : public ButlerVolmer
{
    List<wordList> dependencyNames_;
    List<scalarList> referenceConcentrations_, anodicOrders_, cathodicOrders_;

    List<scalarList> speciesNu_, diffusivity_, chargeNumber_;
    word responseMethod_;
    labelList responseStatus_; // 0=none/startup, 1=Taylor, 2=undefined response
    scalarList responseStepScale_;

    List<scalarField> potentialDiagonal_, potentialSlope_, potentialRhs_;
    List<List<scalarField>> potentialAnchor_, potentialDcDj_;
    label potentialResponseStatus_, potentialSolveIterations_;
    scalar potentialCurrentDefect_, potentialRowResidual_, potentialStepScale_;
    scalarList predictedPotentialMax_;

    virtual void adjustCurrentCandidates
    (
        const List<scalarField>& baseCurrent, List<scalarField>& candidates,
        scalarList& potentialIncrement, bool& historyValid
    );

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

    virtual void report() const;
};

} // namespace Foam
#endif

