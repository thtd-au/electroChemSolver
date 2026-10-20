#ifndef ButlerVolmer_H
#define ButlerVolmer_H

#include "currentManager.H"
#include "scalarList.H"
#include "labelList.H"
#include "scalarField.H"
#include "OFstream.H"

namespace Foam
{

// Shared equipotential-electrode BV controller; ordinary kinetics use fixed
// i0/alpha/Eeq, and derived kinetics supply logarithmic concentration factors.
// The imposed magnitude defines total current through a reference area.
class ButlerVolmer : public currentModel
{
protected:
    scalar appliedCurrentDensity_;
    word referencePatchName_;
    wordList electrodeNames_;
    labelList patchIDs_;
    // Prescribed Eeq (ordinary BV) or Eref (concentration-dependent BV).
    scalarList nuE_, i0_, alphaA_, alphaC_, Ebase_, phiS_;
    scalarList targetCurrent_, area_, actualCurrent_;
    scalarList bvResidual_, boundaryResidual_;
    labelList rootIterations_;
    // Reporting only: proposed redistribution and applied changes for the
    // current transported-and-tested state.
    scalarList currentChange_, appliedCurrentChange_, metalPotentialChange_;
    scalar T_, F_, R_;
    bool referenceState_;
    scalar absTol_, relTol_, rootAbsTol_, rootRelTol_;
    label rootMaxIter_;
    scalar currentRelaxation_;
    word relaxationMethod_;
    // Aitken starts at currentRelaxation_ and is bounded by [relaxationMin_, 1].
    scalar relaxationMin_, activeRelaxation_;
    List<scalarField> previousResidual_;
    label relaxationTimeIndex_;
    bool haveResidual_;
    bool initialized_, evaluated_, converged_;

    autoPtr<OFstream> electrodePotentialLog_, cellVoltageLog_;

    scalar currentTolerance(const scalar value) const;
    // Secondary constructor defers correct() until derived parameters exist.
    ButlerVolmer
    (
        fvMesh&, const dictionary&, const volScalarField&,
        const bool initialise, const bool referenceState
    );
    virtual void logConcentrationFactors
    (
        const label electrode, scalarField& logCa, scalarField& logCc
    ) const;

    // Extension point for concentration-specific interface preconditioning.
    // Return a metal-potential increment paired with the candidate change.
    // Ordinary BV leaves both the candidate and history untouched.
    virtual scalar adjustCurrentCandidate
    (
        const label, const scalarField&, scalarField&,
        const scalarField&, const scalarField&, bool&
    ) { return 0; }

    // Optional batch correction after every frozen-state electrode root exists.
    virtual void adjustCurrentCandidates
    (
        const List<scalarField>&, List<scalarField>&, scalarList&, bool&
    ) {}

public:
    TypeName("ButlerVolmer");

    ButlerVolmer
    (
        fvMesh& mesh,
        const dictionary& currentDict,
        const volScalarField& phiE
    );

    virtual void correct();
    virtual bool requiresPotentialCoupling() const { return true; }
    virtual void evaluateConstraint();
    virtual bool converged() const;
    virtual void report() const;
    virtual void writeDiagnostics();
};

} // namespace Foam

#endif

