#ifndef ButlerVolmer_H
#define ButlerVolmer_H

#include "currentManager.H"
#include "scalarList.H"
#include "labelList.H"
#include "scalarField.H"

namespace Foam
{

// Fixed i0, alpha_a, alpha_c and Eeq per equipotential electrode patch.
// The imposed magnitude defines total current through a reference area.
class ButlerVolmer : public currentModel
{
    scalar appliedCurrentDensity_;
    word referencePatchName_;
    wordList electrodeNames_;
    labelList patchIDs_;
    scalarList nuE_, i0_, alphaA_, alphaC_, Eeq_, phiS_;
    scalarList targetCurrent_, area_, actualCurrent_;
    scalarList bvResidual_, boundaryResidual_;
    labelList rootIterations_;
    scalar T_, F_, R_;
    scalar absTol_, relTol_, rootAbsTol_, rootRelTol_;
    label rootMaxIter_;
    scalar currentRelaxation_;
    word relaxationMethod_;
    scalar relaxationMin_, relaxationMax_, activeRelaxation_;
    List<scalarField> previousResidual_;
    label relaxationTimeIndex_;
    bool haveResidual_;
    bool initialized_, evaluated_, converged_;

    scalar currentTolerance(const scalar value) const;

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
};

} // namespace Foam

#endif

