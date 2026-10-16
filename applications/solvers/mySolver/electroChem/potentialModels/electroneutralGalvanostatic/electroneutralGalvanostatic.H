#ifndef electroneutralGalvanostatic_H
#define electroneutralGalvanostatic_H

#include "potentialManager.H"
#include "surfaceFields.H"

namespace Foam
{

class electroneutralGalvanostatic : public potentialModel
{
    const dictionary& potentialDict_;
    scalar referencePotential_;
    word referencePatchName_;
    label referencePatchID_;
    label referenceCell_;
    scalar minConductivity_;

    scalar absTol_;
    scalar relTol_;
    scalar relaxation_;

    scalarList z_;
    scalarList D_;

    scalarField chargeResidual_;
    scalarField chargeScale_;
    bool evaluated_;
    bool converged_;

    volScalarField kappa_;
    volVectorField jDiff_;
    autoPtr<surfaceScalarField> kappaFace_;
    autoPtr<surfaceScalarField> jDiffFlux_;

    void updateCurrentFields();
    void selectReferenceCell();
    void setCorrectionBoundaryConditions(volScalarField& deltaPhi) const;

public:
    TypeName("electroneutralGalvanostatic");

    electroneutralGalvanostatic
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

    virtual void evaluateConstraint();
    virtual bool converged() const;
    virtual void correct();
    virtual void report() const;
    virtual const scalarField& chargeResidual() const { return chargeResidual_; }
    virtual const scalarField& chargeScale() const { return chargeScale_; }
    virtual void refreshDiagnosticCurrent() { updateCurrentFields(); }
    virtual const volScalarField& conductivityField() const { return kappa_; }
    virtual const volVectorField& diffusionCurrentField() const { return jDiff_; }
};

} // End namespace Foam

#endif

