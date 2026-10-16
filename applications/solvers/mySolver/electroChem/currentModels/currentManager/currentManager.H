#ifndef currentManager_H
#define currentManager_H

#include "fvMesh.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "surfaceFields.H"
#include "volFields.H"

namespace Foam
{

class currentModel
{
protected:
    fvMesh& mesh_;
    const dictionary& currentDict_;
    const volScalarField& phiE_;

    // Registered shared state [A/m2]. Boundary values are j dot n, positive
    // out of the electrolyte. Internal faces are unused and remain zero.
    surfaceScalarField electrodeCurrentDensity_;

public:
    TypeName("currentModel");

    declareRunTimeSelectionTable
    (
        autoPtr, currentModel, dictionary,
        (
            fvMesh& mesh,
            const dictionary& currentDict,
            const volScalarField& phiE
        ),
        (mesh, currentDict, phiE)
    );

    currentModel
    (
        fvMesh& mesh,
        const dictionary& currentDict,
        const volScalarField& phiE
    );
    virtual ~currentModel() = default;

    static autoPtr<currentModel> New
    (
        fvMesh& mesh,
        const dictionary& electroChemistry,
        const volScalarField& phiE
    );

    virtual void correct() = 0;

    const surfaceScalarField& electrodeCurrentDensity() const
    {
        return electrodeCurrentDensity_;
    }
};


class currentManager
{
    autoPtr<currentModel> model_;

public:
    currentManager
    (
        fvMesh& mesh,
        const dictionary& electroChemistry,
        const volScalarField& phiE
    );

    void correct() { model_->correct(); }

    const surfaceScalarField& electrodeCurrentDensity() const
    {
        return model_->electrodeCurrentDensity();
    }
};

} // End namespace Foam

#endif

