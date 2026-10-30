#ifndef constantCurrent_H
#define constantCurrent_H

#include "currentManager.H"
#include "HashTable.H"

namespace Foam
{

class constantCurrent : public currentModel
{
    scalar appliedCurrentDensity_;

    // Read metadata from the live patch's generic write interface once per
    // electrode name; do not depend on a concrete custom-BC C++ class.
    HashTable<scalar> electronStoichCoeffs_;
    scalar electronStoichCoeff(const fvPatchScalarField& patchField);

public:
    TypeName("constantCurrent");

    constantCurrent
    (
        fvMesh& mesh,
        const dictionary& currentDict,
        const volScalarField& phiE
    );

    virtual void correct();
};

} // End namespace Foam

#endif

