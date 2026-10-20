#ifndef concentrationButlerVolmer_H
#define concentrationButlerVolmer_H

#include "ButlerVolmer.H"

namespace Foam
{

// Reference-state kinetics with arbitrary configured species on each electrode.
// Only kinetics live here. The galvanostatic controller, current relaxation,
// exact residual checks and diagnostics are inherited from ButlerVolmer.
// This baseline has no transport-response prediction or Taylor correction.
class concentrationButlerVolmer : public ButlerVolmer
{
    List<wordList> dependencyNames_;
    List<scalarList> referenceConcentrations_, anodicOrders_, cathodicOrders_;

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
};

} // namespace Foam
#endif

