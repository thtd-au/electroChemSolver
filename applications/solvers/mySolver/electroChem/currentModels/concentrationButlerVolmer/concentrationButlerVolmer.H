#ifndef concentrationButlerVolmer_H
#define concentrationButlerVolmer_H

#include "ButlerVolmer.H"

namespace Foam
{

// Reference-state kinetics with arbitrary configured species on each electrode.
// Controller, Aitken coupling and output infrastructure are shared with BV.
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

