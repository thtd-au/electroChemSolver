#ifndef multiplicativeMonod_H
#define multiplicativeMonod_H

#include "reactionManager.H"

namespace Foam
{

class multiplicativeMonod : public reactionModel
{
    labelList substrateIndices_;
    scalarList K_;
    scalarList CThreshold_;
    scalar qMax_;
    scalar X_;

public:
    TypeName("multiplicativeMonod");

    multiplicativeMonod
    (
        const word& name,
        const dictionary& dict,
        const wordList& species
    );

    categoryType category() const override { return kinetic; }

    scalar kineticRate
    (
        const label celli,
        const scalarField& concentrations,
        const scalarField& nu
    ) const override;
};

} // End namespace Foam

#endif

