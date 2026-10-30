#ifndef reversibleKineticReaction_H
#define reversibleKineticReaction_H

#include "reactionManager.H"

namespace Foam
{

class reversibleKineticReaction : public reactionModel
{
    scalar kf_;
    scalar kb_;

public:
    TypeName("reversibleKineticReaction");

    reversibleKineticReaction
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

