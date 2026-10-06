#ifndef fastEquilibriumReaction_H
#define fastEquilibriumReaction_H

#include "reactionManager.H"

namespace Foam
{

class fastEquilibriumReaction : public reactionModel
{
    scalar K_;

public:
    TypeName("fastEquilibriumReaction");

    fastEquilibriumReaction
    (
        const word& name,
        const dictionary& dict,
        const wordList& species
    );

    categoryType category() const override { return equilibrium; }

    void equilibriumResidualAndJacobian
    (
        const label celli,
        const scalarField& concentrations,
        const scalarField& nu,
        scalar& residual,
        scalarField& dFdc,
        scalar& productTerm,
        scalar& reactantTerm
    ) const override;
};

} // End namespace Foam

#endif

