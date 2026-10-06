#include "reversibleKineticReaction.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(reversibleKineticReaction, 0);
addToRunTimeSelectionTable(reactionModel, reversibleKineticReaction, dictionary);

reversibleKineticReaction::reversibleKineticReaction
(
    const word& name,
    const dictionary& dict,
    const wordList&
)
:
    reactionModel(name),
    kf_(dict.get<scalar>("kf")),
    kb_(dict.get<scalar>("kb"))
{
    if (kf_ < 0 || kb_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << ": kf and kb must be non-negative"
            << exit(FatalIOError);
    }
}

scalar reversibleKineticReaction::kineticRate
(
    const label,
    const scalarField& concentrations,
    const scalarField& nu
) const
{
    scalar forward = kf_;
    scalar backward = kb_;

    forAll(nu, speciesI)
    {
        if (nu[speciesI] < 0)
        {
            forward *= std::pow
            (
                concentrations[speciesI],
                -nu[speciesI]
            );
        }
        else if (nu[speciesI] > 0)
        {
            backward *= std::pow
            (
                concentrations[speciesI],
                nu[speciesI]
            );
        }
    }

    return forward - backward;
}

} // End namespace Foam

