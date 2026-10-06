#include "multiplicativeMonod.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(multiplicativeMonod, 0);
addToRunTimeSelectionTable(reactionModel, multiplicativeMonod, dictionary);

multiplicativeMonod::multiplicativeMonod
(
    const word& name,
    const dictionary& dict,
    const wordList& species
)
:
    reactionModel(name),
    substrateIndices_(),
    K_(dict.lookup("K")),
    qMax_(dict.get<scalar>("qMax")),
    X_(dict.get<scalar>("X"))
{
    const wordList substrates(dict.lookup("substrates"));

    if (substrates.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_ << " requires at least one substrate"
            << exit(FatalIOError);
    }

    if (K_.size() != substrates.size())
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_ << ": substrates and K must have "
            << "the same length"
            << exit(FatalIOError);
    }

    if (qMax_ < 0 || X_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_ << ": qMax and X must be non-negative"
            << exit(FatalIOError);
    }

    substrateIndices_.setSize(substrates.size());

    forAll(substrates, i)
    {
        const label speciesI = species.find(substrates[i]);
        if (speciesI < 0)
        {
            FatalIOErrorInFunction(dict)
                << "Reaction " << name_ << " has unknown substrate "
                << substrates[i] << exit(FatalIOError);
        }
        if (K_[i] <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "Reaction " << name_ << ": K[" << i
                << "] must be positive" << exit(FatalIOError);
        }
        substrateIndices_[i] = speciesI;
    }
}

scalar multiplicativeMonod::kineticRate
(
    const label,
    const scalarField& concentrations,
    const scalarField&
) const
{
    scalar rate = qMax_*X_;

    forAll(substrateIndices_, i)
    {
        const scalar ci = concentrations[substrateIndices_[i]];
        rate *= ci/(K_[i] + ci);
    }

    return rate;
}

} // End namespace Foam

