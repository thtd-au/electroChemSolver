#include "multiplicativeMonod.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

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
    forms_(),
    K_(),
    cThreshold_(),
    value_(dict.get<scalar>("value")),
    X_(dict.get<scalar>("X"))
{
    // Reject the old list-based syntax rather than silently ignoring it.
    if
    (
        dict.found("qMax") || dict.found("k1")
     || dict.found("substrates") || dict.found("K")
     || dict.found("C_threshold") || dict.found("c_threshold")
    )
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << ": use rateBasis/value and per-species form, K and c_threshold "
            << "entries in the species dictionary; old top-level kinetic "
            << "entries are not supported" << exit(FatalIOError);
    }

    const word rateBasis(dict.get<word>("rateBasis"));
    if (rateBasis != "maximumSpecificRate" && rateBasis != "firstOrderRateConstant")
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_ << ": unknown rateBasis " << rateBasis
            << "; expected maximumSpecificRate or firstOrderRateConstant"
            << exit(FatalIOError);
    }

    if (!std::isfinite(value_) || value_ < 0 || !std::isfinite(X_) || X_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_ << ": value and X must be finite "
            << "and non-negative" << exit(FatalIOError);
    }

    const dictionary& speciesDict = dict.subDict("species");
    const wordList substrates(speciesDict.toc());
    if (substrates.empty())
    {
        FatalIOErrorInFunction(speciesDict)
            << "Reaction " << name_ << " requires at least one species term"
            << exit(FatalIOError);
    }

    substrateIndices_.setSize(substrates.size());
    forms_.setSize(substrates.size());
    K_.setSize(substrates.size(), 0.0);
    cThreshold_.setSize(substrates.size(), 0.0);
    label nFirstOrder = 0;

    forAll(substrates, i)
    {
        const word& substrate = substrates[i];
        const label speciesI = species.find(substrate);
        if (speciesI < 0)
        {
            FatalIOErrorInFunction(speciesDict)
                << "Reaction " << name_ << " has unknown species "
                << substrate << exit(FatalIOError);
        }
        substrateIndices_[i] = speciesI;

        const dictionary& term = speciesDict.subDict(substrate);
        const word form(term.get<word>("form"));
        if (form == "firstOrder")
        {
            forms_[i] = firstOrder;
            ++nFirstOrder;
            if (term.found("K"))
            {
                FatalIOErrorInFunction(term)
                    << "Reaction " << name_ << ", species " << substrate
                    << ": K is only allowed for form monod"
                    << exit(FatalIOError);
            }
        }
        else if (form == "monod")
        {
            forms_[i] = monod;
            K_[i] = term.get<scalar>("K");
            if (!std::isfinite(K_[i]) || K_[i] <= 0)
            {
                FatalIOErrorInFunction(term)
                    << "Reaction " << name_ << ", species " << substrate
                    << ": K must be finite and positive (mol/m3)"
                    << exit(FatalIOError);
            }
        }
        else
        {
            FatalIOErrorInFunction(term)
                << "Reaction " << name_ << ", species " << substrate
                << ": unknown form " << form << "; expected firstOrder or monod"
                << exit(FatalIOError);
        }

        if (term.found("C_threshold"))
        {
            FatalIOErrorInFunction(term)
                << "Reaction " << name_ << ", species " << substrate
                << ": use lowercase c_threshold instead of C_threshold"
                << exit(FatalIOError);
        }
        cThreshold_[i] = term.lookupOrDefault<scalar>("c_threshold", 0.0);
        if (!std::isfinite(cThreshold_[i]) || cThreshold_[i] < 0)
        {
            FatalIOErrorInFunction(term)
                << "Reaction " << name_ << ", species " << substrate
                << ": c_threshold must be finite and non-negative (mol/m3)"
                << exit(FatalIOError);
        }
    }

    // X*value has units mol/(m3 s) for qMax, or 1/s for k1.
    // Hence qMax permits no dimensional factors, and k1 requires exactly one.
    const label expectedFirstOrder = rateBasis == "firstOrderRateConstant" ? 1 : 0;
    if (nFirstOrder != expectedFirstOrder)
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_ << ": rateBasis " << rateBasis
            << " requires " << expectedFirstOrder << " firstOrder species term(s), "
            << "but found " << nFirstOrder
            << ". All remaining species terms must use form monod."
            << exit(FatalIOError);
    }
}

scalar multiplicativeMonod::kineticRate
(
    const label,
    const scalarField& concentrations,
    const scalarField&
) const
{
    scalar rate = value_*X_;

    forAll(substrateIndices_, i)
    {
        const scalar ci = concentrations[substrateIndices_[i]];
        if (ci <= cThreshold_[i])
        {
            return 0.0;
        }

        const scalar available = ci - cThreshold_[i];
        if (forms_[i] == firstOrder)
        {
            rate *= available;
        }
        else
        {
            // The denominator uses the unshifted concentration.
            rate *= available/(ci + K_[i]);
        }
    }

    return rate;
}

} // End namespace Foam

