#include "fastEquilibriumReaction.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(fastEquilibriumReaction, 0);
addToRunTimeSelectionTable(reactionModel, fastEquilibriumReaction, dictionary);

fastEquilibriumReaction::fastEquilibriumReaction
(
    const word& name,
    const dictionary& dict,
    const wordList&
)
:
    reactionModel(name),
    K_(dict.get<scalar>("K"))
{
    if (K_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_ << ": K must be positive"
            << exit(FatalIOError);
    }
}

void fastEquilibriumReaction::equilibriumResidualAndJacobian
(
    const label,
    const scalarField& c,
    const scalarField& nu,
    scalar& residual,
    scalarField& dFdc,
    scalar& productTerm,
    scalar& reactantTerm
) const
{
    productTerm = 1.0;
    reactantTerm = 1.0;

    forAll(nu, i)
    {
        if (nu[i] > 0)
        {
            productTerm *= std::pow(c[i], nu[i]);
        }
        else if (nu[i] < 0)
        {
            reactantTerm *= std::pow(c[i], -nu[i]);
        }
    }

    residual = productTerm - K_*reactantTerm;
    dFdc.setSize(c.size(), 0.0);

    // Product-form derivatives: intentionally avoids P/c and R/c.
    forAll(nu, i)
    {
        if (nu[i] > 0)
        {
            scalar derivative = nu[i]*std::pow(c[i], nu[i] - 1.0);
            forAll(nu, j)
            {
                if (j != i && nu[j] > 0)
                {
                    derivative *= std::pow(c[j], nu[j]);
                }
            }
            dFdc[i] = derivative;
        }
        else if (nu[i] < 0)
        {
            const scalar exponent = -nu[i];
            scalar derivative = -K_*exponent*std::pow(c[i], exponent - 1.0);
            forAll(nu, j)
            {
                if (j != i && nu[j] < 0)
                {
                    derivative *= std::pow(c[j], -nu[j]);
                }
            }
            dFdc[i] = derivative;
        }
    }
}

} // End namespace Foam

