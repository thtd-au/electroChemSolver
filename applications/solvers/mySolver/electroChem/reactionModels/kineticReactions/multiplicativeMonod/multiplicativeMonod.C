#include "multiplicativeMonod.H"
#include "IOobject.H"
#include "dimensionedScalar.H"

namespace Foam
{

multiplicativeMonod::multiplicativeMonod
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    stoichiometry_(),
    substrates_(dict.lookup("substrates")),
    K_(dict.lookup("K")),
    qMax_(readScalar(dict.lookup("qMax"))),
    X_(readScalar(dict.lookup("X")))
{
    readStoichiometry(dict);

    if (substrates_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << " must contain at least one substrate."
            << exit(FatalIOError);
    }

    if (substrates_.size() != K_.size())
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << ": substrates and K must have the same size." << nl
            << "substrates = " << substrates_ << nl
            << "K = " << K_ << nl
            << exit(FatalIOError);
    }

    if (qMax_ < 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << ": qMax must be non-negative."
            << exit(FatalIOError);
    }

    if (X_ < 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << ": X must be non-negative."
            << exit(FatalIOError);
    }

    forAll(substrates_, substratei)
    {
        const word& substrateName = substrates_[substratei];
        const word substrateFieldName("c_" + substrateName);

        if (!mesh_.foundObject<volScalarField>(substrateFieldName))
        {
            FatalIOErrorInFunction(dict)
                << "Reaction " << name_
                << " refers to missing substrate species "
                << substrateName << nl
                << "Expected concentration field "
                << substrateFieldName
                << exit(FatalIOError);
        }

        if (K_[substratei] <= 0.0)
        {
            FatalIOErrorInFunction(dict)
                << "Reaction " << name_
                << ": all K values must be positive."
                << exit(FatalIOError);
        }
    }
}


void multiplicativeMonod::readStoichiometry
(
    const dictionary& dict
)
{
    if (!dict.found("stoichiometry"))
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << " requires a stoichiometry dictionary."
            << exit(FatalIOError);
    }

    const dictionary& stoichDict = dict.subDict("stoichiometry");

    forAllConstIter(dictionary, stoichDict, iter)
    {
        if (!iter().isStream())
        {
            continue;
        }

        const word speciesName(iter().keyword());
        const scalar coefficient =
            readScalar(stoichDict.lookup(speciesName));

        stoichiometry_.insert(speciesName, coefficient);
    }
}


scalar multiplicativeMonod::stoichiometricCoefficient
(
    const word& speciesName
) const
{
    const HashTable<scalar>::const_iterator iter =
        stoichiometry_.find(speciesName);

    if (iter == stoichiometry_.end())
    {
        return 0.0;
    }

    return iter();
}


tmp<volScalarField> multiplicativeMonod::zeroSource
(
    const word& speciesName
) const
{
    const word concentrationFieldName("c_" + speciesName);

    if (!mesh_.foundObject<volScalarField>(concentrationFieldName))
    {
        FatalErrorInFunction
            << "No concentration field named "
            << concentrationFieldName
            << " exists for species " << speciesName
            << exit(FatalError);
    }

    const volScalarField& concentration =
        mesh_.lookupObject<volScalarField>(concentrationFieldName);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                name_ + "_source_" + speciesName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                concentration.dimensions()/dimTime,
                0.0
            )
        )
    );
}


tmp<volScalarField> multiplicativeMonod::rate() const
{
    const volScalarField& firstSubstrate =
        mesh_.lookupObject<volScalarField>
        (
            "c_" + substrates_[0]
        );

    tmp<volScalarField> tRate
    (
        new volScalarField
        (
            IOobject
            (
                name_ + "_rate",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "rate",
                firstSubstrate.dimensions()/dimTime,
                qMax_*X_
            )
        )
    );

    volScalarField& reactionRate = tRate.ref();

    forAll(substrates_, substratei)
    {
        const volScalarField& concentration =
            mesh_.lookupObject<volScalarField>
            (
                "c_" + substrates_[substratei]
            );

        const dimensionedScalar Ki
        (
            "K_" + substrates_[substratei],
            concentration.dimensions(),
            K_[substratei]
        );

        reactionRate *= concentration/(Ki + concentration);
    }

    return tRate;
}


tmp<volScalarField> multiplicativeMonod::source
(
    const word& speciesName
) const
{
    const scalar nu = stoichiometricCoefficient(speciesName);

    if (mag(nu) < SMALL)
    {
        return zeroSource(speciesName);
    }

    return nu*rate();
}

} // End namespace Foam
