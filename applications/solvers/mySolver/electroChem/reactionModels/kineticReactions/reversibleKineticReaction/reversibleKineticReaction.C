#include "reversibleKineticReaction.H"
#include "IOobject.H"
#include "dimensionedScalar.H"

namespace Foam
{

reversibleKineticReaction::reversibleKineticReaction
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    reactants_(dict.lookup("reactants")),
    products_(dict.lookup("products")),
    kf_(dict.get<scalar>("kf")),
    kb_(dict.get<scalar>("kb")),
    stoichiometry_()
{
    if (reactants_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << " must contain at least one reactant."
            << exit(FatalIOError);
    }

    if (kf_ < 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << ": kf must be non-negative."
            << exit(FatalIOError);
    }

    if (kb_ < 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << ": kb must be non-negative."
            << exit(FatalIOError);
    }

    validateSpecies(reactants_, "reactants", dict);
    validateSpecies(products_, "products", dict);

    buildStoichiometry();

    Info<< "    reversible kinetic reaction " << name_ << ":" << nl
    << "        reactants: " << reactants_ << nl
    << "        products:  " << products_ << nl;

    if (products_.empty())
    {
        Info<< "        backward term is zero-order: kb" << nl;
    }
}


void reversibleKineticReaction::validateSpecies
(
    const wordList& speciesNames,
    const word& listName,
    const dictionary& dict
) const
{
    forAll(speciesNames, i)
    {
        const word& speciesName = speciesNames[i];
        const word fieldName("c_" + speciesName);

        if (!mesh_.foundObject<volScalarField>(fieldName))
        {
            FatalIOErrorInFunction(dict)
                << "Reaction " << name_
                << " refers to unknown species "
                << speciesName << " in " << listName << "." << nl
                << "Expected concentration field "
                << fieldName
                << exit(FatalIOError);
        }
    }
}


void reversibleKineticReaction::buildStoichiometry()
{
    // Reactants have negative stoichiometric coefficients.
    forAll(reactants_, i)
    {
        const word& speciesName = reactants_[i];

        HashTable<scalar>::iterator iter =
            stoichiometry_.find(speciesName);

        if (iter == stoichiometry_.end())
        {
            stoichiometry_.insert(speciesName, -1.0);
        }
        else
        {
            iter() -= 1.0;
        }
    }

    // Products have positive stoichiometric coefficients.
    forAll(products_, i)
    {
        const word& speciesName = products_[i];

        HashTable<scalar>::iterator iter =
            stoichiometry_.find(speciesName);

        if (iter == stoichiometry_.end())
        {
            stoichiometry_.insert(speciesName, 1.0);
        }
        else
        {
            iter() += 1.0;
        }
    }
}


scalar reversibleKineticReaction::stoichiometricCoefficient
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


tmp<volScalarField> reversibleKineticReaction::zeroSource
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
        mesh_.lookupObject<volScalarField>
        (
            concentrationFieldName
        );

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


tmp<volScalarField> reversibleKineticReaction::rate() const
{
    const volScalarField& referenceConcentration =
        mesh_.lookupObject<volScalarField>
        (
            "c_" + reactants_[0]
        );

    const dimensionSet concentrationDimensions =
        referenceConcentration.dimensions();

    const dimensionSet rateDimensions =
        concentrationDimensions/dimTime;

    // /*
    //  * kf dimensions:
    //  *
    //  *     [kf] [C]^Nr = [C]/[t]
    //  *
    //  * kb dimensions:
    //  *
    //  *     [kb] [C]^Np = [C]/[t]
    //  */
    
    dimensionSet kfDimensions(rateDimensions);
    dimensionSet kbDimensions(rateDimensions);

    forAll(reactants_, i)
    {
        kfDimensions /= concentrationDimensions;
    }

    forAll(products_, i)
    {
        kbDimensions /= concentrationDimensions;
    }

    const dimensionedScalar kf
    (
        "kf_" + name_,
        kfDimensions,
        kf_
    );

    const dimensionedScalar kb
    (
        "kb_" + name_,
        kbDimensions,
        kb_
    );

    tmp<volScalarField> tForwardRate
    (
        new volScalarField
        (
            IOobject
            (
                name_ + "_forwardRate",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            kf
        )
    );

    tmp<volScalarField> tBackwardRate
    (
        new volScalarField
        (
            IOobject
            (
                name_ + "_backwardRate",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            kb
        )
    );

    volScalarField& forwardRate = tForwardRate.ref();
    volScalarField& backwardRate = tBackwardRate.ref();

    forAll(reactants_, i)
    {
        const volScalarField& concentration =
            mesh_.lookupObject<volScalarField>
            (
                "c_" + reactants_[i]
            );

        if (concentration.dimensions() != concentrationDimensions)
        {
            FatalErrorInFunction
                << "Reaction " << name_
                << ": concentration field c_" << reactants_[i]
                << " has dimensions " << concentration.dimensions()
                << ", expected " << concentrationDimensions
                << exit(FatalError);
        }

        forwardRate *= max
        (
            concentration,
            dimensionedScalar
            (
                "zero",
                concentrationDimensions,
                0.0
            )
        );
    }

    forAll(products_, i)
    {
        const volScalarField& concentration =
            mesh_.lookupObject<volScalarField>
            (
                "c_" + products_[i]
            );

        if (concentration.dimensions() != concentrationDimensions)
        {
            FatalErrorInFunction
                << "Reaction " << name_
                << ": concentration field c_" << products_[i]
                << " has dimensions " << concentration.dimensions()
                << ", expected " << concentrationDimensions
                << exit(FatalError);
        }

        backwardRate *= max
        (
            concentration,
            dimensionedScalar
            (
                "zero",
                concentrationDimensions,
                0.0
            )
        );
    }

    return tForwardRate - tBackwardRate;
}


tmp<volScalarField> reversibleKineticReaction::source
(
    const word& speciesName
) const
{
    const scalar nu =
        stoichiometricCoefficient(speciesName);

    if (mag(nu) < SMALL)
    {
        return zeroSource(speciesName);
    }

    return nu*rate();
}

} // End namespace Foam

