#include "reactionManager.H"
#include "IOobject.H"
#include "dimensionedScalar.H"

namespace Foam
{

reactionManager::reactionManager
(
    fvMesh& mesh,
    const dictionary& reactionsDict
)
:
    mesh_(mesh),
    reactionsDict_(reactionsDict),
    kineticReactions_(),
    equilibriumReactions_()
{
    readKineticReactions();
    readEquilibriumReactions();

    Info<< "Reaction manager loaded:" << nl
        << "    kinetic reactions: "
        << kineticReactions_.size() << nl
        << "    equilibrium reactions: "
        << equilibriumReactions_.size() << nl;
}


void reactionManager::readKineticReactions()
{
    if (!reactionsDict_.found("kineticReactions"))
    {
        Info<< "No kineticReactions dictionary found." << nl;

        kineticReactions_.setSize(0);
        return;
    }

    const dictionary& kineticDict =
        reactionsDict_.subDict("kineticReactions");

    const wordList reactionNames =
        kineticDict.toc();

    kineticReactions_.setSize(reactionNames.size());

    forAll(reactionNames, reactioni)
    {
        const word& reactionName =
            reactionNames[reactioni];

        if (!kineticDict.isDict(reactionName))
        {
            FatalIOErrorInFunction(kineticDict)
                << "Entry " << reactionName
                << " in kineticReactions is not a dictionary."
                << exit(FatalIOError);
        }

        const dictionary& reactionDict =
            kineticDict.subDict(reactionName);

        const word modelType
        (
            reactionDict.lookup("type")
        );

        if (modelType != "multiplicativeMonod")
        {
            FatalIOErrorInFunction(reactionDict)
                << "Unsupported kinetic reaction type "
                << modelType << nl
                << "Currently supported type: multiplicativeMonod"
                << exit(FatalIOError);
        }

        kineticReactions_.set
        (
            reactioni,
            new multiplicativeMonod
            (
                reactionName,
                mesh_,
                reactionDict
            )
        );
    }
}


void reactionManager::readEquilibriumReactions()
{
    if (!reactionsDict_.found("equilibriumReactions"))
    {
        Info<< "No equilibriumReactions dictionary found." << nl;

        equilibriumReactions_.setSize(0);
        return;
    }

    const dictionary& equilibriumDict =
        reactionsDict_.subDict("equilibriumReactions");

    const wordList reactionNames =
        equilibriumDict.toc();

    equilibriumReactions_.setSize(reactionNames.size());

    forAll(reactionNames, reactioni)
    {
        const word& reactionName =
            reactionNames[reactioni];

        if (!equilibriumDict.isDict(reactionName))
        {
            FatalIOErrorInFunction(equilibriumDict)
                << "Entry " << reactionName
                << " in equilibriumReactions is not a dictionary."
                << exit(FatalIOError);
        }

        const dictionary& reactionDict =
            equilibriumDict.subDict(reactionName);

        const word modelType
        (
            reactionDict.lookup("type")
        );

        if (modelType != "waterEquilibrium")
        {
            FatalIOErrorInFunction(reactionDict)
                << "Unsupported equilibrium reaction type "
                << modelType << nl
                << "Currently supported type: waterEquilibrium"
                << exit(FatalIOError);
        }

        equilibriumReactions_.set
        (
            reactioni,
            new waterEquilibrium
            (
                reactionName,
                mesh_,
                reactionDict
            )
        );
    }
}


tmp<volScalarField> reactionManager::zeroSource
(
    const word& speciesName
) const
{
    const word concentrationFieldName("c_" + speciesName);

    if (!mesh_.foundObject<volScalarField>(concentrationFieldName))
    {
        FatalErrorInFunction
            << "Cannot construct reaction source for species "
            << speciesName << nl
            << "Expected concentration field "
            << concentrationFieldName
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
                "kineticSource_" + speciesName,
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


tmp<volScalarField> reactionManager::kineticSource
(
    const word& speciesName
) const
{
    tmp<volScalarField> tCombinedSource =
        zeroSource(speciesName);

    volScalarField& combinedSource =
        tCombinedSource.ref();

    forAll(kineticReactions_, reactioni)
    {
        combinedSource +=
            kineticReactions_[reactioni].source(speciesName);
    }

    return tCombinedSource;
}


void reactionManager::applyEquilibriumReactions()
{
    forAll(equilibriumReactions_, reactioni)
    {
        equilibriumReactions_[reactioni].apply();
    }
}

} // End namespace Foam

