#ifndef reversibleKineticReaction_H
#define reversibleKineticReaction_H

#include "fvMesh.H"
#include "dictionary.H"
#include "HashTable.H"
#include "wordList.H"
#include "volFields.H"

namespace Foam
{

class reversibleKineticReaction
{
private:

    word name_;
    const fvMesh& mesh_;

    wordList reactants_;
    wordList products_;

    // Dictionary values. Their physical dimensions are inferred from
    // the number of reactants/products when constructing the rate.
    scalar kf_;
    scalar kb_;

    HashTable<scalar> stoichiometry_;

    void validateSpecies
    (
        const wordList& speciesNames,
        const word& listName,
        const dictionary& dict
    ) const;

    void buildStoichiometry();

    scalar stoichiometricCoefficient
    (
        const word& speciesName
    ) const;

    tmp<volScalarField> zeroSource
    (
        const word& speciesName
    ) const;

public:

    reversibleKineticReaction
    (
        const word& name,
        const fvMesh& mesh,
        const dictionary& dict
    );

    tmp<volScalarField> rate() const;

    tmp<volScalarField> source
    (
        const word& speciesName
    ) const;
};

} // End namespace Foam

#endif

