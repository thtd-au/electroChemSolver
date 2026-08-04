#ifndef multiplicativeMonod_H
#define multiplicativeMonod_H

#include "fvMesh.H"
#include "dictionary.H"
#include "HashTable.H"
#include "wordList.H"
#include "scalarList.H"
#include "volFields.H"

namespace Foam
{

class multiplicativeMonod
{
private:

    word name_;
    const fvMesh& mesh_;

    HashTable<scalar> stoichiometry_;

    wordList substrates_;
    scalarList K_;

    scalar qMax_;
    scalar X_;

    void readStoichiometry(const dictionary& dict);

    scalar stoichiometricCoefficient
    (
        const word& speciesName
    ) const;

    tmp<volScalarField> zeroSource
    (
        const word& speciesName
    ) const;

public:

    multiplicativeMonod
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
