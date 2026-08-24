#ifndef waterEquilibrium_H
#define waterEquilibrium_H

#include "fvMesh.H"
#include "dictionary.H"
#include "Switch.H"
#include "volFields.H"

namespace Foam
{

class waterEquilibrium
{
private:

    word name_;
    fvMesh& mesh_;

    word HplusName_;
    word OHminusName_;

    scalar Kw_;
    Switch diagA_;

    volScalarField& Hplus_;
    volScalarField& OHminus_;

public:

    waterEquilibrium
    (
        const word& name,
        fvMesh& mesh,
        const dictionary& dict
    );

    void apply();
};

} // End namespace Foam

#endif

