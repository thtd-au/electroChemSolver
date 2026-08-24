#ifndef waterEquilibrium_H
#define waterEquilibrium_H

#include "fvMesh.H"
#include "dictionary.H"
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

    const word& name() const
    {
        return name_;
    }

    const word& HplusName() const
    {
        return HplusName_;
    }

    const word& OHminusName() const
    {
        return OHminusName_;
    }
    
};

} // End namespace Foam

#endif

