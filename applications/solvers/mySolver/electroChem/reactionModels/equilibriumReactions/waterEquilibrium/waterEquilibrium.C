#include "waterEquilibrium.H"

namespace Foam
{

waterEquilibrium::waterEquilibrium
(
    const word& name,
    fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    HplusName_(dict.lookup("Hplus")),
    OHminusName_(dict.lookup("OHminus")),
    Kw_(readScalar(dict.lookup("Kw"))),
    Hplus_
    (
        mesh_.lookupObjectRef<volScalarField>
        (
            "c_" + HplusName_
        )
    ),
    OHminus_
    (
        mesh_.lookupObjectRef<volScalarField>
        (
            "c_" + OHminusName_
        )
    )
{
    if (HplusName_ == OHminusName_)
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << ": Hplus and OHminus must refer to different species."
            << exit(FatalIOError);
    }

    if (Hplus_.dimensions() != OHminus_.dimensions())
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << ": Hplus and OHminus fields must have identical dimensions."
            << exit(FatalIOError);
    }

    if (Kw_ <= 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "Reaction " << name_
            << ": Kw must be greater than zero."
            << exit(FatalIOError);
    }
}


void waterEquilibrium::apply()
{
    scalarField& cH = Hplus_.primitiveFieldRef();
    scalarField& cOH = OHminus_.primitiveFieldRef();

    forAll(cH, celli)
    {
        const scalar A = cH[celli] - cOH[celli];
        const scalar root = sqrt(sqr(A) + 4.0*Kw_);

        // /*
        //  * Stable evaluation:
        //  * - For A >= 0, Hplus is the larger concentration.
        //  * - For A < 0, OHminus is the larger concentration.
        //  * The smaller concentration is recovered from Kw/c.
        //  */
        if (A >= 0.0)
        {
            const scalar h = 0.5*(A + root);

            cH[celli] = h;
            cOH[celli] = Kw_/h;
        }
        else
        {
            const scalar oh = 0.5*(-A + root);

            cOH[celli] = oh;
            cH[celli] = Kw_/oh;
        }
    }

    // Internal cells are projected; configured patch conditions remain authoritative.
    Hplus_.correctBoundaryConditions();
    OHminus_.correctBoundaryConditions();
}

} // End namespace Foam
