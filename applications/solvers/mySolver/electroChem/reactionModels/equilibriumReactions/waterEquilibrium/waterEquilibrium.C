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
    diagA_(false),
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

    const dictionary& controlDict = mesh_.time().controlDict();

    if (controlDict.found("diagnostics"))
    {
        const dictionary& diagnosticsDict =
            controlDict.subDict("diagnostics");

        if (diagnosticsDict.found("waterEquilibrium"))
        {
            const dictionary& waterDiagDict =
                diagnosticsDict.subDict("waterEquilibrium");

            diagA_ =
                waterDiagDict.getOrDefault<Switch>("diagA", false);
        }
    }
}


void waterEquilibrium::apply()
{
    scalarField& cH = Hplus_.primitiveFieldRef();
    scalarField& cOH = OHminus_.primitiveFieldRef();

    // Capture the transported ion difference before any concentration clipping.
    // This allows a local over-consumption of H+ or OH- to be supplied by the
    // instantaneous water equilibrium while conserving A = [H+] - [OH-].
    scalarField ABefore(cH - cOH);

    forAll(cH, celli)
    {
        // const scalar A = cH[celli] - cOH[celli];
        const scalar A = ABefore[celli];
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

    if (diagA_)
    {
        const scalarField AAfter(cH - cOH);
        const scalarField AResidual(AAfter - ABefore);

        const scalar totalVolume = gSum(mesh_.V());

        const scalar ABeforeMin = gMin(ABefore);
        const scalar ABeforeMean =
            gSum(ABefore*mesh_.V())/(totalVolume + VSMALL);
        const scalar ABeforeMax = gMax(ABefore);
        const scalar ABeforeIntegral = gSum(ABefore*mesh_.V());

        const scalar AAfterMin = gMin(AAfter);
        const scalar AAfterMean =
            gSum(AAfter*mesh_.V())/(totalVolume + VSMALL);
        const scalar AAfterMax = gMax(AAfter);
        const scalar AAfterIntegral = gSum(AAfter*mesh_.V());

        const scalar maxAResidual = gMax(mag(AResidual));
        const scalar integralAResidual = gSum(AResidual*mesh_.V());

        Info<< "Water-equilibrium A diagnostics (" << name_ << "):" << nl
            << "  before min/mean/max = "
            << ABeforeMin << " / "
            << ABeforeMean << " / "
            << ABeforeMax << nl
            << "  after  min/mean/max = "
            << AAfterMin << " / "
            << AAfterMean << " / "
            << AAfterMax << nl
            << "  integral before/after = "
            << ABeforeIntegral << " / "
            << AAfterIntegral << nl
            << "  max|A_after-A_before| = "
            << maxAResidual << nl
            << "  integral(A_after-A_before) = "
            << integralAResidual << nl
            << endl;
    }

    // Internal cells are projected; configured patch conditions remain authoritative.
    Hplus_.correctBoundaryConditions();
    OHminus_.correctBoundaryConditions();
}

} // End namespace Foam

