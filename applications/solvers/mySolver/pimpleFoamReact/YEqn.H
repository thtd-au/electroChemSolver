volScalarField Yt
(
    IOobject
    (
        "Yt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// diagnostics 
volScalarField YTotal
(
    IOobject
    (
        "YTotal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

label inertIndex = species.find(inertSpecie);

const dictionary& diffusivityDict = reactions.subDict("diffusivity");

// Reset reaction source fields
forAll(Y, i)
{
    RY[i] = dimensionedScalar("zero", dimless/dimTime, 0.0);
    RYWaterEq[i] = dimensionedScalar("zero", dimless/dimTime, 0.0); 
}

// Fill RY from constant/reactions
#include "updateReactionSources.H"

forAll(Y, i)
{
    if (i != inertIndex)
    {
        volScalarField& Yi = Y[i];

        dimensionedScalar Di
        (
            "D_" + species[i],
            dimArea/dimTime,
            diffusivityDict.get<scalar>(species[i])
        );

        fvScalarMatrix YiEqn
        (
            fvm::ddt(Yi)
          + fvm::div(phi, Yi)
          - fvm::laplacian(Di, Yi)
         ==
            RY[i]
          + fvOptions(Yi)
        );

        YiEqn.relax();
        fvOptions.constrain(YiEqn);
        YiEqn.solve(Yi.name());
        fvOptions.correct(Yi);

        Yi.max(0.0);
        // Yi.min(1.0); // can hide conservation errors and is therefore not active 
        Yi.correctBoundaryConditions();

        // Yt += Yi;
    }
}

// Fast water dissociation equilibrium is an algebraic correction applied
// after transport and before the inert-species closure.
#include "waterEquilibrium.H"

// Accumulate non-inert species after all algebraic corrections.
forAll(Y, i)
{
    if (i != inertIndex)
    {
        Yt += Y[i];
    }
}

if (inertIndex >= 0)
{
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
    Y[inertIndex].correctBoundaryConditions();
}

if (printYFields)
{
    forAll(Y, i)
    {
        Info<< Y[i].name()
            << " min/max = "
            << min(Y[i]).value() << " / "
            << max(Y[i]).value() << endl;
    
            YTotal += Y[i];
    }
}

if (printDiagnostics)
{
    Info<< "sum non-inert mass fractions Yt min/max = "
        << min(Yt).value() << " / "
        << max(Yt).value() << endl;

    if (max(Yt).value() > 1.0)
    {
        WarningInFunction
            << "Sum of non-inert mass fractions exceeds 1. "
            << "Inert species will become negative before clipping."
            << endl;
    }

    Info<< "total mass fraction YTotal min/max = "
        << min(YTotal).value() << " / "
        << max(YTotal).value() << endl;

    Info<< "mass fraction closure error |YTotal - 1| max = "
        << max(mag(YTotal - scalar(1))).value() << endl;
}

