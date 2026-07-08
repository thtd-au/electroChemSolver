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

const bool hasChargeNumber = reactions.found("chargeNumber");

dimensionedScalar T
(
    "T",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    reactions.getOrDefault<scalar>("temperature", 298.15)
);

dimensionedScalar Rgas
(
    "R",
    dimensionSet(1, 2, -2, -1, -1, 0, 0),
    8.31446261815324
);

dimensionedScalar Fdim
(
    "FaradayConstant",
    dimensionSet(0, 0, 1, 0, -1, 1, 0),    
    reactions.getOrDefault<scalar>("FaradayConstant", 96485.3329)
);

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

        scalar zi = 0.0;

        if (hasChargeNumber)
        {
            const dictionary& zDict = reactions.subDict("chargeNumber");
            zi = zDict.getOrDefault<scalar>(species[i], 0.0);
        }

        surfaceScalarField migrationFlux
        (
            IOobject
            (
                "migrationFlux_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(Di*zi*Fdim/(Rgas*T)*Yi)
        *fvc::snGrad(phiE)
        *mesh.magSf()
        );

        // Transport diagnostics 
        const bool doTransportDiagnostics =
        (
            transportDiagnosticsTotal
        && transportDiagnosticsSpecies.size()
        && transportDiagnosticsSpecies.find(species[i]) >= 0
        );

        if (doTransportDiagnostics)
        {
            volScalarField advTerm
            (
                IOobject
                (
                    "adv_" + Yi.name(),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    transportDiagnosticsWriteFields 
                    ? IOobject::NO_WRITE 
                    : IOobject::NO_WRITE
                ),
                -fvc::div(phi, Yi)
            );

            volScalarField diffTerm
            (
                IOobject
                (
                    "diff_" + Yi.name(), 
                    runTime.timeName(), 
                    mesh,
                    IOobject::NO_READ,
                    transportDiagnosticsWriteFields 
                    ? IOobject::NO_WRITE 
                    : IOobject::NO_WRITE
                ),
                fvc::laplacian(Di, Yi)
            );

            volScalarField migTerm
            (
                IOobject
                (
                    "mig_" + Yi.name(), 
                    runTime.timeName(), 
                    mesh,
                    IOobject::NO_READ,
                    transportDiagnosticsWriteFields 
                    ? IOobject::NO_WRITE 
                    : IOobject::NO_WRITE
                ),
                fvc::div(migrationFlux)
            );

            volScalarField rxnTerm
            (
                IOobject
                (
                    "rxn_" + Yi.name(), 
                    runTime.timeName(), 
                    mesh,
                    IOobject::NO_READ,
                    transportDiagnosticsWriteFields 
                    ? IOobject::NO_WRITE 
                    : IOobject::NO_WRITE
                ),
                RY[i]
            );

            if (transportDiagnosticsWriteFields)
            {
                advYDiag[i] = advTerm;
                diffYDiag[i] = diffTerm;
                migYDiag[i] = migTerm;
                rxnYDiag[i] = rxnTerm;

                advYDiag[i].correctBoundaryConditions();
                diffYDiag[i].correctBoundaryConditions();
                migYDiag[i].correctBoundaryConditions();
                rxnYDiag[i].correctBoundaryConditions();
            }

            const scalar advMax = max(mag(advTerm)).value();
            const scalar diffMax = max(mag(diffTerm)).value();
            const scalar migMax = max(mag(migTerm)).value();
            const scalar rxnMax = max(mag(rxnTerm)).value();

            const scalar advInt = gSum(mag(advTerm.primitiveField())*mesh.V());
            const scalar diffInt = gSum(mag(diffTerm.primitiveField())*mesh.V());
            const scalar migInt = gSum(mag(migTerm.primitiveField())*mesh.V());
            const scalar rxnInt = gSum(mag(rxnTerm.primitiveField())*mesh.V());

            Info<< "Transport diagnostics for " << species[i] << nl
                << "  max|adv diff mig rxn| [1/s] = "
                << advMax << " "
                << diffMax << " "
                << migMax << " "
                << rxnMax << nl
                << "  int|adv diff mig rxn| [m3/s] = "
                << advInt << " "
                << diffInt << " "
                << migInt << " "
                << rxnInt << nl;

            if (transportDiagnosticsFilePtr.valid() && Pstream::master())
            {
                OFstream& os = transportDiagnosticsFilePtr();

                os  << runTime.value() << ' '
                    << species[i] << ' '
                    << advMax << ' '
                    << diffMax << ' '
                    << migMax << ' '
                    << rxnMax << ' '
                    << advInt << ' '
                    << diffInt << ' '
                    << migInt << ' '
                    << rxnInt << nl;
            }
        }

        fvScalarMatrix YiEqn
        (
            fvm::ddt(Yi)
          + fvm::div(phi, Yi)
         ==
            fvm::laplacian(Di, Yi)
          + fvc::div(migrationFlux)
          + RY[i]
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

