// cEqn.H
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

// Reset concentration source terms
forAll(Rc, i)
{
    Rc[i] = dimensionedScalar("zero", dimMoles/dimVolume/dimTime, 0.0);
}


// Assemble reaction source terms in mol/m3/s
if (reactions.found("reactions"))
{
    const dictionary& rxnDicts = reactions.subDict("reactions");

    forAllConstIter(dictionary, rxnDicts, iter)
    {
        if (!iter().isDict()) continue;

        const dictionary& rd = iter().dict();
        const word type(rd.get<word>("type"));

        if (type == "multiplicativeMonod")
        {
            const dictionary& nuDict = rd.subDict("stoichiometry");

            wordList substrates(rd.lookup("substrates"));
            scalarList K(rd.lookup("K"));

            dimensionedScalar qMax
            (
                "qMax",
                dimMoles/dimTime,
                rd.get<scalar>("qMax")
            );

            dimensionedScalar X
            (
                "X",
                dimless/dimVolume,
                rd.get<scalar>("X")
            );

            volScalarField rate
            (
                IOobject
                (
                    "r_" + iter().keyword(),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dimMoles/dimVolume/dimTime,
                    0.0
                )
            );

            rate = X*qMax;

            forAll(substrates, si)
            {
                label subI = species.find(substrates[si]);

                if (subI < 0)
                {
                    FatalIOErrorInFunction(rd)
                        << "Unknown substrate " << substrates[si]
                        << exit(FatalIOError);
                }

                dimensionedScalar Ki
                (
                    "K_" + substrates[si],
                    dimMoles/dimVolume,
                    K[si]
                );

                rate *= c[subI]/(c[subI] + Ki);
            }

            forAllConstIter(dictionary, nuDict, nIter)
            {
                const word specieName = nIter().keyword();
                const label speciei = species.find(specieName);

                if (speciei < 0)
                {
                    FatalIOErrorInFunction(rd)
                        << "Unknown stoichiometric species " << specieName
                        << exit(FatalIOError);
                }

                const scalar nu = readScalar(nIter().stream());

                Rc[speciei] += nu*rate;
            }
        }
        else if (type == "waterEquilibrium")
        {
            // Leave equilibrium rework for later.
        }
        else
        {
            FatalIOErrorInFunction(rd)
                << "Unknown reaction type " << type
                << exit(FatalIOError);
        }
    }
}

// Solve concentration transport
forAll(c, i)
{
    if (i != inertIndex)
    {
        volScalarField& ci = c[i];

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
            fvc::interpolate(Di*zi*Fdim/(Rgas*T)*ci)
          * fvc::snGrad(phiE)
          * mesh.magSf()
        );

        fvScalarMatrix ciEqn
        (
            fvm::ddt(ci)
          + fvm::div(phi, ci)
         ==
            fvm::laplacian(Di, ci)
          + fvc::div(migrationFlux)
          + Rc[i]
          + fvOptions(ci)
        );

        ciEqn.relax();
        fvOptions.constrain(ciEqn);
        ciEqn.solve(ci.name());
        fvOptions.correct(ci);

        ci.max(0.0);
        ci.correctBoundaryConditions();
    }
}

