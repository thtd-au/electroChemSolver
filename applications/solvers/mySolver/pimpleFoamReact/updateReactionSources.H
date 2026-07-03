forAll(RY, i)
{
    RY[i] = dimensionedScalar("zero", dimless/dimTime, 0.0);
}

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
                dimMoles/dimTime,   // adjust if your qMax is not mol/cell/s
                rd.get<scalar>("qMax")
            );

            dimensionedScalar X
            (
                "X",
                dimless/dimVolume,  // cells treated as dimensionless count
                rd.get<scalar>("X")
            );

            volScalarField rate
            (
                IOobject("r_" + iter().keyword(), runTime.timeName(), mesh,
                         IOobject::NO_READ, IOobject::AUTO_WRITE),
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
                const label i = species.find(specieName);

                if (i < 0)
                {
                    FatalIOErrorInFunction(rd)
                        << "Unknown stoichiometric species " << specieName
                        << exit(FatalIOError);
                }

                const scalar nu = readScalar(nIter().stream());

                dimensionedScalar Mi
                (
                    "M_" + specieName,
                    dimMass/dimMoles,
                    molarMassDict.get<scalar>(specieName)
                );

                RY[i] += (Mi/rho)*nu*rate;
            }
        }
        else if (type == "waterEquilibrium")
        {
            // Water equilibrium is an algebraic correction applied after
            // species transport in waterEquilibrium.H.
            // No volumetric source term is assembled here.
        }
        else
        {
            FatalIOErrorInFunction(rd)
                << "Unknown reaction type " << type << exit(FatalIOError);
        }
    }
}

