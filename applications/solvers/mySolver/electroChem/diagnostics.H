// diagnostics.H
// Unified diagnostics for concentration-primary electroChem solver.

{
    const scalar dt = runTime.deltaTValue();

    // pH is a post-convergence diagnostic derived from the accepted H+ field.
    // The floor protects log10 only and never modifies c_Hplus.
    if (pHActive)
    {
        volScalarField& pH = pHPtr();
        const volScalarField& cH = c[HplusIndex];
        const scalar cStandard = 1000.0; // 1 mol/L = 1000 mol/m3
        const scalar cLogFloor = VSMALL*cStandard;

        forAll(pH.primitiveFieldRef(), celli)
        {
            pH.primitiveFieldRef()[celli] =
                -Foam::log10(max(cH[celli], cLogFloor)/cStandard);
        }

        forAll(pH.boundaryFieldRef(), patchi)
        {
            scalarField& pHp = pH.boundaryFieldRef()[patchi];
            const fvPatchScalarField& cHp = cH.boundaryField()[patchi];

            forAll(pHp, facei)
            {
                pHp[facei] =
                    -Foam::log10(max(cHp[facei], cLogFloor)/cStandard);
            }
        }

        if (printPH)
        {
            Info<< "pH: min=" << gMin(pH.primitiveField())
                << " mean="
                << gSum(pH.primitiveField()*mesh.V())/gSum(mesh.V())
                << " max=" << gMax(pH.primitiveField()) << nl;
        }
    }

    const dictionary& diffusivityDict = reactions.subDict("diffusivity");
    const bool hasChargeNumber = reactions.found("chargeNumber");

    const dictionary* zDictPtr = nullptr;

    if (hasChargeNumber)
    {
        zDictPtr = &reactions.subDict("chargeNumber");
    }

    // Diagnostic log cadence is based on physical simulation time and is
    // independent of controlDict writeInterval/writeControl.
    static bool diagLogScheduleInitialised(false);
    static scalar nextDiagLogTime(0.0);

    if (!diagLogScheduleInitialised)
    {
        nextDiagLogTime =
            runTime.startTime().value() + diagLogInterval;

        diagLogScheduleInitialised = true;
    }

    const scalar diagTimeTolerance =
        10.0*SMALL*max(mag(runTime.value()), scalar(1.0));

    bool diagLogTime(false);

    if (runTime.value() + diagTimeTolerance >= nextDiagLogTime)
    {
        diagLogTime = true;

        // If a variable/large timestep jumps over more than one requested
        // diagnostic time, write once at the current solver time and advance
        // the schedule until the next target lies in the future.
        do
        {
            nextDiagLogTime += diagLogInterval;
        }
        while
        (
            nextDiagLogTime
         <= runTime.value() + diagTimeTolerance
        );
    }

    static bool diagnosticsInitialised(false);

    // Cumulative signed additions to the domain [mol].
    // Positive = net addition to the computational domain.
    static List<scalar> N0;
    static List<scalar> cumulativeAdvection;
    static List<scalar> cumulativeDiffusion;
    static List<scalar> cumulativeMigration;
    static List<scalar> cumulativeReaction;
    static List<scalar> cumulativeKineticReaction;
    static List<scalar> cumulativeEquilibriumReaction;

    static autoPtr<OFstream> transportLogPtr;
    static autoPtr<OFstream> conservationLogPtr;
    static autoPtr<OFstream> saturationLogPtr;
    static autoPtr<OFstream> DebyeLogPtr;

    // Patch roles for conservation:
    // 0 = ignored, 1 = flow patch, 2 = electrode patch.
    static List<label> conservationPatchRole;

    if (!diagnosticsInitialised)
    {
        N0.setSize(c.size(), 0.0);
        cumulativeAdvection.setSize(c.size(), 0.0);
        cumulativeDiffusion.setSize(c.size(), 0.0);
        cumulativeMigration.setSize(c.size(), 0.0);
        cumulativeReaction.setSize(c.size(), 0.0);
        cumulativeKineticReaction.setSize(c.size(), 0.0);
        cumulativeEquilibriumReaction.setSize(c.size(), 0.0);

        conservationPatchRole.setSize(mesh.boundary().size(), 0);

        forAll(c, i)
        {
            N0[i] = initialSpeciesInventory[i];
        }

        forAll(diagFlowPatches, i)
        {
            const label patchi =
                mesh.boundaryMesh().findPatchID(diagFlowPatches[i]);

            if (patchi < 0)
            {
                FatalErrorInFunction
                    << "Unknown diagnostics/conservation flow patch "
                    << diagFlowPatches[i]
                    << exit(FatalError);
            }

            if (conservationPatchRole[patchi] != 0)
            {
                FatalErrorInFunction
                    << "Patch " << diagFlowPatches[i]
                    << " is listed more than once in diagnostics/conservation"
                    << exit(FatalError);
            }

            conservationPatchRole[patchi] = 1;
        }

        forAll(diagElectrodePatches, i)
        {
            const label patchi =
                mesh.boundaryMesh().findPatchID(diagElectrodePatches[i]);

            if (patchi < 0)
            {
                FatalErrorInFunction
                    << "Unknown diagnostics/conservation electrode patch "
                    << diagElectrodePatches[i]
                    << exit(FatalError);
            }

            if (conservationPatchRole[patchi] != 0)
            {
                FatalErrorInFunction
                    << "Patch " << diagElectrodePatches[i]
                    << " is listed more than once in diagnostics/conservation"
                    << exit(FatalError);
            }

            conservationPatchRole[patchi] = 2;
        }

        if (Pstream::master())
        {
            fileName diagDir
            (
                runTime.globalPath()/"postProcessing"/"diagnostics"
            );

            mkDir(diagDir);

            if (diagTransportLog)
            {
                transportLogPtr.reset
                (
                    new OFstream(diagDir/"transport.dat")
                );

                transportLogPtr()
                    << "# time species "
                    << "advMax diffMax migMax rxnKineticMax rxnEquilibriumMax rxnTotalMax "
                    << "advInt diffInt migInt rxnKineticInt rxnEquilibriumInt rxnTotalInt"
                    << nl;
            }

            if (diagConservationLog)
            {
                conservationLogPtr.reset
                (
                    new OFstream(diagDir/"conservation.dat")
                );

                if (diagConservationMode == "breakdown")
                {
                    conservationLogPtr()
                        << "# time species N0 Nt dN "
                        << "transport advection diffusion migration "
                        << "reactionKinetic reactionEquilibrium reactionTotal residual relativeResidual"
                        << nl;
                }
                else
                {
                    conservationLogPtr()
                        << "# time species residual relativeResidual"
                        << nl;
                }
            }

            if (diagDebyeLog)
            {
                DebyeLogPtr.reset
                (
                    new OFstream(diagDir/"Debye.dat")
                );

                DebyeLogPtr()
                    << "# time maxE_D"
                    << nl;
            }

            if (diagSatLog)
            {
                saturationLogPtr.reset
                (
                    new OFstream(diagDir/"saturation.dat")
                );

                saturationLogPtr()
                    << "# time species "
                    << "Smin Smax pEqMin pEqMax "
                    << "STotMin STotMax pEqTotMin pEqTotMax"
                    << nl;
            }
        }

        diagnosticsInitialised = true;
    }

// ############################################################ Species state ############################################################

    // Report the accepted post-convergence concentration state.
    if (printCFields)
    {
        Info<< "Concentration fields:" << nl;

        forAll(c, i)
        {
            const scalar cMin =
                gMin(c[i].primitiveField());

            const scalar cMean =
                gSum(c[i].primitiveField()*mesh.V())
               /gSum(mesh.V());

            const scalar cMax =
                gMax(c[i].primitiveField());

            Info<< "  " << c[i].name()
                << " min/mean/max = "
                << cMin << " / "
                << cMean << " / "
                << cMax
                << " mol/m3"
                << nl;
        }

        Info<< nl;
    }

    if (printYFields)
    {
        Info<< "Mass-fraction fields:" << nl;

        forAll(Y, i)
        {
            const scalar YMin = gMin(Y[i].primitiveField());
            const scalar YMax = gMax(Y[i].primitiveField());
            const scalar YMean = gSum(Y[i].primitiveField()*mesh.V())/gSum(mesh.V());

            Info<< "  " << Y[i].name()
                << " min/mean/max = "
                << YMin << " / "
                << YMean << " / "
                << YMax
                << nl;
        }

        Info<< endl;
    }

    // ############################################################ Transport ############################################################ 

    if (diagTransportActive)
    {
        if (diagTransportPrint)
        {
            Info<< "Transport diagnostics:" << nl;
        }

        forAll(c, i)
        {
            if (i == inertIndex)
            {
                continue;
            }

            const word& specieName = species[i];

            if
            (
                diagTransportSpecies.size()
             && diagTransportSpecies.find(specieName) < 0
            )
            {
                continue;
            }

            const volScalarField& ci = c[i];

            dimensionedScalar Di
            (
                "D_" + specieName,
                dimArea/dimTime,
                diffusivityDict.get<scalar>(specieName)
            );

            scalar zi = 0.0;

            if (zDictPtr)
            {
                zi =
                    zDictPtr->getOrDefault<scalar>
                    (
                        specieName,
                        0.0
                    );
            }

            surfaceScalarField migrationFlux
            (
                IOobject
                (
                    "diagMigrationFlux_" + specieName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMoles/dimTime, 0.0)
            );
            if (potentialMgr.active())
            {
                migrationFlux = fvc::interpolate(Di*zi*Fdim/(Rgas*T)*ci)
                    * (fvc::snGrad(phiE)*mesh.magSf());
            }


            volScalarField advTerm
            (
                IOobject
                (
                    "diagAdv_" + specieName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                -fvc::div(phi, ci)
            );

            volScalarField diffTerm
            (
                IOobject
                (
                    "diagDiff_" + specieName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::laplacian(Di, ci)
            );

            volScalarField migTerm
            (
                IOobject
                (
                    "diagMig_" + specieName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::div(migrationFlux)
            );

            const volScalarField& totalReaction = Rc[i];

            const scalar advMax = max(mag(advTerm)).value();
            const scalar diffMax = max(mag(diffTerm)).value();
            const scalar migMax = max(mag(migTerm)).value();
            const scalar rxnMax = max(mag(totalReaction)).value();
            const scalar kineticMax = gMax(mag(RcKinetic[i]));
            const scalar equilibriumMax = gMax(mag(RcEquilibrium[i]));

            const scalar advInt =
                gSum(mag(advTerm.primitiveField())*mesh.V());
            const scalar diffInt =
                gSum(mag(diffTerm.primitiveField())*mesh.V());
            const scalar migInt =
                gSum(mag(migTerm.primitiveField())*mesh.V());
            const scalar rxnInt =
                gSum(mag(totalReaction.primitiveField())*mesh.V());

            const scalar kineticInt = gSum(mag(RcKinetic[i])*mesh.V());
            const scalar equilibriumInt = gSum(mag(RcEquilibrium[i])*mesh.V());

            if (diagTransportPrint)
            {
                Info<< "  " << specieName
                    << " max|adv diff mig rxnKinetic rxnEquilibrium rxnTotal| = "
                    << advMax << " "
                    << diffMax << " "
                    << migMax << " "
                    << kineticMax << " "
                    << equilibriumMax << " "
                    << rxnMax
                    << " [mol/m3/s]"
                    << nl;
            }

            if
            (
                diagLogTime
             && diagTransportLog
             && transportLogPtr.valid()
            )
            {
                OFstream& os = transportLogPtr();

                os  << runTime.value() << ' '
                    << specieName << ' '
                    << advMax << ' '
                    << diffMax << ' '
                    << migMax << ' '
                    << kineticMax << ' '
                    << equilibriumMax << ' '
                    << rxnMax << ' '
                    << advInt << ' '
                    << diffInt << ' '
                    << migInt << ' '
                    << kineticInt << ' '
                    << equilibriumInt << ' '
                    << rxnInt << nl;
            }
        }

        if (diagTransportPrint)
        {
            Info<< endl;
        }
    }

// ############################################################ Conservation ############################################################

    if (diagConservationActive)
    {
        if (diagConservationPrint)
        {
            Info<< "Species conservation diagnostics:" << nl;
        }

        forAll(c, i)
        {
            if (i == inertIndex)
            {
                continue;
            }

            const word& specieName = species[i];
            const volScalarField& ci = c[i];

            dimensionedScalar Di
            (
                "D_" + specieName,
                dimArea/dimTime,
                diffusivityDict.get<scalar>(specieName)
            );

            scalar zi = 0.0;

            if (zDictPtr)
            {
                zi =
                    zDictPtr->getOrDefault<scalar>
                    (
                        specieName,
                        0.0
                    );
            }

            // Outward boundary rates [mol/s].  These are converted below
            // into signed additions to the domain by changing their sign.
            scalar advOutRate = 0.0;
            scalar diffOutRate = 0.0;
            scalar migOutRate = 0.0;

            forAll(mesh.boundary(), patchi)
            {
                if (mesh.boundaryMesh()[patchi].coupled())
                {
                    continue;
                }

                if (conservationPatchRole[patchi] == 0)
                {
                    continue;
                }

                const fvPatchScalarField& cp =
                    ci.boundaryField()[patchi];

                const fvsPatchScalarField& phip =
                    phi.boundaryField()[patchi];

                const fvsPatchScalarField& magSfp =
                    mesh.magSf().boundaryField()[patchi];

                tmp<scalarField> tsnGradCp = cp.snGrad();
                const scalarField& snGradCp = tsnGradCp();

                const fvPatchScalarField& phiEp =
                    phiE.boundaryField()[patchi];
                scalarField snGradPhiE(phiEp.size(), 0.0);
                if (potentialMgr.active()) snGradPhiE = phiEp.snGrad();

                forAll(cp, facei)
                {
                    // Positive values are outward from the domain.
                    advOutRate +=
                        cp[facei]*phip[facei];

                    diffOutRate +=
                       -Di.value()
                       *snGradCp[facei]
                       *magSfp[facei];

                    migOutRate +=
                       -Di.value()
                       *zi
                       *Fdim.value()
                       /(Rgas.value()*T.value())
                       *cp[facei]
                       *snGradPhiE[facei]
                       *magSfp[facei];
                }
            }

            reduce(advOutRate, sumOp<scalar>());
            reduce(diffOutRate, sumOp<scalar>());
            reduce(migOutRate, sumOp<scalar>());

            // Positive cumulative transport means net addition to the domain.
            cumulativeAdvection[i] += -advOutRate*dt;
            cumulativeDiffusion[i] += -diffOutRate*dt;
            cumulativeMigration[i] += -migOutRate*dt;

            const scalar reactionRate =
                gSum(Rc[i].primitiveField()*mesh.V());

            cumulativeReaction[i] += reactionRate*dt;
            if (diagConservationMode == "breakdown")
            {
                cumulativeKineticReaction[i] += gSum(RcKinetic[i]*mesh.V())*dt;
                cumulativeEquilibriumReaction[i] += gSum(RcEquilibrium[i]*mesh.V())*dt;
            }

            const scalar Nt =
                gSum(ci.primitiveField()*mesh.V());

            const scalar dN = Nt - N0[i];

            const scalar transport =
                cumulativeAdvection[i]
              + cumulativeDiffusion[i]
              + cumulativeMigration[i];

            const scalar reaction = cumulativeReaction[i];

            const scalar residual =
                dN - transport - reaction;

            const scalar scale =
                max
                (
                    max
                    (
                        max(mag(dN), mag(cumulativeAdvection[i])),
                        max
                        (
                            mag(cumulativeDiffusion[i]),
                            mag(cumulativeMigration[i])
                        )
                    ),
                    mag(cumulativeReaction[i])
                );

            const scalar relativeResidual =
                residual/(scale + SMALL);

            if (diagConservationPrint)
            {
                if (diagConservationMode == "breakdown")
                {
                    Info<< "  " << specieName << nl
                        << "    dN          = " << dN << " mol" << nl
                        << "    transport   = " << transport << " mol" << nl
                        << "      advection = "
                        << cumulativeAdvection[i] << " mol" << nl
                        << "      diffusion = "
                        << cumulativeDiffusion[i] << " mol" << nl
                        << "      migration = "
                        << cumulativeMigration[i] << " mol" << nl
                        << "    reaction total = " << reaction << " mol" << nl
                        << "      kinetic      = " << cumulativeKineticReaction[i] << " mol" << nl
                        << "      equilibrium  = " << cumulativeEquilibriumReaction[i] << " mol" << nl
                        << "    residual    = " << residual << " mol" << nl
                        << "    rel         = " << relativeResidual << nl;
                }
                else
                {
                    Info<< "  " << specieName
                        << " residual=" << residual
                        << " rel=" << relativeResidual
                        << nl;
                }
            }

            if
            (
                diagLogTime
             && diagConservationLog
             && conservationLogPtr.valid()
            )
            {
                OFstream& os = conservationLogPtr();

                if (diagConservationMode == "breakdown")
                {
                    os  << runTime.value() << ' '
                        << specieName << ' '
                        << N0[i] << ' '
                        << Nt << ' '
                        << dN << ' '
                        << transport << ' '
                        << cumulativeAdvection[i] << ' '
                        << cumulativeDiffusion[i] << ' '
                        << cumulativeMigration[i] << ' '
                        << cumulativeKineticReaction[i] << ' '
                        << cumulativeEquilibriumReaction[i] << ' '
                        << reaction << ' '
                        << residual << ' '
                        << relativeResidual
                        << nl;
                }
                else
                {
                    os  << runTime.value() << ' '
                        << specieName << ' '
                        << residual << ' '
                        << relativeResidual
                        << nl;
                }
            }
        }

        if (diagConservationPrint)
        {
            Info<< endl;
        }
    }

// ############################################################ Debye resolution ############################################################

if (diagDebyeActive)
{
    lambdaD =
        dimensionedScalar
        (
            "zero",
            dimLength,
            0.0
        );

    DebyeNumber =
        dimensionedScalar
        (
            "zero",
            dimless,
            0.0
        );
        
    ionicSumField =
        dimensionedScalar
        (
            "zero",
            dimMoles/dimVolume,
            0.0
        );

    // Debye sum is sum(z_i^2*c_i); the potential charge scale instead
    // uses sum(|z_i|*c_i) with an absolute-tolerance floor.
    forAll(c, i)
    {
        scalar zi = 0.0;

        if (zDictPtr)
        {
            zi =
                zDictPtr->getOrDefault<scalar>
                (
                    species[i],
                    0.0
                );
        }

        if (mag(zi) > SMALL)
        {
            ionicSumField += sqr(zi)*c[i];
        }
    }

    vector directionVector(vector::zero);
    directionVector[diagDebyeDirection] = 1.0;

    scalarField directionalProjectedArea
    (
        mesh.nCells(),
        0.0
    );

    forAll(mesh.Sf(), facei)
    {
        const scalar projectedArea =
            mag(mesh.Sf()[facei] & directionVector);

        const label owner =
            mesh.faceOwner()[facei];

        directionalProjectedArea[owner] +=
            projectedArea;

        if (facei < mesh.nInternalFaces())
        {
            const label neighbour =
                mesh.faceNeighbour()[facei];

            directionalProjectedArea[neighbour] +=
                projectedArea;
        }
    }

    scalarField directionalCellLength
    (
        mesh.nCells(),
        GREAT
    );

    forAll(directionalCellLength, celli)
    {
        // For an orthogonal structured cell, opposite faces both
        // contribute to projectedArea, so divide their sum by two.
        const scalar meanDirectionalFaceArea =
            0.5*directionalProjectedArea[celli];

        if (meanDirectionalFaceArea > SMALL)
        {
            directionalCellLength[celli] = mesh.V()[celli]/meanDirectionalFaceArea;
        }
    }

    scalar localMaxDebyeNumber = 0.0;

    forAll(ionicSumField, celli)
    {
        const scalar ionicSum =
            ionicSumField[celli];

        const scalar cellLength =
            directionalCellLength[celli];

        if
        (
            ionicSum > diagDebyeMinIonicSum
         && cellLength < GREAT
        )
        {
            const scalar lambdaSquared = epsilon.value()*Rgas.value()*T.value()/(Foam::sqr(Fdim.value())*ionicSum);

            lambdaD[celli] = Foam::sqrt(Foam::scalar(lambdaSquared));

            if (lambdaD[celli] > SMALL)
            {
                DebyeNumber[celli] = lambdaD[celli]/cellLength;
                localMaxDebyeNumber = max(localMaxDebyeNumber,DebyeNumber[celli]);
            }
        }
        else
        {
            lambdaD[celli] = 0.0;
            DebyeNumber[celli] = 0.0;
        }
    }
    ionicSumField.correctBoundaryConditions();
    lambdaD.correctBoundaryConditions();
    DebyeNumber.correctBoundaryConditions();

    scalar maxDebyeNumber =
        localMaxDebyeNumber;

    reduce
    (
        maxDebyeNumber,
        maxOp<scalar>()
    );

    if (diagDebyePrint)
    {
        Info<< "Debye Number diagnostics:" << nl
            << "  max(E_D) = "
            << maxDebyeNumber
            << nl
            << endl;
    }

    if (diagLogTime && diagDebyeLog && DebyeLogPtr.valid())
    {
        DebyeLogPtr()
            << runTime.value() << ' '
            << maxDebyeNumber
            << nl;
    }
}

// ############################################################ Saturation ############################################################

    if (diagSaturationActive && henrySpecies.size())
    {
        volScalarField pAbs
        (
            IOobject
            (
                "pAbs_tmp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pAbsRef + rho*p
        );

        pEq_tot = dimensionedScalar("zero", dimPressure, 0.0);

        forAll(henrySpecies, hi)
        {
            const label speciei = henrySpeciesIndex[hi];
            const word& specieName = henrySpecies[hi];

            dimensionedScalar Hi
            (
                "H_" + specieName,
                dimMoles/dimVolume/dimPressure,
                henryConstants[hi]
            );

            pEqHenry[hi] = c[speciei]/Hi;
            saturation[hi] = pEqHenry[hi]/max(pAbs, pSmall);

            pEq_tot += pEqHenry[hi];
        }

        S_tot = pEq_tot/max(pAbs, pSmall);

        // These global reductions must run on all MPI ranks.
        const scalar STotMin = min(S_tot).value();
        const scalar STotMax = max(S_tot).value();
        const scalar pEqTotMin = min(pEq_tot).value();
        const scalar pEqTotMax = max(pEq_tot).value();

        if
        (
            diagSatPrintSi
        || diagSatPrintSTot
        || diagSatPrintPEqi
        || diagSatPrintPEqTot
        )
        {
            Info<< "Saturation diagnostics:" << nl;
        }

        forAll(henrySpecies, hi)
        {
            // Also calculate these on all MPI ranks.
            const scalar SiMin = min(saturation[hi]).value();
            const scalar SiMax = max(saturation[hi]).value();
            const scalar pEqiMin = min(pEqHenry[hi]).value();
            const scalar pEqiMax = max(pEqHenry[hi]).value();

            if (diagSatPrintSi)
            {
                Info<< "  S_" << henrySpecies[hi]
                    << " min/max = "
                    << SiMin << " / " << SiMax
                    << nl;
            }

            if (diagSatPrintPEqi)
            {
                Info<< "  pEq_" << henrySpecies[hi]
                    << " min/max = "
                    << pEqiMin << " / " << pEqiMax
                    << " Pa"
                    << nl;
            }

            // No collective operations inside this master-only block.
            if (diagSatLog && diagLogTime && saturationLogPtr.valid())
            {
                saturationLogPtr()
                    << runTime.value() << ' '
                    << henrySpecies[hi] << ' '
                    << SiMin << ' '
                    << SiMax << ' '
                    << pEqiMin << ' '
                    << pEqiMax << ' '
                    << STotMin << ' '
                    << STotMax << ' '
                    << pEqTotMin << ' '
                    << pEqTotMax
                    << nl;
            }
        }

        if (diagSatPrintSTot)
        {
            Info<< "  S_tot min/max = "
                << STotMin << " / " << STotMax
                << nl;
        }

        if (diagSatPrintPEqTot)
        {
            Info<< "  pEq_tot min/max = "
                << pEqTotMin << " / " << pEqTotMax
                << " Pa"
                << nl;
        }

        if
        (
            diagSatPrintSi
        || diagSatPrintSTot
        || diagSatPrintPEqi
        || diagSatPrintPEqTot
        )
        {
            Info<< endl;
        }
    }

// ############################################################ Inert species closure ############################################################

    if (diagInertClosurePrint)
    {
        volScalarField YSum
        (
            IOobject
            (
                "YSum_diag",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );

        forAll(Y, i)
        {
            YSum += Y[i];
        }

        const scalar maxResidual =
            max(mag(YSum - scalar(1))).value();

        const scalar maxNonInertSum =
            max(Yt).value();

        const scalar maxNonInertExcess =
            max(maxNonInertSum - 1.0, 0.0);

        Info<< "Inert closure diagnostics:" << nl
            << "  max|sum(Y)-1| = "
            << maxResidual << nl
            << "  max sum(non-inert Y) = "
            << maxNonInertSum << nl
            << "  max non-inert excess above one = "
            << maxNonInertExcess << nl
            << endl;
    }
}

