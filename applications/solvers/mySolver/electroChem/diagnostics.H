// diagnostics.H
// Unified diagnostics for concentration-primary electroChem solver.

{
    const scalar dt = runTime.deltaTValue();
    const bool diagWriteTime = runTime.writeTime();

    const dictionary& diffusivityDict = reactions.subDict("diffusivity");
    const bool hasChargeNumber = reactions.found("chargeNumber");

    const dictionary* zDictPtr = nullptr;

    if (hasChargeNumber)
    {
        zDictPtr = &reactions.subDict("chargeNumber");
    }

    static bool diagnosticsInitialised(false);

    static List<scalar> N0;
    static List<scalar> cumulativeFlowIn;
    static List<scalar> cumulativeFlowOut;
    static List<scalar> cumulativeReaction;

    static autoPtr<OFstream> transportLogPtr;
    static autoPtr<OFstream> conservationLogPtr;
    static autoPtr<OFstream> chargeLogPtr;
    static autoPtr<OFstream> saturationLogPtr;
    static autoPtr<OFstream> DebyeLogPtr;

    // Patch roles for conservation:
    // 0 = ignored, 1 = flow patch, 2 = electrode patch
    static List<label> conservationPatchRole;

    if (!diagnosticsInitialised)
    {
        N0.setSize(c.size(), 0.0);
        cumulativeFlowIn.setSize(c.size(), 0.0);
        cumulativeFlowOut.setSize(c.size(), 0.0);
        cumulativeReaction.setSize(c.size(), 0.0);

        conservationPatchRole.setSize(mesh.boundary().size(), 0);

        forAll(c, i)
        {
            N0[i] = gSum(c[i].primitiveField()*mesh.V());
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

            conservationPatchRole[patchi] = 2;
        }

        if (Pstream::master())
        {
            fileName diagDir(runTime.globalPath()/"postProcessing"/"diagnostics");
            mkDir(diagDir);

            if (diagTransportLog)
            {
                transportLogPtr.reset
                (
                    new OFstream(diagDir/"transport.dat")
                );

                transportLogPtr()
                    << "# time species "
                    << "advMax diffMax migMax rxnMax "
                    << "advInt diffInt migInt rxnInt"
                    << nl;
            }

            if (diagConservationLog)
            {
                conservationLogPtr.reset
                (
                    new OFstream(diagDir/"conservation.dat")
                );

                conservationLogPtr()
                    << "# time species "
                    << "N0 Nt dN "
                    << "flowIn flowOut reaction "
                    << "residual relativeResidual"
                    << nl;
            }

            if (diagChargeLog)
            {
                chargeLogPtr.reset
                (
                    new OFstream(diagDir/"charge.dat")
                );

                OFstream& os = chargeLogPtr();

                os << "# time";

                if
                (
                    diagChargeOutput == "molarCharge"
                || diagChargeOutput == "both"
                )
                {
                    os  << " molarMin"
                        << " molarMean"
                        << " molarMax"
                        << " totalMolarCharge";
                }

                if
                (
                    diagChargeOutput == "chargeDensity"
                || diagChargeOutput == "both"
                )
                {
                    os  << " chargeMin"
                        << " chargeMean"
                        << " chargeMax"
                        << " totalCharge";
                }

                os << nl;
            }

            if (diagDebyeLog)
            {
                DebyeLogPtr.reset
                (
                    new OFstream
                    (
                        diagDir/"Debye.dat"
                    )
                );

                DebyeLogPtr()
                    << "# time maxE_D"
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
                // print only, no file needed
            }

            if
            (
                diagSatWriteSi
             || diagSatWriteSTot
             || diagSatWritePEqi
             || diagSatWritePEqTot
            )
            {
                // fields written by AUTO_WRITE, no log needed
            }

            if
            (
                diagSaturationActive
             && (
                    diagSatPrintSi
                 || diagSatPrintSTot
                 || diagSatPrintPEqi
                 || diagSatPrintPEqTot
                )
            )
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

// ############################################################ print species  ############################################################

    if (printCFields || printCFieldsBeforeClipping)
    {
        Info<< "Concentration fields:" << nl;

        if (printCFieldsBeforeClipping)
        {
            Info<< "  Pre-clipping:" << nl;

            forAll(c, i)
            {
                Info<< "    " << c[i].name()
                    << " min/mean/max = "
                    << cMinBeforeClipping[i] << " / "
                    << cMeanBeforeClipping[i] << " / "
                    << cMaxBeforeClipping[i]
                    << " mol/m3"
                    << nl;
            }
        }

        if (printCFields)
        {
            if (printCFieldsBeforeClipping)
            {
                Info<< "  After clipping:" << nl;
            }

            forAll(c, i)
            {
                const scalar cMin =
                    gMin(c[i].primitiveField());

                const scalar cMean =
                    gSum(c[i].primitiveField()*mesh.V())
                /gSum(mesh.V());

                const scalar cMax =
                    gMax(c[i].primitiveField());

                Info<<
                    (
                        printCFieldsBeforeClipping
                    ? "    "
                    : "  "
                    )
                    << c[i].name()
                    << " min/mean/max = "
                    << cMin << " / "
                    << cMean << " / "
                    << cMax
                    << " mol/m3"
                    << nl;
            }
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
                zi = zDictPtr->getOrDefault<scalar>(specieName, 0.0);
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
                fvc::interpolate(Di*zi*Fdim/(Rgas*T)*ci)
              * (fvc::interpolate(gradPhiE) & mesh.Sf())
            );

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

            const scalar advMax = max(mag(advTerm)).value();
            const scalar diffMax = max(mag(diffTerm)).value();
            const scalar migMax = max(mag(migTerm)).value();
            const scalar rxnMax = max(mag(Rc[i])).value();

            const scalar advInt =
                gSum(mag(advTerm.primitiveField())*mesh.V());

            const scalar diffInt =
                gSum(mag(diffTerm.primitiveField())*mesh.V());

            const scalar migInt =
                gSum(mag(migTerm.primitiveField())*mesh.V());

            const scalar rxnInt =
                gSum(mag(Rc[i].primitiveField())*mesh.V());

            if (diagTransportPrint)
            {
                Info<< "  " << specieName
                    << " max|adv diff mig rxn| = "
                    << advMax << " "
                    << diffMax << " "
                    << migMax << " "
                    << rxnMax
                    << " [mol/m3/s]"
                    << nl;
            }

            if (diagWriteTime && diagTransportLog && transportLogPtr.valid())
            {
                transportLogPtr()
                    << runTime.value() << ' '
                    << specieName << ' '
                    << advMax << ' '
                    << diffMax << ' '
                    << migMax << ' '
                    << rxnMax << ' '
                    << advInt << ' '
                    << diffInt << ' '
                    << migInt << ' '
                    << rxnInt
                    << nl;
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
                zi = zDictPtr->getOrDefault<scalar>(specieName, 0.0);
            }

            scalar flowInRate = 0.0;
            scalar flowOutRate = 0.0;
            scalar electrodeGenRate = 0.0;

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

                const fvPatchVectorField& gradPhiEp =
                    gradPhiE.boundaryField()[patchi];

                const fvsPatchVectorField& Sfp =
                    mesh.Sf().boundaryField()[patchi];

                scalar patchFlowIn = 0.0;
                scalar patchFlowOut = 0.0;
                scalar patchElectrodeGen = 0.0;

                forAll(cp, facei)
                {
                    const scalar adv =
                        cp[facei]*phip[facei];

                    const scalar diff =
                       -Di.value()
                       *snGradCp[facei]
                       *magSfp[facei];

                    const scalar mig =
                       -Di.value()
                       *zi
                       *Fdim.value()
                       /(Rgas.value()*T.value())
                       *cp[facei]
                       *(gradPhiEp[facei] & Sfp[facei]);

                    const scalar flow = adv + diff + mig;

                    patchFlowOut += max(flow, 0.0);
                    patchFlowIn  += max(-flow, 0.0);

                    if (conservationPatchRole[patchi] == 2)
                    {
                        // Generation into domain from non-advective flux:
                        // inward diff+mig is positive generation.
                        const scalar nonAdvectiveOut = diff + mig;
                        patchElectrodeGen += -nonAdvectiveOut;
                    }
                }

                if (conservationPatchRole[patchi] == 1)
                {
                    flowInRate  += patchFlowIn;
                    flowOutRate += patchFlowOut;
                }
                else if (conservationPatchRole[patchi] == 2)
                {
                    electrodeGenRate += patchElectrodeGen;
                }
            }

            reduce(flowInRate, sumOp<scalar>());
            reduce(flowOutRate, sumOp<scalar>());
            reduce(electrodeGenRate, sumOp<scalar>());

            const scalar reactionRate =
                gSum(Rc[i].primitiveField()*mesh.V());

            const scalar totalGenerationRate =
                reactionRate + electrodeGenRate;

            cumulativeFlowIn[i] += flowInRate*dt;
            cumulativeFlowOut[i] += flowOutRate*dt;
            cumulativeReaction[i] += totalGenerationRate*dt;

            const scalar Nt =
                gSum(ci.primitiveField()*mesh.V());

            const scalar dN = Nt - N0[i];

            const scalar residual =
                dN
              + cumulativeFlowOut[i]
              - cumulativeFlowIn[i]
              - cumulativeReaction[i];

            const scalar scale =
                max
                (
                    max(mag(dN), mag(cumulativeReaction[i])),
                    max(mag(cumulativeFlowIn[i]), mag(cumulativeFlowOut[i]))
                );

            const scalar relativeResidual = residual/(scale + SMALL);

            if (diagConservationPrint)
            {
                if (diagConservationMode == "full")
                {
                    Info<< "  " << specieName
                        << " dN=" << dN
                        << " flowIn=" << cumulativeFlowIn[i]
                        << " flowOut=" << cumulativeFlowOut[i]
                        << " generation=" << cumulativeReaction[i]
                        << " residual=" << residual
                        << " rel=" << relativeResidual
                        << nl;
                }
                else
                {
                    Info<< "  " << specieName
                        << " residual=" << residual
                        << " rel=" << relativeResidual
                        << nl;
                }
            }

            if (diagWriteTime && diagConservationLog && conservationLogPtr.valid())
            {
                OFstream& os = conservationLogPtr();

                if (diagConservationMode == "full")
                {
                    os  << runTime.value() << ' '
                        << specieName << ' '
                        << N0[i] << ' '
                        << Nt << ' '
                        << dN << ' '
                        << cumulativeFlowIn[i] << ' '
                        << cumulativeFlowOut[i] << ' '
                        << cumulativeReaction[i] << ' '
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

// ############################################################ Charge ############################################################

    if (diagChargeActive)
    {
        molarCharge =
            dimensionedScalar
            (
                "zero",
                dimMoles/dimVolume,
                0.0
            );

        forAll(c, i)
        {
            const word& specieName = species[i];

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

            if (mag(zi) > SMALL)
            {
                molarCharge += zi*c[i];
            }
        }

        chargeDensity = Fdim*molarCharge;

        const scalar molarChargeMin =
            min(molarCharge).value();

        const scalar molarChargeMax =
            max(molarCharge).value();

        const scalar molarChargeMean =
            gSum(molarCharge.primitiveField()*mesh.V())
        /gSum(mesh.V());

        const scalar chargeDensityMin =
            min(chargeDensity).value();

        const scalar chargeDensityMax =
            max(chargeDensity).value();

        const scalar chargeDensityMean =
            gSum(chargeDensity.primitiveField()*mesh.V())
        /gSum(mesh.V());

        const scalar totalMolarCharge =
            gSum(molarCharge.primitiveField()*mesh.V());

        const scalar totalCharge =
            gSum(chargeDensity.primitiveField()*mesh.V());

        if (diagChargePrint)
        {
            Info<< "Charge diagnostics:" << nl;

            if
            (
                diagChargeOutput == "molarCharge"
            || diagChargeOutput == "both"
            )
            {
                Info<< "  molarCharge min/mean/max = "
                    << molarChargeMin << " / "
                    << molarChargeMean << " / "
                    << molarChargeMax
                    << " mol/m3"
                    << nl;

                Info<< "  total molar charge = "
                    << totalMolarCharge
                    << " mol"
                    << nl;
            }

            if
            (
                diagChargeOutput == "chargeDensity"
            || diagChargeOutput == "both"
            )
            {
                Info<< "  chargeDensity min/mean/max = "
                    << chargeDensityMin << " / "
                    << chargeDensityMean << " / "
                    << chargeDensityMax
                    << " C/m3"
                    << nl;

                Info<< "  total charge = "
                    << totalCharge
                    << " C"
                    << nl;
            }

            Info<< endl;
        }

        if (diagWriteTime && diagChargeLog && chargeLogPtr.valid())
        {
            OFstream& os = chargeLogPtr();

            os << runTime.value();

            if
            (
                diagChargeOutput == "molarCharge"
            || diagChargeOutput == "both"
            )
            {
                os  << ' ' << molarChargeMin
                    << ' ' << molarChargeMean
                    << ' ' << molarChargeMax
                    << ' ' << totalMolarCharge;
            }

            if
            (
                diagChargeOutput == "chargeDensity"
            || diagChargeOutput == "both"
            )
            {
                os  << ' ' << chargeDensityMin
                    << ' ' << chargeDensityMean
                    << ' ' << chargeDensityMax
                    << ' ' << totalCharge;
            }

            os << nl;
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

    if (diagWriteTime && diagDebyeLog && DebyeLogPtr.valid())
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
            if (diagWriteTime && saturationLogPtr.valid())
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

