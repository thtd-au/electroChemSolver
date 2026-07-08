// Species conservation monitor

if (speciesConservationActive)
{
    const scalar dt = runTime.deltaTValue();

    const dimensionedScalar rho
    (
        "rho",
        dimDensity,
        laminarTransport
    );

    const dictionary& diffusivityDict = reactions.subDict("diffusivity");
    const dictionary& molarMassDict = reactions.subDict("molarMass");

    if (speciesConservationLogToInfo)
    {
        Info<< "Species conservation:" << nl;
    }

    forAll(Y, i)
    {
        if (i == inertIndex)
        {
            continue;
        }

        const volScalarField& Yi = Y[i];
        const scalar Di = diffusivityDict.get<scalar>(species[i]);
        const scalar Mi = molarMassDict.get<scalar>(species[i]);
        const scalar rhoByMi = rho.value()/Mi;
        const scalar rhoDiByMi = rho.value()*Di/Mi;

        scalar advInRate = 0.0;
        scalar advOutRate = 0.0;
        scalar diffInRate = 0.0;
        scalar diffOutRate = 0.0;
        scalar flowInRate = 0.0;
        scalar flowOutRate = 0.0;
        // Total species generation [mol/s]
        // Includes:
        //   - electrode boundary generation
        //   - volumetric reaction generation
        scalar genRate = 0.0;

        List<scalar> patchAdvIn(mesh.boundary().size(), 0.0);
        List<scalar> patchAdvOut(mesh.boundary().size(), 0.0);
        List<scalar> patchDiffIn(mesh.boundary().size(), 0.0);
        List<scalar> patchDiffOut(mesh.boundary().size(), 0.0);
        List<scalar> patchFlowIn(mesh.boundary().size(), 0.0);
        List<scalar> patchFlowOut(mesh.boundary().size(), 0.0);
        List<scalar> patchGen(mesh.boundary().size(), 0.0);

        forAll(mesh.boundary(), patchi)
        {
            // Skip processor/coupled patches. They are decomposition interfaces,
            // not physical domain boundaries, and coupledFvPatchField::snGrad()
            // is not implemented.
            if (mesh.boundaryMesh()[patchi].coupled())
            {
                continue;
            }

            if
            (
                !speciesConservationAllPatchDiagnostics
             && !speciesConservationPrintPatchBreakdown
             && speciesConservationPatchRole[patchi] == 0
            )
            {
                continue;
            }

            const fvPatchScalarField& Yp = Yi.boundaryField()[patchi];
            const fvsPatchScalarField& phip = phi.boundaryField()[patchi];
            const fvsPatchScalarField& magSfp = mesh.magSf().boundaryField()[patchi];

            tmp<scalarField> tsnGradYp = Yp.snGrad();
            const scalarField& snGradYp = tsnGradYp();

            scalar pAdvIn = 0.0;
            scalar pAdvOut = 0.0;
            scalar pDiffIn = 0.0;
            scalar pDiffOut = 0.0;
            scalar pFlowIn = 0.0;
            scalar pFlowOut = 0.0;
            scalar pGen = 0.0;

            forAll(Yp, facei)
            {
                // Molar advective rate through the face [mol/s].
                // phi is positive outward by OpenFOAM boundary convention.
                const scalar adv = rhoByMi*Yp[facei]*phip[facei];

                // Physical molar diffusive rate through the face [mol/s].
                // Positive is outward: -D grad(c).n A.
                const scalar diff = -rhoDiByMi*snGradYp[facei]*magSfp[facei];

                // Combined physical boundary transport rate [mol/s].
                // This is the quantity used in the conservation residual.
                const scalar flow = adv + diff;

                // Mechanism split, diagnostics only.
                pAdvOut += max(adv, 0.0);
                pAdvIn  += max(-adv, 0.0);

                pDiffOut += max(diff, 0.0);
                pDiffIn  += max(-diff, 0.0);

                // Net transport split, used for conservation accounting.
                pFlowOut += max(flow, 0.0);
                pFlowIn  += max(-flow, 0.0);

                // For electrode patches, positive generation into the domain is
                // rho*D/M*grad(Y).n A = s*j/(nF) A with the implemented BC sign.
                if (speciesConservationPatchRole[patchi] == 2)
                {
                    pGen += rhoDiByMi*snGradYp[facei]*magSfp[facei];
                }
            }

            patchAdvIn[patchi] = pAdvIn;
            patchAdvOut[patchi] = pAdvOut;
            patchDiffIn[patchi] = pDiffIn;
            patchDiffOut[patchi] = pDiffOut;
            patchFlowIn[patchi] = pFlowIn;
            patchFlowOut[patchi] = pFlowOut;
            patchGen[patchi] = pGen;

            if (speciesConservationPatchRole[patchi] == 1)
            {
                advInRate += pAdvIn;
                advOutRate += pAdvOut;
                diffInRate += pDiffIn;
                diffOutRate += pDiffOut;

                flowInRate += pFlowIn;
                flowOutRate += pFlowOut;
            }
            else if (speciesConservationPatchRole[patchi] == 2)
            {
                genRate += pGen;
            }
        }

        reduce(advInRate, sumOp<scalar>());
        reduce(advOutRate, sumOp<scalar>());
        reduce(diffInRate, sumOp<scalar>());
        reduce(diffOutRate, sumOp<scalar>());
        reduce(flowInRate, sumOp<scalar>());
        reduce(flowOutRate, sumOp<scalar>());
        reduce(genRate, sumOp<scalar>());

        // Add volumetric reaction source to generation [mol/s]
        const scalar rxnGenRate = rhoByMi*gSum(RY[i].primitiveField()*mesh.V());
        genRate += rxnGenRate;

        // Add pseudo source from water equilibrium
        scalar waterEqGenRate = 0.0;

        if (runWaterEquilibrium)
        {
            waterEqGenRate = rhoByMi*gSum(RYWaterEq[i].primitiveField()*mesh.V());

            genRate += waterEqGenRate;
        }

        speciesConservationAdvIn[i] += advInRate*dt;
        speciesConservationAdvOut[i] += advOutRate*dt;
        speciesConservationDiffIn[i] += diffInRate*dt;
        speciesConservationDiffOut[i] += diffOutRate*dt;

        speciesConservationFlowIn[i] += flowInRate*dt;
        speciesConservationFlowOut[i] += flowOutRate*dt;

        speciesConservationGen[i] += genRate*dt;
        speciesConservationWaterEqGen[i] += waterEqGenRate*dt;

        const scalar Nt = rhoByMi*gSum(Yi.primitiveField()*mesh.V());
        const scalar dN = Nt - speciesConservationN0[i];

        const scalar residual =
            dN
          + speciesConservationFlowOut[i]
          - speciesConservationFlowIn[i]
          - speciesConservationGen[i];

        const scalar scale = max
        (
            max(mag(dN), mag(speciesConservationGen[i])),
            max(mag(speciesConservationFlowIn[i]), mag(speciesConservationFlowOut[i]))
        );

        const scalar relativeResidual = residual/(scale + SMALL);

        if (speciesConservationLogToInfo)
        {
            if (speciesConservationMode == "full")
            {
                Info<< "  " << species[i] << nl
                    << "    inventory:" << nl
                    << "      N0=" << speciesConservationN0[i]
                    << " Nt=" << Nt
                    << " dN=" << dN << nl

                    << "    rates [mol/s]:" << nl
                    << "      flowIn=" << flowInRate
                    << " flowOut=" << flowOutRate
                    << " gen=" << genRate << nl
                    << " waterEqGenRate=" << waterEqGenRate << nl

                    << "      advIn=" << advInRate
                    << " advOut=" << advOutRate
                    << " diffIn=" << diffInRate
                    << " diffOut=" << diffOutRate << nl

                    << "    totals [mol]:" << nl
                    << "      flowIn=" << speciesConservationFlowIn[i]
                    << " flowOut=" << speciesConservationFlowOut[i]
                    << " gen=" << speciesConservationGen[i] << nl
                    << "      advIn=" << speciesConservationAdvIn[i]
                    << " advOut=" << speciesConservationAdvOut[i]
                    << " diffIn=" << speciesConservationDiffIn[i]
                    << " diffOut=" << speciesConservationDiffOut[i] << nl

                    << "    residual:" << nl
                    << "      R=" << residual
                    << " relR=" << relativeResidual << nl;
            }
            else
            {
                Info<< "  " << species[i]
                    << " dN=" << dN
                    << " flowIn=" << speciesConservationFlowIn[i]
                    << " flowOut=" << speciesConservationFlowOut[i]
                    << " gen=" << speciesConservationGen[i]
                    << " waterEqGen=" << speciesConservationWaterEqGen[i]
                    << " R=" << residual
                    << " relR=" << relativeResidual << nl;
            }

            if (speciesConservationAllPatchDiagnostics)
            {
                forAll(mesh.boundary(), patchi)
                {
                    if (speciesConservationPatchRole[patchi] == 0)
                    {
                        const scalar leakRate =
                            patchFlowOut[patchi] - patchFlowIn[patchi];

                        if (mag(leakRate) > speciesConservationWarningTolerance)
                        {
                            WarningInFunction
                                << "species " << species[i]
                                << ": unaccounted patch "
                                << mesh.boundaryMesh()[patchi].name()
                                << " has net transport rate " << leakRate
                                << " mol/s" << nl;
                        }
                    }
                }
            }

            if (speciesConservationPrintPatchBreakdown)
            {
                Info<< "    patch diagnostics:" << nl;

                forAll(mesh.boundary(), patchi)
                {
                    Info<< "      " << mesh.boundaryMesh()[patchi].name()
                        << " flowIn=" << patchFlowIn[patchi]
                        << " flowOut=" << patchFlowOut[patchi]
                        << " advIn=" << patchAdvIn[patchi]
                        << " advOut=" << patchAdvOut[patchi]
                        << " diffIn=" << patchDiffIn[patchi]
                        << " diffOut=" << patchDiffOut[patchi]
                        << " gen=" << patchGen[patchi]
                        << nl;
                }
            }
        }

        if (speciesConservationWriteToFile && speciesConservationFilePtr.valid())
        {
            OFstream& os = speciesConservationFilePtr();

            if (speciesConservationMode == "full")
            {
                os  << runTime.value() << ' ' << species[i] << ' '
                    << speciesConservationN0[i] << ' ' << Nt << ' ' << dN << ' '
                    << speciesConservationFlowIn[i] << ' '
                    << speciesConservationFlowOut[i] << ' '
                    << speciesConservationAdvIn[i] << ' '
                    << speciesConservationAdvOut[i] << ' '
                    << speciesConservationDiffIn[i] << ' '
                    << speciesConservationDiffOut[i] << ' '
                    << speciesConservationGen[i] << ' '
                    << speciesConservationWaterEqGen[i] << ' '
                    << residual << ' ' << relativeResidual << nl;
            }
            else
            {
                os  << runTime.value() << ' ' << species[i] << ' '
                    << dN << ' '
                    << speciesConservationFlowIn[i] << ' '
                    << speciesConservationFlowOut[i] << ' '
                    << speciesConservationGen[i] << ' '
                    << speciesConservationWaterEqGen[i] << ' '
                    << residual << ' ' << relativeResidual << nl;
            }
        }
    }

    // if (speciesConservationLogToInfo)
    // {
    //     Info<< endl;
    // }
}

