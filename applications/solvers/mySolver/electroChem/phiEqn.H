// phiEqn.H
// ENPhi galvanostatic potential solve.
// Diagnostic version: compares the current enforced by the potential equation
// against the charge-weighted species transport current on the SAME final
// (equilibrium-projected) concentration state.

{
    const dictionary& diffusivityDict = reactions.subDict("diffusivity");
    const bool hasChargeNumber = reactions.found("chargeNumber");

    const dictionary* zDictPtr = nullptr;

    if (hasChargeNumber)
    {
        zDictPtr = &reactions.subDict("chargeNumber");
    }

    // Keep the legacy cell-centred diagnostic fields up to date.
    jDiff =
        dimensionedVector
        (
            "zero",
            dimCurrentDensity,
            vector::zero
        );

    kappa =
        dimensionedScalar
        (
            "zero",
            dimConductivity,
            0.0
        );

    // Authoritative face fields used by the ENPhi potential discretisation.
    jDiffFlux =
        dimensionedScalar
        (
            "zero",
            dimCurrentDensity*dimArea,
            0.0
        );

    kappaFace =
        dimensionedScalar
        (
            "zero",
            dimConductivity,
            0.0
        );

    // --------------------------------------------------------------------- //
    // Build potential-equation coefficients from the CURRENT
    // equilibrium-projected concentration state.
    // --------------------------------------------------------------------- //

    forAll(c, i)
    {
        scalar zi = 0.0;

        if (zDictPtr)
        {
            zi = zDictPtr->getOrDefault<scalar>(species[i], 0.0);
        }

        if (mag(zi) <= SMALL)
        {
            continue;
        }

        dimensionedScalar Di
        (
            "D_" + species[i],
            dimArea/dimTime,
            diffusivityDict.get<scalar>(species[i])
        );

        if (Di.value() <= SMALL)
        {
            FatalIOErrorInFunction(reactions)
                << "Diffusivity for charged species " << species[i]
                << " must be positive"
                << exit(FatalIOError);
        }

        // Legacy cell-centred diagnostics.
        jDiff -= Fdim*zi*Di*fvc::grad(c[i]);
        kappa += sqr(Fdim)*sqr(zi)*Di*c[i]/(Rgas*T);

        // Exact face diffusive-current representation used by the potential
        // equation.
        jDiffFlux -=
            Fdim*zi*Di
           *fvc::snGrad(c[i])
           *mesh.magSf();

        // Exact face migration coefficient.  This is algebraically identical
        // to F*z_i times the concentration coefficient used in the species
        // migration flux.
        kappaFace +=
            Fdim*zi
           *fvc::interpolate
            (
                Di*zi*Fdim/(Rgas*T)*c[i]
            );
    }

    jDiff.correctBoundaryConditions();
    kappa.correctBoundaryConditions();

    const scalar kappaMin = gMin(kappaFace.primitiveField());

    if (kappaMin < minConductivity.value())
    {
        FatalErrorInFunction
            << "ENPhi potential update encountered face conductivity below "
            << "minConductivity." << nl
            << "min(kappaFace) = " << kappaMin << " S/m" << nl
            << "minConductivity = " << minConductivity.value() << " S/m"
            << exit(FatalError);
    }

    // --------------------------------------------------------------------- //
    // Potential equation
    // --------------------------------------------------------------------- //
    //
    //     div(jDiffFlux - kappaFace*grad(phiE).Sf) = 0
    //
    fvScalarMatrix phiEEqn
    (
        fvm::laplacian
        (
            kappaFace,
            phiE,
            "laplacian(kappa,phiE)"
        )
     == fvc::div(jDiffFlux)
    );

    phiEEqn.setReference(phiERefCell, phiERefValue);
    phiEEqn.relax();
    phiEEqn.solve(phiE.name());
    phiE.correctBoundaryConditions();

    E = -fvc::grad(phiE);
    E.correctBoundaryConditions();

    // --------------------------------------------------------------------- //
    // Current enforced by phiEqn
    // --------------------------------------------------------------------- //

    surfaceScalarField ionicCurrentFlux
    (
        IOobject
        (
            "ionicCurrentFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        jDiffFlux
      - kappaFace
       *fvc::snGrad(phiE)
       *mesh.magSf()
    );

    volScalarField currentResidual
    (
        IOobject
        (
            "currentResidual",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(ionicCurrentFlux)
    );

    enPhiCurrentLocalError =
        gSum(mag(currentResidual.primitiveField())*mesh.V());

    enPhiCurrentGlobalError =
        gSum(currentResidual.primitiveField()*mesh.V());


    // ===================================================================== //
    // SAME-STATE SPECIES-CURRENT DIAGNOSTIC
    // ===================================================================== //
    // Reconstruct the complete charge-weighted species current using the
    // SAME concentration fields and the newly solved phiE used above.
    //
    //   I_adv  = F sum_i z_i (phi*c_i,f)
    //   I_diff = F sum_i z_i (-D_i snGrad(c_i)|Sf|)
    //   I_mig  = F sum_i z_i (-migrationFlux_i)
    //
    // and
    //
    //   I_species = I_adv + I_diff + I_mig
    //
    // The potential equation enforces only I_diff + I_mig under the bulk-EN
    // assumption.  Therefore this diagnostic tells us directly whether the
    // remaining discrepancy is convection or something else.
    // ===================================================================== //

    surfaceScalarField speciesAdvCurrentFlux
    (
        IOobject
        (
            "speciesAdvCurrentFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimCurrent,
            0.0
        )
    );

    surfaceScalarField speciesDiffCurrentFlux
    (
        IOobject
        (
            "speciesDiffCurrentFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimCurrent,
            0.0
        )
    );

    surfaceScalarField speciesMigCurrentFlux
    (
        IOobject
        (
            "speciesMigCurrentFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimCurrent,
            0.0
        )
    );

    forAll(c, i)
    {
        scalar zi = 0.0;

        if (zDictPtr)
        {
            zi = zDictPtr->getOrDefault<scalar>(species[i], 0.0);
        }

        if (mag(zi) <= SMALL)
        {
            continue;
        }

        dimensionedScalar Di
        (
            "D_" + species[i],
            dimArea/dimTime,
            diffusivityDict.get<scalar>(species[i])
        );

        surfaceScalarField speciesMigrationFlux
        (
            IOobject
            (
                "phiDiagMigrationFlux_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvc::interpolate
            (
                Di*zi*Fdim/(Rgas*T)*c[i]
            )
           *fvc::snGrad(phiE)
           *mesh.magSf()
        );

        speciesAdvCurrentFlux +=
            Fdim*zi*fvc::interpolate(c[i])*phi;

        speciesDiffCurrentFlux +=
            Fdim*zi
           *(
                -Di*fvc::snGrad(c[i])*mesh.magSf()
            );

        speciesMigCurrentFlux +=
            Fdim*zi*(-speciesMigrationFlux);
    }

    surfaceScalarField speciesNoAdvCurrentFlux
    (
        IOobject
        (
            "speciesNoAdvCurrentFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        speciesDiffCurrentFlux + speciesMigCurrentFlux
    );

    surfaceScalarField speciesTotalCurrentFlux
    (
        IOobject
        (
            "speciesTotalCurrentFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        speciesAdvCurrentFlux + speciesNoAdvCurrentFlux
    );

    surfaceScalarField noAdvMinusPhiCurrentFlux
    (
        IOobject
        (
            "noAdvMinusPhiCurrentFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        speciesNoAdvCurrentFlux - ionicCurrentFlux
    );

    surfaceScalarField totalMinusPhiCurrentFlux
    (
        IOobject
        (
            "totalMinusPhiCurrentFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        speciesTotalCurrentFlux - ionicCurrentFlux
    );

    volScalarField divSpeciesAdvCurrent
    (
        IOobject
        (
            "divSpeciesAdvCurrent",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(speciesAdvCurrentFlux)
    );

    volScalarField divSpeciesDiffCurrent
    (
        IOobject
        (
            "divSpeciesDiffCurrent",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(speciesDiffCurrentFlux)
    );

    volScalarField divSpeciesMigCurrent
    (
        IOobject
        (
            "divSpeciesMigCurrent",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(speciesMigCurrentFlux)
    );

    volScalarField divSpeciesNoAdvCurrent
    (
        IOobject
        (
            "divSpeciesNoAdvCurrent",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(speciesNoAdvCurrentFlux)
    );

    volScalarField divSpeciesTotalCurrent
    (
        IOobject
        (
            "divSpeciesTotalCurrent",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(speciesTotalCurrentFlux)
    );

    volScalarField divNoAdvMinusPhi
    (
        IOobject
        (
            "divNoAdvMinusPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(noAdvMinusPhiCurrentFlux)
    );

    volScalarField divTotalMinusPhi
    (
        IOobject
        (
            "divTotalMinusPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(totalMinusPhiCurrentFlux)
    );

    Info<< "CURRENT-COMPARE same-state summary:" << nl
        << "  max|div(phi current)|       = "
        << gMax(mag(currentResidual.primitiveField()))
        << " A/m3" << nl
        << "  max|div(species adv)|       = "
        << gMax(mag(divSpeciesAdvCurrent.primitiveField()))
        << " A/m3" << nl
        << "  max|div(species diff)|      = "
        << gMax(mag(divSpeciesDiffCurrent.primitiveField()))
        << " A/m3" << nl
        << "  max|div(species mig)|       = "
        << gMax(mag(divSpeciesMigCurrent.primitiveField()))
        << " A/m3" << nl
        << "  max|div(species no-adv)|    = "
        << gMax(mag(divSpeciesNoAdvCurrent.primitiveField()))
        << " A/m3" << nl
        << "  max|div(species total)|     = "
        << gMax(mag(divSpeciesTotalCurrent.primitiveField()))
        << " A/m3" << nl
        << "  max|div(noAdv - phi)|       = "
        << gMax(mag(divNoAdvMinusPhi.primitiveField()))
        << " A/m3" << nl
        << "  max|div(total - phi)|       = "
        << gMax(mag(divTotalMinusPhi.primitiveField()))
        << " A/m3" << nl;

    // Identify the cell with the largest species-total vs phi-current
    // divergence difference and report the complete decomposition there.
    {
        const scalarField& d = divTotalMinusPhi.primitiveField();

        scalar localMax = -GREAT;
        label localCell = -1;

        forAll(d, celli)
        {
            if (mag(d[celli]) > localMax)
            {
                localMax = mag(d[celli]);
                localCell = celli;
            }
        }

        if (localCell >= 0)
        {
            Info<< "  same-state worst local cell = " << localCell
                << ", centre = " << mesh.C()[localCell] << nl
                << "    div(phi)     = " << currentResidual[localCell]
                << " A/m3" << nl
                << "    div(adv)     = " << divSpeciesAdvCurrent[localCell]
                << " A/m3" << nl
                << "    div(diff)    = " << divSpeciesDiffCurrent[localCell]
                << " A/m3" << nl
                << "    div(mig)     = " << divSpeciesMigCurrent[localCell]
                << " A/m3" << nl
                << "    div(noAdv)   = " << divSpeciesNoAdvCurrent[localCell]
                << " A/m3" << nl
                << "    div(total)   = " << divSpeciesTotalCurrent[localCell]
                << " A/m3" << nl
                << "    noAdv-phi    = " << divNoAdvMinusPhi[localCell]
                << " A/m3" << nl
                << "    total-phi    = " << divTotalMinusPhi[localCell]
                << " A/m3" << nl;
        }
    }

    // Patch-integrated current comparison.
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];

        if (p.size())
        {
            const scalar Iphi =
                gSum(ionicCurrentFlux.boundaryField()[patchi]);

            const scalar Iadv =
                gSum(speciesAdvCurrentFlux.boundaryField()[patchi]);

            const scalar Idiff =
                gSum(speciesDiffCurrentFlux.boundaryField()[patchi]);

            const scalar Imig =
                gSum(speciesMigCurrentFlux.boundaryField()[patchi]);

            const scalar Itotal =
                gSum(speciesTotalCurrentFlux.boundaryField()[patchi]);

            Info<< "  CURRENT-COMPARE patch " << p.name()
                << ": phi = " << Iphi
                << " A, adv = " << Iadv
                << " A, diff = " << Idiff
                << " A, mig = " << Imig
                << " A, total = " << Itotal
                << " A, total-phi = " << Itotal - Iphi
                << " A" << nl;
        }
    }

    // --------------------------------------------------------------------- //
    // Compare against the ACTUAL cEqn fvMatrix transport flux from before
    // equilibrium projection, when that field exists.  The initial phiEqn.H
    // call occurs before any cEqn.H call, so this must be registry-guarded.
    // --------------------------------------------------------------------- //

    if (mesh.foundObject<surfaceScalarField>("chargeMatrixTransportFlux"))
    {
        const surfaceScalarField& chargeMatrixTransportFlux =
            mesh.lookupObject<surfaceScalarField>
            (
                "chargeMatrixTransportFlux"
            );

        surfaceScalarField cEqnMatrixCurrentFlux
        (
            IOobject
            (
                "cEqnMatrixCurrentFlux",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Fdim*chargeMatrixTransportFlux
        );

        volScalarField divCEqnMatrixCurrent
        (
            IOobject
            (
                "divCEqnMatrixCurrent",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::div(cEqnMatrixCurrentFlux)
        );

        surfaceScalarField cEqnMatrixMinusPhiFlux
        (
            IOobject
            (
                "cEqnMatrixMinusPhiFlux",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            cEqnMatrixCurrentFlux - ionicCurrentFlux
        );

        volScalarField divCEqnMatrixMinusPhi
        (
            IOobject
            (
                "divCEqnMatrixMinusPhi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::div(cEqnMatrixMinusPhiFlux)
        );

        Info<< "CURRENT-COMPARE cEqn-matrix vs phi:" << nl
            << "  max|div(cEqn matrix current)| = "
            << gMax(mag(divCEqnMatrixCurrent.primitiveField()))
            << " A/m3" << nl
            << "  max|div(cEqn matrix - phi)|   = "
            << gMax(mag(divCEqnMatrixMinusPhi.primitiveField()))
            << " A/m3" << nl;
    }
}

