#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "reactionManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    reactionManager reactionMgr(mesh, reactions);
    // Establish instantaneous equilibria before calculating initial conductivity and potential. 
    // The model assumes this projected initial state is already bulk electroneutral.
    reactionMgr.applyEquilibriumReactions();
    {
        const dictionary* zDictPtr = nullptr;

        if (reactions.found("chargeNumber"))
        {
            zDictPtr = &reactions.subDict("chargeNumber");
        }

        scalar localENAbs = 0.0;
        scalar localENRel = 0.0;

        forAll(mesh.V(), celli)
        {
            scalar chargeSum = 0.0;
            scalar chargeScale = enPhiENScale;

            forAll(c, i)
            {
                const scalar zi = zDictPtr
                    ? zDictPtr->getOrDefault<scalar>(species[i], 0.0)
                    : 0.0;

                chargeSum += zi*c[i][celli];
                chargeScale += mag(zi*c[i][celli]);
            }

            localENAbs = max(localENAbs, mag(chargeSum));
            localENRel = max(localENRel, mag(chargeSum)/chargeScale);
        }

        reduce(localENAbs, maxOp<scalar>());
        reduce(localENRel, maxOp<scalar>());

        if
        (
            localENAbs > enPhiENAbsTol
         && localENRel > enPhiENRelTol
        )
        {
            FatalErrorInFunction
                << "Initial equilibrium-projected concentrations are not "
                << "electroneutral." << nl
                << "max |sum(z_i c_i)| = " << localENAbs << " mol/m3" << nl
                << "max relative EN residual = " << localENRel << nl
                << "ENPhi requires an electroneutral initial state."
                << exit(FatalError);
        }
    }
    // Initialise the potential needed by the first Nernst-Planck transport solve.
    #include "phiEqn.H" 

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        // --- Electrochemistry ENPhi Picard corrector loop
        bool enPhiConverged = false;

        // Controlled Req-feedback test:
        // start each physical timestep from zero algebraic equilibrium source,
        // then build the source cumulatively over the ENPhi iterations.
        forAll(Req, i)
        {
            Req[i] =
                dimensionedScalar
                (
                    "zero",
                    dimMoles/dimVolume/dimTime,
                    0.0
                );
        }

        for (label enPhiIter = 1; enPhiIter <= enPhiMaxIter; ++enPhiIter)
        {
            Info<< "ENPhi: iteration " << enPhiIter << nl;

            List<scalarField> cPicardOld(c.size());

            forAll(c, i)
            {
                cPicardOld[i] = c[i].primitiveField();
            }

            const scalarField phiEPicardOld(phiE.primitiveField());

            #include "cEqn.H"
            #include "phiEqn.H"

            scalar cMaxAbsCorrection = 0.0;
            scalar cMaxRelCorrection = 0.0;

            forAll(c, i)
            {
                if (i == inertIndex)
                {
                    continue;
                }

                const scalarField& cNew = c[i].primitiveField();
                const scalarField& cOld = cPicardOld[i];

                forAll(cNew, celli)
                {
                    const scalar dc = mag(cNew[celli] - cOld[celli]);
                    const scalar scale = max
                    (
                        max(mag(cNew[celli]), mag(cOld[celli])),
                        enPhiCScale
                    );

                    cMaxAbsCorrection = max(cMaxAbsCorrection, dc);
                    cMaxRelCorrection = max(cMaxRelCorrection, dc/scale);
                }
            }

            reduce(cMaxAbsCorrection, maxOp<scalar>());
            reduce(cMaxRelCorrection, maxOp<scalar>());

            scalar phiMaxAbsCorrection = 0.0;
            scalar phiMaxRelCorrection = 0.0;

            const scalarField& phiENew = phiE.primitiveField();

            forAll(phiENew, celli)
            {
                const scalar dPhi =
                    mag(phiENew[celli] - phiEPicardOld[celli]);

                const scalar scale = max
                (
                    max
                    (
                        mag(phiENew[celli]),
                        mag(phiEPicardOld[celli])
                    ),
                    enPhiPhiScale
                );

                phiMaxAbsCorrection = max(phiMaxAbsCorrection, dPhi);
                phiMaxRelCorrection = max(phiMaxRelCorrection, dPhi/scale);
            }

            reduce(phiMaxAbsCorrection, maxOp<scalar>());
            reduce(phiMaxRelCorrection, maxOp<scalar>());

            const dictionary* zDictPtr = nullptr;

            if (reactions.found("chargeNumber"))
            {
                zDictPtr = &reactions.subDict("chargeNumber");
            }

            scalar enMaxAbs = 0.0;
            scalar enMaxRel = 0.0;

            forAll(mesh.V(), celli)
            {
                scalar chargeSum = 0.0;
                scalar chargeScale = enPhiENScale;

                forAll(c, i)
                {
                    const scalar zi = zDictPtr
                        ? zDictPtr->getOrDefault<scalar>(species[i], 0.0)
                        : 0.0;

                    chargeSum += zi*c[i][celli];
                    chargeScale += mag(zi*c[i][celli]);
                }

                enMaxAbs = max(enMaxAbs, mag(chargeSum));
                enMaxRel = max(enMaxRel, mag(chargeSum)/chargeScale);
            }

            reduce(enMaxAbs, maxOp<scalar>());
            reduce(enMaxRel, maxOp<scalar>());

            // ------------------------------------------------------------- //
            // FINAL-STATE DISCRETE CHARGE-BALANCE DIAGNOSTIC
            //
            // Evaluate, on the same final c and phiE state used by the ENPhi
            // convergence test,
            //
            //   (q^{n+1} - q^n)/dt
            // + div(sum_i z_i N_i)
            // - sum_i z_i Rc_i
            // - sum_i z_i Req_i
            //
            // where
            //
            //   q = sum_i z_i c_i
            //
            // and N_i is reconstructed with the same finite-volume operators
            // used by the species equations:
            //
            //   N_i.Sf =
            //       interpolate(c_i)*phi
            //     - D_i*snGrad(c_i)*|Sf|
            //     - interpolate(D_i*z_i*F/(R*T)*c_i)
            //       *snGrad(phiE)*|Sf|.
            //
            // For a fully consistent converged state this complete balance
            // should be close to the linear-solver residual level.
            // ------------------------------------------------------------- //

            volScalarField qFinalState
            (
                IOobject
                (
                    "qFinalState",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dimMoles/dimVolume,
                    0.0
                )
            );

            volScalarField qOldState
            (
                IOobject
                (
                    "qOldState",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dimMoles/dimVolume,
                    0.0
                )
            );

            volScalarField chargeRcFinalState
            (
                IOobject
                (
                    "chargeRcFinalState",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dimMoles/dimVolume/dimTime,
                    0.0
                )
            );

            volScalarField chargeReqFinalState
            (
                IOobject
                (
                    "chargeReqFinalState",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dimMoles/dimVolume/dimTime,
                    0.0
                )
            );

            surfaceScalarField chargeMolarFluxFinalState
            (
                IOobject
                (
                    "chargeMolarFluxFinalState",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dimMoles/dimTime,
                    0.0
                )
            );

            forAll(c, i)
            {
                const scalar zi = zDictPtr
                    ? zDictPtr->getOrDefault<scalar>(species[i], 0.0)
                    : 0.0;

                if (mag(zi) <= SMALL)
                {
                    continue;
                }

                qFinalState += zi*c[i];
                qOldState += zi*c[i].oldTime();

                chargeRcFinalState += zi*Rc[i];
                chargeReqFinalState += zi*Req[i];

                dimensionedScalar Di
                (
                    "D_" + species[i],
                    dimArea/dimTime,
                    diffusivityDict.get<scalar>(species[i])
                );

                surfaceScalarField migrationFluxFinalState
                (
                    IOobject
                    (
                        "migrationFluxFinalState_" + species[i],
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

                chargeMolarFluxFinalState +=
                    zi*
                    (
                        fvc::interpolate(c[i])*phi
                      - Di*fvc::snGrad(c[i])*mesh.magSf()
                      - migrationFluxFinalState
                    );
            }

            volScalarField qDdtFinalState
            (
                IOobject
                (
                    "qDdtFinalState",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                (qFinalState - qOldState)/runTime.deltaT()
            );

            volScalarField divChargeMolarFluxFinalState
            (
                IOobject
                (
                    "divChargeMolarFluxFinalState",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                fvc::div(chargeMolarFluxFinalState)
            );

            volScalarField fullChargeBalanceFinalState
            (
                IOobject
                (
                    "fullChargeBalanceFinalState",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                qDdtFinalState
              + divChargeMolarFluxFinalState
              - chargeRcFinalState
              - chargeReqFinalState
            );

            label localWorstENCell = -1;
            scalar localWorstEN = -GREAT;

            forAll(qFinalState.primitiveField(), celli)
            {
                const scalar magQ = mag(qFinalState[celli]);

                if (magQ > localWorstEN)
                {
                    localWorstEN = magQ;
                    localWorstENCell = celli;
                }
            }

            Info<< "FINAL-STATE CHARGE BALANCE summary:" << nl
                << "  max|qNow|                       = "
                << gMax(mag(qFinalState.primitiveField()))
                << " mol/m3" << nl
                << "  max|qOld|                       = "
                << gMax(mag(qOldState.primitiveField()))
                << " mol/m3" << nl
                << "  max|dq/dt|                      = "
                << gMax(mag(qDdtFinalState.primitiveField()))
                << " mol/m3/s" << nl
                << "  max|div(sum z_i N_i)|           = "
                << gMax(mag(divChargeMolarFluxFinalState.primitiveField()))
                << " mol/m3/s" << nl
                << "  max|sum z_i Rc_i|               = "
                << gMax(mag(chargeRcFinalState.primitiveField()))
                << " mol/m3/s" << nl
                << "  max|sum z_i Req_i|              = "
                << gMax(mag(chargeReqFinalState.primitiveField()))
                << " mol/m3/s" << nl
                << "  max|FULL charge balance|        = "
                << gMax(mag(fullChargeBalanceFinalState.primitiveField()))
                << " mol/m3/s" << nl;

            if (localWorstENCell >= 0)
            {
                const label celli = localWorstENCell;

                Info<< "FINAL-STATE CHARGE BALANCE worst-EN local cell = "
                    << celli
                    << ", centre = " << mesh.C()[celli] << nl
                    << "  qNow                           = "
                    << qFinalState[celli] << " mol/m3" << nl
                    << "  qOld                           = "
                    << qOldState[celli] << " mol/m3" << nl
                    << "  (qNow-qOld)/dt                 = "
                    << qDdtFinalState[celli] << " mol/m3/s" << nl
                    << "  div(sum z_i N_i)               = "
                    << divChargeMolarFluxFinalState[celli]
                    << " mol/m3/s" << nl
                    << "  sum(z_i Rc_i)                  = "
                    << chargeRcFinalState[celli] << " mol/m3/s" << nl
                    << "  sum(z_i Req_i)                 = "
                    << chargeReqFinalState[celli] << " mol/m3/s" << nl
                    << "  FULL balance                   = "
                    << fullChargeBalanceFinalState[celli]
                    << " mol/m3/s" << nl;
            }

            Info<< "electroneutrality errors : max molar = "
                << enMaxAbs << " mol/m3, max relative = "
                << enMaxRel << nl
                << "current conservation errors : sum local = "
                << enPhiCurrentLocalError << " A, global = "
                << enPhiCurrentGlobalError << " A" << nl
                << "ENPhi corrections : c max abs = "
                << cMaxAbsCorrection << " mol/m3, c max rel = "
                << cMaxRelCorrection << ", phiE max abs = "
                << phiMaxAbsCorrection << " V, phiE max rel = "
                << phiMaxRelCorrection << nl
                << "equilibrium projection correction : max abs = "
                << eqProjectionMaxAbsCorrection << " mol/m3" << nl
                << "Req update : max abs = "
                << reqUpdateMaxAbs << " mol/m3/s, max |Req| = "
                << reqMagnitudeMaxAbs << " mol/m3/s, max charge source = "
                << reqChargeMaxAbs << " mol/m3/s" << nl;

            const bool cConverged =
                cMaxAbsCorrection <= enPhiCAbsTol
             || cMaxRelCorrection <= enPhiCRelTol;

            const bool phiConverged =
                phiMaxAbsCorrection <= enPhiPhiAbsTol
             || phiMaxRelCorrection <= enPhiPhiRelTol;

            const bool enConverged =
                enMaxAbs <= enPhiENAbsTol
             || enMaxRel <= enPhiENRelTol;

            // In the cumulative Req formulation, a proper fixed point also
            // requires the equilibrium projection to become negligible.
            const bool equilibriumProjectionConverged =
                eqProjectionMaxAbsCorrection <= enPhiCAbsTol;

            if
            (
                enPhiIter >= enPhiMinIter
             && cConverged
             && phiConverged
             && enConverged
             && equilibriumProjectionConverged
            )
            {
                enPhiConverged = true;
                Info<< "ENPhi: converged in " << enPhiIter
                    << " iteration(s)" << nl << endl;
                break;
            }

            Info<< endl;
        }

        if (!enPhiConverged)
        {
            FatalErrorInFunction
                << "ENPhi Picard corrector failed to converge in "
                << enPhiMaxIter << " iteration(s)."
                << exit(FatalError);
        }

        #include "YEqn.H"
        #include "diagnostics.H"

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

