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
#include "potentialManager.H"
#include "currentManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient reactive electrochemical transport with independently selectable"
        " potential and fast-equilibrium models on a moving mesh."
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
    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const auto updatePotentialOutput = [&]()
    // Read-only projection of the potential model's authoritative constraint.
    {
        const scalarField& charge = potentialMgr.chargeResidual();
        const scalarField& scale = potentialMgr.chargeScale();
        forAll(charge, celli)
        {
            E[celli] = charge[celli];
            relE[celli] = mag(charge[celli])/max(scale[celli], SMALL);
        }
        E.correctBoundaryConditions();
        relE.correctBoundaryConditions();
        if (writeGradPhiE || printGradPhiE)
        {
            gradPhiE = fvc::grad(phiE);
        }
        if (writeKappa || writeJDiff || printKappa || printJDiff)
        {
            potentialMgr.refreshDiagnosticCurrent();
            if (writeKappa || printKappa)
                kappaOutput = potentialMgr.conductivityField();
            if (writeJDiff || printJDiff)
                jDiffOutput = potentialMgr.diffusionCurrentField();
        }
    };

    // Initial diagnostic snapshot. Never alter the initial transport/potential state.
    potentialMgr.evaluateConstraint();
    updatePotentialOutput();
    if (runTime.value() == 0.0)
    {
        if (writeE) E.write();
        if (writeRelE) relE.write();
        if (writeGradPhiE) gradPhiE.write();
        if (writeKappa) kappaOutput.write();
        if (writeJDiff) jDiffOutput.write();
        // phiE already exists at t=0 as the solver's MUST_READ input.
    }

    Info<< "\nStarting time loop\n" << endl;

    label picardTimeSteps = 0;
    label picardMaxUsed = 0;
    scalar picardTotalUsed = 0.0;
    label electroneutralityMaxUsed = 0;
    scalar electroneutralityTotalUsed = 0.0;
    label equilibriumMaxUsed = 0;
    scalar equilibriumTotalUsed = 0.0;
    label currentMaxUsed = 0;
    scalar currentTotalUsed = 0.0;

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

        // --- Coupled electrochemical Picard loop
        const dictionary& picardDict = mesh.solutionDict().subDict("Picard");
        const label maxPicardIter = picardDict.get<label>("maxIter");
        const dictionary& transportPicard = picardDict.subDict("transport");
        const scalar negativeConcentrationTolerance = transportPicard.get<scalar>("negativeConcentrationTolerance");

        if (negativeConcentrationTolerance < 0.0)
        {
            FatalIOErrorInFunction(transportPicard)
                << "negativeConcentrationTolerance must be >= 0"
                << exit(FatalIOError);
        }

        if (maxPicardIter <= 0)
        {
            FatalIOErrorInFunction(picardDict)
                << "Picard/maxIter must be positive, found "
                << maxPicardIter
                << exit(FatalIOError);
        }

        // Fast-equilibrium pseudo-sources are timestep-local iterative state.
        // phiE is deliberately NOT reset: it starts from the previously
        // converged physical-time value read/retained by the solver.
        reactionMgr.readOutputControls();
        reactionMgr.resetEquilibriumSources();

        // Retained across cEqn.H include scope and all Picard sweeps.
        scalarList cPreMin(c.size(), 0.0);
        scalarList cPreMean(c.size(), 0.0);
        scalarList cPreMax(c.size(), 0.0);

        bool electroChemConverged = false;
        label testedPicardStates = 0;
        label electroneutralityLastFailed = 0;
        label equilibriumLastFailed = 0;
        label currentLastFailed = 0;

        for (label picardI = 0; picardI < maxPicardIter; ++picardI)
        {
            testedPicardStates = picardI + 1;

            Info<< "Picard state " << testedPicardStates
                << '/' << maxPicardIter << nl;

            // Solve electrode metal potentials at frozen phiE, then freeze
            // the shared face currents for the complete species sweep.
            currentMgr.correct();

            // One common kinetic state for the complete species sweep.
            reactionMgr.updateKineticRates();

            // Transport all active species using the current lambda and phiE.
            // cEqn also supplies S_i = V/aP_i to reactionMgr.
            #include "cEqn.H"

            // Both constraints are evaluated on this SAME transported state.
            reactionMgr.evaluateEquilibrium();
            potentialMgr.evaluateConstraint();
            if (currentMgr.requiresPotentialCoupling())
                potentialMgr.refreshDiagnosticCurrent();
            currentMgr.evaluateConstraint();

            const bool currentOK = currentMgr.converged();
            const bool equilibriumOK = reactionMgr.equilibriumConverged();
            const bool electroneutralityOK = potentialMgr.converged();
            if (!electroneutralityOK) electroneutralityLastFailed = testedPicardStates;
            if (!equilibriumOK) equilibriumLastFailed = testedPicardStates;
            if (!currentOK) currentLastFailed = testedPicardStates;

            const bool finalState = picardI == maxPicardIter - 1;
            const bool constraintsOK = equilibriumOK && electroneutralityOK && currentOK;

            // Keep each constraint report adjacent to its applied correction.
            currentMgr.report();
            potentialMgr.report();
            if ((!electroneutralityOK || !currentOK) && !finalState)
                potentialMgr.correct(!currentOK);

            reactionMgr.reportEquilibrium();
            if (!equilibriumOK && !finalState)
                reactionMgr.correctEquilibrium();
            else
                reactionMgr.reportNoEquilibriumCorrection
                (equilibriumOK ? "equilibrium converged" : "maxIter reached");

            if (constraintsOK)
            {
                electroChemConverged = true;
                break;
            }
            // Never apply a correction without a following tested transport.
            if (finalState) break;

            // Any lambda/phiE correction is intentionally followed by the
            // next complete cEqn transport before either constraint is tested.
        }

        if (!electroChemConverged)
        {
            Info<< nl
                << "Final electrochemical constraint state after "
                << testedPicardStates << " transported states:" << nl;
            reactionMgr.reportEquilibrium();
            potentialMgr.report();
            currentMgr.report();

            FatalErrorInFunction
                << "Electrochemical Picard loop failed to converge after "
                << testedPicardStates << " transported-and-tested states at time "
                << runTime.timeName() << nl
                << "No untested final lambda or phiE correction was applied."
                << exit(FatalError);
        }

        ++picardTimeSteps;
        picardMaxUsed = max(picardMaxUsed, testedPicardStates);
        picardTotalUsed += testedPicardStates;
        // First state in the final uninterrupted sequence of passing states.
        // A constraint passing throughout the timestep uses one tested state.
        if (potentialMgr.active())
        {
            const label used = electroneutralityLastFailed + 1;
            electroneutralityMaxUsed = max(electroneutralityMaxUsed, used);
            electroneutralityTotalUsed += used;
        }
        if (reactionMgr.nEquilibriumReactions())
        {
            const label used = equilibriumLastFailed + 1;
            equilibriumMaxUsed = max(equilibriumMaxUsed, used);
            equilibriumTotalUsed += used;
        }

        if (currentMgr.requiresPotentialCoupling())
        {
            const label used = currentLastFailed + 1;
            currentMaxUsed = max(currentMaxUsed, used);
            currentTotalUsed += used;
        }

        updatePotentialOutput();
        // Accepted-state internal-cell statistics; means are volume weighted.
        const auto printScalarStats = [&](const word& name,
            const scalarField& values, const char* units)
        {
            Info<< "  " << name << " min/mean/max [" << units << "]: "
                << gMin(values) << " / "
                << gSum(values*mesh.V())/gSum(mesh.V()) << " / "
                << gMax(values) << nl;
        };
        const auto printVectorStats = [&](const word& name,
            const vectorField& values, const boolList& directions, const char* units)
        {
            Info<< "  " << name << " components:" << nl;
            for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
            {
                if (!directions[cmpt]) continue;
                const word componentName(name + "_" + word(cmpt == 0 ? "x" : cmpt == 1 ? "y" : "z"));
                const scalarField componentValues(values.component(cmpt));
                printScalarStats(componentName, componentValues, units);
            }
        };
        bool printElectrodePotentialStats = false;
        if (runTime.controlDict().found("potential"))
        {
            const dictionary& d = runTime.controlDict().subDict("potential");
            if (d.found("print") && currentMgr.requiresPotentialCoupling())
            {
                const dictionary& p = d.subDict("print");
                printElectrodePotentialStats =
                    p.getOrDefault<Switch>("phiS", false)
                 || p.getOrDefault<Switch>("Eeq", false)
                 || p.getOrDefault<Switch>("eta", false)
                 || p.getOrDefault<Switch>("cellVoltage", false)
                 || p.getOrDefault<Switch>("voltageBudget", false)
                 || p.getOrDefault<Switch>("Eref", false)
                 || p.getOrDefault<Switch>("etaRef", false);
            }
        }
        if (printPhiE || printGradPhiE || printKappa || printJDiff
         || printElectrodePotentialStats)
        {
            Info<< "Potential field statistics:" << nl;
            if (printPhiE) printScalarStats("phiE", phiE.primitiveField(), "V");
            if (printGradPhiE) printVectorStats("gradPhiE", gradPhiE.primitiveField(), printGradPhiEDirections, "V/m");
            if (printKappa) printScalarStats("kappa", kappaOutput.primitiveField(), "S/m");
            if (printJDiff) printVectorStats("jDiff", jDiffOutput.primitiveField(), printJDiffDirections, "A/m2");
        }

        currentMgr.writeDiagnostics();

        // Optional summary of the final transport state before clipping.
        if (printCFieldsBeforeClipping)
        {
            Info<< "Concentrations before clipping (final Picard sweep):" << nl;
            forAll(c, i)
            {
                if (i == inertIndex) continue; // inert species is reconstructed later
                Info<< "  " << c[i].name() << " min/mean/max = "
                    << cPreMin[i] << " / " << cPreMean[i]
                    << " / " << cPreMax[i] << " mol/m3" << nl;
            }
        }
        #include "YEqn.H"
        #include "diagnostics.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "Picard iteration summary:" << nl
        << "    completed timesteps: " << picardTimeSteps << nl
        << "    max iterations used: " << picardMaxUsed << nl;
    if (picardTimeSteps)
        Info<< "    mean iterations used: "
            << picardTotalUsed/scalar(picardTimeSteps) << nl;
    else Info<< "    mean iterations used: n/a (no completed timesteps)" << nl;

    const auto reportConstraintIterations = [&](const char* name,
        const bool active, const label maxUsed, const scalar totalUsed)
    {
        Info<< name << ':' << nl;
        if (!active) Info<< "    inactive" << nl;
        else if (!picardTimeSteps) Info<< "    no completed timesteps" << nl;
        else Info<< "    max iterations used: " << maxUsed << nl
            << "    mean iterations used: " << totalUsed/scalar(picardTimeSteps) << nl;
    };
    reportConstraintIterations("electroneutrality", potentialMgr.active(),
        electroneutralityMaxUsed, electroneutralityTotalUsed);
    reportConstraintIterations("equilibrium", reactionMgr.nEquilibriumReactions() > 0,
        equilibriumMaxUsed, equilibriumTotalUsed);
    reportConstraintIterations("ButlerVolmer", currentMgr.requiresPotentialCoupling(),
        currentMaxUsed, currentTotalUsed);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

