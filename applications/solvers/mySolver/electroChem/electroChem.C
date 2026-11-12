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
#include "OSspecific.H"
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>

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
    const bool profilingEnabled = profilingPrint || profilingLog;
    #include "createUfIfPresent.H"
    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const auto updatePotentialOutput = [&](const bool fieldWriteTime)
    // Read-only diagnostics; do not touch output-only fields between writes.
    {
        if (fieldWriteTime && (writeE || writeRelE))
        {
            const scalarField& charge = potentialMgr.chargeResidual();
            const scalarField& scale = potentialMgr.chargeScale();
            forAll(charge, celli)
            {
                if (writeE) E[celli] = charge[celli];
                if (writeRelE)
                    relE[celli] = mag(charge[celli])/max(scale[celli], SMALL);
            }
            if (writeE) E.correctBoundaryConditions();
            if (writeRelE) relE.correctBoundaryConditions();
        }
        if (printGradPhiE || (fieldWriteTime && writeGradPhiE))
        {
            gradPhiE = fvc::grad(phiE);
        }
        if (printKappa || printJDiff
         || (fieldWriteTime && (writeKappa || writeJDiff)))
        {
            potentialMgr.refreshDiagnosticCurrent();
            if (printKappa || (fieldWriteTime && writeKappa))
                kappaOutput = potentialMgr.conductivityField();
            if (printJDiff || (fieldWriteTime && writeJDiff))
                jDiffOutput = potentialMgr.diffusionCurrentField();
        }
    };

    // Initial diagnostic snapshot. Never alter the initial transport/potential state.
    potentialMgr.evaluateConstraint();
    // The first snapshot is explicitly written at t=0 when requested.
    updatePotentialOutput(runTime.value() == 0.0);
    if (runTime.value() == 0.0)
    {
        if (writeE) E.write();
        if (writeRelE) relE.write();
        if (writeGradPhiE) gradPhiE.write();
        if (writeKappa) kappaOutput.write();
        if (writeJDiff) jDiffOutput.write();
        // phiE already exists at t=0 as the solver's MUST_READ input.
    }

    // Profiling uses wall-clock elapsed time, not CPU time. Each rank measures
    // its own work; global maxima are reduced before master writes the log.
    // Profiling and Picard log switches are initialized in createFields.H.
    using ProfileClock = std::chrono::steady_clock;
    const auto profileNow = []() -> ProfileClock::time_point
    {
        return ProfileClock::now();
    };
    const auto profileSeconds = [](const ProfileClock::time_point& t) -> scalar
    {
        return std::chrono::duration_cast<std::chrono::duration<scalar>>
        (ProfileClock::now() - t).count();
    };
    // Timers are inclusive at the call-site; `other` is the uninstrumented
    // remainder of the outer iteration. Values are local until reduced.
    const label nProfileSections = 7;
    const char* profileNames[] =
    {"current", "kinetics", "transport", "constraints", "potential",
     "equilibrium", "reports"};
    scalarList profileTotals(nProfileSections, 0.0);
    labelList profileCalls(nProfileSections, 0);
    scalar profileOuterTotal = 0.0;
    label profileOuterCalls = 0;
    const auto profileAppend = [&](const fileName& path,
        const std::string& header, const std::string& line)
    {
        if (!Pstream::master()) return;
        std::ifstream previous(path.c_str());
        const bool empty = !previous.good()
            || previous.peek() == std::ifstream::traits_type::eof();
        previous.close();
        std::ofstream output(path.c_str(), std::ios::app);
        if (!output)
            FatalErrorInFunction << "Cannot open profiling file " << path
                << exit(FatalError);
        if (empty) output << header << '\n';
        output << line << '\n';
        output.flush();
        if (!output)
            FatalErrorInFunction << "Cannot write profiling file " << path
                << exit(FatalError);
    };
    if (profilingLog && Pstream::master())
        mkDir(runTime.globalPath()/"postProcessing"/"profiling");

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
        const auto profileStepStart = profileNow();
        scalarList profileStepSections(nProfileSections, 0.0);
        labelList profileStepCalls(nProfileSections, 0);
        scalar profileStepOuter = 0.0;


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

        // Initialize fast-equilibrium pseudo-sources from the selected guess.
        // Only fully accepted states are eligible for lastAccepted history.
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
            const auto profileOuterStart = profileNow();
            scalarList profileIterationSections(nProfileSections, 0.0);
            const auto profileRecord = [&](const label section,
                const ProfileClock::time_point& start)
            {
                if (profilingEnabled)
                {
                    const scalar dt = profileSeconds(start);
                    profileTotals[section] += dt;
                    profileStepSections[section] += dt;
                    profileIterationSections[section] += dt;
                    ++profileCalls[section];
                    ++profileStepCalls[section];
                }
            };


            Info<< "Picard state " << testedPicardStates
                << '/' << maxPicardIter << nl;

            // Solve electrode metal potentials at frozen phiE, then freeze
            // the shared face currents for the complete species sweep.
            const auto profileCurrentStart = profileNow();
            currentMgr.correct();
            profileRecord(0, profileCurrentStart);

            // One common kinetic state for the complete species sweep.
            const auto profileKineticStart = profileNow();
            reactionMgr.updateKineticRates();
            profileRecord(1, profileKineticStart);

            // A retained lambda is only an initial guess. Its first complete
            // transport must be feasible before any constraints are evaluated.
            const bool safeguardInitialGuess = picardI == 0
                && reactionMgr.usingAcceptedEquilibriumGuess();
            PtrList<volScalarField> transportStart
                (safeguardInitialGuess ? c.size() : 0);
            forAll(transportStart, i)
                transportStart.set
                (
                    i,
                    new volScalarField
                    (
                        IOobject(c[i].name(), runTime.timeName(), mesh,
                            IOobject::NO_READ, IOobject::NO_WRITE, false),
                        c[i]
                    )
                );

            scalar initialGuessScale = 1.0;
            label initialGuessBacktracks = 0;
            for (;;)
            {
                const bool permitInitialGuessRetry = safeguardInitialGuess
                    && initialGuessScale > 0;
                bool transportGuessRejected = false;

                // cEqn also supplies S_i = V/aP_i to reactionMgr. Current,
                // potential and kinetic rates stay frozen across retries.
                const auto profileTransportStart = profileNow();
                #include "cEqn.H"
                profileRecord(2, profileTransportStart);
                if (!transportGuessRejected) break;

                // Roll back EVERY species, including boundary values. Keep
                // physical oldTime fields intact; this is not a new timestep.
                forAll(transportStart, i)
                {
                    c[i].primitiveFieldRef() = transportStart[i].primitiveField();
                    forAll(c[i].boundaryField(), patchi)
                        c[i].boundaryFieldRef()[patchi]
                            == transportStart[i].boundaryField()[patchi];
                }
                // At most eight halvings, then retry the original zero guess.
                // If zero also fails, cEqn retains its original fatal check.
                initialGuessScale = ++initialGuessBacktracks <= 8
                    ? scalar(0.5)*initialGuessScale : scalar(0);
                reactionMgr.scaleInitialEquilibriumGuess(initialGuessScale);
                Info<< "    Retrying initial transport with lambda = "
                    << initialGuessScale << " * lastAccepted" << nl;
            }

            // Both constraints are evaluated on this SAME transported state.
            const auto profileConstraintStart = profileNow();
            reactionMgr.evaluateEquilibrium();
            potentialMgr.evaluateConstraint();
            if (currentMgr.requiresPotentialCoupling())
                potentialMgr.refreshConstraintCurrent();
            currentMgr.evaluateConstraint();
            profileRecord(3, profileConstraintStart);

            const bool currentOK = currentMgr.converged();
            const bool equilibriumOK = reactionMgr.equilibriumConverged();
            const bool electroneutralityOK = potentialMgr.converged();
            if (!electroneutralityOK) electroneutralityLastFailed = testedPicardStates;
            if (!equilibriumOK) equilibriumLastFailed = testedPicardStates;
            if (!currentOK) currentLastFailed = testedPicardStates;

            const bool finalState = picardI == maxPicardIter - 1;
            const bool constraintsOK = equilibriumOK && electroneutralityOK && currentOK;

            // Keep each constraint report adjacent to its applied correction.
            const auto profileReportStart = profileNow();
            currentMgr.report();
            potentialMgr.report();
            profileRecord(6, profileReportStart);
            if ((!electroneutralityOK || !currentOK) && !finalState)
            {
                const auto profilePotentialStart = profileNow();
                potentialMgr.correct(!currentOK);
                profileRecord(4, profilePotentialStart);
            }

            reactionMgr.reportEquilibrium();
            if (!equilibriumOK && !finalState)
            {
                const auto profileEquilibriumStart = profileNow();
                reactionMgr.correctEquilibrium();
                profileRecord(5, profileEquilibriumStart);
            }
            else
                reactionMgr.reportNoEquilibriumCorrection
                (equilibriumOK ? "equilibrium converged" : "maxIter reached");

            if (profilingEnabled)
            {
                const scalar outerDt = profileSeconds(profileOuterStart);
                profileOuterTotal += outerDt;
                profileStepOuter += outerDt;
                ++profileOuterCalls;
                scalar other = outerDt;
                forAll(profileIterationSections, si)
                    other -= profileIterationSections[si];
                // MPI collective: slowest rank determines elapsed work cost.
                scalarList times(profileIterationSections);
                forAll(times, si) reduce(times[si], maxOp<scalar>());
                scalar maxOuter = outerDt;
                scalar maxOther = other;
                reduce(maxOuter, maxOp<scalar>());
                reduce(maxOther, maxOp<scalar>());
                if (profilingLog && Pstream::master())
                {
                    std::ostringstream row;
                    row << std::setprecision(12) << runTime.value() << ' '
                        << testedPicardStates << ' ' << maxOuter;
                    forAll(times, si) row << ' ' << times[si];
                    row << ' ' << maxOther;
                    profileAppend
                    (
                        runTime.globalPath()/"postProcessing"/"profiling"/"iterations.dat",
                        "# time picard outerWall_s current_s kinetics_s transport_s constraints_s potential_s equilibrium_s reports_s other_s",
                        row.str()
                    );
                }
            }
            if (constraintsOK)
            {
                reactionMgr.acceptEquilibriumSources();
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

        updatePotentialOutput(runTime.writeTime());
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

        // Optional accepted-state statistics at the normal field-write times.
        // Use the global case path and master rank only, including in parallel.
        if (Pstream::master() && runTime.writeTime()
         && picardLog)
        {
            const fileName picardDir
                (runTime.globalPath()/"postProcessing"/"picard");
            mkDir(picardDir);
            const fileName picardPath(picardDir/"picard.dat");
            // Append across restarts/toggle changes; header only for empty files.
            std::ifstream existing(picardPath.c_str());
            const bool writeHeader =
                existing.peek() == std::ifstream::traits_type::eof();
            existing.close();
            std::ofstream picardData(picardPath.c_str(), std::ios::app);
            if (!picardData)
                FatalErrorInFunction << "Cannot open Picard log "
                    << picardPath << exit(FatalError);
            if (writeHeader)
                picardData << "# time EqPicardIter ENPicardIter "
                    << "BVPicardIter totalPicardIter runTime\n";
            // Match the summary: first sweep in the final uninterrupted passing
            // sequence. Inactive constraints use zero. Retry solves are excluded.
            picardData << runTime.value() << ' '
                << (reactionMgr.nEquilibriumReactions() ? equilibriumLastFailed + 1 : 0)
                << ' ' << (potentialMgr.active() ? electroneutralityLastFailed + 1 : 0)
                << ' ' << (currentMgr.requiresPotentialCoupling() ? currentLastFailed + 1 : 0)
                << ' ' << testedPicardStates
                << ' ' << runTime.elapsedCpuTime() << '\n';
            picardData.flush();
            if (!picardData)
                FatalErrorInFunction << "Cannot write Picard log "
                    << picardPath << exit(FatalError);
        }

        if (profilingEnabled)
        {
            scalar stepWall = profileSeconds(profileStepStart);
            scalar maxOuter = profileStepOuter;
            reduce(stepWall, maxOp<scalar>());
            reduce(maxOuter, maxOp<scalar>());
            scalarList maxSections(profileStepSections);
            forAll(maxSections, si)
                reduce(maxSections[si], maxOp<scalar>());
            if (Pstream::master())
            {
                std::ostringstream row;
                row << std::setprecision(12) << runTime.value() << ' '
                    << testedPicardStates << ' ' << stepWall << ' '
                    << maxOuter;
                forAll(maxSections, si) row << ' ' << maxSections[si];
                if (profilingLog) profileAppend
                (
                    runTime.globalPath()/"postProcessing"/"profiling"/"timesteps.dat",
                    "# time picard stepWall_s outerWall_s current_s kinetics_s transport_s constraints_s potential_s equilibrium_s reports_s",
                    row.str()
                );
                if (profilingPrint)
                {
                    Info<< "Profiling: timestep wall=" << stepWall
                        << " s, Picard wall=" << maxOuter << " s" << nl;
                    forAll(maxSections, si)
                        Info<< "    " << profileNames[si] << ": "
                            << maxSections[si] << " s (local calls "
                            << profileStepCalls[si] << ") " << nl;
                }
            }
        }
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

    if (profilingPrint)
    {
        Info<< "Profiling cumulative local work (seconds, calls):" << nl;
        forAll(profileTotals, si)
            Info<< "    " << profileNames[si] << ": " << profileTotals[si]
                << " s / " << profileCalls[si] << " calls" << nl;
        Info<< "    Picard: " << profileOuterTotal << " s / "
            << profileOuterCalls << " iterations" << nl;
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

