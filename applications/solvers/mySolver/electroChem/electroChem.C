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
        reactionMgr.resetEquilibriumSources();

        // Retained across cEqn.H include scope and all Picard sweeps.
        scalarList cPreMin(c.size(), 0.0);
        scalarList cPreMean(c.size(), 0.0);
        scalarList cPreMax(c.size(), 0.0);

        bool electroChemConverged = false;
        label testedPicardStates = 0;

        for (label picardI = 0; picardI < maxPicardIter; ++picardI)
        {
            testedPicardStates = picardI + 1;

            Info<< "Electrochemical Picard state " << testedPicardStates
                << '/' << maxPicardIter << nl;

            // One common kinetic state for the complete species sweep.
            reactionMgr.updateKineticRates();

            // Transport all active species using the current lambda and phiE.
            // cEqn also supplies S_i = V/aP_i to reactionMgr.
            #include "cEqn.H"

            // Both constraints are evaluated on this SAME transported state.
            reactionMgr.evaluateEquilibrium();
            potentialMgr.evaluateConstraint();

            const bool equilibriumOK = reactionMgr.equilibriumConverged();
            const bool electroneutralityOK = potentialMgr.converged();

            reactionMgr.reportEquilibrium();
            potentialMgr.report();

            if (equilibriumOK && electroneutralityOK)
            {
                electroChemConverged = true;
                break;
            }

            // maxIter counts transported-and-tested states.  Never calculate
            // a correction after the final state because its physical effect
            // could not be transported and tested in this timestep.
            if (picardI == maxPicardIter - 1)
            {
                break;
            }

            if (!equilibriumOK)
            {
                reactionMgr.correctEquilibrium();
            }

            if (!electroneutralityOK)
            {
                potentialMgr.correct();
            }

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

            FatalErrorInFunction
                << "Electrochemical Picard loop failed to converge after "
                << testedPicardStates << " transported-and-tested states at time "
                << runTime.timeName() << nl
                << "No untested final lambda or phiE correction was applied."
                << exit(FatalError);
        }

        updatePotentialOutput();
        if (printE || printRelE)
        {
            if (printE)
                Info<< "E min/max/maxAbs [mol/m3]: "
                    << gMin(E.primitiveField()) << " / "
                    << gMax(E.primitiveField()) << " / "
                    << gMax(mag(E.primitiveField())) << nl;
            if (printRelE)
                Info<< "relE max [-]: " << gMax(relE.primitiveField()) << nl;
        }

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

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

