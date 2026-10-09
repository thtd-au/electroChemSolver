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
        const dictionary diagnosticMode = picardDict.found("diagnosticMode")
            ? picardDict.subDict("diagnosticMode") : dictionary();
        const bool forceFixedIterations =
            diagnosticMode.getOrDefault<bool>("forceFixedIterations", false);
        const bool writeChargeFields1I =
            diagnosticMode.getOrDefault<bool>("writeChargeFields", false);
        const bool runFullDiagnostics1I =
            diagnosticMode.getOrDefault<bool>("runFullDiagnostics", !forceFixedIterations);
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
        scalarList cPreMin1I(c.size(), 0.0);
        scalarList cPreMean1I(c.size(), 0.0);
        scalarList cPreMax1I(c.size(), 0.0);

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

            if (!forceFixedIterations && equilibriumOK && electroneutralityOK)
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

        if (!electroChemConverged && !forceFixedIterations)
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

        if (forceFixedIterations)
        {
            Info<< "TEST1I FORCED STATE (NOT CONVERGED/VALIDATED): time="
                << runTime.timeName() << ", sweeps=" << testedPicardStates
                << ", equilibriumOK=" << reactionMgr.equilibriumConverged()
                << ", electroneutralityOK=" << potentialMgr.converged()
                << nl;
        }

        // Optional spatial charge diagnostics, evaluated on the final transported state.
        if (writeChargeFields1I)
        {
            const dictionary* zDict1I = reactions.found("chargeNumber")
                ? &reactions.subDict("chargeNumber") : nullptr;
            molarCharge = dimensionedScalar
            (
                "zero", molarCharge.dimensions(), 0.0
            );
            volScalarField relativeCharge1I
            (
                IOobject("relativeCharge1I", runTime.timeName(), mesh,
                    IOobject::NO_READ, IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            );
            scalarField ionicScale1I(mesh.nCells(), 0.0);
            forAll(c, i)
            {
                const scalar zi = zDict1I
                    ? zDict1I->getOrDefault<scalar>(species[i], 0.0) : 0.0;
                if (mag(zi) <= SMALL) continue;
                molarCharge += zi*c[i];
                forAll(ionicScale1I, celli)
                {
                    ionicScale1I[celli] += mag(zi)*c[i][celli];
                }
            }
            scalar maxAbs1I = 0.0;
            scalar maxRel1I = 0.0;
            scalar integratedAbs1I = 0.0;
            scalar integratedSigned1I = 0.0;
            forAll(ionicScale1I, celli)
            {
                const scalar E = molarCharge[celli];
                const scalar rel = mag(E)/max(ionicScale1I[celli], SMALL);
                relativeCharge1I[celli] = rel;
                maxAbs1I = max(maxAbs1I, mag(E));
                maxRel1I = max(maxRel1I, rel);
                integratedAbs1I += mag(E)*mesh.V()[celli];
                integratedSigned1I += E*mesh.V()[celli];
            }
            reduce(maxAbs1I, maxOp<scalar>());
            reduce(maxRel1I, maxOp<scalar>());
            reduce(integratedAbs1I, sumOp<scalar>());
            reduce(integratedSigned1I, sumOp<scalar>());
            Info<< "TEST1I CHARGE: maxAbs=" << maxAbs1I
                << " mol/m3, maxRel=" << maxRel1I
                << ", integralAbs=" << integratedAbs1I
                << " mol, integralSigned=" << integratedSigned1I << " mol" << nl;
            if (runTime.writeTime())
            {
                molarCharge.write();
                relativeCharge1I.write();
            }
        }

        // TEST1I: configurable species summaries for the FINAL transported state.
        // Print even in forced mode, without enabling the full diagnostics.H.
        if (printCFieldsBeforeClipping)
        {
            Info<< "TEST1I concentrations before clipping (final Picard sweep):" << nl;
            forAll(c, i)
            {
                if (i == inertIndex) continue; // inert species is reconstructed later
                Info<< "  " << c[i].name() << " min/mean/max = "
                    << cPreMin1I[i] << " / " << cPreMean1I[i]
                    << " / " << cPreMax1I[i] << " mol/m3" << nl;
            }
        }
        if (printCFields)
        {
            Info<< "TEST1I concentrations after clipping (final Picard sweep):" << nl;
            const scalar volume1I = gSum(mesh.V());
            forAll(c, i)
            {
                if (i == inertIndex) continue;
                Info<< "  " << c[i].name() << " min/mean/max = "
                    << gMin(c[i].primitiveField()) << " / "
                    << gSum(c[i].primitiveField()*mesh.V())/volume1I << " / "
                    << gMax(c[i].primitiveField()) << " mol/m3" << nl;
            }
        }

        #include "YEqn.H"

        if (printYFields)
        {
            Info<< "TEST1I mass fractions (final Picard sweep):" << nl;
            const scalar volume1I = gSum(mesh.V());
            forAll(Y, i)
            {
                Info<< "  " << Y[i].name() << " min/mean/max = "
                    << gMin(Y[i].primitiveField()) << " / "
                    << gSum(Y[i].primitiveField()*mesh.V())/volume1I << " / "
                    << gMax(Y[i].primitiveField()) << nl;
            }
        }

        if (runFullDiagnostics1I)
        {
            #include "diagnostics.H"
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

