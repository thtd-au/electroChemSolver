// cEqn.H
//
// Diagnostic version:
//   - keeps the original species transport formulation
//   - restores the pre-equilibrium concentration snapshot
//   - retains matrix-residual / fvOptions checks
//   - adds an unrelaxed fvMatrix::flux() vs manual face-flux comparison
//
// For this diagnostic run, use the ORIGINAL implicit
// galvanostaticMolarFlux boundary condition.
//
// Also set in system/fvSchemes:
//
// fluxRequired
// {
//     default         yes;
// }
//

label inertIndex = species.find(inertSpecie);

const dictionary& diffusivityDict = reactions.subDict("diffusivity");
const bool hasChargeNumber = reactions.found("chargeNumber");

const dictionary* cEqnZDictPtr = nullptr;

if (hasChargeNumber)
{
    cEqnZDictPtr = &reactions.subDict("chargeNumber");
}


// Face-integrated potential-gradient quantity used by migration:
//
//     snGrad(phiE) * |Sf|
//
// dimensions: V*m
surfaceScalarField gradPhiEFlux
(
    IOobject
    (
        "gradPhiEFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::snGrad(phiE)*mesh.magSf()
);


// ------------------------------------------------------------------------- //
// Kinetic source fields
// ------------------------------------------------------------------------- //

forAll(Rc, i)
{
    tmp<volScalarField> tSource
    (
        reactionMgr.kineticSource(species[i])
    );

    Rc[i] = tSource();
}


const scalar domainVolume = gSum(mesh.V());


// Initialise concentration diagnostics for all species, including inert.
// Active species are overwritten after their transport solves.
forAll(c, i)
{
    cMinBeforeClipping[i] =
        gMin(c[i].primitiveField());

    cMeanBeforeClipping[i] =
        gSum(c[i].primitiveField()*mesh.V())/domainVolume;

    cMaxBeforeClipping[i] =
        gMax(c[i].primitiveField());
}


// ------------------------------------------------------------------------- //
// Charge-weighted diagnostics for the ACTUAL species matrices
// ------------------------------------------------------------------------- //

// Charge-weighted algebraic residual of the relaxed/constrained matrices.
// Units after division by V: mol/m3/s.
volScalarField chargeMatrixResidual
(
    IOobject
    (
        "chargeMatrixResidual",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimMoles/dimVolume/dimTime,
        0.0
    )
);


// Charge-weighted transport flux extracted from the UNRELAXED physical
// fvMatrix, with the explicit migration flux added using the equation sign.
//
// Species equation:
//
//   ddt(c) + div(phi,c)
//     = laplacian(D,c) + div(migrationFlux) + sources
//
// so the physical outward Nernst-Planck face flux is
//
//   matrixSpatialFlux - migrationFlux
//
surfaceScalarField chargeMatrixTransportFlux
(
    IOobject
    (
        "chargeMatrixTransportFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimMoles/dimTime,
        0.0
    )
);


// Manual reconstruction used for comparison:
//
//   phi*c_f - D*snGrad(c)*|Sf| - migrationFlux
//
surfaceScalarField chargeManualTransportFlux
(
    IOobject
    (
        "chargeManualTransportFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimMoles/dimTime,
        0.0
    )
);


// ------------------------------------------------------------------------- //
// Solve concentration transport
// ------------------------------------------------------------------------- //

forAll(c, i)
{
    if (i != inertIndex)
    {
        volScalarField& ci = c[i];

        dimensionedScalar Di
        (
            "D_" + species[i],
            dimArea/dimTime,
            diffusivityDict.get<scalar>(species[i])
        );

        scalar zi = 0.0;

        if (cEqnZDictPtr)
        {
            zi =
                cEqnZDictPtr->getOrDefault<scalar>
                (
                    species[i],
                    0.0
                );
        }


        // Explicit electro-migration contribution.
        //
        // NOTE:
        // This diagnostic version intentionally uses the normal ci field,
        // i.e. the original solver formulation.  Do not use the temporary
        // Picard-lagged BC experiment for this run.
        surfaceScalarField migrationFlux
        (
            IOobject
            (
                "migrationFlux_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate
            (
                Di*zi*Fdim/(Rgas*T)*ci
            )
           *gradPhiEFlux
        );


        fvScalarMatrix ciEqn
        (
            fvm::ddt(ci)
          + fvm::div(phi, ci)
         ==
            fvm::laplacian(Di, ci)
          + fvc::div(migrationFlux)
          + Rc[i]
          + Req[i]   // lagged/cumulative algebraic equilibrium source
          + fvOptions(ci)
        );


        // Preserve the equation BEFORE under-relaxation and fvOptions
        // constraints.  This lets fvMatrix::flux() represent the original
        // discretised transport operator rather than the relaxed matrix.
        fvScalarMatrix ciEqnPhysical(ciEqn);


        ciEqn.relax();
        fvOptions.constrain(ciEqn);

        ciEqn.solve(ci.name());


        // ----------------------------------------------------------------- //
        // Matrix diagnostics immediately after solve
        // ----------------------------------------------------------------- //

        if (mag(zi) > SMALL)
        {
            // Residual of the matrix that was actually solved
            // (relaxed + constrained).
            tmp<scalarField> tMatrixResidual =
                ciEqn.residual();

            const scalarField& matrixResidual =
                tMatrixResidual();

            chargeMatrixResidual.primitiveFieldRef() +=
                zi*matrixResidual/mesh.V();

            Info<< "  matrix residual "
                << species[i]
                << ": max = "
                << gMax(mag(matrixResidual/mesh.V()))
                << " mol/m3/s"
                << nl;


            // ------------------------------------------------------------- //
            // Exact fvMatrix transport-flux diagnostic
            //
            // ciEqnPhysical.flux() represents the spatial face flux of the
            // unrelaxed implicit matrix:
            //
            //   convection - diffusion
            //
            // Explicit migration is on the RHS of ciEqn, so subtract it to
            // reconstruct the complete outward species transport flux.
            // ------------------------------------------------------------- //

            tmp<surfaceScalarField> tImplicitMatrixFlux =
                ciEqnPhysical.flux();

            const surfaceScalarField& implicitMatrixFlux =
                tImplicitMatrixFlux();

            surfaceScalarField matrixTotalFlux
            (
                IOobject
                (
                    "matrixTotalFlux_" + species[i],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                implicitMatrixFlux - migrationFlux
            );


            // Manual reconstruction evaluated from the just-solved ci.
            surfaceScalarField manualTotalFlux
            (
                IOobject
                (
                    "manualTotalFlux_" + species[i],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                fvc::interpolate(ci)*phi
              - Di*fvc::snGrad(ci)*mesh.magSf()
              - migrationFlux
            );


            chargeMatrixTransportFlux +=
                zi*matrixTotalFlux;

            chargeManualTransportFlux +=
                zi*manualTotalFlux;


            surfaceScalarField fluxDifference
            (
                IOobject
                (
                    "matrixMinusManualFlux_" + species[i],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                matrixTotalFlux - manualTotalFlux
            );


            Info<< "  FV-FLUX "
                << species[i]
                << ": max|matrix-manual| = "
                << gMax(mag(fluxDifference.primitiveField()))
                << " mol/s"
                << nl;


            // Patch-integrated species fluxes.
            // These are especially important at cathode/anode.
            forAll(mesh.boundary(), patchi)
            {
                const fvPatch& p =
                    mesh.boundary()[patchi];

                if (p.size())
                {
                    const scalar matrixPatch =
                        gSum
                        (
                            matrixTotalFlux.boundaryField()[patchi]
                        );

                    const scalar manualPatch =
                        gSum
                        (
                            manualTotalFlux.boundaryField()[patchi]
                        );

                    Info<< "    patch "
                        << p.name()
                        << ": matrix = "
                        << matrixPatch
                        << " mol/s, manual = "
                        << manualPatch
                        << " mol/s, diff = "
                        << matrixPatch - manualPatch
                        << " mol/s"
                        << nl;
                }
            }
        }


        // Check whether fvOptions.correct() modifies the solved field.
        const scalarField cBeforeFvOptionsCorrect
        (
            ci.primitiveField()
        );

        fvOptions.correct(ci);

        scalar fvOptionsCorrectMaxChange = 0.0;

        forAll(ci.primitiveField(), celli)
        {
            fvOptionsCorrectMaxChange =
                max
                (
                    fvOptionsCorrectMaxChange,
                    mag
                    (
                        ci[celli]
                      - cBeforeFvOptionsCorrect[celli]
                    )
                );
        }

        reduce
        (
            fvOptionsCorrectMaxChange,
            maxOp<scalar>()
        );

        if (mag(zi) > SMALL)
        {
            Info<< "  fvOptions.correct change "
                << species[i]
                << ": max = "
                << fvOptionsCorrectMaxChange
                << " mol/m3"
                << nl;
        }


        // Save raw transport solution before positivity clipping.
        if (printCFieldsBeforeClipping)
        {
            cMinBeforeClipping[i] =
                gMin(ci.primitiveField());

            cMeanBeforeClipping[i] =
                gSum(ci.primitiveField()*mesh.V())
               /domainVolume;

            cMaxBeforeClipping[i] =
                gMax(ci.primitiveField());
        }


        // Positivity clipping is intentionally postponed until AFTER the
        // fast-equilibrium projection.
        // ci.max(0.0);

        ci.correctBoundaryConditions();
    }
}


// ------------------------------------------------------------------------- //
// Charge-weighted fvMatrix flux summary
// ------------------------------------------------------------------------- //

volScalarField divChargeMatrixTransportFlux
(
    IOobject
    (
        "divChargeMatrixTransportFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(chargeMatrixTransportFlux)
);


volScalarField divChargeManualTransportFlux
(
    IOobject
    (
        "divChargeManualTransportFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(chargeManualTransportFlux)
);


volScalarField divChargeFluxDifference
(
    IOobject
    (
        "divChargeFluxDifference",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    divChargeMatrixTransportFlux
  - divChargeManualTransportFlux
);


Info<< "FV-FLUX charge summary:" << nl
    << "  max|div(matrix charge transport flux)| = "
    << gMax
       (
           mag
           (
               divChargeMatrixTransportFlux.primitiveField()
           )
       )
    << " mol/m3/s" << nl
    << "  max|div(manual charge transport flux)| = "
    << gMax
       (
           mag
           (
               divChargeManualTransportFlux.primitiveField()
           )
       )
    << " mol/m3/s" << nl
    << "  max|div(matrix-manual)| = "
    << gMax
       (
           mag
           (
               divChargeFluxDifference.primitiveField()
           )
       )
    << " mol/m3/s"
    << nl;


// Report the local worst cell for the matrix/manual divergence difference.
// This is local-rank information; in the current serial debugging case that
// is the actual mesh cell label.
{
    const scalarField& d =
        divChargeFluxDifference.primitiveField();

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
        Info<< "  local worst cell = "
            << localCell
            << ", centre = "
            << mesh.C()[localCell]
            << ", matrix div = "
            << divChargeMatrixTransportFlux[localCell]
            << ", manual div = "
            << divChargeManualTransportFlux[localCell]
            << ", difference = "
            << divChargeFluxDifference[localCell]
            << " mol/m3/s"
            << nl;
    }
}


Info<< "  charge-weighted solved-matrix residual: max = "
    << gMax
       (
           mag
           (
               chargeMatrixResidual.primitiveField()
           )
       )
    << " mol/m3/s"
    << nl;


// ======================================================================== //
// EQUILIBRIUM-PROJECTION TIMING DIAGNOSTIC
// ======================================================================== //
//
// We now know that the species fvMatrix face flux agrees with the manual
// reconstruction.  This block measures what the fast-equilibrium projection
// itself does to:
//
//   1) local charge q = sum(z_i c_i)
//   2) advective, diffusive and migrative charge-current divergence
//   3) the concentrations of each charged species at the most relevant cells
//
// Everything in this block uses the SAME phiE: the potential from the
// beginning of the current ENPhi iteration.  Thus PRE -> POST differences
// isolate the equilibrium projection rather than a potential update.
//
// PRE state  = immediately after all transport solves
// POST state = immediately after applyEquilibriumReactions()
//
// The later phiEqn.H diagnostic can then compare the final cEqn state against
// the newly updated potential.
// ======================================================================== //


// ------------------------------------------------------------------------- //
// Store concentrations immediately BEFORE fast-equilibrium projection
// ------------------------------------------------------------------------- //

List<scalarField> cBeforeEquilibrium(c.size());

forAll(c, i)
{
    cBeforeEquilibrium[i] =
        c[i].primitiveField();
}


// ------------------------------------------------------------------------- //
// PRE-equilibrium charge and current decomposition
// ------------------------------------------------------------------------- //

volScalarField qPreEquilibrium
(
    IOobject
    (
        "qPreEquilibrium",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimMoles/dimVolume,
        0.0
    )
);

surfaceScalarField preEqAdvCurrentFlux
(
    IOobject
    (
        "preEqAdvCurrentFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimCurrent, 0.0)
);

surfaceScalarField preEqDiffCurrentFlux
(
    IOobject
    (
        "preEqDiffCurrentFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimCurrent, 0.0)
);

surfaceScalarField preEqMigCurrentFlux
(
    IOobject
    (
        "preEqMigCurrentFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimCurrent, 0.0)
);

forAll(c, i)
{
    scalar zi = 0.0;

    if (cEqnZDictPtr)
    {
        zi =
            cEqnZDictPtr->getOrDefault<scalar>
            (
                species[i],
                0.0
            );
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

    qPreEquilibrium += zi*c[i];

    surfaceScalarField preEqMigrationFlux
    (
        IOobject
        (
            "preEqMigrationFlux_" + species[i],
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

    preEqAdvCurrentFlux +=
        Fdim*zi*fvc::interpolate(c[i])*phi;

    preEqDiffCurrentFlux +=
        Fdim*zi
       *(
            -Di*fvc::snGrad(c[i])*mesh.magSf()
        );

    preEqMigCurrentFlux +=
        Fdim*zi*(-preEqMigrationFlux);
}

surfaceScalarField preEqNoAdvCurrentFlux
(
    IOobject
    (
        "preEqNoAdvCurrentFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    preEqDiffCurrentFlux + preEqMigCurrentFlux
);

surfaceScalarField preEqTotalCurrentFlux
(
    IOobject
    (
        "preEqTotalCurrentFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    preEqAdvCurrentFlux + preEqNoAdvCurrentFlux
);

volScalarField divPreEqAdvCurrent
(
    IOobject
    (
        "divPreEqAdvCurrent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(preEqAdvCurrentFlux)
);

volScalarField divPreEqDiffCurrent
(
    IOobject
    (
        "divPreEqDiffCurrent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(preEqDiffCurrentFlux)
);

volScalarField divPreEqMigCurrent
(
    IOobject
    (
        "divPreEqMigCurrent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(preEqMigCurrentFlux)
);

volScalarField divPreEqNoAdvCurrent
(
    IOobject
    (
        "divPreEqNoAdvCurrent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(preEqNoAdvCurrentFlux)
);

volScalarField divPreEqTotalCurrent
(
    IOobject
    (
        "divPreEqTotalCurrent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(preEqTotalCurrentFlux)
);

Info<< "EQ-PROJECTION PRE summary:" << nl
    << "  max|q|                 = "
    << gMax(mag(qPreEquilibrium.primitiveField()))
    << " mol/m3" << nl
    << "  max|div(Iadv)|         = "
    << gMax(mag(divPreEqAdvCurrent.primitiveField()))
    << " A/m3" << nl
    << "  max|div(Idiff)|        = "
    << gMax(mag(divPreEqDiffCurrent.primitiveField()))
    << " A/m3" << nl
    << "  max|div(Imig)|         = "
    << gMax(mag(divPreEqMigCurrent.primitiveField()))
    << " A/m3" << nl
    << "  max|div(Idiff+Imig)|   = "
    << gMax(mag(divPreEqNoAdvCurrent.primitiveField()))
    << " A/m3" << nl
    << "  max|div(Itotal)|       = "
    << gMax(mag(divPreEqTotalCurrent.primitiveField()))
    << " A/m3" << nl;


// Optional water-equilibrium diagnostic before equilibrium projection.
if
(
    diagWaterEquilibriumA
 && reactionMgr.hasWaterEquilibrium()
)
{
    waterAPreclip.setSize
    (
        reactionMgr.nWaterEquilibriumReactions()
    );

    for
    (
        label reactioni = 0;
        reactioni
      < reactionMgr.nWaterEquilibriumReactions();
        ++reactioni
    )
    {
        const waterEquilibrium& waterEq =
            reactionMgr.waterEquilibriumReaction(reactioni);

        const volScalarField& cH =
            mesh.lookupObject<volScalarField>
            (
                "c_" + waterEq.HplusName()
            );

        const volScalarField& cOH =
            mesh.lookupObject<volScalarField>
            (
                "c_" + waterEq.OHminusName()
            );

        waterAPreclip[reactioni] =
            cH.primitiveField()
          - cOH.primitiveField();
    }
}


// ------------------------------------------------------------------------- //
// Apply all fast-equilibrium projections
// ------------------------------------------------------------------------- //

reactionMgr.applyEquilibriumReactions();


// ------------------------------------------------------------------------- //
// POST-equilibrium charge and current decomposition
//
// NOTE: this is deliberately evaluated immediately after the equilibrium
// projection, before the later positivity/BC correction.  Hence PRE -> POST
// isolates the projection operation itself.
// ------------------------------------------------------------------------- //

volScalarField qPostEquilibrium
(
    IOobject
    (
        "qPostEquilibrium",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimMoles/dimVolume,
        0.0
    )
);

surfaceScalarField postEqAdvCurrentFlux
(
    IOobject
    (
        "postEqAdvCurrentFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimCurrent, 0.0)
);

surfaceScalarField postEqDiffCurrentFlux
(
    IOobject
    (
        "postEqDiffCurrentFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimCurrent, 0.0)
);

surfaceScalarField postEqMigCurrentFlux
(
    IOobject
    (
        "postEqMigCurrentFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimCurrent, 0.0)
);

forAll(c, i)
{
    scalar zi = 0.0;

    if (cEqnZDictPtr)
    {
        zi =
            cEqnZDictPtr->getOrDefault<scalar>
            (
                species[i],
                0.0
            );
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

    qPostEquilibrium += zi*c[i];

    surfaceScalarField postEqMigrationFlux
    (
        IOobject
        (
            "postEqMigrationFlux_" + species[i],
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

    postEqAdvCurrentFlux +=
        Fdim*zi*fvc::interpolate(c[i])*phi;

    postEqDiffCurrentFlux +=
        Fdim*zi
       *(
            -Di*fvc::snGrad(c[i])*mesh.magSf()
        );

    postEqMigCurrentFlux +=
        Fdim*zi*(-postEqMigrationFlux);
}

surfaceScalarField postEqNoAdvCurrentFlux
(
    IOobject
    (
        "postEqNoAdvCurrentFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    postEqDiffCurrentFlux + postEqMigCurrentFlux
);

surfaceScalarField postEqTotalCurrentFlux
(
    IOobject
    (
        "postEqTotalCurrentFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    postEqAdvCurrentFlux + postEqNoAdvCurrentFlux
);

volScalarField divPostEqAdvCurrent
(
    IOobject
    (
        "divPostEqAdvCurrent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(postEqAdvCurrentFlux)
);

volScalarField divPostEqDiffCurrent
(
    IOobject
    (
        "divPostEqDiffCurrent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(postEqDiffCurrentFlux)
);

volScalarField divPostEqMigCurrent
(
    IOobject
    (
        "divPostEqMigCurrent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(postEqMigCurrentFlux)
);

volScalarField divPostEqNoAdvCurrent
(
    IOobject
    (
        "divPostEqNoAdvCurrent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(postEqNoAdvCurrentFlux)
);

volScalarField divPostEqTotalCurrent
(
    IOobject
    (
        "divPostEqTotalCurrent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(postEqTotalCurrentFlux)
);

volScalarField deltaQEquilibrium
(
    IOobject
    (
        "deltaQEquilibrium",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    qPostEquilibrium - qPreEquilibrium
);

volScalarField deltaDivNoAdvCurrentEquilibrium
(
    IOobject
    (
        "deltaDivNoAdvCurrentEquilibrium",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    divPostEqNoAdvCurrent - divPreEqNoAdvCurrent
);

volScalarField deltaDivTotalCurrentEquilibrium
(
    IOobject
    (
        "deltaDivTotalCurrentEquilibrium",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    divPostEqTotalCurrent - divPreEqTotalCurrent
);

Info<< "EQ-PROJECTION POST summary:" << nl
    << "  max|q|                 = "
    << gMax(mag(qPostEquilibrium.primitiveField()))
    << " mol/m3" << nl
    << "  max|delta q|           = "
    << gMax(mag(deltaQEquilibrium.primitiveField()))
    << " mol/m3" << nl
    << "  max|div(Iadv)|         = "
    << gMax(mag(divPostEqAdvCurrent.primitiveField()))
    << " A/m3" << nl
    << "  max|div(Idiff)|        = "
    << gMax(mag(divPostEqDiffCurrent.primitiveField()))
    << " A/m3" << nl
    << "  max|div(Imig)|         = "
    << gMax(mag(divPostEqMigCurrent.primitiveField()))
    << " A/m3" << nl
    << "  max|div(Idiff+Imig)|   = "
    << gMax(mag(divPostEqNoAdvCurrent.primitiveField()))
    << " A/m3" << nl
    << "  max|div(Itotal)|       = "
    << gMax(mag(divPostEqTotalCurrent.primitiveField()))
    << " A/m3" << nl
    << "  max|delta div(noAdv)|  = "
    << gMax(mag(deltaDivNoAdvCurrentEquilibrium.primitiveField()))
    << " A/m3" << nl
    << "  max|delta div(total)|  = "
    << gMax(mag(deltaDivTotalCurrentEquilibrium.primitiveField()))
    << " A/m3" << nl;


// ------------------------------------------------------------------------- //
// POST-equilibrium PHYSICAL TRANSPORT residual diagnostic
//
// Re-assemble the UNRELAXED species transport equation on the POST-equilibrium
// concentrations, but DO NOT solve it.  This measures whether the projected
// state still satisfies the transport equation:
//
//   ddt(c_i) + div(phi,c_i)
// = laplacian(D_i,c_i) + div(migration_i) + Rc_i + fvOptions(c_i)
//
// A small residual here would mean the equilibrium projection preserves the
// transport solution.  A large residual, especially at the known electrode
// cell, proves that the projection moves the state off the transport equation.
//
// The migration term uses the same phiE as the transport solve in this ENPhi
// iteration, so this is a clean transport-vs-projection test.
// ------------------------------------------------------------------------- //

volScalarField chargePostEqPhysicalTransportResidual
(
    IOobject
    (
        "chargePostEqPhysicalTransportResidual",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimMoles/dimVolume/dimTime,
        0.0
    )
);

Info<< "POST-EQ PHYSICAL TRANSPORT residuals:" << nl;

forAll(c, i)
{
    if (i == inertIndex)
    {
        continue;
    }

    volScalarField& ci = c[i];

    dimensionedScalar Di
    (
        "D_" + species[i],
        dimArea/dimTime,
        diffusivityDict.get<scalar>(species[i])
    );

    scalar zi = 0.0;

    if (cEqnZDictPtr)
    {
        zi =
            cEqnZDictPtr->getOrDefault<scalar>
            (
                species[i],
                0.0
            );
    }

    surfaceScalarField postEqMigrationFluxForResidual
    (
        IOobject
        (
            "postEqMigrationFluxForResidual_" + species[i],
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fvc::interpolate
        (
            Di*zi*Fdim/(Rgas*T)*ci
        )
       *gradPhiEFlux
    );

    fvScalarMatrix ciEqnPostEqPhysical
    (
        fvm::ddt(ci)
      + fvm::div(phi, ci)
     ==
        fvm::laplacian(Di, ci)
      + fvc::div(postEqMigrationFluxForResidual)
      + Rc[i]
      + Req[i]
      + fvOptions(ci)
    );

    // Match constraints that define the physical discretized equation,
    // but deliberately DO NOT under-relax and DO NOT solve.
    fvOptions.constrain(ciEqnPostEqPhysical);

    tmp<scalarField> tPostEqPhysicalResidual =
        ciEqnPostEqPhysical.residual();

    const scalarField& postEqPhysicalResidual =
        tPostEqPhysicalResidual();

    const scalarField postEqPhysicalResidualPerV =
        postEqPhysicalResidual/mesh.V();

    Info<< "  " << species[i]
        << ": max|R_transport(postEq)| = "
        << gMax(mag(postEqPhysicalResidualPerV))
        << " mol/m3/s" << nl;

    if (mag(zi) > SMALL)
    {
        chargePostEqPhysicalTransportResidual.primitiveFieldRef() +=
            zi*postEqPhysicalResidualPerV;
    }
}

Info<< "  charge-weighted max|sum(z_i R_i,transport postEq)| = "
    << gMax
       (
           mag
           (
               chargePostEqPhysicalTransportResidual.primitiveField()
           )
       )
    << " mol/m3/s" << nl;


// Identify the worst local post-equilibrium transport-residual cell.
label localWorstPostEqTransportResidualCell = -1;
scalar localWorstPostEqTransportResidual = -GREAT;

forAll
(
    chargePostEqPhysicalTransportResidual.primitiveField(),
    celli
)
{
    const scalar magResidual =
        mag(chargePostEqPhysicalTransportResidual[celli]);

    if (magResidual > localWorstPostEqTransportResidual)
    {
        localWorstPostEqTransportResidual = magResidual;
        localWorstPostEqTransportResidualCell = celli;
    }
}

if (localWorstPostEqTransportResidualCell >= 0)
{
    const label celli =
        localWorstPostEqTransportResidualCell;

    Info<< "POST-EQ PHYSICAL TRANSPORT worst charge-residual local cell = "
        << celli
        << ", centre = " << mesh.C()[celli] << nl
        << "  sum(z_i R_i,transport postEq) = "
        << chargePostEqPhysicalTransportResidual[celli]
        << " mol/m3/s" << nl
        << "  q PRE                           = "
        << qPreEquilibrium[celli] << " mol/m3" << nl
        << "  q POST                          = "
        << qPostEquilibrium[celli] << " mol/m3" << nl
        << "  delta q                         = "
        << deltaQEquilibrium[celli] << " mol/m3" << nl
        << "  PRE div(Idiff+Imig)             = "
        << divPreEqNoAdvCurrent[celli] << " A/m3" << nl
        << "  POST div(Idiff+Imig)            = "
        << divPostEqNoAdvCurrent[celli] << " A/m3" << nl;
}


// ------------------------------------------------------------------------- //
// Report local cells of interest and charged-species PRE -> POST changes
// ------------------------------------------------------------------------- //

label localWorstQCell = -1;
scalar localWorstQ = -GREAT;

forAll(qPostEquilibrium.primitiveField(), celli)
{
    const scalar magQ =
        mag(qPostEquilibrium[celli]);

    if (magQ > localWorstQ)
    {
        localWorstQ = magQ;
        localWorstQCell = celli;
    }
}

label localWorstCurrentChangeCell = -1;
scalar localWorstCurrentChange = -GREAT;

forAll(deltaDivNoAdvCurrentEquilibrium.primitiveField(), celli)
{
    const scalar magChange =
        mag(deltaDivNoAdvCurrentEquilibrium[celli]);

    if (magChange > localWorstCurrentChange)
    {
        localWorstCurrentChange = magChange;
        localWorstCurrentChangeCell = celli;
    }
}


if (localWorstQCell >= 0)
{
    const label celli = localWorstQCell;

    Info<< "EQ-PROJECTION worst POST-charge local cell = "
        << celli
        << ", centre = " << mesh.C()[celli] << nl
        << "  q PRE                 = "
        << qPreEquilibrium[celli] << " mol/m3" << nl
        << "  q POST                = "
        << qPostEquilibrium[celli] << " mol/m3" << nl
        << "  delta q               = "
        << deltaQEquilibrium[celli] << " mol/m3" << nl
        << "  PRE div(Idiff+Imig)   = "
        << divPreEqNoAdvCurrent[celli] << " A/m3" << nl
        << "  POST div(Idiff+Imig)  = "
        << divPostEqNoAdvCurrent[celli] << " A/m3" << nl
        << "  delta div(noAdv)      = "
        << deltaDivNoAdvCurrentEquilibrium[celli] << " A/m3" << nl
        << "  PRE div(Itotal)       = "
        << divPreEqTotalCurrent[celli] << " A/m3" << nl
        << "  POST div(Itotal)      = "
        << divPostEqTotalCurrent[celli] << " A/m3" << nl;

    Info<< "  charged species PRE -> POST at this cell:" << nl;

    forAll(c, i)
    {
        scalar zi = 0.0;

        if (cEqnZDictPtr)
        {
            zi =
                cEqnZDictPtr->getOrDefault<scalar>
                (
                    species[i],
                    0.0
                );
        }

        if (mag(zi) <= SMALL)
        {
            continue;
        }

        const scalar cPre =
            cBeforeEquilibrium[i][celli];

        const scalar cPost =
            c[i][celli];

        Info<< "    " << species[i]
            << " z=" << zi
            << ": pre=" << cPre
            << ", post=" << cPost
            << ", delta=" << cPost - cPre
            << " mol/m3"
            << nl;
    }
}


if
(
    localWorstCurrentChangeCell >= 0
 && localWorstCurrentChangeCell != localWorstQCell
)
{
    const label celli =
        localWorstCurrentChangeCell;

    Info<< "EQ-PROJECTION worst current-change local cell = "
        << celli
        << ", centre = " << mesh.C()[celli] << nl
        << "  q PRE                 = "
        << qPreEquilibrium[celli] << " mol/m3" << nl
        << "  q POST                = "
        << qPostEquilibrium[celli] << " mol/m3" << nl
        << "  delta q               = "
        << deltaQEquilibrium[celli] << " mol/m3" << nl
        << "  PRE div(Idiff+Imig)   = "
        << divPreEqNoAdvCurrent[celli] << " A/m3" << nl
        << "  POST div(Idiff+Imig)  = "
        << divPostEqNoAdvCurrent[celli] << " A/m3" << nl
        << "  delta div(noAdv)      = "
        << deltaDivNoAdvCurrentEquilibrium[celli] << " A/m3" << nl;

    Info<< "  charged species PRE -> POST at this cell:" << nl;

    forAll(c, i)
    {
        scalar zi = 0.0;

        if (cEqnZDictPtr)
        {
            zi =
                cEqnZDictPtr->getOrDefault<scalar>
                (
                    species[i],
                    0.0
                );
        }

        if (mag(zi) <= SMALL)
        {
            continue;
        }

        const scalar cPre =
            cBeforeEquilibrium[i][celli];

        const scalar cPost =
            c[i][celli];

        Info<< "    " << species[i]
            << " z=" << zi
            << ": pre=" << cPre
            << ", post=" << cPost
            << ", delta=" << cPost - cPre
            << " mol/m3"
            << nl;
    }
}


const scalar invDeltaT =
    1.0/runTime.deltaTValue();

// ------------------------------------------------------------------------- //
// LAGGED / CUMULATIVE EQUILIBRIUM-SOURCE FEEDBACK TEST
//
// IMPORTANT:
//   deltaCeq/dt is the *missing source correction* for the transport solve
//   that used the current lagged Req.  It is therefore accumulated:
//
//       Req^(k+1) = Req^k + omegaReq * (cPostEq - cPreEq)/dt
//
//   It is NOT used as an absolute replacement for Req.
//
// At convergence:
//   - Req may remain finite (it represents the algebraic fast-reaction source)
//   - the projection increment cPostEq-cPreEq should tend to zero
//   - Req^(k+1)-Req^k should tend to zero
// ------------------------------------------------------------------------- //

const scalar reqRelax = 0.35;

scalar eqProjectionMaxAbsCorrection = 0.0;
scalar reqUpdateMaxAbs = 0.0;
scalar reqMagnitudeMaxAbs = 0.0;

volScalarField reqChargeSource
(
    IOobject
    (
        "reqChargeSource",
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

forAll(c, i)
{
    const scalarField deltaCeq
    (
        c[i].primitiveField()
      - cBeforeEquilibrium[i]
    );

    eqProjectionMaxAbsCorrection =
        max
        (
            eqProjectionMaxAbsCorrection,
            gMax(mag(deltaCeq))
        );

    scalarField& ReqInternal =
        Req[i].primitiveFieldRef();

    const scalarField reqIncrement
    (
        reqRelax*deltaCeq*invDeltaT
    );

    reqUpdateMaxAbs =
        max
        (
            reqUpdateMaxAbs,
            gMax(mag(reqIncrement))
        );

    ReqInternal += reqIncrement;

    reqMagnitudeMaxAbs =
        max
        (
            reqMagnitudeMaxAbs,
            gMax(mag(ReqInternal))
        );

    // Equilibrium source bookkeeping is internal-cell only.
    forAll(Req[i].boundaryField(), patchi)
    {
        Req[i].boundaryFieldRef()[patchi] == 0.0;
    }

    scalar ziReq = 0.0;

    if (cEqnZDictPtr)
    {
        ziReq =
            cEqnZDictPtr->getOrDefault<scalar>
            (
                species[i],
                0.0
            );
    }

    if (mag(ziReq) > SMALL)
    {
        reqChargeSource.primitiveFieldRef() +=
            ziReq*ReqInternal;
    }
}

reduce(eqProjectionMaxAbsCorrection, maxOp<scalar>());
reduce(reqUpdateMaxAbs, maxOp<scalar>());
reduce(reqMagnitudeMaxAbs, maxOp<scalar>());

const scalar reqChargeMaxAbs =
    gMax(mag(reqChargeSource.primitiveField()));

Info<< "REQ-FEEDBACK summary:" << nl
    << "  relaxation omega                  = " << reqRelax << nl
    << "  max|projection delta c|           = "
    << eqProjectionMaxAbsCorrection << " mol/m3" << nl
    << "  max|Req update|                   = "
    << reqUpdateMaxAbs << " mol/m3/s" << nl
    << "  max|Req|                          = "
    << reqMagnitudeMaxAbs << " mol/m3/s" << nl
    << "  max|sum(z_i Req_i)|               = "
    << reqChargeMaxAbs << " mol/m3/s" << nl;


// ------------------------------------------------------------------------- //
// Positivity enforcement only AFTER equilibrium projection
// ------------------------------------------------------------------------- //

forAll(c, i)
{
    if (i != inertIndex)
    {
        c[i].max(0.0);
        c[i].correctBoundaryConditions();
    }
}


// ------------------------------------------------------------------------- //
// pH update
// ------------------------------------------------------------------------- //

if (pHActive)
{
    volScalarField& pH =
        pHPtr();

    const volScalarField& cHplus =
        c[HplusIndex];

    const dimensionedScalar standardConcentration
    (
        "standardConcentration",
        cHplus.dimensions(),
        1000.0
    );

    const dimensionedScalar concentrationFloor
    (
        "concentrationFloor",
        cHplus.dimensions(),
        SMALL
    );

    pH =
        -log10
        (
            max
            (
                cHplus,
                concentrationFloor
            )
           /standardConcentration
        );

    pH.correctBoundaryConditions();
}

