// enPhiChargeDiagnostics.H
//
// Temporary diagnostic for the ENPhi electroneutrality plateau.
//
// Charge-weighted species balance:
//
//   q = sum_i z_i c_i
//
//   ddt(q) + sum_i z_i div(phi,c_i)
//   = sum_i z_i laplacian(D_i,c_i)
//   + sum_i z_i div(migrationFlux_i)
//   + sum_i z_i Rc_i
//   + sum_i z_i Req_i
//
// where q has units mol/m3 and the balance terms have units mol/m3/s.

{
    volScalarField qNow
    (
        IOobject
        (
            "diagEN_qNow",
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

    volScalarField qDdt
    (
        IOobject
        (
            "diagEN_qDdt",
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

    volScalarField qAdv
    (
        IOobject
        (
            "diagEN_qAdv",
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

    volScalarField qDiff
    (
        IOobject
        (
            "diagEN_qDiff",
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

    volScalarField qMig
    (
        IOobject
        (
            "diagEN_qMig",
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

    volScalarField qKin
    (
        IOobject
        (
            "diagEN_qKin",
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

    volScalarField qEq
    (
        IOobject
        (
            "diagEN_qEq",
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

    // Old-time molar charge, evaluated directly from oldTime().
    scalarField qOld(mesh.nCells(), 0.0);

    const dictionary& diagDiffusivityDict =
        reactions.subDict("diffusivity");

    const dictionary* diagZDictPtr = nullptr;

    if (reactions.found("chargeNumber"))
    {
        diagZDictPtr =
            &reactions.subDict("chargeNumber");
    }

    forAll(c, i)
    {
        const scalar zi =
            diagZDictPtr
          ? diagZDictPtr->getOrDefault<scalar>
            (
                species[i],
                0.0
            )
          : 0.0;

        if (mag(zi) <= SMALL)
        {
            continue;
        }

        const volScalarField& ci = c[i];

        dimensionedScalar Di
        (
            "D_" + species[i],
            dimArea/dimTime,
            diagDiffusivityDict.get<scalar>(species[i])
        );

        qNow += zi*ci;

        qOld +=
            zi*ci.oldTime().primitiveField();

        qDdt +=
            zi*fvc::ddt(ci);

        qAdv +=
            zi*fvc::div(phi, ci);

        qDiff +=
            zi*fvc::laplacian(Di, ci);

        surfaceScalarField diagMigrationFlux
        (
            IOobject
            (
                "diagENMigrationFlux_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate
            (
                Di*zi*Fdim/(Rgas*T)*ci
            )
           *fvc::snGrad(phiE)
           *mesh.magSf()
        );

        qMig +=
            zi*fvc::div(diagMigrationFlux);

        qKin +=
            zi*Rc[i];

        qEq +=
            zi*Req[i];
    }


    // Electro-diffusive part.
    volScalarField qElectro
    (
        IOobject
        (
            "diagEN_qElectro",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        qDiff + qMig
    );

    // Total chemical charge source.
    volScalarField qReaction
    (
        IOobject
        (
            "diagEN_qReaction",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        qKin + qEq
    );

    // Discrete diagnostic residual:
    //
    // ddt + adv - diff - mig - kin - eq = 0
    volScalarField qBalanceResidual
    (
        IOobject
        (
            "diagEN_qBalanceResidual",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        qDdt
      + qAdv
      - qDiff
      - qMig
      - qKin
      - qEq
    );


    // ------------------------------------------------------------
    // Find cell containing maximum internal |sum(z_i c_i)|
    // ------------------------------------------------------------

    scalar localMaxCharge = -GREAT;
    label localMaxCell = -1;

    const scalarField& qNowInternal =
        qNow.primitiveField();

    forAll(qNowInternal, celli)
    {
        const scalar absQ =
            mag(qNowInternal[celli]);

        if (absQ > localMaxCharge)
        {
            localMaxCharge = absQ;
            localMaxCell = celli;
        }
    }

    const scalar globalMaxCharge =
        returnReduce
        (
            localMaxCharge,
            maxOp<scalar>()
        );


    // ------------------------------------------------------------
    // Global magnitudes
    // ------------------------------------------------------------

    const scalar maxQDdt =
        gMax(mag(qDdt.primitiveField()));

    const scalar maxQAdv =
        gMax(mag(qAdv.primitiveField()));

    const scalar maxQDiff =
        gMax(mag(qDiff.primitiveField()));

    const scalar maxQMig =
        gMax(mag(qMig.primitiveField()));

    const scalar maxQElectro =
        gMax(mag(qElectro.primitiveField()));

    const scalar maxQKin =
        gMax(mag(qKin.primitiveField()));

    const scalar maxQEq =
        gMax(mag(qEq.primitiveField()));

    const scalar maxQReaction =
        gMax(mag(qReaction.primitiveField()));

    const scalar maxQBalanceResidual =
        gMax(mag(qBalanceResidual.primitiveField()));


    Info<< "EN charge-balance diagnostic:" << nl
        << "  max|q|                 = "
        << globalMaxCharge
        << " mol/m3" << nl

        << "  max|sum z ddt(c)|       = "
        << maxQDdt
        << " mol/m3/s" << nl

        << "  max|sum z div(phi,c)|   = "
        << maxQAdv
        << " mol/m3/s" << nl

        << "  max|sum z diffusion|    = "
        << maxQDiff
        << " mol/m3/s" << nl

        << "  max|sum z migration|    = "
        << maxQMig
        << " mol/m3/s" << nl

        << "  max|diff+mig|           = "
        << maxQElectro
        << " mol/m3/s" << nl

        << "  max|sum z kinetic R|    = "
        << maxQKin
        << " mol/m3/s" << nl

        << "  max|sum z equilibrium R|= "
        << maxQEq
        << " mol/m3/s" << nl

        << "  max|total charge rxn|   = "
        << maxQReaction
        << " mol/m3/s" << nl

        << "  max|charge eq residual| = "
        << maxQBalanceResidual
        << " mol/m3/s" << nl;


    // Print detailed data on the processor that owns the global max.
    //
    // In your present serial case this is simply the max-charge cell.
    if
    (
        localMaxCell >= 0
     && mag(localMaxCharge - globalMaxCharge)
        <= SMALL*max(globalMaxCharge, scalar(1))
    )
    {
        const vector& x =
            mesh.C()[localMaxCell];

        Info<< "  max-charge cell:" << nl
            << "    processor = "
            << Pstream::myProcNo() << nl
            << "    cell      = "
            << localMaxCell << nl
            << "    centre    = "
            << x << nl
            << "    q(old)    = "
            << qOld[localMaxCell]
            << " mol/m3" << nl
            << "    q(now)    = "
            << qNow[localMaxCell]
            << " mol/m3" << nl

            << "    qDdt      = "
            << qDdt[localMaxCell] << nl
            << "    qAdv      = "
            << qAdv[localMaxCell] << nl
            << "    qDiff     = "
            << qDiff[localMaxCell] << nl
            << "    qMig      = "
            << qMig[localMaxCell] << nl
            << "    qElectro  = "
            << qElectro[localMaxCell] << nl
            << "    qKin      = "
            << qKin[localMaxCell] << nl
            << "    qEq       = "
            << qEq[localMaxCell] << nl
            << "    qReaction = "
            << qReaction[localMaxCell] << nl
            << "    residual  = "
            << qBalanceResidual[localMaxCell]
            << " mol/m3/s" << nl;

        Info<< "    charged-species contributions:"
            << nl;

        forAll(c, i)
        {
            const scalar zi =
                diagZDictPtr
              ? diagZDictPtr->getOrDefault<scalar>
                (
                    species[i],
                    0.0
                )
              : 0.0;

            if (mag(zi) <= SMALL)
            {
                continue;
            }

            Info<< "      "
                << species[i]
                << ": z=" << zi
                << ", c="
                << c[i][localMaxCell]
                << ", z*c="
                << zi*c[i][localMaxCell]
                << " mol/m3"
                << nl;
        }
    }

    Info<< endl;
}

