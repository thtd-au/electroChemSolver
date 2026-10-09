// cEqn.H -- PotentialReaction reactive Nernst-Planck transport

const dictionary& diffusivityDict = reactions.subDict("diffusivity");
const bool hasChargeNumber = reactions.found("chargeNumber");

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

// Kinetic rates are frozen for the complete species sweep and are updated
// exactly once by electroChem.C before this include is entered.

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
    if (hasChargeNumber)
    {
        zi = reactions.subDict("chargeNumber")
            .getOrDefault<scalar>(species[i], 0.0);
    }

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
        fvc::interpolate(Di*zi*Fdim/(Rgas*T)*ci)*gradPhiEFlux
    );

    tmp<volScalarField> tReactionSource = reactionMgr.source(i);

    // Store the total reaction source for post-convergence diagnostics.
    // No kinetic/equilibrium split is inferred from the transported state.
    Rc[i] = tReactionSource();

    fvScalarMatrix ciEqn
    (
        fvm::ddt(ci)
      + fvm::div(phi, ci)
     ==
        fvm::laplacian(Di, ci)
      + fvc::div(migrationFlux)
      + tReactionSource()
      + fvOptions(ci)
    );

    ciEqn.relax();
    fvOptions.constrain(ciEqn);

    // S_i = V/aP_i must come from the final constrained transport matrix.
    // Exact v2512 fvScalarMatrix::diag() semantics remain compile-verification item.
    reactionMgr.setTransportResponse(i, ciEqn);

    // TEST1J: OH- depletion audit. Runs BEFORE the positivity safeguard,
    // including the sweep which first produces a prohibited negative value.
    scalar minBeforeSolve1J = GREAT;
    if (species[i] == "OHminus")
    {
        minBeforeSolve1J = gMin(ci.primitiveField());
    }

    ciEqn.solve(ci.name());

    if (species[i] == "OHminus")
    {
        const scalar minAfterSolve1J = gMin(ci.primitiveField());
        if (minAfterSolve1J < 1.5e-4)
        {
            label worst1J = -1;
            scalar localMin1J = GREAT;
            forAll(ci.primitiveField(), celli)
            {
                if (ci[celli] < localMin1J)
                {
                    localMin1J = ci[celli];
                    worst1J = celli;
                }
            }
            Info<< "TEST1J OHminus: time=" << runTime.timeName()
                << " Picard=" << (picardI+1)
                << " minBeforeSolve=" << minBeforeSolve1J
                << " minAfterSolve=" << minAfterSolve1J << " mol/m3" << nl;

            if (worst1J >= 0 && localMin1J <= minAfterSolve1J + SMALL)
            {
                // Evaluate each term with the operators used in cEqn.H.
                // These are post-solve field evaluations, NOT an exact
                // decomposition of the implicit fvMatrix residual.
                const volScalarField advDiv1J(fvc::div(phi, ci));
                const volScalarField diffDiv1J(fvc::laplacian(Di, ci));
                const volScalarField migDiv1J(fvc::div(migrationFlux));
                const volScalarField ddt1J(fvc::ddt(ci));
                const scalarField matrixResidual1J(ciEqn.residual());
                const scalarField& V1J = mesh.V();
                const vector gradPhi1J = fvc::grad(phiE)()[worst1J];

                Info<< "TEST1J minimum cell=" << worst1J
                    << " centre=" << mesh.C()[worst1J]
                    << " volume=" << V1J[worst1J]
                    << " cOldTime=" << ci.oldTime()[worst1J]
                    << " cNow=" << ci[worst1J]
                    << " Hplus=";
                forAll(c, jj)
                {
                    if (species[jj] == "Hplus")
                    {
                        Info<< c[jj][worst1J];
                    }
                }
                Info<< " phiE=" << phiE[worst1J]
                    << " gradPhiE=" << gradPhi1J << nl;
                Info<< "TEST1J terms [mol/m3/s]: ddt=" << ddt1J[worst1J]
                    << " advDiv=" << advDiv1J[worst1J]
                    << " diffDiv=" << diffDiv1J[worst1J]
                    << " migDiv=" << migDiv1J[worst1J]
                    << " reaction=" << Rc[i][worst1J]
                    << " matrixResidual/V="
                    << matrixResidual1J[worst1J]/V1J[worst1J]
                    << nl;
                Info<< "TEST1J boundary patches touching min cell:";
                const cell& faces1J = mesh.cells()[worst1J];
                forAll(faces1J, fi)
                {
                    const label facei = faces1J[fi];
                    if (facei >= mesh.nInternalFaces())
                    {
                        const label patchi = mesh.boundaryMesh().whichPatch(facei);
                        if (patchi >= 0)
                        {
                            Info<< " " << mesh.boundaryMesh()[patchi].name();
                        }
                    }
                }
                Info<< nl;

                // TEST1K: face-resolved audit of the same minimum cell.
                // Signed values below are OUTWARD from this cell.
                // Matrix flux is the implicit convection+diffusion flux;
                // migration is explicit and has physical sign -migrationFlux.
                const surfaceScalarField matrixFace1K(ciEqn.flux());
                const surfaceScalarField advFace1K(phi*fvc::interpolate(ci));
                const surfaceScalarField diffFace1K
                (
                    -Di*fvc::snGrad(ci)*mesh.magSf()
                );
                const surfaceScalarField migFace1K(-migrationFlux);
                scalar sumMatrix1K=0, sumAdv1K=0, sumDiff1K=0, sumMig1K=0;
                Info<< "TEST1K face audit cell=" << worst1J
                    << " time=" << runTime.timeName() << nl;
                forAll(faces1J, fi)
                {
                    const label facei = faces1J[fi];
                    scalar orientation = 1.0;
                    scalar mat=0, adv=0, diff=0, mig=0, flow=0;
                    if (facei < mesh.nInternalFaces())
                    {
                        const label own = mesh.owner()[facei];
                        const label nei = mesh.neighbour()[facei];
                        orientation = own == worst1J ? 1.0 : -1.0;
                        mat = matrixFace1K[facei];
                        adv = advFace1K[facei];
                        diff = diffFace1K[facei];
                        mig = migFace1K[facei];
                        flow = phi[facei];
                        const label other = own == worst1J ? nei : own;
                        Info<< "TEST1K internal face=" << facei
                            << " neighbor=" << other
                            << " neighborC=" << ci[other];
                    }
                    else
                    {
                        const label patchi = mesh.boundaryMesh().whichPatch(facei);
                        const label localFace = facei-mesh.boundaryMesh()[patchi].start();
                        mat = matrixFace1K.boundaryField()[patchi][localFace];
                        adv = advFace1K.boundaryField()[patchi][localFace];
                        diff = diffFace1K.boundaryField()[patchi][localFace];
                        mig = migFace1K.boundaryField()[patchi][localFace];
                        flow = phi.boundaryField()[patchi][localFace];
                        Info<< "TEST1K boundary face=" << facei
                            << " patch=" << mesh.boundaryMesh()[patchi].name()
                            << " cBC=" << ci.boundaryField()[patchi].type()
                            << " faceC=" << ci.boundaryField()[patchi][localFace];
                    }
                    sumMatrix1K += orientation*mat;
                    sumAdv1K += orientation*adv;
                    sumDiff1K += orientation*diff;
                    sumMig1K += orientation*mig;
                    Info<< " phiOut=" << orientation*flow
                        << " matrixOut=" << orientation*mat
                        << " advOut=" << orientation*adv
                        << " diffOut=" << orientation*diff
                        << " migOut=" << orientation*mig << nl;
                }
                Info<< "TEST1K outward sums [mol/s]: matrix=" << sumMatrix1K
                    << " adv=" << sumAdv1K << " diff=" << sumDiff1K
                    << " mig=" << sumMig1K
                    << " totalMatrixPlusMig=" << sumMatrix1K+sumMig1K
                    << " /V=" << (sumMatrix1K+sumMig1K)/V1J[worst1J]
                    << " mol/m3/s" << nl;
            }
        }
    }

    fvOptions.correct(ci);
    ci.correctBoundaryConditions();

    // cMinAfterTransport[i] = gMin(ci.primitiveField());
    // cMeanAfterTransport[i] =
    //     gSum(ci.primitiveField()*mesh.V())/domainVolume;
    // cMaxAfterTransport[i] = gMax(ci.primitiveField());

    // if (cMinAfterTransport[i] < 0.0)
    // {
    //     FatalErrorInFunction
    //         << "Negative concentration after transport for " << species[i]
    //         << ": min(c)=" << cMinAfterTransport[i] << nl
    //         << "PotentialReaction does not clip concentrations because that "
    //         << "would alter the transported physical state."
    //         << exit(FatalError);
    // }

    if (printCFieldsBeforeClipping)
    {
        cPreMin1I[i] = gMin(ci.primitiveField());
        cPreMean1I[i] = gSum(ci.primitiveField()*mesh.V())/gSum(mesh.V());
        cPreMax1I[i] = gMax(ci.primitiveField());
    }

    const scalar minCi = gMin(c[i]);

    if (minCi < -negativeConcentrationTolerance)
    {
        FatalErrorInFunction
            << "Negative concentration after transport for "
            << species[i] << ": min(c)=" << minCi << " mol/m3" << nl
            << "negativeConcentrationTolerance="
            << negativeConcentrationTolerance << " mol/m3" << nl
            << "The undershoot exceeds the permitted numerical tolerance."
            << exit(FatalError);
    }

    if (minCi < 0.0)
    {
        Info<< "Clipping numerical concentration undershoot for "
            << species[i]
            << ": min(c)=" << minCi << " mol/m3"
            << " -> 0"
            << nl;

        c[i].max(0.0);
    }
}

// No fast-equilibrium concentration projection is permitted here.
// lambda corrections are calculated by reactionManager after this complete
// transported state has been evaluated, and affect concentrations only through
// the next cEqn solve.

