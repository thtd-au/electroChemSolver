// phiEqn.H
// Algebraic simple-current potential update.
// gradPhiE is the authoritative electric potential gradient.
// phiE is reconstructed only for reference/output.

{
    const dictionary& diffusivityDict = reactions.subDict("diffusivity");
    const bool hasChargeNumber = reactions.found("chargeNumber");

    const dictionary* zDictPtr = nullptr;

    if (hasChargeNumber)
    {
        zDictPtr = &reactions.subDict("chargeNumber");
    }

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

        jDiff -= Fdim*zi*Di*fvc::grad(c[i]);
        kappa += sqr(Fdim)*sqr(zi)*Di*c[i]/(Rgas*T);
    }

    const scalar kappaMin = gMin(kappa.primitiveField());

    if (kappaMin < minConductivity.value())
    {
        FatalErrorInFunction
            << "simpleCurrent potential update encountered conductivity "
            << "below minConductivity." << nl
            << "min(kappa) = " << kappaMin << " S/m" << nl
            << "minConductivity = " << minConductivity.value() << " S/m"
            << exit(FatalError);
    }

    scalarField gradPhiNormal(mesh.nCells(), 0.0);

    forAll(gradPhiNormal, celli)
    {
        const scalar jDiffNormal =
            jDiff[celli] & electrodeNormal;

        gradPhiNormal[celli] =
            (jDiffNormal - jApplied.value())/kappa[celli];

        gradPhiE[celli] =
            gradPhiNormal[celli]*electrodeNormal;
    }

    // The potential-gradient boundary values are taken directly from the
    // adjacent internal cells. The gradient is the primary resolved field.
    forAll(gradPhiE.boundaryField(), patchi)
    {
        gradPhiE.boundaryFieldRef()[patchi] ==
            gradPhiE.boundaryField()[patchi].patchInternalField();

        jDiff.boundaryFieldRef()[patchi] ==
            jDiff.boundaryField()[patchi].patchInternalField();

        kappa.boundaryFieldRef()[patchi] ==
            kappa.boundaryField()[patchi].patchInternalField();
    }

    // Reconstruct phiE column-by-column from the reference patch.
    // This first implementation assumes a serial structured mesh with
    // cell-centre connections aligned with electrodeNormal.
    scalarField& phiCells = phiE.primitiveFieldRef();
    boolList resolved(mesh.nCells(), false);

    const fvPatch& refPatch = mesh.boundary()[phiEReferencePatch];
    const labelUList& refFaceCells = refPatch.faceCells();
    const vectorField& refCf = mesh.Cf().boundaryField()[phiEReferencePatch];

    label nResolved = 0;

    forAll(refFaceCells, facei)
    {
        const label celli = refFaceCells[facei];
        const scalar ds =
            (mesh.C()[celli] - refCf[facei]) & electrodeNormal;

        const scalar reconstructedPhi =
            phiERef.value() + gradPhiNormal[celli]*ds;

        if (!resolved[celli])
        {
            phiCells[celli] = reconstructedPhi;
            resolved[celli] = true;
            ++nResolved;
        }
        else
        {
            const scalar consistencyScale =
                max(mag(phiCells[celli]), scalar(1.0));

            if
            (
                mag(phiCells[celli] - reconstructedPhi)
              > 1e-10*consistencyScale
            )
            {
                FatalErrorInFunction
                    << "Multiple faces on referencePatch reconstruct "
                    << "different potentials for cell " << celli << '.' << nl
                    << "This is inconsistent with the current structured "
                    << "1-D reconstruction assumption."
                    << exit(FatalError);
            }
        }
    }

    while (nResolved < mesh.nCells())
    {
        const label nResolvedBefore = nResolved;

        for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
        {
            const label owner = mesh.faceOwner()[facei];
            const label neighbour = mesh.faceNeighbour()[facei];

            if (resolved[owner] == resolved[neighbour])
            {
                continue;
            }

            const vector delta =
                mesh.C()[neighbour] - mesh.C()[owner];

            const scalar ds = delta & electrodeNormal;
            const vector transverseDelta =
                delta - ds*electrodeNormal;

            if
            (
                mag(ds) <= SMALL
             || mag(transverseDelta)
              > directionTolerance*max(mag(delta), SMALL)
            )
            {
                continue;
            }

            const scalar meanGradient =
                0.5*(gradPhiNormal[owner] + gradPhiNormal[neighbour]);

            if (resolved[owner])
            {
                phiCells[neighbour] =
                    phiCells[owner] + meanGradient*ds;

                resolved[neighbour] = true;
                ++nResolved;
            }
            else
            {
                phiCells[owner] =
                    phiCells[neighbour] - meanGradient*ds;

                resolved[owner] = true;
                ++nResolved;
            }
        }

        if (nResolved == nResolvedBefore)
        {
            label firstUnresolved = -1;

            forAll(resolved, celli)
            {
                if (!resolved[celli])
                {
                    firstUnresolved = celli;
                    break;
                }
            }

            FatalErrorInFunction
                << "Could not reconstruct phiE through all cells from "
                << "referencePatch " << phiEReferencePatchName << '.' << nl
                << "First unresolved cell = " << firstUnresolved << nl
                << "The simpleCurrent implementation requires a structured "
                << "mesh whose cell-centre columns are aligned with "
                << "electrodeNormal " << electrodeNormal << '.'
                << exit(FatalError);
        }
    }

    // Reconstruct boundary potential values from adjacent cells.
    // The reference patch itself is fixed to referencePotential.
    forAll(phiE.boundaryField(), patchi)
    {
        fvPatchScalarField& phiPatch =
            phiE.boundaryFieldRef()[patchi];

        const labelUList& faceCells =
            mesh.boundary()[patchi].faceCells();

        const vectorField& patchCf =
            mesh.Cf().boundaryField()[patchi];

        if (patchi == phiEReferencePatch)
        {
            phiPatch == phiERef.value();
            continue;
        }

        scalarField patchValues(phiPatch.size(), 0.0);

        forAll(patchValues, facei)
        {
            const label celli = faceCells[facei];
            const scalar ds =
                (patchCf[facei] - mesh.C()[celli]) & electrodeNormal;

            patchValues[facei] =
                phiCells[celli] + gradPhiNormal[celli]*ds;
        }

        phiPatch == patchValues;
    }

    Info<< "Potential fields:" << nl;

    const scalar domainVolumePotential = gSum(mesh.V());

    if (printPhiE)
    {
        const scalar phiMin = gMin(phiE.primitiveField());
        const scalar phiMean =
            gSum(phiE.primitiveField()*mesh.V())
            /domainVolumePotential;
        const scalar phiMax = gMax(phiE.primitiveField());

        Info<< "  phiE min/mean/max = "
            << phiMin << " / " << phiMean << " / " << phiMax
            << " V" << nl;
    }

    if (printGradPhiE)
    {
        const scalar gradMin = gMin(gradPhiNormal);
        const scalar gradMean =
            gSum(gradPhiNormal*mesh.V())
            /domainVolumePotential;
        const scalar gradMax = gMax(gradPhiNormal);

        Info<< "  gradPhiE(normal) min/mean/max = "
            << gradMin << " / " << gradMean << " / " << gradMax
            << " V/m" << nl;
    }

    if (printKappa)
    {
        const scalar kappaMean =
            gSum(kappa.primitiveField()*mesh.V())
            /domainVolumePotential;

        Info<< "  kappa min/mean/max = "
            << gMin(kappa.primitiveField())  << " / "
            << kappaMean << " / "
            << gMax(kappa.primitiveField()) 
            << " S/m" << nl;
    }

    if (printJDiff)
    {
        scalarField jDiffNormal(mesh.nCells(), 0.0);

        forAll(jDiffNormal, celli)
        {
            jDiffNormal[celli] =
                jDiff[celli] & electrodeNormal;
        }

        const scalar jDiffMean =
            gSum(jDiffNormal*mesh.V())
            /domainVolumePotential;

        Info<< "  jDiff(normal) min/mean/max = "
            << gMin(jDiffNormal) << " / "
            << jDiffMean << " / "
            << gMax(jDiffNormal)
            << " A/m2" << nl;
    }

    Info<< endl;
}

