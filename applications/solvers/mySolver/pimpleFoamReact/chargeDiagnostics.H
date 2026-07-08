if (chargeDiagnosticsActive)
{
    const dictionary& molarMassDict = reactions.subDict("molarMass");

    dimensionedScalar Fdim
    (
        "FaradayConstant",
        dimCurrent*dimTime/dimMoles,
        reactions.getOrDefault<scalar>("FaradayConstant", 96485.3329)
    );

    c_charge =
        dimensionedScalar("zero", dimMoles/dimVolume, 0.0);

    List<scalar> speciesChargeIntegral(Y.size(), 0.0);

    forAll(Y, i)
    {
        if (mag(chargeNumber[i]) <= SMALL)
        {
            continue;
        }

        dimensionedScalar Mi
        (
            "M_" + species[i],
            dimMass/dimMoles,
            molarMassDict.get<scalar>(species[i])
        );

        volScalarField ci
        (
            IOobject
            (
                "cChargeTmp_" + species[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rho*Y[i]/Mi
        );

        c_charge += chargeNumber[i]*ci;

        speciesChargeIntegral[i] =
            chargeNumber[i]*gSum(ci.primitiveField()*mesh.V());
    }

    c_charge.correctBoundaryConditions();

    rho_e = Fdim*c_charge;
    rho_e.correctBoundaryConditions();

    const scalar Vtot = gSum(mesh.V());

    const scalar Cmin = min(c_charge).value();
    const scalar Cmax = max(c_charge).value();
    const scalar Cmean =
        gSum(c_charge.primitiveField()*mesh.V())/Vtot;
    const scalar Cint =
        gSum(c_charge.primitiveField()*mesh.V());

    const scalar rhoEmin = min(rho_e).value();
    const scalar rhoEmax = max(rho_e).value();
    const scalar rhoEmean =
        gSum(rho_e.primitiveField()*mesh.V())/Vtot;
    const scalar rhoEint =
        gSum(rho_e.primitiveField()*mesh.V());

    if (chargeDiagnosticsLogToInfo)
    {
        Info<< "Charge diagnostics:" << nl;

        if
        (
            chargeDiagnosticsOutput == "molarCharge"
         || chargeDiagnosticsOutput == "both"
        )
        {
            Info<< "  c_charge [mol/m3] min/mean/max/int = "
                << Cmin << " / "
                << Cmean << " / "
                << Cmax << " / "
                << Cint << nl;
        }

        if
        (
            chargeDiagnosticsOutput == "chargeDensity"
         || chargeDiagnosticsOutput == "both"
        )
        {
            Info<< "  rho_e [C/m3] min/mean/max/int = "
                << rhoEmin << " / "
                << rhoEmean << " / "
                << rhoEmax << " / "
                << rhoEint << nl;
        }

        if
        (
            chargeDiagnosticsMode == "full"
        )
        {
            Info<< "  species charge inventory [mol charge]:" << nl;

            forAll(species, i)
            {
                if (mag(chargeNumber[i]) > SMALL)
                {
                    Info<< "    " << species[i]
                        << " z=" << chargeNumber[i]
                        << " integral=" << speciesChargeIntegral[i]
                        << nl;
                }
            }
        }

        if (mag(Cint) > chargeDiagnosticsWarningTolerance)
        {
            WarningInFunction
                << "Domain-integrated molar charge = "
                << Cint
                << " exceeds warningTolerance = "
                << chargeDiagnosticsWarningTolerance
                << nl;
        }

        Info<< endl;
    }

    if (chargeDiagnosticsWriteToFile && chargeDiagnosticsFilePtr.valid())
    {
        OFstream& os = chargeDiagnosticsFilePtr();

        os  << runTime.value() << ' '
            << Cmin << ' '
            << Cmean << ' '
            << Cmax << ' '
            << Cint << ' '
            << rhoEmin << ' '
            << rhoEmean << ' '
            << rhoEmax << ' '
            << rhoEint << nl;
    }
}

