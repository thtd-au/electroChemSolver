#include "pHMolarConcentration.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "mathematicalConstants.H"

namespace Foam
{

defineTypeNameAndDebug(pHMolarConcentration, 0);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    pHMolarConcentration,
    dictionary
);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    pHMolarConcentration,
    patchMapper
);


pHMolarConcentration::pHMolarConcentration
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    equilibriumReaction_(word::null),
    pH_(7.0)
{}


pHMolarConcentration::pHMolarConcentration
(
    const pHMolarConcentration& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    equilibriumReaction_(ptf.equilibriumReaction_),
    pH_(ptf.pH_)
{}


pHMolarConcentration::pHMolarConcentration
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    equilibriumReaction_(dict.get<word>("equilibriumReaction")),
    pH_(dict.get<scalar>("pH"))
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=
        (
            patchInternalField()
        );
    }
}


pHMolarConcentration::pHMolarConcentration
(
    const pHMolarConcentration& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    equilibriumReaction_(ptf.equilibriumReaction_),
    pH_(ptf.pH_)
{}


pHMolarConcentration::pHMolarConcentration
(
    const pHMolarConcentration& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    equilibriumReaction_(ptf.equilibriumReaction_),
    pH_(ptf.pH_)
{}


void pHMolarConcentration::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (!mesh.foundObject<IOdictionary>("reactions"))
    {
        FatalErrorInFunction
            << "pHMolarConcentration requires the registered "
            << "constant/reactions dictionary." << nl
            << "Patch: " << patch().name() << nl
            << "Field: " << internalField().name()
            << exit(FatalError);
    }

    const IOdictionary& reactions =
        mesh.lookupObject<IOdictionary>("reactions");

    if (!reactions.found("equilibriumReactions"))
    {
        FatalIOErrorInFunction(reactions)
            << "pHMolarConcentration requires an equilibriumReactions "
            << "dictionary in constant/reactions." << nl
            << "Requested equilibriumReaction: " << equilibriumReaction_
            << exit(FatalIOError);
    }

    const dictionary& equilibriumReactions =
        reactions.subDict("equilibriumReactions");

    if (!equilibriumReactions.found(equilibriumReaction_))
    {
        FatalIOErrorInFunction(reactions)
            << "Unknown equilibriumReaction '" << equilibriumReaction_
            << "' for pHMolarConcentration." << nl
            << "Available equilibrium reactions are: "
            << equilibriumReactions.toc()
            << exit(FatalIOError);
    }

    const dictionary& equilibriumDict =
        equilibriumReactions.subDict(equilibriumReaction_);

    if (!equilibriumDict.found("type"))
    {
        FatalIOErrorInFunction(equilibriumDict)
            << "Equilibrium reaction '" << equilibriumReaction_
            << "' has no type entry."
            << exit(FatalIOError);
    }

    const word equilibriumType = equilibriumDict.get<word>("type");

    if
    (
        equilibriumType != "waterEquilibrium"
     && equilibriumType != "waterCarbonEquilibrium"
    )
    {
        FatalIOErrorInFunction(equilibriumDict)
            << "Equilibrium reaction '" << equilibriumReaction_
            << "' has unsupported type '" << equilibriumType << "'." << nl
            << "pHMolarConcentration currently supports: " << nl
            << "    waterEquilibrium" << nl
            << "    waterCarbonEquilibrium"
            << exit(FatalIOError);
    }

    if
    (
        !equilibriumDict.found("Hplus")
     || !equilibriumDict.found("OHminus")
     || !equilibriumDict.found("Kw")
    )
    {
        FatalIOErrorInFunction(equilibriumDict)
            << "Equilibrium reaction '" << equilibriumReaction_
            << "' must define Hplus, OHminus and Kw for "
            << "pHMolarConcentration."
            << exit(FatalIOError);
    }

    const word HplusName = equilibriumDict.get<word>("Hplus");
    const word OHminusName = equilibriumDict.get<word>("OHminus");
    const scalar Kw = equilibriumDict.get<scalar>("Kw");

    if (Kw <= 0.0)
    {
        FatalIOErrorInFunction(equilibriumDict)
            << "Kw must be greater than zero for equilibrium reaction '"
            << equilibriumReaction_ << "'."
            << exit(FatalIOError);
    }

    const word fieldName = internalField().name();
    const word HplusField("c_" + HplusName);
    const word OHminusField("c_" + OHminusName);

    const bool isHplus =
        fieldName == HplusName || fieldName == HplusField;

    const bool isOHminus =
        fieldName == OHminusName || fieldName == OHminusField;

    if (!isHplus && !isOHminus)
    {
        FatalErrorInFunction
            << "pHMolarConcentration is applied to field '"
            << fieldName << "', which is not a water ion of equilibrium "
            << "reaction '" << equilibriumReaction_ << "'." << nl
            << "That equilibrium reaction defines:" << nl
            << "    Hplus   " << HplusName << nl
            << "    OHminus " << OHminusName << nl
            << "Valid concentration fields for the current solver are:" << nl
            << "    " << HplusField << nl
            << "    " << OHminusField
            << exit(FatalError);
    }

    // The solver uses mol/m3, whereas pH is defined relative to 1 mol/L.
    // Therefore c_H+ [mol/m3] = 1000 * 10^(-pH).
    const scalar cHplus = 1000.0*pow(10.0, -pH_);

    if (cHplus <= 0.0)
    {
        FatalErrorInFunction
            << "pHMolarConcentration produced a non-positive Hplus "
            << "concentration for pH = " << pH_ << nl
            << "Patch: " << patch().name()
            << exit(FatalError);
    }

    const scalar concentration =
        isHplus ? cHplus : Kw/cHplus;

    operator==(scalarField(patch().size(), concentration));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void pHMolarConcentration::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("equilibriumReaction", equilibriumReaction_);
    os.writeEntry("pH", pH_);
    writeEntry("value", os);
}

} // End namespace Foam

