#ifndef potentialManager_H
#define potentialManager_H

#include "fvMesh.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "PtrList.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "error.H"
#include "boolList.H"
#include "token.H"
#include "fvPatchFields.H"
#include "PstreamReduceOps.H"
#include <cmath>
#include <limits>

namespace Foam
{

// Shared parallel operations used by the potential and current models.
namespace electroChemParallel
{

// Constraint patches are interfaces/geometry, not electrochemical boundaries.
inline bool preserveConstraint(const fvPatchScalarField& p)
{
    return p.coupled() || fvPatch::constraintType(p.patch().type());
}

inline word constraintType(const fvPatchScalarField& p)
{
    return p.coupled() ? p.type() : p.patch().type();
}

// Collective: every rank must call, including ranks with no reference faces.
// Select the lexicographically first face centre on the named physical patch.
// This keeps the physical reference independent of processor numbering for a
// static mesh with unique face centres. Only its owner returns a local cell.
inline label referenceCell(const fvMesh& mesh, const word& patchName)
{
    const label patchi = mesh.boundaryMesh().findPatchID(patchName);
    bool valid = patchi >= 0;
    if (valid) valid = !mesh.boundary()[patchi].coupled();
    reduce(valid, andOp<bool>());
    if (!valid)
    {
        FatalErrorInFunction
            << "referencePatch " << patchName
            << " must exist as a physical patch on every rank"
            << exit(FatalError);
    }

    const vectorField& centres = mesh.boundary()[patchi].Cf();
    const labelUList& cells = mesh.boundary()[patchi].faceCells();
    const scalar infinity = std::numeric_limits<scalar>::infinity();
    vector selected(infinity, infinity, infinity);
    label selectedFace = -1;
    valid = true;
    forAll(centres, facei)
    {
        const vector& centre = centres[facei];
        if (!std::isfinite(centre.x()) || !std::isfinite(centre.y())
         || !std::isfinite(centre.z()))
        {
            valid = false;
            continue;
        }
        if (selectedFace < 0 || centre.x() < selected.x()
         || (centre.x() == selected.x() && centre.y() < selected.y())
         || (centre.x() == selected.x() && centre.y() == selected.y()
             && centre.z() < selected.z()))
        {
            selected = centre;
            selectedFace = facei;
        }
    }
    reduce(valid, andOp<bool>());
    if (!valid)
    {
        FatalErrorInFunction << "Non-finite reference face centre on "
            << patchName << exit(FatalError);
    }

    scalar x = selected.x();
    reduce(x, minOp<scalar>());
    scalar y = selected.x() == x ? selected.y() : infinity;
    reduce(y, minOp<scalar>());
    scalar z = selected.x() == x && selected.y() == y
        ? selected.z() : infinity;
    reduce(z, minOp<scalar>());
    label owner = selectedFace >= 0 && selected.x() == x
        && selected.y() == y && selected.z() == z
        ? Pstream::myProcNo() : Pstream::nProcs();
    reduce(owner, minOp<label>());
    if (owner == Pstream::nProcs())
    {
        FatalErrorInFunction << "referencePatch " << patchName
            << " has no faces globally" << exit(FatalError);
    }
    return Pstream::myProcNo() == owner ? cells[selectedFace] : -1;
}

// Collective, with exactly one valid local reference index across all ranks.
inline scalar referenceValue(const volScalarField& field, const label celli)
{
    scalar value = celli >= 0 ? field[celli] : scalar(0);
    reduce(value, sumOp<scalar>());
    if (!std::isfinite(value))
    {
        FatalErrorInFunction << "Non-finite reference value for "
            << field.name() << exit(FatalError);
    }
    return value;
}

} // namespace electroChemParallel

// Shared output policy for the driver and both electrochemical managers.
inline label picardOutputLevel(const dictionary& controls)
{
    const word level = controls.lookupOrDefault<word>("picardOutput", "picard");
    if (level == "minimal") return 0;
    if (level == "picard") return 1;
    if (level == "full") return 2;
    if (level == "debug") return 3;
    FatalIOErrorInFunction(controls)
        << "picardOutput must be minimal, picard, full, or debug"
        << exit(FatalIOError);
    return 0;
}

// true/false or an unparenthesized sequence of component indices.
inline boolList potentialPrintDirections(const dictionary& controls, const word& key)
{
    boolList selected(3, false);
    if (!controls.found(key)) return selected;
    const ITstream& input = controls.lookup(key);
    if (input.size() == 1 && input[0].isWord())
    {
        const word value(input[0].wordToken());
        if (value == "true") return boolList(3, true);
        if (value == "false") return selected;
    }
    if (input.empty())
        FatalIOErrorInFunction(controls)
            << key << " requires true, false, or indices 0=x, 1=y, 2=z"
            << exit(FatalIOError);
    forAll(input, i)
    {
        if (!input[i].isLabel() || input[i].labelToken() < 0 || input[i].labelToken() > 2)
            FatalIOErrorInFunction(controls)
                << key << " requires true, false, or indices 0=x, 1=y, 2=z"
                << exit(FatalIOError);
        selected[input[i].labelToken()] = true;
    }
    return selected;
}

class potentialModel
{
protected:
    fvMesh& mesh_;
    const dictionary& electroChemistry_;
    const dictionary& reactions_;
    const wordList& species_;
    PtrList<volScalarField>& c_;
    volScalarField& phiE_;
    const dimensionedScalar& T_;
    const dimensionedScalar& Rgas_;
    const dimensionedScalar& F_;

public:
    TypeName("potentialModel");

    declareRunTimeSelectionTable
    (
        autoPtr, potentialModel, dictionary,
        (
            fvMesh& mesh,
            const dictionary& electroChemistry,
            const dictionary& reactions,
            const wordList& species,
            PtrList<volScalarField>& c,
            volScalarField& phiE,
            const dimensionedScalar& T,
            const dimensionedScalar& Rgas,
            const dimensionedScalar& F
        ),
        (mesh, electroChemistry, reactions, species, c, phiE, T, Rgas, F)
    );

    potentialModel
    (
        fvMesh& mesh,
        const dictionary& electroChemistry,
        const dictionary& reactions,
        const wordList& species,
        PtrList<volScalarField>& c,
        volScalarField& phiE,
        const dimensionedScalar& T,
        const dimensionedScalar& Rgas,
        const dimensionedScalar& F
    );

    virtual ~potentialModel() = default;

    static autoPtr<potentialModel> New
    (
        fvMesh& mesh,
        const dictionary& electroChemistry,
        const dictionary& reactions,
        const wordList& species,
        PtrList<volScalarField>& c,
        volScalarField& phiE,
        const dimensionedScalar& T,
        const dimensionedScalar& Rgas,
        const dimensionedScalar& F
    );

    virtual void evaluateConstraint() = 0;
    virtual bool converged() const = 0;
    virtual void correct(const bool force = false) = 0;
    virtual void report() const = 0;
    virtual const scalarField& chargeResidual() const = 0;
    virtual const scalarField& chargeScale() const = 0;
    virtual void refreshDiagnosticCurrent() = 0;
    virtual const volScalarField& conductivityField() const = 0;
    virtual const volVectorField& diffusionCurrentField() const = 0;
};


class potentialManager
{
    autoPtr<potentialModel> model_;

public:
    potentialManager
    (
        fvMesh& mesh,
        const dictionary& electroChemistry,
        const dictionary& reactions,
        const wordList& species,
        PtrList<volScalarField>& c,
        volScalarField& phiE,
        const dimensionedScalar& T,
        const dimensionedScalar& Rgas,
        const dimensionedScalar& F
    );

    bool active() const { return model_->type() != "none"; }

    void evaluateConstraint() { model_->evaluateConstraint(); }
    bool converged() const { return model_->converged(); }
    void correct(const bool force = false) { model_->correct(force); }
    void report() const { model_->report(); }
    const scalarField& chargeResidual() const { return model_->chargeResidual(); }
    const scalarField& chargeScale() const { return model_->chargeScale(); }
    void refreshDiagnosticCurrent() { model_->refreshDiagnosticCurrent(); }
    const volScalarField& conductivityField() const { return model_->conductivityField(); }
    const volVectorField& diffusionCurrentField() const { return model_->diffusionCurrentField(); }
};

} // End namespace Foam

#endif

