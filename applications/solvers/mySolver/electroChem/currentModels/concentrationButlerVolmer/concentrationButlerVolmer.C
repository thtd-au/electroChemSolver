#include "concentrationButlerVolmer.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"
#include "dimensionSets.H"
#include "PstreamReduceOps.H"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Foam
{
namespace
{
namespace electroChemConcentration
{
// Returns log((c/cref)^order). Exact zero factors extinguish a branch;
// no positive concentration floor or artificial reaction is introduced.
template<class Real>
Real logFactor(Real concentration, Real reference, Real order)
{
    if (!std::isfinite(concentration) || concentration < Real(0))
        throw std::runtime_error("BV requires finite non-negative surface concentrations");
    if (!std::isfinite(reference) || reference <= Real(0) || !std::isfinite(order))
        throw std::runtime_error("Invalid concentration reference/order");
    if (order == Real(0)) return Real(0);
    if (concentration == Real(0))
    {
        if (order < Real(0))
            throw std::runtime_error("Negative kinetic order is singular at zero concentration");
        return -std::numeric_limits<Real>::infinity();
    }
    const Real value=order*(std::log(concentration)-std::log(reference));
    if (!std::isfinite(value))
        throw std::runtime_error("BV logarithmic concentration factor exceeds scalar range");
    return value;
}
} // namespace electroChemConcentration
}

defineTypeNameAndDebug(concentrationButlerVolmer, 0);
addToRunTimeSelectionTable(currentModel, concentrationButlerVolmer, dictionary);

concentrationButlerVolmer::concentrationButlerVolmer
(
    fvMesh& mesh, const dictionary& currentDict, const volScalarField& phiE
)
:
    ButlerVolmer(mesh, currentDict, phiE, false, true),
    dependencyNames_(electrodeNames_.size()),
    referenceConcentrations_(electrodeNames_.size()),
    anodicOrders_(electrodeNames_.size()),
    cathodicOrders_(electrodeNames_.size())
{
    const IOdictionary& reactions=mesh.lookupObject<IOdictionary>("reactions");
    const wordList species(reactions.lookup("species"));
    const word inert(reactions.get<word>("inertSpecie"));
    forAll(electrodeNames_, e)
    {
        const dictionary& electrode=currentDict.subDict("electrodes").subDict(electrodeNames_[e]);
        // An empty dictionary is an explicit concentration-independent test.
        const dictionary& dependencies=electrode.subDict("concentrationDependence");
        dependencyNames_[e]=dependencies.sortedToc();
        referenceConcentrations_[e].setSize(dependencyNames_[e].size());
        anodicOrders_[e].setSize(dependencyNames_[e].size());
        cathodicOrders_[e].setSize(dependencyNames_[e].size());
        forAll(dependencyNames_[e], i)
        {
            const word& specie=dependencyNames_[e][i];
            const dictionary& entry=dependencies.subDict(specie);
            if (species.find(specie)<0)
                FatalIOErrorInFunction(entry) << "Unknown concentration species " << specie
                    << "; available: " << species << exit(FatalIOError);
            referenceConcentrations_[e][i]=entry.get<scalar>("c_ref");
            anodicOrders_[e][i]=entry.get<scalar>("anodicOrder");
            cathodicOrders_[e][i]=entry.get<scalar>("cathodicOrder");
            if (!std::isfinite(referenceConcentrations_[e][i])
             || referenceConcentrations_[e][i]<=0
             || !std::isfinite(anodicOrders_[e][i]) || !std::isfinite(cathodicOrders_[e][i]))
                FatalIOErrorInFunction(entry) << "Require positive finite c_ref and finite branch orders"
                    << exit(FatalIOError);
            if (specie==inert && (anodicOrders_[e][i]!=0 || cathodicOrders_[e][i]!=0))
                FatalIOErrorInFunction(entry)
                    << "Inert closure species " << specie << " is updated after Picard convergence. "
                    << "Use unit solvent activity (omit its dependency or set both orders zero)."
                    << exit(FatalIOError);
            const volScalarField& field=mesh.lookupObject<volScalarField>("c_"+specie);
            if (field.dimensions()!=dimMoles/dimVolume)
                FatalIOErrorInFunction(entry) << field.name() << " must have concentration dimensions"
                    << exit(FatalIOError);
        }
    }
    // All derived state now exists; virtual factors participate in this root.
    correct();
}

void concentrationButlerVolmer::logConcentrationFactors
(
    const label e, scalarField& logCa, scalarField& logCc
) const
{
    ButlerVolmer::logConcentrationFactors(e,logCa,logCc);
    bool valid=true;
    word failedSpecies;
    forAll(dependencyNames_[e], i)
    {
        const volScalarField& field=mesh_.lookupObject<volScalarField>("c_"+dependencyNames_[e][i]);
        const fvPatchScalarField& patch=field.boundaryField()[patchIDs_[e]];
        // Use stored Picard face values: evaluating a BC here would change
        // the frozen state and could create a current/concentration recursion.
        try
        {
            forAll(patch, f)
            {
                logCa[f]+=electroChemConcentration::logFactor
                    (patch[f],referenceConcentrations_[e][i],anodicOrders_[e][i]);
                logCc[f]+=electroChemConcentration::logFactor
                    (patch[f],referenceConcentrations_[e][i],cathodicOrders_[e][i]);
            }
        }
        catch (const std::exception&)
        {
            valid=false;
            failedSpecies=dependencyNames_[e][i];
        }
    }
    // Reject invalid data collectively before the controller's global sums.
    reduce(valid,andOp<bool>());
    if (!valid)
        FatalErrorInFunction << "Invalid surface concentration or singular order on electrode "
            << electrodeNames_[e] << "; local species=" << failedSpecies << nl
            << "Require non-negative finite concentrations; negative orders need positive concentrations."
            << exit(FatalError);
}

} // namespace Foam

