#include "reactionManager.H"
#include "IOobject.H"
#include "dimensionedScalar.H"
#include "DynamicList.H"
#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(reactionModel, 0);
defineRunTimeSelectionTable(reactionModel, dictionary);


autoPtr<reactionModel> reactionModel::New
(
    const word& modelType,
    const word& name,
    const dictionary& dict,
    const wordList& species
)
{
    Info<< "    Selecting reaction model " << modelType
        << " for " << name << nl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "reaction model",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<reactionModel>(ctorPtr(name, dict, species));
}


scalar reactionModel::kineticRate
(
    const label,
    const scalarField&,
    const scalarField&
) const
{
    FatalErrorInFunction
        << "Reaction " << name_ << " is not a kinetic model"
        << exit(FatalError);

    return 0.0;
}


void reactionModel::equilibriumResidualAndJacobian
(
    const label,
    const scalarField&,
    const scalarField&,
    scalar&,
    scalarField&,
    scalar&,
    scalar&
) const
{
    FatalErrorInFunction
        << "Reaction " << name_ << " is not an equilibrium model"
        << exit(FatalError);
}


reactionManager::reactionRecord::reactionRecord()
:
    name(),
    type(),
    speciesIndices(),
    stoichCoeffs(),
    model(nullptr)
{}


reactionManager::reactionManager
(
    fvMesh& mesh,
    const dictionary& reactionsDict,
    const wordList& species,
    PtrList<volScalarField>& c,
    const label inertIndex
)
:
    mesh_(mesh), reactionsDict_(reactionsDict), species_(species), c_(c),
    inertIndex_(inertIndex), reactions_(), kineticReactionIndices_(),
    equilibriumReactionIndices_(), kineticRates_(), lambda_(),
    transportResponse_(), residual_(), residualScale_(),
    equilibriumAbsTol_(), equilibriumRelTol_(), equilibriumRelaxation_(1.0),
    equilibriumGram_(0), equilibriumN_(), equilibriumEvaluated_(false),
    equilibriumConverged_(false)
{
    readReactions();
    readEquilibriumControls();
    buildAndValidateEquilibriumStoichiometry();

    transportResponse_.setSize(species_.size());
    forAll(transportResponse_, i)
    {
        transportResponse_[i].setSize(mesh_.nCells(), 0.0);
    }

    Info<< "Reaction manager loaded:" << nl
        << "    reactions:             " << reactions_.size() << nl
        << "    kinetic reactions:     " << nKineticReactions() << nl
        << "    fast equilibria:       " << nEquilibriumReactions() << nl;
}


label reactionManager::speciesIndex(const word& name, const dictionary& dict) const
{
    const label i = species_.find(name);
    if (i < 0)
    {
        FatalIOErrorInFunction(dict) << "Unknown species " << name << nl
            << "Available species are " << species_ << exit(FatalIOError);
    }
    return i;
}


void reactionManager::validateStoichiometry
(
    const word& reactionName,
    const dictionary& reactionDict,
    reactionRecord& record
) const
{
    if (!reactionDict.found("stoichiometry"))
    {
        FatalIOErrorInFunction(reactionDict) << "Reaction " << reactionName
            << " is missing stoichiometry" << exit(FatalIOError);
    }

    const dictionary& stoichDict = reactionDict.subDict("stoichiometry");
    const wordList names = stoichDict.toc();
    if (names.empty())
    {
        FatalIOErrorInFunction(stoichDict) << "Reaction " << reactionName
            << " has empty stoichiometry" << exit(FatalIOError);
    }

    record.speciesIndices.setSize(names.size());
    record.stoichCoeffs.setSize(names.size());

    forAll(names, j)
    {
        const label speciesI = speciesIndex(names[j], stoichDict);
        const scalar nu = stoichDict.get<scalar>(names[j]);
        if (mag(nu) <= SMALL)
        {
            FatalIOErrorInFunction(stoichDict) << "Reaction " << reactionName
                << " lists zero stoichiometric coefficient for " << names[j]
                << exit(FatalIOError);
        }
        if (speciesI == inertIndex_)
        {
            FatalIOErrorInFunction(stoichDict) << "Reaction " << reactionName
                << " explicitly includes inertSpecie " << names[j] << nl
                << "The v1 PotentialReaction formulation omits the inert "
                << "closure species from active reaction stoichiometry."
                << exit(FatalIOError);
        }
        record.speciesIndices[j] = speciesI;
        record.stoichCoeffs[j] = nu;
    }
}


void reactionManager::readReactions()
{
    if (!reactionsDict_.found("reactions"))
    {
        FatalIOErrorInFunction(reactionsDict_) << "Missing reactions subdictionary"
            << exit(FatalIOError);
    }

    const dictionary& allReactions = reactionsDict_.subDict("reactions");
    const wordList names = allReactions.toc();
    reactions_.setSize(names.size());
    DynamicList<label> kineticIds, equilibriumIds;

    forAll(names, reactionI)
    {
        const word& reactionName = names[reactionI];
        if (!allReactions.isDict(reactionName))
        {
            FatalIOErrorInFunction(allReactions) << "Entry " << reactionName
                << " in reactions is not a dictionary" << exit(FatalIOError);
        }

        const dictionary& reactionDict = allReactions.subDict(reactionName);
        reactionRecord& record = reactions_[reactionI];
        record.name = reactionName;
        record.type = reactionDict.get<word>("type");
        validateStoichiometry(reactionName, reactionDict, record);
        record.model = reactionModel::New(record.type, reactionName, reactionDict, species_);

        if (record.model->category() == reactionModel::kinetic) kineticIds.append(reactionI);
        else equilibriumIds.append(reactionI);
    }

    kineticReactionIndices_.transfer(kineticIds);
    equilibriumReactionIndices_.transfer(equilibriumIds);

    kineticRates_.setSize(nKineticReactions());
    forAll(kineticRates_, j) kineticRates_[j].setSize(mesh_.nCells(), 0.0);

    lambda_.setSize(nEquilibriumReactions());
    residual_.setSize(nEquilibriumReactions());
    residualScale_.setSize(nEquilibriumReactions());
    forAll(lambda_, r)
    {
        lambda_[r].setSize(mesh_.nCells(), 0.0);
        residual_[r].setSize(mesh_.nCells(), 0.0);
        residualScale_[r].setSize(mesh_.nCells(), 1.0);
    }
}


void reactionManager::readEquilibriumControls()
{
    if (!mesh_.solutionDict().found("Picard"))
    {
        FatalIOErrorInFunction(mesh_.solutionDict())
            << "Missing Picard dictionary in fvSolution" << exit(FatalIOError);
    }
    const dictionary& picard = mesh_.solutionDict().subDict("Picard");
    if (!picard.found("equilibrium"))
    {
        FatalIOErrorInFunction(picard) << "Missing Picard/equilibrium dictionary"
            << exit(FatalIOError);
    }
    const dictionary& eqControls = picard.subDict("equilibrium");
    equilibriumRelaxation_ = eqControls.get<scalar>("relaxation");
    if (equilibriumRelaxation_ <= 0 || equilibriumRelaxation_ > 1)
    {
        FatalIOErrorInFunction(eqControls) << "equilibrium relaxation must satisfy 0 < relaxation <= 1"
            << exit(FatalIOError);
    }

    equilibriumAbsTol_.setSize(nEquilibriumReactions());
    equilibriumRelTol_.setSize(nEquilibriumReactions());

    forAll(equilibriumReactionIndices_, r)
    {
        const word& name = reactions_[equilibriumReactionIndices_[r]].name;
        if (!eqControls.found(name) || !eqControls.isDict(name))
        {
            FatalIOErrorInFunction(eqControls) << "Missing Picard/equilibrium/" << name
                << " tolerance dictionary" << exit(FatalIOError);
        }
        const dictionary& d = eqControls.subDict(name);
        equilibriumAbsTol_[r] = d.get<scalar>("absTol");
        equilibriumRelTol_[r] = d.get<scalar>("relTol");
        if (equilibriumAbsTol_[r] <= 0 || equilibriumRelTol_[r] <= 0)
        {
            FatalIOErrorInFunction(d) << "absTol and relTol must be positive"
                << exit(FatalIOError);
        }
    }

    const wordList entries = eqControls.toc();
    forAll(entries, i)
    {
        if (entries[i] == "relaxation" || !eqControls.isDict(entries[i])) continue;
        bool found = false;
        forAll(equilibriumReactionIndices_, r)
            found = found || reactions_[equilibriumReactionIndices_[r]].name == entries[i];
        if (!found)
        {
            FatalIOErrorInFunction(eqControls) << "Picard/equilibrium entry " << entries[i]
                << " does not name a fastEquilibriumReaction" << exit(FatalIOError);
        }
    }
}


scalar reactionManager::reactionStoich(const reactionRecord& record, const label speciesI) const
{
    forAll(record.speciesIndices, j)
        if (record.speciesIndices[j] == speciesI) return record.stoichCoeffs[j];
    return 0.0;
}


bool reactionManager::solveLocalSystem
(
    scalarSquareMatrix A,
    scalarField b,
    scalarField& x,
    scalar& minRelativePivot
) const
{
    const label n = b.size();
    x.setSize(n, 0.0);
    minRelativePivot = GREAT;
    if (!n) return true;

    scalarField rowScale(n, 0.0);
    for (label i=0; i<n; ++i)
        for (label j=0; j<n; ++j) rowScale[i] = max(rowScale[i], mag(A(i,j)));

    const scalar pivotTol = 100*SMALL;
    for (label k=0; k<n; ++k)
    {
        label pivot = k;
        scalar best = -1;
        for (label i=k; i<n; ++i)
        {
            const scalar rel = rowScale[i] > VSMALL ? mag(A(i,k))/rowScale[i] : 0.0;
            if (rel > best) { best = rel; pivot = i; }
        }
        minRelativePivot = min(minRelativePivot, max(best, scalar(0)));
        if (best <= pivotTol) return false;

        if (pivot != k)
        {
            for (label j=0; j<n; ++j) Swap(A(k,j), A(pivot,j));
            Swap(b[k], b[pivot]);
            Swap(rowScale[k], rowScale[pivot]);
        }

        for (label i=k+1; i<n; ++i)
        {
            const scalar f = A(i,k)/A(k,k);
            A(i,k) = 0;
            for (label j=k+1; j<n; ++j) A(i,j) -= f*A(k,j);
            b[i] -= f*b[k];
        }
    }

    for (label ii=n; ii-- > 0; )
    {
        scalar rhs = b[ii];
        for (label j=ii+1; j<n; ++j) rhs -= A(ii,j)*x[j];
        x[ii] = rhs/A(ii,ii);
        if (!std::isfinite(x[ii])) return false;
    }
    return true;
}


void reactionManager::buildAndValidateEquilibriumStoichiometry()
{
    const label ne = nEquilibriumReactions();
    equilibriumN_.setSize(species_.size());
    forAll(equilibriumN_, i) equilibriumN_[i].setSize(ne, 0.0);

    forAll(equilibriumReactionIndices_, r)
    {
        const reactionRecord& rec = reactions_[equilibriumReactionIndices_[r]];
        forAll(rec.speciesIndices, j)
            equilibriumN_[rec.speciesIndices[j]][r] = rec.stoichCoeffs[j];
    }

    if (!ne) return;

    // N has full column rank iff N^T N is nonsingular. This is only a
    // dimensionless startup structural check, so scaled pivoting is adequate.
    equilibriumGram_.setSize(ne);
    equilibriumGram_ = scalar(0);
    for (label r=0; r<ne; ++r)
        for (label s=0; s<ne; ++s)
            forAll(equilibriumN_, i)
                equilibriumGram_(r,s) += equilibriumN_[i][r]*equilibriumN_[i][s];

    scalarField rhs(ne, 0.0), dummy;
    rhs[0] = 1.0;
    scalar minPivot = GREAT;
    if (!solveLocalSystem(equilibriumGram_, rhs, dummy, minPivot))
    {
        FatalErrorInFunction << "Fast-equilibrium stoichiometric columns are structurally dependent. "
            << "rank(N) < number of fast equilibria. Reactions: ";
        forAll(equilibriumReactionIndices_, r)
            FatalError << reactions_[equilibriumReactionIndices_[r]].name << ' ';
        FatalError << exit(FatalError);
    }
}


void reactionManager::resetEquilibriumSources()
{
    forAll(lambda_, r) lambda_[r] = 0.0;
    equilibriumEvaluated_ = false;
    equilibriumConverged_ = false;
}


void reactionManager::updateKineticRates()
{
    scalarField concentrations(species_.size(), 0.0), nu(species_.size(), 0.0);
    forAll(kineticReactionIndices_, kineticI)
    {
        const reactionRecord& record = reactions_[kineticReactionIndices_[kineticI]];
        nu = 0.0;
        forAll(record.speciesIndices, j) nu[record.speciesIndices[j]] = record.stoichCoeffs[j];
        scalarField& rate = kineticRates_[kineticI];
        forAll(rate, celli)
        {
            forAll(c_, i) concentrations[i] = c_[i][celli];
            rate[celli] = record.model->kineticRate(celli, concentrations, nu);
            if (!std::isfinite(rate[celli]))
                FatalErrorInFunction << "Non-finite kinetic rate in reaction " << record.name
                    << " at cell " << celli << exit(FatalError);
        }
    }
}


tmp<volScalarField> reactionManager::source(const label speciesI) const
{
    if (speciesI < 0 || speciesI >= species_.size())
        FatalErrorInFunction << "Invalid species index " << speciesI << exit(FatalError);

    tmp<volScalarField> tSource
    (
        new volScalarField
        (
            IOobject("R_" + species_[speciesI], mesh_.time().timeName(), mesh_,
                IOobject::NO_READ, IOobject::NO_WRITE, false),
            mesh_, dimensionedScalar("zero", dimMoles/dimVolume/dimTime, 0.0)
        )
    );
    scalarField& src = tSource.ref().primitiveFieldRef();

    forAll(kineticReactionIndices_, k)
    {
        const scalar nu = reactionStoich(reactions_[kineticReactionIndices_[k]], speciesI);
        if (mag(nu) > SMALL) src += nu*kineticRates_[k];
    }
    forAll(equilibriumReactionIndices_, r)
    {
        const scalar nu = equilibriumN_[speciesI][r];
        if (mag(nu) > SMALL) src += nu*lambda_[r];
    }
    return tSource;
}


void reactionManager::setTransportResponse(const label speciesI, const fvScalarMatrix& cEqn)
{
    if (speciesI < 0 || speciesI >= species_.size() || speciesI == inertIndex_)
        FatalErrorInFunction << "Invalid active species index " << speciesI << exit(FatalError);

    const scalarField& aP = cEqn.diag();
    scalarField& S = transportResponse_[speciesI];
    const scalarField& V = mesh_.V();
    forAll(S, celli)
    {
        if (!std::isfinite(aP[celli]) || mag(aP[celli]) <= VSMALL)
            FatalErrorInFunction << "Invalid transport diagonal for " << species_[speciesI]
                << " at cell " << celli << ": aP=" << aP[celli] << exit(FatalError);
        S[celli] = V[celli]/aP[celli];
        if (!std::isfinite(S[celli]) || S[celli] <= 0)
            FatalErrorInFunction << "Invalid transport response S=V/aP for " << species_[speciesI]
                << " at cell " << celli << ": S=" << S[celli] << exit(FatalError);
    }
}


bool reactionManager::cellEquilibriumConverged(const label celli) const
{
    forAll(residual_, r)
    {
        const scalar absR = mag(residual_[r][celli]);
        const scalar relR = absR/residualScale_[r][celli];
        if (!(absR <= equilibriumAbsTol_[r] || relR <= equilibriumRelTol_[r])) return false;
    }
    return true;
}


void reactionManager::evaluateEquilibrium()
{
    scalarField concentrations(species_.size(), 0.0), nu(species_.size(), 0.0), dFdc(species_.size(), 0.0);
    equilibriumConverged_ = true;

    forAll(equilibriumReactionIndices_, r)
    {
        const reactionRecord& rec = reactions_[equilibriumReactionIndices_[r]];
        nu = 0.0;
        forAll(rec.speciesIndices, j) nu[rec.speciesIndices[j]] = rec.stoichCoeffs[j];

        forAll(mesh_.V(), celli)
        {
            forAll(c_, i) concentrations[i] = c_[i][celli];
            scalar F = 0, P = 0, KR = 0;
            rec.model->equilibriumResidualAndJacobian(celli, concentrations, nu, F, dFdc, P, KR);
            if (!std::isfinite(F) || !std::isfinite(P) || !std::isfinite(KR))
                FatalErrorInFunction << "Non-finite equilibrium residual in " << rec.name
                    << " at cell " << celli << exit(FatalError);
            residual_[r][celli] = F;
            residualScale_[r][celli] = max(max(mag(P), mag(KR)), equilibriumAbsTol_[r]);
        }
    }

    forAll(mesh_.V(), celli)
        if (!cellEquilibriumConverged(celli)) equilibriumConverged_ = false;

    reduce(equilibriumConverged_, andOp<bool>());
    equilibriumEvaluated_ = true;
}


bool reactionManager::equilibriumConverged() const
{
    if (!equilibriumEvaluated_)
        FatalErrorInFunction << "equilibriumConverged() called before evaluateEquilibrium()"
            << exit(FatalError);
    return equilibriumConverged_;
}


void reactionManager::correctEquilibrium()
{
    if (!equilibriumEvaluated_)
        FatalErrorInFunction << "correctEquilibrium() called before evaluateEquilibrium()"
            << exit(FatalError);

    const label ne = nEquilibriumReactions();
    if (!ne || equilibriumConverged_) return;

    scalarField concentrations(species_.size(), 0.0), nu(species_.size(), 0.0), dFdc(species_.size(), 0.0);
    scalarSquareMatrix J(ne);
    scalarField rhs(ne), delta(ne);

    forAll(mesh_.V(), celli)
    {
        if (cellEquilibriumConverged(celli)) continue;
        J = scalar(0);
        rhs = 0.0;

        forAll(equilibriumReactionIndices_, r)
        {
            const reactionRecord& rec = reactions_[equilibriumReactionIndices_[r]];
            nu = 0.0;
            forAll(rec.speciesIndices, j) nu[rec.speciesIndices[j]] = rec.stoichCoeffs[j];
            forAll(c_, i) concentrations[i] = c_[i][celli];

            scalar F=0, P=0, KR=0;
            rec.model->equilibriumResidualAndJacobian(celli, concentrations, nu, F, dFdc, P, KR);
            const scalar q = residualScale_[r][celli];
            rhs[r] = -F/q;

            for (label s=0; s<ne; ++s)
            {
                scalar value = 0.0;
                forAll(species_, i)
                    value += dFdc[i]*transportResponse_[i][celli]*equilibriumN_[i][s];
                J(r,s) = value/q;
            }
        }

        scalar minPivot = GREAT;
        if (!solveLocalSystem(J, rhs, delta, minPivot))
        {
            FatalErrorInFunction << "Singular local equilibrium response system at cell " << celli
                << " centre " << mesh_.C()[celli] << nl << "Scaled J_lambda=" << nl << J
                << nl << "scaled RHS=" << rhs << nl << "minimum relative pivot=" << minPivot
                << exit(FatalError);
        }

        forAll(lambda_, r)
        {
            lambda_[r][celli] += equilibriumRelaxation_*delta[r];
            if (!std::isfinite(lambda_[r][celli]))
                FatalErrorInFunction << "Non-finite lambda for "
                    << reactions_[equilibriumReactionIndices_[r]].name << " at cell " << celli
                    << exit(FatalError);
        }
    }

    // The corrected lambda has not yet been transported/tested.
    equilibriumEvaluated_ = false;
    equilibriumConverged_ = false;
}


void reactionManager::reportEquilibrium() const
{
    if (!equilibriumEvaluated_)
    {
        Info<< "Equilibrium residuals: no transported state currently evaluated" << nl;
        return;
    }

    Info<< "Equilibrium residuals:" << nl;
    forAll(equilibriumReactionIndices_, r)
    {
        scalar maxAbs = 0, maxRel = 0;
        forAll(residual_[r], celli)
        {
            maxAbs = max(maxAbs, mag(residual_[r][celli]));
            maxRel = max(maxRel, mag(residual_[r][celli])/residualScale_[r][celli]);
        }
        reduce(maxAbs, maxOp<scalar>());
        reduce(maxRel, maxOp<scalar>());
        Info<< "    " << reactions_[equilibriumReactionIndices_[r]].name
            << ": max|F|=" << maxAbs << ", maxRel=" << maxRel
            << ", absTol=" << equilibriumAbsTol_[r]
            << ", relTol=" << equilibriumRelTol_[r] << nl;
    }
}

} // End namespace Foam

