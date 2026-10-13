#include "reactionManager.H"
#include "IOobject.H"
#include "dimensionedScalar.H"
#include "DynamicList.H"
#include "token.H"
#include "equilibriumKrylov.H"
#include "PstreamReduceOps.H"
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
    equilibriumFBL_(false), equilibriumFractionToBoundary_(0.99),
    equilibriumResponse_("diagonal"), responseSolver_(), transportMatrices_(),
    coupledRestart_(20), coupledMaxIter_(100),
    coupledRelTol_(1e-6), coupledAbsTol_(1e-12),
    coupledMaxBacktracks_(30), coupledReduction_(0.5),
    coupledSufficientDecrease_(1e-4),
    equilibriumDiagnostics_(false), diagnosticStepValid_(false),
    diagnosticBefore_(), diagnosticDeltaC_(), diagnosticResponse_(), diagnosticKineticBefore_(),
    diagnosticDeltaLambda_(), diagnosticAlpha_(), diagnosticLimiter_(),
    diagnosticTransported_(),
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
    if (!nEquilibriumReactions()) return;
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
    equilibriumResponse_ =
        eqControls.lookupOrDefault<word>("response", "diagonal");
    if (equilibriumResponse_ != "diagonal"
     && equilibriumResponse_ != "transportMatrix"
     && equilibriumResponse_ != "coupledNewton")
    {
        FatalIOErrorInFunction(eqControls)
            << "response must be diagonal, transportMatrix, or coupledNewton"
            << exit(FatalIOError);
    }
    if (equilibriumResponse_ != "diagonal")
    {
        responseSolver_ = eqControls.subDict("responseSolver");
        transportMatrices_.setSize(species_.size());
    }
    if (equilibriumResponse_ == "coupledNewton")
    {
        const dictionary& k = eqControls.subDict("coupledSolver");
        coupledRestart_ = k.lookupOrDefault<label>("restart", 20);
        coupledMaxIter_ = k.lookupOrDefault<label>("maxIter", 100);
        coupledRelTol_ = k.lookupOrDefault<scalar>("relTol", 1e-6);
        coupledAbsTol_ = k.lookupOrDefault<scalar>("absTol", 1e-12);
        coupledMaxBacktracks_ = k.lookupOrDefault<label>("maxBacktracks", 30);
        coupledReduction_ = k.lookupOrDefault<scalar>("reduction", 0.5);
        coupledSufficientDecrease_ =
            k.lookupOrDefault<scalar>("sufficientDecrease", 1e-4);
        if (coupledRestart_ < 1 || coupledMaxIter_ < 1
         || !std::isfinite(coupledRelTol_) || coupledRelTol_ <= 0
         || coupledRelTol_ >= 1
         || !std::isfinite(coupledAbsTol_) || coupledAbsTol_ <= 0)
            FatalIOErrorInFunction(k)
                << "Require restart,maxIter >= 1, 0 < relTol < 1, absTol > 0"
                << exit(FatalIOError);
        if (coupledMaxBacktracks_ < 0
         || !std::isfinite(coupledReduction_)
         || coupledReduction_ <= 0 || coupledReduction_ >= 1
         || !std::isfinite(coupledSufficientDecrease_)
         || coupledSufficientDecrease_ <= 0
         || coupledSufficientDecrease_ >= 1)
            FatalIOErrorInFunction(k)
                << "Require maxBacktracks>=0, 0<reduction,sufficientDecrease<1"
                << exit(FatalIOError);
    }
    equilibriumDiagnostics_ =
        eqControls.lookupOrDefault<bool>("diagnostics", false);
    // Inspect the token so both a scalar and the FBL selector are accepted.
    ITstream& relaxationStream = eqControls.lookup("relaxation");
    const token relaxationToken(relaxationStream);
    if (relaxationToken.isNumber())
    {
        equilibriumRelaxation_ = relaxationToken.number();
    }
    else if (relaxationToken.isWord() && relaxationToken.wordToken() == "FBL")
    {
        equilibriumFBL_ = true;
        equilibriumRelaxation_ =
            eqControls.lookupOrDefault<scalar>("maxRelaxation", 1.0);
        equilibriumFractionToBoundary_ =
            eqControls.lookupOrDefault<scalar>("fractionToBoundary", 0.99);
    }
    else
    {
        FatalIOErrorInFunction(eqControls)
            << "equilibrium relaxation must be a scalar or FBL"
            << exit(FatalIOError);
    }
    if (!relaxationStream.eof())
    {
        FatalIOErrorInFunction(eqControls)
            << "equilibrium relaxation must contain exactly one value"
            << exit(FatalIOError);
    }
    if (!std::isfinite(equilibriumRelaxation_)
     || equilibriumRelaxation_ <= 0 || equilibriumRelaxation_ > 1)
    {
        FatalIOErrorInFunction(eqControls)
            << "equilibrium relaxation (or FBL maxRelaxation) must satisfy 0 < value <= 1"
            << exit(FatalIOError);
    }
    if (!std::isfinite(equilibriumFractionToBoundary_)
     || equilibriumFractionToBoundary_ <= 0
     || equilibriumFractionToBoundary_ >= 1)
    {
        FatalIOErrorInFunction(eqControls)
            << "fractionToBoundary must satisfy 0 < value < 1"
            << exit(FatalIOError);
    }
    if (!equilibriumFBL_
     && (eqControls.found("maxRelaxation") || eqControls.found("fractionToBoundary")))
    {
        FatalIOErrorInFunction(eqControls)
            << "maxRelaxation and fractionToBoundary require relaxation FBL"
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
        if (entries[i] == "relaxation"
         || entries[i] == "responseSolver"
         || entries[i] == "coupledSolver"
         || entries[i] == "maxRelaxation"
         || entries[i] == "fractionToBoundary"
         || !eqControls.isDict(entries[i])) continue;
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
    diagnosticStepValid_ = false;
    if (!nEquilibriumReactions()) return;
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


tmp<volScalarField> reactionManager::source
(
    const label speciesI,
    scalarField* kineticSource,
    scalarField* equilibriumSource
) const
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
    if (kineticSource) kineticSource->setSize(mesh_.nCells(), 0.0);
    if (equilibriumSource) equilibriumSource->setSize(mesh_.nCells(), 0.0);
    if (kineticSource) *kineticSource = 0.0;
    if (equilibriumSource) *equilibriumSource = 0.0;

    forAll(kineticReactionIndices_, k)
    {
        const scalar nu = reactionStoich(reactions_[kineticReactionIndices_[k]], speciesI);
        if (mag(nu) > SMALL)
        {
            src += nu*kineticRates_[k];
            if (kineticSource) *kineticSource += nu*kineticRates_[k];
        }
    }
    forAll(equilibriumReactionIndices_, r)
    {
        const scalar nu = equilibriumN_[speciesI][r];
        if (mag(nu) > SMALL)
        {
            src += nu*lambda_[r];
            if (equilibriumSource) *equilibriumSource += nu*lambda_[r];
        }
    }
    return tSource;
}


void reactionManager::setTransportResponse(const label speciesI, const fvScalarMatrix& cEqn)
{
    if (!nEquilibriumReactions()) return;
    if (speciesI < 0 || speciesI >= species_.size() || speciesI == inertIndex_)
        FatalErrorInFunction << "Invalid active species index " << speciesI << exit(FatalError);

    if (equilibriumResponse_ != "diagonal")
    {
        bool participating = false;
        forAll(equilibriumN_[speciesI], r)
            participating = participating || equilibriumN_[speciesI][r] != 0;
        if (participating)
            transportMatrices_.set(speciesI, new fvScalarMatrix(cEqn));
    }

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
    if (!nEquilibriumReactions()) return;
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
    if (!nEquilibriumReactions()) return true;
    if (!equilibriumEvaluated_)
        FatalErrorInFunction << "equilibriumConverged() called before evaluateEquilibrium()"
            << exit(FatalError);
    return equilibriumConverged_;
}


void reactionManager::correctEquilibrium()
{
    if (!nEquilibriumReactions()) return;
    if (!equilibriumEvaluated_)
        FatalErrorInFunction << "correctEquilibrium() called before evaluateEquilibrium()"
            << exit(FatalError);

    const label ne = nEquilibriumReactions();
    if (!ne || equilibriumConverged_) return;

    if (equilibriumDiagnostics_)
    {
        diagnosticBefore_.setSize(species_.size());
        diagnosticDeltaC_.setSize(species_.size());
        diagnosticResponse_ = transportResponse_;
        diagnosticKineticBefore_.setSize(species_.size());
        diagnosticDeltaLambda_.setSize(ne);
        diagnosticAlpha_.setSize(mesh_.nCells(), 0.0);
        diagnosticLimiter_.setSize(mesh_.nCells(), -1);
        diagnosticTransported_.setSize(species_.size(), false);
        forAll(species_, i)
        {
            diagnosticBefore_[i] = c_[i].primitiveField();
            diagnosticDeltaC_[i].setSize(mesh_.nCells(), 0.0);
            diagnosticDeltaC_[i] = 0.0;
            diagnosticKineticBefore_[i].setSize(mesh_.nCells(), 0.0);
            diagnosticKineticBefore_[i] = 0.0;
            forAll(kineticReactionIndices_, k)
                diagnosticKineticBefore_[i] +=
                    reactionStoich(reactions_[kineticReactionIndices_[k]], i)
                   *kineticRates_[k];
        }
        forAll(diagnosticDeltaLambda_, r)
        {
            diagnosticDeltaLambda_[r].setSize(mesh_.nCells(), 0.0);
            diagnosticDeltaLambda_[r] = 0.0;
        }
        diagnosticStepValid_ = true;
    }

    if (equilibriumResponse_ == "coupledNewton")
    {
        applyTransportMatrixCorrection(coupledNewtonDirection());
        equilibriumEvaluated_ = false;
        equilibriumConverged_ = false;
        return;
    }

    scalarField concentrations(species_.size(), 0.0), nu(species_.size(), 0.0), dFdc(species_.size(), 0.0);
    scalarSquareMatrix J(ne);
    scalarField rhs(ne), delta(ne);
    List<scalarField> direction(ne);
    if (equilibriumResponse_ == "transportMatrix")
        forAll(direction, r) direction[r].setSize(mesh_.nCells(), 0.0);

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

        if (equilibriumResponse_ == "transportMatrix")
        {
            // Keep the existing local Newton direction, but predict its
            // concentration response with the whole transport matrix below.
            forAll(delta, r) direction[r][celli] = delta[r];
            continue;
        }

        scalar alpha = equilibriumRelaxation_;
        label limitingSpecies = -1;
        if (equilibriumFBL_)
        {
            // One factor for the coupled step: deltaC = S*N*deltaLambda.
            // Limit only species moving towards the c=0 boundary.
            forAll(species_, i)
            {
                if (i == inertIndex_) continue;
                scalar deltaC = 0.0;
                forAll(delta, r)
                    deltaC += equilibriumN_[i][r]*delta[r];
                deltaC *= transportResponse_[i][celli];
                const scalar ci = c_[i][celli];
                if (!std::isfinite(deltaC) || !std::isfinite(ci) || ci < 0)
                {
                    FatalErrorInFunction
                        << "Invalid FBL concentration response for " << species_[i]
                        << " at cell " << celli << ": c=" << ci
                        << ", deltaC=" << deltaC << exit(FatalError);
                }
                if (deltaC < 0)
                {
                    // Compute the ratio only if the current cap would violate
                    // the safety margin; this also avoids a huge c/(-deltaC).
                    const scalar available = equilibriumFractionToBoundary_*ci;
                    if (alpha*(-deltaC) > available)
                    {
                        alpha = available/(-deltaC);
                        limitingSpecies = i;
                    }
                }
            }
            if (alpha <= 0)
            {
                FatalErrorInFunction
                    << "FBL equilibrium step is blocked at cell " << celli
                    << ": a zero concentration has a negative predicted increment. "
                    << "No positive fraction-to-boundary step is feasible."
                    << exit(FatalError);
            }
        }

        if (equilibriumDiagnostics_)
        {
            diagnosticAlpha_[celli] = alpha;
            diagnosticLimiter_[celli] = limitingSpecies;
            forAll(delta, r)
                diagnosticDeltaLambda_[r][celli] = alpha*delta[r];
            forAll(species_, i)
            {
                scalar dc = 0.0;
                forAll(delta, r) dc += equilibriumN_[i][r]*delta[r];
                diagnosticDeltaC_[i][celli] =
                    alpha*transportResponse_[i][celli]*dc;
            }
        }

        forAll(lambda_, r)
        {
            lambda_[r][celli] += alpha*delta[r];
            if (!std::isfinite(lambda_[r][celli]))
                FatalErrorInFunction << "Non-finite lambda for "
                    << reactions_[equilibriumReactionIndices_[r]].name << " at cell " << celli
                    << exit(FatalError);
        }
    }

    if (equilibriumResponse_ == "transportMatrix")
        applyTransportMatrixCorrection(direction);

    // The corrected lambda has not yet been transported/tested.
    equilibriumEvaluated_ = false;
    equilibriumConverged_ = false;
}


List<scalarField> reactionManager::coupledNewtonDirection() const
{
    const label ne = nEquilibriumReactions(), nc = mesh_.nCells();
    typedef std::vector<scalar> Vector;
    List<List<scalarField>> derivative(ne);
    List<scalarSquareMatrix> preconditioner(nc);
    Vector b(ne*nc, 0.0);
    scalarField concentrations(species_.size()), nu(species_.size()), dFdc;
    forAll(derivative, r)
    {
        derivative[r].setSize(species_.size());
        forAll(species_, i) derivative[r][i].setSize(nc, 0.0);
    }
    // Include ALL cells, even locally converged ones: neighbour corrections
    // can perturb their equilibrium residuals.
    for (label celli=0; celli<nc; ++celli)
    {
        preconditioner[celli].setSize(ne);
        preconditioner[celli] = scalar(0);
        forAll(c_, i) concentrations[i] = c_[i][celli];
        forAll(equilibriumReactionIndices_, r)
        {
            const reactionRecord& rec = reactions_[equilibriumReactionIndices_[r]];
            nu = 0.0;
            forAll(rec.speciesIndices, j) nu[rec.speciesIndices[j]] = rec.stoichCoeffs[j];
            scalar F=0, P=0, KR=0;
            rec.model->equilibriumResidualAndJacobian
                (celli, concentrations, nu, F, dFdc, P, KR);
            const scalar scale = residualScale_[r][celli];
            b[r*nc+celli] = -F/scale;
            forAll(species_, i)
            {
                if (!std::isfinite(dFdc[i]))
                    FatalErrorInFunction << "Non-finite coupled derivative at cell "
                        << celli << exit(FatalError);
                derivative[r][i][celli] = dFdc[i]/scale;
            }
            for (label s=0; s<ne; ++s)
                forAll(species_, i)
                    preconditioner[celli](r,s) += derivative[r][i][celli]
                        *transportResponse_[i][celli]*equilibriumN_[i][s];
        }
    }
    // Right preconditioning uses the existing local equilibrium Jacobian.
    const auto precondition = [&](const Vector& v)
    {
        Vector z(v.size(), 0.0);
        scalarField rhs(ne), solution;
        for (label celli=0; celli<nc; ++celli)
        {
            for (label r=0; r<ne; ++r) rhs[r] = v[r*nc+celli];
            scalar pivot=GREAT;
            if (!solveLocalSystem(preconditioner[celli], rhs, solution, pivot))
                FatalErrorInFunction
                    << "Singular coupled Newton preconditioner at cell "
                    << celli << exit(FatalError);
            for (label r=0; r<ne; ++r) z[r*nc+celli] = solution[r];
        }
        return z;
    };
    const auto apply = [&](const Vector& v)
    {
        const Vector z = precondition(v);
        List<scalarField> direction(ne);
        for (label r=0; r<ne; ++r)
        {
            direction[r].setSize(nc);
            for (label celli=0; celli<nc; ++celli)
                direction[r][celli] = z[r*nc+celli];
        }
        const List<scalarField> dc = transportMatrixResponse(direction);
        Vector result(v.size(), 0.0);
        for (label r=0; r<ne; ++r)
            forAll(species_, i)
                for (label celli=0; celli<nc; ++celli)
                    result[r*nc+celli] += derivative[r][i][celli]*dc[i][celli];
        return result;
    };
    const auto dot = [&](const Vector& a, const Vector& v)
    {
        scalar sum = 0.0;
        for (size_t i=0; i<a.size(); ++i) sum += a[i]*v[i];
        reduce(sum, sumOp<scalar>());
        return sum;
    };
    Vector y;
    int iterations = 0;
    scalar norm = GREAT;
    const bool ok = equilibriumKrylov::solve
    (
        b, y, apply, dot, coupledRestart_, coupledMaxIter_,
        coupledRelTol_, coupledAbsTol_, iterations, norm
    );
    const scalar initialNorm = std::sqrt(max(dot(b,b), scalar(0)));
    Info<< "Coupled equilibrium Newton: Krylov iterations=" << iterations
        << ", initial scaled norm=" << initialNorm
        << ", final true scaled norm=" << norm
        << ", converged=" << ok << nl;
    if (!ok)
        FatalErrorInFunction << "Coupled equilibrium Newton linear solve failed. "
            << "No reaction correction was applied."
            << exit(FatalError);
    const Vector z = precondition(y);
    List<scalarField> direction(ne);
    for (label r=0; r<ne; ++r)
    {
        direction[r].setSize(nc);
        for (label celli=0; celli<nc; ++celli)
            direction[r][celli] = z[r*nc+celli];
    }
    return direction;
}


List<scalarField> reactionManager::transportMatrixResponse
(
    const List<scalarField>& direction
) const
{
    List<scalarField> response(species_.size());
    forAll(species_, i)
    {
        response[i].setSize(mesh_.nCells(), 0.0);
        bool participating = false;
        forAll(equilibriumN_[i], r)
            participating = participating || equilibriumN_[i][r] != 0;
        if (!participating || i == inertIndex_) continue;
        if (!transportMatrices_.set(i))
            FatalErrorInFunction << "Missing transport matrix for " << species_[i]
                << exit(FatalError);
        const fvScalarMatrix& stored = transportMatrices_[i];

        // Use an unregistered clone with the same concentration name, so
        // species-dependent patch types retain their species identity.
        // The real concentration field is never modified by this solve.
        volScalarField trial
        (
            IOobject(c_[i].name(), mesh_.time().timeName(), mesh_,
                IOobject::NO_READ, IOobject::NO_WRITE, false),
            c_[i]
        );
        fvScalarMatrix responseEqn(trial, stored.dimensions());
        responseEqn.diag() = stored.diag();
        if (!stored.diagonal()) responseEqn.upper() = stored.upper();
        if (stored.asymmetric()) responseEqn.lower() = stored.lower();
        responseEqn.internalCoeffs() = stored.internalCoeffs();
        responseEqn.boundaryCoeffs() = stored.boundaryCoeffs();

        // residual = b - A*c, including boundary/interface contributions.
        // Subtract it to anchor the trial at the current transported state:
        // A*cTrial = A*c + V*N*direction. This avoids treating any baseline
        // linear-solve error as part of the reaction response.
        if (equilibriumResponse_ == "coupledNewton")
        {
            // Direct increment solve: avoid subtracting nearly equal absolute
            // concentrations in repeated Krylov matrix-vector products.
            trial.primitiveFieldRef() = 0.0;
            forAll(trial.boundaryField(), patchi)
            {
                trial.boundaryFieldRef()[patchi] == scalar(0);
                if (!trial.boundaryField()[patchi].coupled())
                    responseEqn.boundaryCoeffs()[patchi] = 0.0;
            }
            responseEqn.source() = 0.0;
        }
        else
        {
            responseEqn.source() = stored.source() - stored.residual()();
        }
        forAll(direction, r)
            responseEqn.source() +=
                mesh_.V()*equilibriumN_[i][r]*direction[r];
        if (equilibriumResponse_ == "coupledNewton")
        {
            scalar maxSource = 0.0;
            forAll(responseEqn.source(), celli)
                maxSource = max(maxSource, mag(responseEqn.source()[celli]));
            reduce(maxSource, maxOp<scalar>());
            if (maxSource == 0) continue;
        }
        if (equilibriumDiagnostics_)
            Info<< "Equilibrium transport response solve: species="
                << species_[i] << nl;
        const auto performance = responseEqn.solve(responseSolver_);
        bool responseOK = performance.converged();
        reduce(responseOK, andOp<bool>());
        if (!responseOK)
            FatalErrorInFunction << "Neighbour response solve failed for "
                << species_[i] << ": final residual=" << performance.finalResidual()
                << ". Check responseSolver and transport matrix."
                << exit(FatalError);
        if (equilibriumResponse_ == "coupledNewton")
            response[i] = trial.primitiveField();
        else
            response[i] = trial.primitiveField() - c_[i].primitiveField();
    }
    return response;
}

void reactionManager::applyTransportMatrixCorrection
(
    const List<scalarField>& direction
)
{
    const List<scalarField> response = transportMatrixResponse(direction);
    scalar alpha = equilibriumRelaxation_;
    label localLimiter = -1, localLimiterCell = -1;
    forAll(species_, i)
    {
        forAll(response[i], celli)
        {
            const scalar ci = c_[i][celli];
            const scalar dc = response[i][celli];
            if (!std::isfinite(ci) || ci < 0 || !std::isfinite(dc))
                FatalErrorInFunction << "Invalid neighbour response for "
                    << species_[i] << " at cell " << celli
                    << ": c=" << ci << ", deltaC=" << dc << exit(FatalError);
            if (equilibriumFBL_ && dc < 0
             && alpha*(-dc) > equilibriumFractionToBoundary_*ci)
            {
                alpha = equilibriumFractionToBoundary_*ci/(-dc);
                localLimiter = i;
                localLimiterCell = celli;
            }
        }
    }

    // A single scalar commutes with A^-1; cell-dependent damping does not.
    const scalar localAlpha = alpha;
    reduce(alpha, minOp<scalar>());
    if (!std::isfinite(alpha) || alpha <= 0)
        FatalErrorInFunction << "No positive global neighbour-response FBL step "
            << "is feasible (a zero species may have a negative response)."
            << exit(FatalError);
    if (equilibriumResponse_ == "coupledNewton")
    {
        scalar baseSquared = 0.0;
        forAll(residual_, r)
            forAll(residual_[r], celli)
                baseSquared += sqr(residual_[r][celli]/residualScale_[r][celli]);
        reduce(baseSquared, sumOp<scalar>());
        scalarField concentrations(species_.size()), nu(species_.size()), dFdc;
        bool accepted = false;
        scalar trialSquared = GREAT;
        label backtracks = 0;
        for (; backtracks <= coupledMaxBacktracks_; ++backtracks)
        {
            trialSquared = 0.0;
            bool valid = true;
            forAll(mesh_.V(), celli)
            {
                forAll(c_, i)
                {
                    concentrations[i] = c_[i][celli] + alpha*response[i][celli];
                    valid = valid && std::isfinite(concentrations[i])
                        && concentrations[i] >= 0;
                }
                if (!valid) break;
                forAll(equilibriumReactionIndices_, r)
                {
                    const reactionRecord& rec = reactions_[equilibriumReactionIndices_[r]];
                    nu = 0.0;
                    forAll(rec.speciesIndices, j)
                        nu[rec.speciesIndices[j]] = rec.stoichCoeffs[j];
                    scalar F=0, P=0, KR=0;
                    rec.model->equilibriumResidualAndJacobian
                        (celli, concentrations, nu, F, dFdc, P, KR);
                    valid = valid && std::isfinite(F);
                    if (!valid) break;
                    trialSquared += sqr(F/residualScale_[r][celli]);
                }
                if (!valid) break;
            }
            reduce(valid, andOp<bool>());
            reduce(trialSquared, sumOp<scalar>());
            if (valid && std::isfinite(trialSquared)
             && trialSquared <= sqr(1-coupledSufficientDecrease_*alpha)*baseSquared)
            {
                accepted = true;
                break;
            }
            if (backtracks < coupledMaxBacktracks_) alpha *= coupledReduction_;
        }
        Info<< "Coupled equilibrium predicted line search: backtracks="
            << backtracks << ", base scaled norm=" << std::sqrt(baseSquared)
            << ", trial scaled norm=" << std::sqrt(trialSquared)
            << ", accepted=" << accepted << nl;
        if (!accepted)
            FatalErrorInFunction << "Coupled predicted line search failed. "
                << "No reaction correction was applied." << exit(FatalError);
    }
    Info<< "Equilibrium neighbour response: global alpha=" << alpha
        << ", cap=" << equilibriumRelaxation_ << nl;
    if (localLimiter >= 0 && localAlpha == alpha)
        Pout<< "Equilibrium global limiter: species=" << species_[localLimiter]
            << " localCell=" << localLimiterCell
            << " centre=" << mesh_.C()[localLimiterCell] << nl;

    forAll(lambda_, r)
    {
        lambda_[r] += alpha*direction[r];
        forAll(lambda_[r], celli)
            if (!std::isfinite(lambda_[r][celli]))
                FatalErrorInFunction << "Non-finite neighbour-corrected lambda"
                    << exit(FatalError);
    }
    if (equilibriumDiagnostics_)
    {
        forAll(diagnosticAlpha_, celli)
        {
            diagnosticAlpha_[celli] = alpha;
            // -2 denotes a global limiter, possibly on another processor.
            diagnosticLimiter_[celli] = -2;
        }
        forAll(direction, r) diagnosticDeltaLambda_[r] = alpha*direction[r];
        forAll(response, i) diagnosticDeltaC_[i] = alpha*response[i];
    }
}


void reactionManager::beginTransportDiagnostics()
{
    if (equilibriumDiagnostics_ && diagnosticStepValid_)
        diagnosticTransported_ = false;
}


void reactionManager::reportTransportDiagnostics(const label speciesI)
{
    if (!equilibriumDiagnostics_) return;
    if (diagnosticStepValid_) diagnosticTransported_[speciesI] = true;
    const scalarField& actual = c_[speciesI].primitiveField();
    label minCell = -1, mismatchCell = -1;
    scalar minInternal = GREAT, maxMismatch = -1;
    forAll(actual, celli)
    {
        if (actual[celli] < minInternal)
        {
            minInternal = actual[celli];
            minCell = celli;
        }
        if (diagnosticStepValid_)
        {
            const scalar error = mag(actual[celli] - diagnosticBefore_[speciesI][celli]
                - diagnosticDeltaC_[speciesI][celli]);
            if (error > maxMismatch) { maxMismatch = error; mismatchCell = celli; }
        }
    }
    scalar minBoundary = GREAT;
    label minPatch = -1, minFace = -1;
    forAll(c_[speciesI].boundaryField(), patchi)
    {
        const fvPatchScalarField& patch = c_[speciesI].boundaryField()[patchi];
        // Processor patches are interfaces, not physical boundaries.
        if (patch.coupled()) continue;
        forAll(patch, facei)
        {
            if (patch[facei] < minBoundary)
            {
                minBoundary = patch[facei];
                minPatch = patchi;
                minFace = facei;
            }
        }
    }
    Pout<< "EQDIAG time=" << mesh_.time().timeName()
        << " species=" << species_[speciesI]
        << " localInternalMin=" << minInternal << " localCell=" << minCell
        << " localPhysicalBoundaryMin=" << minBoundary;
    if (minPatch >= 0)
        Pout<< " patch=" << mesh_.boundary()[minPatch].name()
            << " localFace=" << minFace;
    Pout<< nl;

    // Print the minimum cell and the worst prediction mismatch. All cell
    // indices are processor-local; Pout identifies the rank in parallel.
    for (label sample=0; sample<2; ++sample)
    {
        const label celli = sample == 0 ? minCell : mismatchCell;
        if (celli < 0 || (sample == 1 && celli == minCell)) continue;
        Pout<< "EQDIAG sample=" << (sample == 0 ? "minimum" : "maxMismatch")
            << " species=" << species_[speciesI]
            << " localCell=" << celli << " centre=" << mesh_.C()[celli];
        if (!diagnosticStepValid_)
        {
            Pout<< " noPreviousEquilibriumCorrection" << nl;
            continue;
        }
        Pout<< " alpha=" << diagnosticAlpha_[celli] << " limiter=";
        const label limiter = diagnosticLimiter_[celli];
        if (limiter >= 0) Pout<< species_[limiter];
        else if (limiter == -2) Pout<< "global(seeCorrectionReport)";
        else Pout<< "none/cap";
        Pout<< nl;
        forAll(lambda_, r)
            Pout<< "EQDIAG reaction="
                << reactions_[equilibriumReactionIndices_[r]].name
                << " appliedDeltaLambda=" << diagnosticDeltaLambda_[r][celli]
                << " lambdaUsed=" << lambda_[r][celli] << nl;
        forAll(species_, i)
        {
            bool participating = false;
            forAll(equilibriumN_[i], r)
                participating = participating || equilibriumN_[i][r] != 0;
            if (!participating) continue;
            scalar kineticNow = 0.0;
            forAll(kineticReactionIndices_, k)
                kineticNow += reactionStoich(reactions_[kineticReactionIndices_[k]], i)
                    *kineticRates_[k][celli];
            const scalar before = diagnosticBefore_[i][celli];
            const scalar predictedDC = diagnosticDeltaC_[i][celli];
            Pout<< "EQDIAG species=" << species_[i]
                << " before=" << before
                << " responseUsed=" << diagnosticResponse_[i][celli]
                << " predictedDeltaC=" << predictedDC
                << " predictedC=" << before + predictedDC
                << " kineticBefore=" << diagnosticKineticBefore_[i][celli]
                << " kineticUsed=" << kineticNow;
            if (diagnosticTransported_[i])
                Pout<< " actualC=" << c_[i][celli]
                    << " actualDeltaC=" << c_[i][celli] - before
                    << " predictionError=" << c_[i][celli] - before - predictedDC;
            else Pout<< " actual=pendingTransport";
            Pout<< nl;
        }
    }
}


void reactionManager::reportEquilibrium() const
{
    if (!nEquilibriumReactions()) return;
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

