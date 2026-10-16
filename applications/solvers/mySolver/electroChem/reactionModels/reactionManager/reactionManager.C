#include "reactionManager.H"
#include "potentialManager.H"
#include "IOobject.H"
#include "dimensionedScalar.H"
#include "DynamicList.H"
#include "token.H"
#include "PstreamReduceOps.H"
#include <cmath>

#include <vector>
#include <algorithm>
#include <limits>

namespace
{
namespace equilibriumKrylov
{
// Restarted GMRES. apply is the right-preconditioned operator; dot must
// perform the global reduction when the vectors are distributed.
template<class Real, class Apply, class Dot>
bool solve
(
    const std::vector<Real>& b, std::vector<Real>& x,
    Apply apply, Dot dot, const int restart, const int maxIter,
    const Real relTol, const Real absTol, int& iterations, Real& finalNorm
)
{
    using Vec = std::vector<Real>;
    x.assign(b.size(), Real(0));
    iterations = 0;
    const Real bNorm = std::sqrt(std::max(dot(b,b), Real(0)));
    const Real target = std::max(absTol, relTol*bNorm);
    finalNorm = bNorm;
    if (finalNorm <= target) return true;
    while (iterations < maxIter)
    {
        Vec ax = apply(x), r(b.size());
        for (size_t i=0; i<b.size(); ++i) r[i] = b[i]-ax[i];
        const Real beta = std::sqrt(std::max(dot(r,r), Real(0)));
        finalNorm = beta;
        if (!std::isfinite(beta)) return false;
        if (beta <= target) return true;
        const int m = std::min(restart, maxIter-iterations);
        std::vector<Vec> v(m+1, Vec(b.size(), Real(0)));
        std::vector<Vec> h(m+1, Vec(m, Real(0)));
        Vec cs(m), sn(m), g(m+1, Real(0));
        for (size_t i=0; i<b.size(); ++i) v[0][i] = r[i]/beta;
        g[0] = beta;
        int used = 0;
        bool breakdown = false;
        for (int j=0; j<m; ++j)
        {
            Vec w = apply(v[j]);
            // Two-pass modified Gram-Schmidt improves orthogonality.
            for (int pass=0; pass<2; ++pass)
                for (int k=0; k<=j; ++k)
                {
                    const Real a = dot(w,v[k]);
                    h[k][j] += a;
                    for (size_t i=0; i<w.size(); ++i) w[i] -= a*v[k][i];
                }
            h[j+1][j] = std::sqrt(std::max(dot(w,w), Real(0)));
            if (!std::isfinite(h[j+1][j])) return false;
            breakdown = h[j+1][j] <= std::numeric_limits<Real>::min();
            if (!breakdown)
                for (size_t i=0; i<w.size(); ++i) v[j+1][i] = w[i]/h[j+1][j];
            for (int k=0; k<j; ++k)
            {
                const Real a = cs[k]*h[k][j] + sn[k]*h[k+1][j];
                h[k+1][j] = -sn[k]*h[k][j] + cs[k]*h[k+1][j];
                h[k][j] = a;
            }
            const Real denom = std::hypot(h[j][j], h[j+1][j]);
            if (!std::isfinite(denom) || denom <= std::numeric_limits<Real>::min())
                return false;
            cs[j] = h[j][j]/denom;
            sn[j] = h[j+1][j]/denom;
            h[j][j] = denom;
            h[j+1][j] = 0;
            g[j+1] = -sn[j]*g[j];
            g[j] *= cs[j];
            ++iterations;
            used = j+1;
            if (std::abs(g[j+1]) <= target || breakdown) break;
        }
        Vec y(used);
        for (int j=used; j-- > 0;)
        {
            Real value = g[j];
            for (int k=j+1; k<used; ++k) value -= h[j][k]*y[k];
            y[j] = value/h[j][j];
            if (!std::isfinite(y[j])) return false;
        }
        for (int j=0; j<used; ++j)
            for (size_t i=0; i<x.size(); ++i) x[i] += y[j]*v[j][i];
        // Verify the true residual, not just the Arnoldi estimate.
        ax = apply(x);
        for (size_t i=0; i<r.size(); ++i) r[i] = b[i]-ax[i];
        finalNorm = std::sqrt(std::max(dot(r,r), Real(0)));
        if (!std::isfinite(finalNorm)) return false;
        if (finalNorm <= target) return true;
        if (breakdown) return false;
    }
    return false;
}
}
} // anonymous namespace

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
    responseSolver_(), transportMatrices_(),
    coupledRestart_(20), coupledMaxIter_(100),
    coupledRelTol_(1e-6), coupledAbsTol_(1e-12),
    coupledMaxBacktracks_(30), coupledReduction_(0.5),
    coupledSufficientDecrease_(1e-4),
    outputLevel_(1), equilibriumDiagnostics_(false), diagnosticStepValid_(false),
    diagnosticBefore_(), diagnosticDeltaC_(), diagnosticResponse_(), diagnosticKineticBefore_(),
    diagnosticDeltaLambda_(), diagnosticAlpha_(), diagnosticLimiter_(),
    diagnosticTransported_(),
    equilibriumGram_(0), equilibriumN_(), equilibriumEvaluated_(false),
    equilibriumConverged_(false)
{
    readOutputControls();
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
    if (eqControls.found("response") || eqControls.found("diagnostics"))
        FatalIOErrorInFunction(eqControls)
            << "Remove response and diagnostics: coupled Newton is mandatory. "
            << "Set picardOutput minimal/picard/full/debug in controlDict."
            << exit(FatalIOError);

    responseSolver_.add("solver", word("PBiCGStab"));
    responseSolver_.add("preconditioner", word("DILU"));
    responseSolver_.add("tolerance", scalar(1e-12));
    responseSolver_.add("relTol", scalar(0));
    responseSolver_.add("minIter", label(1));
    responseSolver_.add("maxIter", label(1000));
    if (eqControls.found("responseSolver"))
        responseSolver_.merge(eqControls.subDict("responseSolver"));
    transportMatrices_.setSize(species_.size());

    dictionary k;
    if (eqControls.found("coupledSolver")) k = eqControls.subDict("coupledSolver");
    coupledRestart_ = k.lookupOrDefault<label>("restart", 20);
    coupledMaxIter_ = k.lookupOrDefault<label>("maxIter", 100);
    coupledRelTol_ = k.lookupOrDefault<scalar>("relTol", 1e-6);
    coupledAbsTol_ = k.lookupOrDefault<scalar>("absTol", 1e-12);
    if (k.found("maxBacktracks") || k.found("reduction") || k.found("sufficientDecrease"))
        FatalIOErrorInFunction(k)
            << "Move maxBacktracks, reduction, sufficientDecrease to equilibrium/lineSearch"
            << exit(FatalIOError);
    dictionary search;
    if (eqControls.found("lineSearch")) search = eqControls.subDict("lineSearch");
    coupledMaxBacktracks_ = search.lookupOrDefault<label>("maxBacktracks", 30);
    coupledReduction_ = search.lookupOrDefault<scalar>("reduction", 0.5);
    coupledSufficientDecrease_ = search.lookupOrDefault<scalar>("sufficientDecrease", 1e-4);
    if (coupledRestart_ < 1 || coupledMaxIter_ < 1
     || !std::isfinite(coupledRelTol_) || coupledRelTol_ <= 0 || coupledRelTol_ >= 1
     || !std::isfinite(coupledAbsTol_) || coupledAbsTol_ <= 0)
        FatalIOErrorInFunction(eqControls)
            << "Require restart,maxIter>=1, 0<relTol<1, absTol>0"
            << exit(FatalIOError);
    if (coupledMaxBacktracks_ < 0
     || !std::isfinite(coupledReduction_) || coupledReduction_ <= 0 || coupledReduction_ >= 1
     || !std::isfinite(coupledSufficientDecrease_)
     || coupledSufficientDecrease_ <= 0 || coupledSufficientDecrease_ >= 1)
        FatalIOErrorInFunction(eqControls)
            << "Require maxBacktracks>=0, 0<reduction,sufficientDecrease<1"
            << exit(FatalIOError);
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
         || entries[i] == "lineSearch"
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

    applyTransportMatrixCorrection(coupledNewtonDirection());
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
    if (outputLevel_ >= 2) Info<< "Coupled equilibrium Newton: Krylov iterations=" << iterations
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

        // Homogeneous perturbation of the frozen transport operator.
        trial.primitiveFieldRef() = 0.0;
        forAll(trial.boundaryField(), patchi)
        {
            trial.boundaryFieldRef()[patchi] == scalar(0);
            if (!trial.boundaryField()[patchi].coupled())
                responseEqn.boundaryCoeffs()[patchi] = 0.0;
        }
        responseEqn.source() = 0.0;
        forAll(direction, r)
            responseEqn.source() +=
                mesh_.V()*equilibriumN_[i][r]*direction[r];
        {
            scalar maxSource = 0.0;
            forAll(responseEqn.source(), celli)
                maxSource = max(maxSource, mag(responseEqn.source()[celli]));
            reduce(maxSource, maxOp<scalar>());
            if (maxSource == 0) continue;
        }
        if (outputLevel_ >= 2)
            Info<< "Equilibrium transport response solve: species="
                << species_[i] << nl;
        dictionary controls(responseSolver_);
        controls.set("log", label(outputLevel_ >= 2 ? 1 : 0));
        const auto performance = responseEqn.solve(controls);
        bool responseOK = performance.converged();
        reduce(responseOK, andOp<bool>());
        if (!responseOK)
            FatalErrorInFunction << "Neighbour response solve failed for "
                << species_[i] << ": final residual=" << performance.finalResidual()
                << ". Check responseSolver and transport matrix."
                << exit(FatalError);
        response[i] = trial.primitiveField();
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
        if (outputLevel_ >= 2) Info<< "Coupled equilibrium predicted line search: backtracks="
            << backtracks << ", base scaled norm=" << std::sqrt(baseSquared)
            << ", trial scaled norm=" << std::sqrt(trialSquared)
            << ", accepted=" << accepted << nl;
        if (!accepted)
            FatalErrorInFunction << "Coupled predicted line search failed. "
                << "No reaction correction was applied." << exit(FatalError);
    }
    if (outputLevel_ >= 2) Info<< "Equilibrium neighbour response: global alpha=" << alpha
        << ", cap=" << equilibriumRelaxation_ << nl;
    if (outputLevel_ >= 2 && localLimiter >= 0 && localAlpha == alpha)
        Pout<< "Equilibrium global limiter: species=" << species_[localLimiter]
            << " localCell=" << localLimiterCell
            << " centre=" << mesh_.C()[localLimiterCell] << nl;

    scalar maxAbsDeltaLambda = 0.0;
    forAll(direction, r)
        forAll(direction[r], celli)
            maxAbsDeltaLambda = max(maxAbsDeltaLambda, mag(direction[r][celli]));
    reduce(maxAbsDeltaLambda, maxOp<scalar>());

    Info<< "    max|deltalambda|=" << maxAbsDeltaLambda << " mol/m3/s"
        << ", relaxation=" << alpha
        << ", applied max|alpha*deltalambda|="
        << alpha*maxAbsDeltaLambda << " mol/m3/s" << nl;

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


void reactionManager::readOutputControls()
{
    outputLevel_ = picardOutputLevel(mesh_.time().controlDict());
    equilibriumDiagnostics_ = outputLevel_ == 3;
}

void reactionManager::reportNoEquilibriumCorrection(const char* reason) const
{
    if (nEquilibriumReactions() && !equilibriumConverged_)
        Info<< "    relaxation factor = not applied (" << reason << ')' << nl;
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

    Info<< "    Equilibrium residuals: converged="
        << (equilibriumConverged_ ? "yes" : "no") << nl;
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

