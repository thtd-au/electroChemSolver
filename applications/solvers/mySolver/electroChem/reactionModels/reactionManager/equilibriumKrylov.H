#ifndef equilibriumKrylov_H
#define equilibriumKrylov_H
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

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
#endif
