#ifndef BISICLES_THICKNESS_CONSTRAINT_HPP
#define BISICLES_THICKNESS_CONSTRAINT_HPP

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include "BoxIterator.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "MayDay.H"
#include "Vector.H"
#include "NamespaceHeader.H"

// The PFASST unknown is cell-average thickness on a fixed AMR hierarchy.
// Its implicit term is the reaction enforcing H >= 0, not an ordinary
// single-valued RHS: H = rhs + dtq*g, H >= 0, g >= 0, H*g = 0.
// Only call project on solution states. Corrections, quadrature integrals,
// explicit tendencies, and residuals must retain their signs.
namespace BisiclesThicknessConstraint
{
inline void fail(const char* reason)
{
    const std::string message = std::string("PFASST thickness constraint: ") + reason;
    MayDay::Error(message.c_str());
    std::abort(); // MayDay::Error must not return, including in release builds.
}

typedef Vector<LevelData<FArrayBox>* > Fields;

inline void requireHierarchy(const Fields& fields)
{
    if (fields.size() == 0) fail("empty thickness hierarchy");
    for (int lev = 0; lev < static_cast<int>(fields.size()); ++lev)
    {
        if (fields[lev] == NULL || !fields[lev]->isDefined())
            fail("missing or undefined thickness level");
        if (fields[lev]->nComp() != 1)
            fail("thickness must have exactly one component");
        for (int other = 0; other < lev; ++other)
            if (fields[lev] == fields[other])
                fail("AMR levels must use independent storage");
    }
}

inline void requireSameLayout(const LevelData<FArrayBox>& a,
                              const LevelData<FArrayBox>& b)
{
    // Matching DataIndex ownership/order is required for paired local access.
    // Ghost widths may differ: the constraint only acts on valid unknowns.
    if (!(a.getBoxes() == b.getBoxes()))
        fail("constraint vectors have different box layouts");
}

inline void requireFinite(Real value)
{
    if (!std::isfinite(value)) fail("nonfinite thickness, tendency, or reaction");
}

inline void requireSameHierarchy(const Fields& a, const Fields& b)
{
    if (a.size() != b.size())
        fail("constraint vectors have different AMR level counts");
    for (int lev = 0; lev < static_cast<int>(a.size()); ++lev)
        requireSameLayout(*a[lev], *b[lev]);
}

inline void requireIndependentReaction(const Fields& g, const Fields& a,
                                       const Fields& b)
{
    for (int lev = 0; lev < static_cast<int>(g.size()); ++lev)
        for (int other = 0; other < static_cast<int>(a.size()); ++other)
            if (g[lev] == a[other] || g[lev] == b[other])
                fail("reaction storage must not alias solution, RHS, or tendency");
}

inline void project(const Fields& fields)
{
    requireHierarchy(fields);
    for (int lev = 0; lev < static_cast<int>(fields.size()); ++lev)
    {
        LevelData<FArrayBox>& h = *fields[lev];
        const DisjointBoxLayout& grids = h.disjointBoxLayout();
        for (DataIterator dit = h.dataIterator(); dit.ok(); ++dit)
        {
            FArrayBox& fab = h[dit()];
            for (BoxIterator bit(grids[dit()]); bit.ok(); ++bit)
            {
                Real& value = fab(bit(), 0);
                requireFinite(value);
                if (value <= 0.0) value = 0.0;
            }
        }
        // Only same-level/periodic ghosts are exchanged here. BISICLES fills
        // coarse-fine and physical ghosts when installing the complete state.
        h.exchange();
    }
}

inline void solve(const Fields& yFields, const Fields& rhsFields,
                  const Fields& gFields, Real dtq)
{
    if (!std::isfinite(dtq) || dtq <= 0.0)
        fail("implicit substep dtq must be finite and positive");
    requireHierarchy(yFields);
    requireHierarchy(rhsFields);
    requireHierarchy(gFields);
    requireSameHierarchy(yFields, rhsFields);
    requireSameHierarchy(yFields, gFields);
    requireIndependentReaction(gFields, yFields, rhsFields);
    // An in-place solve is legal on corresponding levels only. Cross-level
    // aliases could overwrite an RHS before that level has been processed.
    for (int lev = 0; lev < static_cast<int>(yFields.size()); ++lev)
        for (int other = 0; other < static_cast<int>(rhsFields.size()); ++other)
            if (lev != other && yFields[lev] == rhsFields[other])
                fail("solution must not alias another AMR level of the RHS");

    for (int lev = 0; lev < static_cast<int>(yFields.size()); ++lev)
    {
        LevelData<FArrayBox>& y = *yFields[lev];
        const LevelData<FArrayBox>& rhs = *rhsFields[lev];
        LevelData<FArrayBox>& g = *gFields[lev];
        const DisjointBoxLayout& grids = y.disjointBoxLayout();
        for (DataIterator dit = y.dataIterator(); dit.ok(); ++dit)
        {
            FArrayBox& yf = y[dit()];
            const FArrayBox& rf = rhs[dit()];
            FArrayBox& gf = g[dit()];
            for (BoxIterator bit(grids[dit()]); bit.ok(); ++bit)
            {
                // Save rhs before writing y, so an in-place solution is safe.
                const Real r = rf(bit(), 0);
                requireFinite(r);
                const Real reaction = r < 0.0 ? -r / dtq : 0.0;
                requireFinite(reaction);
                yf(bit(), 0) = r > 0.0 ? r : 0.0;
                gf(bit(), 0) = reaction;
            }
        }
        y.exchange();
        g.exchange();
    }
    // Do not average down either result: every stored valid cell must still
    // satisfy y = rhs + dtq*g and y*g = 0, including covered coarse cells.
}

inline void refresh(const Fields& yFields, const Fields& fFields,
                    const Fields& gFields, bool atStart)
{
    requireHierarchy(yFields);
    requireHierarchy(fFields);
    requireHierarchy(gFields);
    requireSameHierarchy(yFields, fFields);
    requireSameHierarchy(yFields, gFields);
    requireIndependentReaction(gFields, yFields, fFields);

    for (int lev = 0; lev < static_cast<int>(yFields.size()); ++lev)
    {
        const LevelData<FArrayBox>& y = *yFields[lev];
        const LevelData<FArrayBox>& f = *fFields[lev];
        LevelData<FArrayBox>& g = *gFields[lev];
        const DisjointBoxLayout& grids = y.disjointBoxLayout();
        for (DataIterator dit = y.dataIterator(); dit.ok(); ++dit)
        {
            const FArrayBox& yf = y[dit()];
            const FArrayBox& ff = f[dit()];
            FArrayBox& gf = g[dit()];
            for (BoxIterator bit(grids[dit()]); bit.ok(); ++bit)
            {
                const Real h = yf(bit(), 0);
                const Real raw = ff(bit(), 0);
                requireFinite(h);
                requireFinite(raw);
                if (h < 0.0) fail("project the solution before refreshing its reaction");
                if (h > 0.0)
                    gf(bit(), 0) = 0.0;
                else if (atStart)
                    // At the left endpoint there is no preceding implicit solve.
                    // Choose the tangent reaction so H'=f+g cannot point below 0.
                    gf(bit(), 0) = std::max(-raw, Real(0.0));
                else
                {
                    // At later nodes g is a solved/iterated multiplier. Replacing
                    // it by max(-f,0) would destroy the implicit substep equation.
                    const Real oldReaction = gf(bit(), 0);
                    requireFinite(oldReaction);
                    gf(bit(), 0) = std::max(oldReaction, Real(0.0));
                }
            }
        }
        g.exchange();
    }
}
}

#include "NamespaceFooter.H"
#endif
