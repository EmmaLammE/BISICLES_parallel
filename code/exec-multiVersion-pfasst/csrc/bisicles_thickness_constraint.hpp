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

// The PFASST unknown is cell-average thickness on a fixed level-0 grid.
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

inline LevelData<FArrayBox>& checkedLevel(const Fields& fields)
{
    if (fields.size() != 1)
        fail("only a fixed single AMR level is supported");
    if (fields[0] == NULL || !fields[0]->isDefined())
        fail("missing or undefined thickness level");
    if (fields[0]->nComp() != 1)
        fail("thickness must have exactly one component");
    return *fields[0];
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

inline void project(const Fields& fields)
{
    LevelData<FArrayBox>& h = checkedLevel(fields);
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
    // Refresh inter-box/periodic ghosts. Physical boundary ghosts are not
    // unknowns; BISICLES must fill its geometry/calving boundary conditions.
    // In particular, DomainEdgeCalvingModel removes exterior ghosts, not the
    // positive valid cells adjacent to a domain edge.
    h.exchange();
}

inline void solve(const Fields& yFields, const Fields& rhsFields,
                  const Fields& gFields, Real dtq)
{
    if (!std::isfinite(dtq) || dtq <= 0.0)
        fail("implicit substep dtq must be finite and positive");
    LevelData<FArrayBox>& y = checkedLevel(yFields);
    const LevelData<FArrayBox>& rhs = checkedLevel(rhsFields);
    LevelData<FArrayBox>& g = checkedLevel(gFields);
    requireSameLayout(y, rhs);
    requireSameLayout(y, g);
    if (&g == &y || &g == &rhs)
        fail("reaction storage must not alias solution or implicit RHS");

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

inline void refresh(const Fields& yFields, const Fields& fFields,
                    const Fields& gFields, bool atStart)
{
    const LevelData<FArrayBox>& y = checkedLevel(yFields);
    const LevelData<FArrayBox>& f = checkedLevel(fFields);
    LevelData<FArrayBox>& g = checkedLevel(gFields);
    requireSameLayout(y, f);
    requireSameLayout(y, g);
    if (&g == &y || &g == &f)
        fail("reaction storage must not alias solution or explicit tendency");

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

#include "NamespaceFooter.H"
#endif