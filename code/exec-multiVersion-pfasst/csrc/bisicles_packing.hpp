#ifndef BISICLES_PACKING_HPP
#define BISICLES_PACKING_HPP

#include <algorithm>
#include <limits>
#include <vector>
#include <mpi.h>
#include "Vector.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "MayDay.H"
#include "NamespaceHeader.H"

// Count valid unknowns owned by this spatial rank; ghosts are not sent.
inline int BisiclesLocalVectorSize(const Vector<LevelData<FArrayBox>* >& H)
{
    long long count = 0;
    for (int lev = 0; lev < H.size(); ++lev)
    {
        if (H[lev] == NULL) MayDay::Error("Null PFASST thickness level");
        const LevelData<FArrayBox>& data = *H[lev];
        const DisjointBoxLayout& grids = data.disjointBoxLayout();
        for (DataIterator dit = data.dataIterator(); dit.ok(); ++dit)
        {
            count += static_cast<long long>(grids[dit()].numPts()) * data.nComp();
            if (count > std::numeric_limits<int>::max())
                MayDay::Error("PFASST local vector exceeds C/MPI integer range");
        }
    }
    return static_cast<int>(count);
}

// Corresponding spatial ranks on different time slices must own identical
// boxes in identical order.  Check once, before PFASST communicates vectors.
inline int BisiclesPfasstBufferSize(const Vector<LevelData<FArrayBox>* >& H,
                                    MPI_Comm time_comm)
{
    const int local_count = BisiclesLocalVectorSize(H);
    std::vector<int> layout;
    layout.push_back(H.size());
    layout.push_back(local_count);
    for (int lev = 0; lev < H.size(); ++lev)
    {
        const LevelData<FArrayBox>& data = *H[lev];
        const DisjointBoxLayout& grids = data.disjointBoxLayout();
        layout.push_back(data.nComp());
        int nboxes = 0;
        for (DataIterator dit = data.dataIterator(); dit.ok(); ++dit) ++nboxes;
        layout.push_back(nboxes);
        for (DataIterator dit = data.dataIterator(); dit.ok(); ++dit)
        {
            const Box& box = grids[dit()];
            for (int d = 0; d < SpaceDim; ++d)
            {
                layout.push_back(box.smallEnd(d));
                layout.push_back(box.bigEnd(d));
            }
        }
    }
    if (layout.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        MayDay::Error("PFASST layout metadata exceeds MPI integer range");
    const int length = static_cast<int>(layout.size());
    int min_length, max_length;
    MPI_Allreduce(&length, &min_length, 1, MPI_INT, MPI_MIN, time_comm);
    MPI_Allreduce(&length, &max_length, 1, MPI_INT, MPI_MAX, time_comm);
    if (min_length != max_length)
        MayDay::Error("PFASST time replicas have different local box layouts");

    std::vector<int> lo(length), hi(length);
    MPI_Allreduce(layout.data(), lo.data(), length, MPI_INT, MPI_MIN, time_comm);
    MPI_Allreduce(layout.data(), hi.data(), length, MPI_INT, MPI_MAX, time_comm);
    if (lo != hi)
        MayDay::Error("PFASST time replicas have different local box layouts/order");

    // LibPFASST requires a positive buffer length, even for an empty rank.
    return std::max(1, local_count);
}

#include "NamespaceFooter.H"
#endif