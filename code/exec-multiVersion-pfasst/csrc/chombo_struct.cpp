#include "chombo_struct.hpp"

void ChomboStruct::SetComm(MPI_Comm in_comm)
{
   comm = in_comm;
   MPI_Comm_rank(comm, &myid);
   MPI_Comm_size(comm, &num_procs);
}

void ChomboStruct:: InitGrid(int num_grid_points)
{
   num_of_rows = num_grid_points*num_grid_points; // currently only for 2d
}


int ChomboStruct::GetNumRows(void)
{
   return num_of_rows;
}
