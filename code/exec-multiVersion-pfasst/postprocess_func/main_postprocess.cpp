// #include "compute_spatial_L2norm_pf_ref.H"
// #include "read_hdf5.H"
#include <H5Cpp.h>
#include <iostream>
#include <string>
#include <vector>
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif


const H5std_string FILE_NAME( "ref_stream.L1L2.0512.petsc.000023.3d.hdf5" );
const H5std_string DATASET_NAME( "thickness" );


int main(int argc, char *argv[])
{


	H5::H5File file( FILE_NAME, H5F_ACC_RDONLY );

	H5::DataSet dataset = file.openDataSet( DATASET_NAME );

}