//---------------------------------------------------------------------------
//    $Id: simple_mpi_01.cc 23327 2011-02-11 03:19:07Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// PETScWrappers: document bug with PETSc SparseMatrix and clear_rows()
// until now, also the PETSc-internal SparsityPattern removes the
// rows that are emptied with clear_rows(). This results in errors
// when reusing the matrix later.

#include "../tests.h"

#include <lac/compressed_simple_sparsity_pattern.h>
#include <lac/petsc_sparse_matrix.h>    
#include <lac/petsc_parallel_sparse_matrix.h>    
#include <base/logstream.h>
#include <base/utilities.h>

#include <fstream>
//#include <mpi.h>


void test()
{
  

  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::System::get_n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;


  CompressedSimpleSparsityPattern csp(2*numprocs);
  for (unsigned int i=0;i<numprocs*2;++i)
    csp.add(i,i);
  csp.add(1,0);
  
  PETScWrappers::MPI::SparseMatrix mat;
  std::vector< unsigned int > local_rows(numprocs,2);
  
  mat.reinit(MPI_COMM_WORLD, csp, local_rows, local_rows, myid);
  
  mat.add(2*myid,2*myid,1.0);
  mat.add(2*myid+1,2*myid+1,1.0);
  mat.add(1,0,42.0);
  
  mat.add((2*myid+2)%(2*numprocs),(2*myid+2)%(2*numprocs),0.1);
  
  mat.compress();

  std::vector<unsigned int> rows(1,1);
  mat.clear_rows(rows);
  
//    mat.write_ascii();
  if (myid==0)
    deallog << "2nd try" << std::endl;
  
  mat = 0;
  mat.add(1,0,42.0);
  mat.add(2*myid,2*myid,1.0);
  mat.add(2*myid+1,2*myid+1,1.0);
  
  mat.add((2*myid+2)%(2*numprocs),(2*myid+2)%(2*numprocs),0.1);
  
  mat.compress();
//    mat.write_ascii();

  if (myid==0)
    deallog << "done" << std::endl;

}


int main(int argc, char *argv[])
{
  PetscInitialize(&argc,&argv,0,0);

  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::ofstream logfile(output_file_for_mpi("petsc_01").c_str());
//      deallog.attach(logfile);
//      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

  PetscFinalize();
}
