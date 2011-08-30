//---------------------------------------------------------------------------
//    $Id: simple_mpi_01.cc 23327 2011-02-11 03:19:07Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2009, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// PETScWrappers: document bug with PETSc SparseMatrix. If only one CPU
// does matrix-assembly, it calls compress() inside and the others don't.
// We should implement this like in PETSc::MPI::Vector.

#include "../tests.h"

#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <fstream>
//#include <mpi.h>


void test()
{
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  unsigned int numprocs = Utilities::System::get_n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;

  CompressedSimpleSparsityPattern csp(2);
  csp.add(0,0);
  csp.add(1,1);

  PETScWrappers::MPI::SparseMatrix mat;
  std::vector< unsigned int > local_rows(numprocs, 0);
  local_rows[0]=2;

  mat.reinit(MPI_COMM_WORLD, csp, local_rows, local_rows, myid);

  if (myid == 0 )
    mat.add(0, 0, 1.0);


  mat.compress();

  if (myid==0)
    deallog << "done" << std::endl;
}


int main(int argc, char *argv[])
{
  PetscInitialize(&argc,&argv,0,0);

  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::ofstream logfile(output_file_for_mpi("petsc_02").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

  PetscFinalize();
}
