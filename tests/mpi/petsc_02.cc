// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



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
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;

  CompressedSimpleSparsityPattern csp(2);
  csp.add(0,0);
  csp.add(1,1);

  PETScWrappers::MPI::SparseMatrix mat;
  std::vector<types::global_dof_index> local_rows(numprocs, 0);
  local_rows[0]=2;

  mat.reinit(MPI_COMM_WORLD, csp, local_rows, local_rows, myid);

  if (myid == 0 )
    mat.add(0, 0, 1.0);


  mat.compress(VectorOperation::add);

  if (myid==0)
    deallog << "done" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();


}
