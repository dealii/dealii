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



// PETScWrappers: document bug with PETSc SparseMatrix and clear_rows()
// until now, also the PETSc-internal SparsityPattern removes the
// rows that are emptied with clear_rows(). This results in errors
// when reusing the matrix later.

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

  CompressedSimpleSparsityPattern csp(2*numprocs);
  for (unsigned int i=0; i<numprocs*2; ++i)
    csp.add(i,i);
  csp.add(1,0);

  PETScWrappers::MPI::SparseMatrix mat;
  std::vector<types::global_dof_index> local_rows(numprocs,2);

  mat.reinit(MPI_COMM_WORLD, csp, local_rows, local_rows, myid);

  mat.add(2*myid,2*myid,1.0);
  mat.add(2*myid+1,2*myid+1,1.0);
  mat.add(1,0,42.0);

  mat.add((2*myid+2)%(2*numprocs),(2*myid+2)%(2*numprocs),0.1);

  mat.compress(VectorOperation::add);

  std::vector<types::global_dof_index> rows(1,1);
  mat.clear_rows(rows);

//    mat.write_ascii();
  if (myid==0)
    deallog << "2nd try" << std::endl;

  mat = 0;
  mat.add(1,0,42.0);
  mat.add(2*myid,2*myid,1.0);
  mat.add(2*myid+1,2*myid+1,1.0);

  mat.add((2*myid+2)%(2*numprocs),(2*myid+2)%(2*numprocs),0.1);

  mat.compress(VectorOperation::add);
//    mat.write_ascii();

  if (myid==0)
    deallog << "done" << std::endl;

}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

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
