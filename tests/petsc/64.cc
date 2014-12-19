// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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



// This test should be run on multiple processors. note that this test also
// started to fail with the upgrade to petsc 2.2.1 which required a fix in
// PETScWrappers::MatrixBase::operator=


#include "../tests.h"
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>
#include <vector>


template<typename MatrixType>
void test (MatrixType &m)
{
  m.add(0,0,1);
  m = 0;
  m.compress();

  Assert(fabs(m.frobenius_norm())<1e-15, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        const unsigned int n_dofs=420;
        // check
        // PETScWrappers::SparseMatrix
        PETScWrappers::SparseMatrix
        v1 (n_dofs, n_dofs, 5);
        test (v1);

        // check
        // PETScWrappers::MPI::SparseMatrix
        MPI_Comm mpi_communicator (MPI_COMM_WORLD);
        int n_jobs=1;
        MPI_Comm_size (mpi_communicator, &n_jobs);
        const unsigned int n_mpi_processes=static_cast<unsigned int>(n_jobs);
        Assert(n_dofs%n_mpi_processes==0, ExcInternalError());
        const unsigned int n_local_dofs=n_dofs/n_mpi_processes;
        PETScWrappers::MPI::SparseMatrix
        v2 (mpi_communicator, n_dofs, n_dofs, n_local_dofs, n_local_dofs, 5);
        test (v2);
      }

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
