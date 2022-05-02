// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// This test should be run on multiple processors.


#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <Epetra_Comm.h>
#include <Epetra_Map.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <typename MatrixType>
void
test(MatrixType &m)
{
  m.set(0, 0, 1.);
  m.compress(VectorOperation::insert);
  m = 0;
  m.compress(VectorOperation::insert);

  Assert(fabs(m.frobenius_norm()) < 1e-15, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());


  try
    {
      {
        const unsigned int n_dofs = 420;
        // check
        // TrilinosWrappers::SparseMatrix
        TrilinosWrappers::SparseMatrix v1(n_dofs, n_dofs, 5U);
        test(v1);

        // check
        // TrilinosWrappers::SparseMatrix
        const unsigned int n_jobs =
          Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
        const unsigned int my_id =
          Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
        Assert(n_dofs % n_jobs == 0, ExcInternalError());
        const unsigned int n_local_dofs = n_dofs / n_jobs;
        IndexSet           local_rows(n_dofs);
        local_rows.add_range(n_local_dofs * my_id, n_local_dofs * (my_id + 1));
        TrilinosWrappers::SparseMatrix v2(local_rows, MPI_COMM_WORLD, 5);
        test(v2);
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
