// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
  catch (const std::exception &exc)
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
