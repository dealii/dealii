// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Tests basic stuff of Trilinos sparsity patterns

#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  TrilinosWrappers::SparsityPattern sp;
  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;

  deallog << "Creating entries..." << std::endl;

  sp.reinit(5, 7, 3);
  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;

  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 7; ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        sp.add(i, j);

  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;

  sp.compress();

  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;
  deallog << "Number of entries: " << sp.n_nonzero_elements() << std::endl;
  deallog << "Number of rows: " << sp.n_rows() << std::endl;
  deallog << "Number of columns: " << sp.n_cols() << std::endl;
  deallog << "Local size: " << sp.local_size() << std::endl;
  deallog << "Max row length: " << sp.max_entries_per_row() << std::endl;
  deallog << "SP::row_length(0): " << sp.row_length(0) << std::endl;
  deallog << "Bandwidth: " << sp.bandwidth() << std::endl;
  deallog << "SP::empty(): " << sp.empty() << std::endl;

  sp.compress();
  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;

  deallog << "Clearing..." << std::endl;

  sp.clear();

  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;
  deallog << "Bandwidth: " << sp.bandwidth() << std::endl;
  deallog << "SP::empty(): " << sp.empty() << std::endl;
  deallog << "Number of rows: " << sp.n_rows() << std::endl;

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
      test();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
