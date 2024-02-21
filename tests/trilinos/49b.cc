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



// like 49, but do the test for
//  TrilinosWrappers::MPI::BlockVector
//         ::operator = (dealii::BlockVector<TrilinosScalar>)
// with block vectors instead of plain vectors

#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(TrilinosWrappers::MPI::BlockVector &v)
{
  std::vector<types::global_dof_index> sizes(2, 3);
  dealii::BlockVector<TrilinosScalar>  w(sizes);

  for (unsigned int i = 0; i < w.size(); ++i)
    w(i) = i;

  v = w;


  // make sure we get the expected result
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      AssertThrow(w(i) == i, ExcInternalError());
      AssertThrow(v(i) == i, ExcInternalError());
    }

  // now also check the reverse assignment
  w = v;
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      AssertThrow(w(i) == i, ExcInternalError());
      AssertThrow(v(i) == i, ExcInternalError());
    }


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
        std::vector<IndexSet> sizes(2, complete_index_set(3));

        TrilinosWrappers::MPI::BlockVector v;
        v.reinit(sizes, MPI_COMM_WORLD);
        test(v);
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
