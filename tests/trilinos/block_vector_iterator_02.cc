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



// like _01, except that we use operator[] instead of operator*

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <iostream>

#include "../tests.h"


void
test()
{
  TrilinosWrappers::MPI::BlockVector v;
  v.reinit(2);
  for (unsigned int i = 0; i < v.n_blocks(); ++i)
    v.block(i).reinit(complete_index_set(1), MPI_COMM_WORLD);
  v.collect_sizes();

  v(0) = 1;
  v(1) = 2;

  // first check reading through a const
  // iterator
  {
    TrilinosWrappers::MPI::BlockVector::const_iterator i = v.begin();
    AssertThrow(i[0] == 1, ExcInternalError());
    AssertThrow(i[1] == 2, ExcInternalError());
  }

  // same, but create iterator in a different
  // way
  {
    TrilinosWrappers::MPI::BlockVector::const_iterator i =
      const_cast<const TrilinosWrappers::MPI::BlockVector &>(v).begin();
    AssertThrow(i[0] == 1, ExcInternalError());
    AssertThrow(i[1] == 2, ExcInternalError());
  }

  // read through a read-write iterator
  {
    TrilinosWrappers::MPI::BlockVector::iterator i = v.begin();
    AssertThrow(i[0] == 1, ExcInternalError());
    AssertThrow(i[1] == 2, ExcInternalError());
  }

  // write through a read-write iterator
  {
    TrilinosWrappers::MPI::BlockVector::iterator i = v.begin();
    i[0]                                           = 2;
    i[1]                                           = 3;
  }

  // and read again
  {
    TrilinosWrappers::MPI::BlockVector::iterator i = v.begin();
    AssertThrow(i[0] == 2, ExcInternalError());
    AssertThrow(i[1] == 3, ExcInternalError());
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
        test();
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
