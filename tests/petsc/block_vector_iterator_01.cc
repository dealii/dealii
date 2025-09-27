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



// make sure that block vector iterator allows reading and writing correctly

#include <deal.II/lac/petsc_block_vector.h>

#include <iostream>

#include "../tests.h"


void
test()
{
  PETScWrappers::MPI::BlockVector v(2, MPI_COMM_WORLD, 1, 1);
  v(0) = 1;
  v(1) = 2;

  // first check reading through a const
  // iterator
  {
    PETScWrappers::MPI::BlockVector::const_iterator i = v.begin();
    AssertThrow(*i == 1, ExcInternalError());
    ++i;
    AssertThrow(*i == 2, ExcInternalError());
  }

  // same, but create iterator in a different
  // way
  {
    PETScWrappers::MPI::BlockVector::const_iterator i =
      const_cast<const PETScWrappers::MPI::BlockVector &>(v).begin();
    AssertThrow(*i == 1, ExcInternalError());
    ++i;
    AssertThrow(*i == 2, ExcInternalError());
  }

  // read through a read-write iterator
  {
    PETScWrappers::MPI::BlockVector::iterator i = v.begin();
    AssertThrow(*i == 1, ExcInternalError());
    ++i;
    AssertThrow(*i == 2, ExcInternalError());
  }

  // write through a read-write iterator
  {
    PETScWrappers::MPI::BlockVector::iterator i = v.begin();

    *i = 2;
    ++i;
    *i = 3;
  }

  // and read again
  {
    PETScWrappers::MPI::BlockVector::iterator i = v.begin();
    AssertThrow(*i == 2, ExcInternalError());
    ++i;
    AssertThrow(*i == 3, ExcInternalError());
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
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
