// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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



// make sure that block vector iterator allows reading and writing correctly

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_parallel_block_vector.h>

#include <iostream>

#include "../tests.h"


void
test()
{
  LinearAlgebra::TpetraWrappers::MPI::BlockVector v;
  v.reinit(2);
  for (unsigned int i = 0; i < v.n_blocks(); ++i)
    v.block(i).reinit(complete_index_set(1), MPI_COMM_WORLD);
  v.collect_sizes();

  v(0) = 1;
  v(1) = 2;

  // first check reading through a const
  // iterator
  {
    LinearAlgebra::TpetraWrappers::MPI::BlockVector::const_iterator i =
      v.begin();
    AssertThrow(*i == 1, ExcInternalError());
    ++i;
    AssertThrow(*i == 2, ExcInternalError());
  }

  // same, but create iterator in a different
  // way
  {
    LinearAlgebra::TpetraWrappers::MPI::BlockVector::const_iterator i =
      const_cast<const LinearAlgebra::TpetraWrappers::MPI::BlockVector &>(v)
        .begin();
    AssertThrow(*i == 1, ExcInternalError());
    ++i;
    AssertThrow(*i == 2, ExcInternalError());
  }

  // read through a read-write iterator
  {
    LinearAlgebra::TpetraWrappers::MPI::BlockVector::iterator i = v.begin();
    AssertThrow(*i == 1, ExcInternalError());
    ++i;
    AssertThrow(*i == 2, ExcInternalError());
  }

  // write through a read-write iterator
  {
    LinearAlgebra::TpetraWrappers::MPI::BlockVector::iterator i = v.begin();

    *i = 2;
    ++i;
    *i = 3;
  }

  // and read again
  {
    LinearAlgebra::TpetraWrappers::MPI::BlockVector::iterator i = v.begin();
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
