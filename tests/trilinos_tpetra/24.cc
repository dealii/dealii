// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this test used to check for
// LinearAlgebra::TpetraWrappers::Vector<double>::clear(). However, this
// function has since been removed, so we test for v=0 instead, although that
// may be covered by one of the other tests

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> &v)
{
  // set some entries of the vector
  for (unsigned int i = 0; i < v.size(); ++i)
    if (i % 3 == 0)
      v(i) = i + 1.;
  v.compress(VectorOperation::insert);

  // then clear it again and make sure the
  // vector is really empty
  const unsigned int sz = v.size();
  v                     = 0;
  AssertThrow(v.size() == sz, ExcInternalError());
  AssertThrow(v.l2_norm() == 0, ExcInternalError());

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
        LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> v;
        v.reinit(complete_index_set(100), MPI_COMM_WORLD);
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
