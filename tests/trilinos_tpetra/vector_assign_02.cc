// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this is equivalent to the trilinos_vector_assign_02 test, except that we use
// operator+= instead of operator=. Now, this does not present a problem,
// since the compiler does not automatically generate a version of this
// operator, but simply performs the conversion to TrilinosScalar, i.e. the
// argument to the user-defined operator+=. This is not exciting, but since I
// wrote the test to make sure it works this way, let's keep it then...

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> &v,
     LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> &w)
{
  // set the first vector
  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = i;

  // add elements by reference
  for (unsigned int i = 0; i < v.size(); ++i)
    w(i) += v(i);

  // check that they're equal
  AssertThrow(v == w, ExcInternalError());

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
        LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> w;
        w.reinit(complete_index_set(100), MPI_COMM_WORLD);
        test(v, w);
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
