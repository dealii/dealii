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



// this is equivalent to the petsc_parallel_vector_assign_01 test, except that
// we use operator+= instead of operator=. This is also not exciting...

#include <deal.II/lac/petsc_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(PETScWrappers::MPI::Vector &v, PETScWrappers::MPI::Vector &w)
{
  // set the first vector
  for (unsigned int k = 0; k < v.size(); ++k)
    v(k) = PetscScalar(k, k * k);

  // add elements by reference
  for (unsigned int k = 0; k < v.size(); ++k)
    w(k) += v(k);

  // check that they're equal
  AssertThrow(v == w, ExcInternalError());

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
        PETScWrappers::MPI::Vector v(MPI_COMM_WORLD, 20, 20);
        PETScWrappers::MPI::Vector w(MPI_COMM_WORLD, 20, 20);
        test(v, w);

        // Some output
        deallog << "Complex vectors: " << std::endl;
        deallog << "v: " << std::flush;
        v.print(deallog.get_file_stream(), 0, false, true);
        deallog << "w: " << std::flush;
        w.print(deallog.get_file_stream(), 0, false, true);
        deallog << "OK" << std::endl;
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
