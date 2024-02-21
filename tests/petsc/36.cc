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



// check PETScWrappers::MPI::Vector::operator-=(Vector)

#include <deal.II/lac/petsc_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(PETScWrappers::MPI::Vector &v, PETScWrappers::MPI::Vector &w)
{
  // set only certain elements of each
  // vector
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      v(i) = i;
      if (i % 3 == 0)
        w(i) = i + 1.;
    }

  v.compress(VectorOperation::insert);
  w.compress(VectorOperation::insert);

  v -= w;

  // make sure we get the expected result
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      if (i % 3 == 0)
        {
          AssertThrow(w(i) == i + 1., ExcInternalError());
          AssertThrow(v(i) == i - (i + 1.), ExcInternalError());
        }
      else
        {
          AssertThrow(w(i) == 0, ExcInternalError());
          AssertThrow(v(i) == i, ExcInternalError());
        }
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
        IndexSet indices(100);
        indices.add_range(0, 100);
        PETScWrappers::MPI::Vector v(indices, MPI_COMM_WORLD);
        PETScWrappers::MPI::Vector w(indices, MPI_COMM_WORLD);
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
