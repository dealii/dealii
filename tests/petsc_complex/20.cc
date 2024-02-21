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



// check PETScWrappers::MPI::Vector::operator*(Vector) on two vectors that are
// not orthogonal

#include <deal.II/lac/petsc_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(PETScWrappers::MPI::Vector &v, PETScWrappers::MPI::Vector &w)
{
  // set only certain elements of each
  // vector, and record the expected scalar
  // product
  PetscScalar product = PetscScalar(0);
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      const PetscScalar vi = std::complex<double>(i, (i % 2) * 2.0);
      v(i)                 = vi;
      if (i % 3 == 0)
        {
          const PetscScalar wi = std::complex<double>(5.0 - i, 2.5 * (i % 6));
          w(i)                 = wi;
          product += PetscConj(vi) * wi;
        }
    }

  v.compress(VectorOperation::insert);
  w.compress(VectorOperation::insert);

  // make sure the scalar product is zero
  AssertThrow(v * w == product, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();
  deallog.depth_console(0);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        PETScWrappers::MPI::Vector v(MPI_COMM_WORLD, 100, 100);
        PETScWrappers::MPI::Vector w(MPI_COMM_WORLD, 100, 100);
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
