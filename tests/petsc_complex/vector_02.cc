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



// check assignment of elements in Vector

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  const unsigned int         s = 10;
  PETScWrappers::MPI::Vector v(MPI_COMM_WORLD, s, s);
  for (unsigned int k = 0; k < v.size(); ++k)
    v(k) = k;

  v.compress(VectorOperation::insert);

  PETScWrappers::MPI::Vector v2(MPI_COMM_WORLD, s, s);
  for (int k = 0; (unsigned int)k < v2.size(); ++k)
    v2(k) = PetscScalar(k, -k);

  v2.compress(VectorOperation::insert);

  // we now print the vector v, it is all real. Then print after
  // adding v2 which is complex, and then subtract v2 again to get the
  // original vector back.

  deallog << "before: " << std::endl;
  for (unsigned int k = 0; k < s; ++k)
    deallog << '(' << v(k).real() << ',' << v(k).imag() << "i) ";
  deallog << std::endl;

  v.add(1.0, v2);

  deallog << "after: " << std::endl;
  for (unsigned int k = 0; k < s; ++k)
    deallog << '(' << v(k).real() << ',' << v(k).imag() << "i) ";
  deallog << std::endl;

  v.add(-1.0, v2);

  deallog << "back to original: " << std::endl;
  for (unsigned int k = 0; k < s; ++k)
    deallog << '(' << v(k).real() << ',' << v(k).imag() << "i) ";
  deallog << std::endl;

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
