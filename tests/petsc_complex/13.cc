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



// check PETScWrappers::MPI::Vector::operator() in add-mode

#include <deal.II/lac/petsc_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(PETScWrappers::MPI::Vector &v)
{
  // set only certain elements of the vector. have a bit pattern of
  // where we actually wrote elements to
  std::vector<bool> pattern(v.size(), false);
  for (unsigned int i = 0; i < v.size(); i += 1 + i)
    {
      v(i) += PetscScalar(i, i);
      pattern[i] = true;
    }

  v.compress(VectorOperation::add);

  // check that they are ok, and this time all of them
  for (unsigned int i = 0; i < v.size(); ++i)
    AssertThrow(
      ((pattern[i] == true) && (v(i).real() == i) && (v(i).imag() == i)) ||
        //&& (v(i)==std::complex<double> (1.,1.)) ) ||
        ((pattern[i] == false) && (v(i).real() == 0.) && (v(i).imag() == 0.)),
      ExcInternalError());

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
