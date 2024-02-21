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


// deal.II includes
#include <deal.II/lac/petsc_sparse_matrix.h>

#include <cassert>
#include <complex>
#include <iostream>

#include "../tests.h"

// test read/write access to matrices using explicit cast always.

// sparse matrix elements
void
test_matrix(PETScWrappers::SparseMatrix &m)
{
  deallog << "Check 01 matrix access" << std::endl;

  // fill up a matrix with some numbers
  for (unsigned int k = 0; k < m.m(); ++k)
    for (unsigned int l = 0; l < m.n(); ++l)
      {
        PetscReal el_r = static_cast<double>(k + l);
        PetscReal el_i = static_cast<double>(-1. * (k + l));
        m.set(k, l, std::complex<double>(el_r, el_i));
      }

  m.compress(VectorOperation::insert);

  for (unsigned int k = 0; k < m.m(); ++k)
    for (unsigned int l = 0; l < m.n(); ++l)
      {
        PetscReal el_r = static_cast<double>(k + l);
        PetscReal el_i = static_cast<double>(-1. * (k + l));

        AssertThrow((m(k, l).real() == el_r) && (m(k, l).imag() == el_i),
                    ExcInternalError());
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
        PETScWrappers::SparseMatrix m(5, 5, 5);
        test_matrix(m);

        deallog << "matrix:" << std::endl;
        m.print(deallog.get_file_stream());
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
    }

  deallog.get_file_stream() << std::endl;

  return 0;
}
