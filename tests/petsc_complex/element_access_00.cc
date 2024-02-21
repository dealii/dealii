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
#include <deal.II/lac/petsc_vector.h>

#include <cassert>
#include <complex>
#include <iostream>

#include "../tests.h"

// test dealii::internal::VectorReference::real()
// test dealii::internal::VectorReference::imag()
// on vector

// vector elements
void
test_vector(PETScWrappers::MPI::Vector &v)
{
  deallog << "Check vector access" << std::endl;

  // fill up a vector with some numbers
  for (unsigned int k = 0; k < v.size(); ++k)
    v(k) = std::complex<double>(k, v.size() - k);

  v.compress(VectorOperation::insert);

  // check that is what we get by casting PetscScalar to std::real()
  // and std::imag()
  for (unsigned int k = 0; k < v.size(); ++k)
    AssertThrow((static_cast<std::complex<double>>(v(k)).real() == k) &&
                  (static_cast<std::complex<double>>(v(k)).imag() ==
                   v.size() - k),
                ExcInternalError());

  // check that is what we get by
  // dealii::internal::VectorReference::real() and
  // dealii::internal::VectorReference::imag()
  for (unsigned int k = 0; k < v.size(); ++k)
    AssertThrow((v(k).real() == k) && (v(k).imag() == v.size() - k),
                ExcInternalError());

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
        PETScWrappers::MPI::Vector v(MPI_COMM_WORLD, 5, 5);
        test_vector(v);

        deallog << "vector:" << std::endl;
        v.print(deallog.get_file_stream());
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
