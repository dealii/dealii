// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Tests eigenvalues and eigenvectors of LAPACKFullMatrix with complex numbers

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"

/*
 * Eigenvalues and -vectors of this system are
 * lambda = 3     v = (0.707,  0.707)
 * lambda = 5     v = (0.707, -0.707)
 */
const double mat1[2][2] = {{4., -1.}, {-1., 4.}};

/*
 * Eigenvalues and -vectors of this system are
 * lambda = 1+1i     v = (0.707, -0.707i)
 * lambda = 1-1i     v = (0.707,  0.707i)
 */
const double mat2[2][2] = {{1., -1.}, {1., 1.}};


void
check_matrix(const double *matrix_pointer, const bool make_imaginary)
{
  LAPACKFullMatrix<std::complex<double>> LA(2, 2);
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      {
        LA(i, j).real(matrix_pointer[i * 2 + j]);
        if (make_imaginary)
          LA(i, j).imag(1.0);
      }

  deallog << "Checking matrix " << std::endl;
  LA.print_formatted(deallog.get_file_stream());

  LA.compute_eigenvalues(true, true);
  deallog << "Eigenvalues: ";
  for (unsigned int i = 0; i < LA.m(); ++i)
    {
      std::complex<double> lambda = LA.eigenvalue(i);
      deallog << lambda << " ";
    }
  deallog << std::endl;

  // There are a variety of ways to normalize complex eigenvectors. In this
  // specific example the first entry is always nonzero: normalize twice, first
  // by dividing through by that value (so that it is 1 + 0i) and then again by
  // the magnitude of the vector.
  const auto print_eigenvectors =
    [](const FullMatrix<std::complex<double>> &eigenvectors) {
      for (unsigned int col = 0; col < eigenvectors.n(); ++col)
        {
          const auto first_value = eigenvectors(0, col);
          // we assume this assertion will always pass for all matrices present
          // in this test (obviously in general it does not hold)
          Assert(std::abs(first_value) > 1.0e-10, ExcInternalError());
          // Similarly, this is not Kahan summation, but it is fine for this
          // test with very small matrices
          double norm = 0.0;
          for (unsigned int row = 0; row < eigenvectors.n(); ++row)
            norm += std::norm(eigenvectors(row, col) / first_value);
          norm = std::sqrt(norm);

          for (unsigned int row = 0; row < eigenvectors.n(); ++row)
            deallog << (eigenvectors(row, col) / (first_value * norm)) << "  ";
          deallog << std::endl;
        }
    };
  deallog << "Right eigenvectors" << std::endl;
  print_eigenvectors(LA.get_right_eigenvectors());
  deallog << "Left eigenvectors" << std::endl;
  print_eigenvectors(LA.get_left_eigenvectors());
  deallog << std::endl;
}



int
main()
{
  initlog();

  // Test symmetric system
  check_matrix(&mat1[0][0], false);

  // Test non-symmetric system with complex eigenvalues/eigenvectors
  check_matrix(&mat2[0][0], false);

  // Test with imaginary contribution in matrix
  check_matrix(&mat1[0][0], true);
  check_matrix(&mat2[0][0], true);
}
