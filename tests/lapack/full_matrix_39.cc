// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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

  // LAPACK might produce different signs in the eigenvectors depending on the
  // implementation. To avoid these problems, print eigenvectors normalized by
  // the sign of the first non-zero real part of the eigenvectors
  const auto print_eigenvectors =
    [](const FullMatrix<std::complex<double>> &eigenvectors) {
      for (unsigned int col = 0; col < eigenvectors.n(); ++col)
        {
          double sign = 1.;
          for (unsigned int row = 0; row < eigenvectors.n(); ++row)
            if (std::abs(eigenvectors(row, col).real()) > 1e-12)
              {
                sign = (eigenvectors(row, col).real() > 0.) ? 1.0 : -1.0;
                break;
              }
          for (unsigned int row = 0; row < eigenvectors.n(); ++row)
            deallog << sign * eigenvectors(row, col) << "  ";
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
