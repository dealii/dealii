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


// Tests eigenvectors of LAPACKFullMatrix

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

/*
 * Eigenvalues and -vectors of this system are
 * lambda = 1     v = (1, 0, 0)
 * lambda = 1     v = (0, 1, 0)
 * lambda = 1     v = (-1, 0, 0)
 */
const double mat3[3][3] = {{1., 0, 2.}, {0., 1., 0.}, {0., 0., 1.}};


void
check_matrix(const double *matrix_pointer, const unsigned int size)
{
  FullMatrix<double>       A(size, size, matrix_pointer);
  LAPACKFullMatrix<double> LA(size, size);
  LA = A;

  deallog << "Checking matrix " << std::endl;
  LA.print_formatted(deallog.get_file_stream());

  LA.compute_eigenvalues(true, true);
  deallog << "Eigenvalues: ";
  for (unsigned int i = 0; i < A.m(); ++i)
    {
      std::complex<double> lambda = LA.eigenvalue(i);
      deallog << lambda << " ";
    }
  deallog << std::endl;
  deallog << "Right eigenvectors" << std::endl;
  LA.get_right_eigenvectors().print(deallog.get_file_stream());
  deallog << "Left eigenvectors" << std::endl;
  LA.get_left_eigenvectors().print(deallog.get_file_stream());
  deallog << std::endl;
}



int
main()
{
  initlog();

  // Test symmetric system
  check_matrix(&mat1[0][0], 2);

  // Test non-symmetric system with complex eigenvalues/eigenvectors
  check_matrix(&mat2[0][0], 2);

  // Test case of one eigenvalue with two associated eigenvectors
  check_matrix(&mat3[0][0], 3);
}
