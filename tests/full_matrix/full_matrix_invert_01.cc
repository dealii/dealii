// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test left and right inversion of FullMatrix

#include <limits>

#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
fill_matrix_invertible(FullMatrix<number> &A)
{
  for (unsigned int i = 0; i < A.m(); ++i)
    for (unsigned int j = 0; j < A.n(); ++j)
      {
        A(i, j) = number(i * j);
        if (i == j)
          A(i, i) += i + A.m();
      }
}


template <typename number>
bool
calculate(const FullMatrix<number> A, bool disp = true)
{
  bool retval = true;
  // Different tolerance for different number types
  const number tol = 1000 * std::numeric_limits<number>::epsilon();

  // Test left invert
  if (A.m() >= A.n())
    {
      FullMatrix<number> A_l_inv(A.n(), A.m());
      FullMatrix<number> identity(A.n(), A.n());
      A_l_inv.left_invert(A);
      A_l_inv.mmult(identity, A);

      FullMatrix<double> M(IdentityMatrix(identity.n()));
      M.add(-1, identity);
      if (disp || M.linfty_norm() > tol)
        {
          // deallog << "A matrix" << std::endl;
          // display_matrix(A);
          deallog << "Left inverse" << std::endl;
          display_matrix(A_l_inv);
          // deallog << "Identity = A_l_inv*A" << std::endl;
          // display_matrix(identity);
          retval = false;
        }
    }

  // Test right invert
  if (A.m() <= A.n())
    {
      FullMatrix<number> A_r_inv(A.n(), A.m());
      FullMatrix<number> identity(A.m(), A.m());
      A_r_inv.right_invert(A);
      A.mmult(identity, A_r_inv);

      FullMatrix<double> M(IdentityMatrix(identity.n()));
      M.add(-1, identity);
      if (disp || M.linfty_norm() > tol)
        {
          // deallog << "A matrix"<< std::endl;
          // display_matrix(A);
          deallog << "Right inverse" << std::endl;
          display_matrix(A_r_inv);
          // deallog << "Identity = A*A_r_inv" << std::endl;
          // display_matrix(identity);
          // deallog << std::endl;
          retval = false;
        }
    }

  return retval;
}


template <typename number>
void
check()
{
  unsigned int nFails = 0;
  unsigned int maxDim = 10;
  for (unsigned int i = 1; i <= maxDim; ++i)
    {
      for (unsigned int j = 1; j <= maxDim; ++j)
        {
          FullMatrix<number> A(i, j);
          fill_matrix_invertible(A);
          if (!calculate(A, false))
            {
              nFails++;
            }
        }
    }

  if (nFails > 0)
    deallog << nFails << " out of " << maxDim * maxDim
            << " calls with wrong result" << std::endl;
  else
    deallog << "OK" << std::endl;
}
