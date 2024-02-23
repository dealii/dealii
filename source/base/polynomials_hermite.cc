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

#include <deal.II/base/polynomials_hermite.h>
#include <deal.II/base/utilities.h>

DEAL_II_NAMESPACE_OPEN

namespace Polynomials
{
  namespace
  {
    std::vector<double>
    hermite_poly_coeffs(const unsigned int regularity, const unsigned int index)
    {
      AssertIndexRange(index, 2 * regularity + 2);

      const unsigned int curr_index = index % (regularity + 1);
      const unsigned int side       = (index > regularity) ? 1 : 0;

      // Signed ints are used here to protect against underflow errors
      const int loop_control_1 = static_cast<int>(regularity + 1 - curr_index);
      const int loop_control_2 = (side == 1) ?
                                   static_cast<int>(curr_index + 1) :
                                   static_cast<int>(regularity + 2);

      std::vector<double> poly_coeffs(2 * regularity + 2, 0.0);

      if (side == 1) // right side: g polynomials
        {
          int binomial_1 = (curr_index % 2) ? -1 : 1;

          for (int i = 0; i < loop_control_2; ++i)
            {
              int inv_binomial = 1;

              for (int j = 0; j < loop_control_1; ++j)
                {
                  int binomial_2 = 1;

                  for (int k = 0; k < j + 1; ++k)
                    {
                      poly_coeffs[regularity + i + k + 1] +=
                        binomial_1 * inv_binomial * binomial_2;
                      binomial_2 *= k - j;
                      binomial_2 /= k + 1;
                    }
                  inv_binomial *= regularity + j + 1;
                  inv_binomial /= j + 1;
                }
              // ints used here to protect against underflow errors
              binomial_1 *= -static_cast<int>(curr_index - i);
              binomial_1 /= i + 1;
            }
        }
      else // left side: f polynomials
        {
          int binomial = 1;

          for (int i = 0; i < loop_control_2; ++i)
            {
              int inv_binomial = 1;

              for (int j = 0; j < loop_control_1; ++j)
                {
                  poly_coeffs[curr_index + i + j] += binomial * inv_binomial;
                  inv_binomial *= regularity + j + 1;
                  inv_binomial /= j + 1;
                }
              // Protection needed here against underflow errors
              binomial *= -static_cast<int>(regularity + 1 - i);
              binomial /= i + 1;
            }
        }

      // rescale coefficients by a factor of 4^curr_index to account for reduced
      // L2-norms
      double precond_factor = Utilities::pow(4, curr_index);
      for (auto &it : poly_coeffs)
        it *= precond_factor;

      return poly_coeffs;
    }
  } // namespace



  PolynomialsHermite::PolynomialsHermite(const unsigned int regularity,
                                         const unsigned int index)
    : Polynomial<double>(hermite_poly_coeffs(regularity, index))
    , degree(2 * regularity + 1)
    , regularity(regularity)
    , side_index(index % (regularity + 1))
    , side((index >= regularity + 1) ? 1 : 0)
  {
    AssertIndexRange(index, 2 * (regularity + 1));
  }



  std::vector<Polynomial<double>>
  PolynomialsHermite::generate_complete_basis(const unsigned int regularity)
  {
    std::vector<Polynomial<double>> polys;
    const unsigned int              sz = 2 * regularity + 2;
    polys.reserve(sz);

    for (unsigned int i = 0; i < sz; ++i)
      polys.emplace_back(PolynomialsHermite(regularity, i));

    return polys;
  }
} // namespace Polynomials

DEAL_II_NAMESPACE_CLOSE
