// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/polynomials_bernstein.h>

#include <boost/math/special_functions/binomial.hpp>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace
{
  template <typename number>
  std::vector<number>
  get_bernstein_coefficients(const unsigned int k, const unsigned int n)
  {
    Assert(n > 0,
           ExcMessage("Bernstein polynomial needs to be of degree > 0."));
    AssertIndexRange(k, n + 1);
    std::vector<number> coeff(n + 1, number(0.0));
    for (unsigned int i = k; i < n + 1; ++i)
      {
        coeff[i] = ((i - k) % 2 == 0 ? 1 : -1) *
                   boost::math::binomial_coefficient<number>(n, i) *
                   boost::math::binomial_coefficient<number>(i, k);
      }
    return coeff;
  }
} // namespace

template <typename number>
PolynomialsBernstein<number>::PolynomialsBernstein(const unsigned int index,
                                                   const unsigned int degree)
  : Polynomials::Polynomial<number>(
      get_bernstein_coefficients<number>(index, degree))
{}


template class dealii::PolynomialsBernstein<double>;


DEAL_II_NAMESPACE_CLOSE
