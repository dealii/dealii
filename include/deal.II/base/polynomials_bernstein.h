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

#ifndef dealii_polynomials_bernstein_h
#define dealii_polynomials_bernstein_h


#include <deal.II/base/config.h>

#include <deal.II/base/polynomial.h>

#include <fstream>
#include <iostream>


DEAL_II_NAMESPACE_OPEN

/**
 * This class implements Bernstein basis polynomials of desire degree as
 * described in
 * http://www.idav.ucdavis.edu/education/CAGDNotes/Bernstein-Polynomials.pdf in
 * the paragraph "Converting from the Bernstein Basis to the Power Basis".
 *
 * They are used to create the Bernstein finite element FE_Bernstein.
 *
 * @ingroup Polynomials
 */
template <typename number>
class PolynomialsBernstein : public Polynomials::Polynomial<number>
{
public:
  /**
   * Construct the @p index -th Bernstein Polynomial of degree @p degree.
   *
   * @f{align*}{
   * B_{\text{index}, \text{degree}} (t)
   *   &= \text{binom}(\text{degree}, \text{index})
   *      \cdot t^{\text{index}}
   *      \cdot (1 - t)^{\text{degree} - \text{index}} \\
   *   &= \sum_{i = \text{index}}^\text{degree}
   *      \cdot (-1)^{i - \text{index}}
   *      \cdot \text{binom}(\text{degree}, i)
   *      \cdot \text{binom}(i, \text{index})
   *      \cdot t^i
   * @f}
   *
   * @param index
   * @param degree
   */
  PolynomialsBernstein(const unsigned int index, const unsigned int degree);
};


template <typename number>
std::vector<Polynomials::Polynomial<number>>
generate_complete_bernstein_basis(const unsigned int degree)
{
  std::vector<Polynomials::Polynomial<number>> v;
  v.reserve(degree + 1);
  for (unsigned int i = 0; i < degree + 1; ++i)
    v.push_back(PolynomialsBernstein<number>(i, degree));
  return v;
}

DEAL_II_NAMESPACE_CLOSE

#endif
