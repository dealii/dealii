// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_polynomials_integrated_legendre_sz_h
#define dealii_polynomials_integrated_legendre_sz_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>


DEAL_II_NAMESPACE_OPEN

/**
 * Class implementing the integrated Legendre polynomials described in the PhD
 * thesis of Sabine Zaglmayr.
 *
 * This class was written based upon the existing deal.II Legendre class as a
 * base, but with the coefficients adjusted so that the recursive formula is for
 * the integrated Legendre polynomials described in the PhD thesis of Sabine
 * Zaglmayr. The polynomials can be generated recursively from:
 *
 * - $L_{0}(x) = -1$ (added so that it can be generated recursively from 0)
 * - $L_{1}(x) = x$
 * - $L_{2}(x) = \frac{(x^2 - 1)}{2}$
 * - $(n+1)L_{n+1} = (2n-1)L_{n} - (n-2)L_{n-1}$.
 *
 * However, it is also possible to generate them directly from the Legendre
 * polynomials:
 *
 * $L_{n} = \frac{l_{n} - l_{n-2}}{2n-1)}$
 *
 */
class IntegratedLegendreSZ : public Polynomials::Polynomial<double>
{
public:
  /**
   * Constructor generating the coefficients of the polynomials at degree p.
   */
  IntegratedLegendreSZ(const unsigned int p);

  /**
   * Return the complete set of Integrated Legendre polynomials up to the
   * given degree.
   */
  static std::vector<Polynomials::Polynomial<double>>
  generate_complete_basis(const unsigned int degree);

private:
  /**
   * Main function to compute the co-efficients of the polynomial at degree p.
   */
  static const std::vector<double>
  get_coefficients(const unsigned int k);
};

DEAL_II_NAMESPACE_CLOSE

#endif
