// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__polynomials_integrated_legendre_sz_h
#define dealii__polynomials_integrated_legendre_sz_h

#include <deal.II/base/config.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/point.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/thread_management.h>


DEAL_II_NAMESPACE_OPEN

/**
 * Class implementing the integrated Legendre polynomials described in the PhD thesis of Sabine Zaglmayer.
 *
 * This class was written based upon the existing deal.II Legendre class as a base, but with the coefficents adjusted
 * so that the recursive formula is for the integrated Legendre polynomials described in the PhD thesis of
 * Sabine Zaglmayer. The polynomials can be generated recursively from:
 *
 * - $L_{0}(x) = -1$ (added so that it can be generated recursively from 0)
 * - $L_{1}(x) = x$
 * - $L_{2}(x) = \frac{(x^2 - 1)}{2}$
 * - $(n+1)L_{n+1} = (2n-1)L_{n} - (n-2)L_{n-1}$.
 *
 * However, it is also possible to generate them directly from the Legendre polynomials:
 *
 * $L_{n} = \frac{l_{n} - l_{n-2}}{2n-1)}$
 *
 */
class IntegratedLegendreSZ : public Polynomials::Polynomial<double>
{
public:
  /**
   * Constructor generating the coefficient of the polynomials up to degree p.
   */
  IntegratedLegendreSZ (const unsigned int p);


  /**
   * Returns the complete set of Integrated Legendre polynomials up to the given degree.
   */
  static std::vector<Polynomials::Polynomial<double>> generate_complete_basis (const unsigned int degree);


private:
  /**
   * Lock that guarantees that at most one thread is changing and accessing the recursive_coefficients array.
   */
  static Threads::Mutex coefficients_lock;


  /**
   * Vector with already computed coefficients. For each degree of the
   * polynomial, we keep one pointer to the list of coefficients; we do so
   * rather than keeping a vector of vectors in order to simplify
   * programming multithread-safe. In order to avoid memory leak, we use a
   * shared_ptr in order to correctly free the memory of the vectors when
   * the global destructor is called.
   */
  static std::vector<std::shared_ptr<const std::vector<double>>> recursive_coefficients;


  /**
   * Main function to compute the co-efficients of the polyonial at degree p.
   */
  static void compute_coefficients (const unsigned int p);


  /**
   * Get coefficients for constructor.
   */
  static const std::vector<double> &get_coefficients (const unsigned int k);
};

DEAL_II_NAMESPACE_CLOSE

#endif
