// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2017 by the deal.II authors
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

#ifndef dealii_polynomial_space_h
#define dealii_polynomial_space_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Representation of the space of polynomials of degree at most n in higher
 * dimensions.
 *
 * Given a vector of <i>n</i> one-dimensional polynomials <i>P<sub>0</sub></i>
 * to <i>P<sub>n</sub></i>, where <i>P<sub>i</sub></i> has degree <i>i</i>,
 * this class generates all dim-dimensional polynomials of the form <i>
 * P<sub>ijk</sub>(x,y,z) =
 * P<sub>i</sub>(x)P<sub>j</sub>(y)P<sub>k</sub>(z)</i>, where the sum of
 * <i>i</i>, <i>j</i> and <i>k</i> is less than or equal <i>n</i>.
 *
 * The output_indices() function prints the ordering of the polynomials, i.e.
 * for each dim-dimensional polynomial in the polynomial space it gives the
 * indices i,j,k of the one-dimensional polynomials in x,y and z direction.
 * The ordering of the dim-dimensional polynomials can be changed by using the
 * set_numbering() function.
 *
 * The standard ordering of polynomials is that indices for the first space
 * dimension vary fastest and the last space dimension is slowest. In
 * particular, if we take for simplicity the vector of monomials
 * <i>x<sup>0</sup>, x<sup>1</sup>, x<sup>2</sup>,..., x<sup>n</sup></i>, we
 * get
 *
 * <dl> <dt> 1D <dd> <i> x<sup>0</sup>, x<sup>1</sup>,...,x<sup>n</sup></i>
 * <dt> 2D: <dd> <i> x<sup>0</sup>y<sup>0</sup>,
 * x<sup>1</sup>y<sup>0</sup>,..., x<sup>n</sup>y<sup>0</sup>,
 * <br>
 * x<sup>0</sup>y<sup>1</sup>, x<sup>1</sup>y<sup>1</sup>,...,
 * x<sup>n-1</sup>y<sup>1</sup>,
 * <br>
 * x<sup>0</sup>y<sup>2</sup>,... x<sup>n-2</sup>y<sup>2</sup>,
 * <br>
 * ...
 * <br>
 * x<sup>0</sup>y<sup>n-1</sup>, x<sup>1</sup>y<sup>n-1</sup>,
 * <br>
 * x<sup>0</sup>y<sup>n</sup> </i> <dt> 3D: <dd> <i>
 * x<sup>0</sup>y<sup>0</sup>z<sup>0</sup>,...,
 * x<sup>n</sup>y<sup>0</sup>z<sup>0</sup>,
 * <br>
 * x<sup>0</sup>y<sup>1</sup>z<sup>0</sup>,...,
 * x<sup>n-1</sup>y<sup>1</sup>z<sup>0</sup>,
 * <br>
 * ...
 * <br>
 * x<sup>0</sup>y<sup>n</sup>z<sup>0</sup>,
 * <br>
 * x<sup>0</sup>y<sup>0</sup>z<sup>1</sup>,...
 * x<sup>n-1</sup>y<sup>0</sup>z<sup>1</sup>,
 * <br>
 * ...
 * <br>
 * x<sup>0</sup>y<sup>n-1</sup>z<sup>1</sup>,
 * <br>
 * x<sup>0</sup>y<sup>0</sup>z<sup>2</sup>,...
 * x<sup>n-2</sup>y<sup>0</sup>z<sup>2</sup>,
 * <br>
 * ...
 * <br>
 * x<sup>0</sup>y<sup>0</sup>z<sup>n</sup> </i> </dl>
 *
 * @ingroup Polynomials
 * @author Guido Kanschat, Wolfgang Bangerth, Ralf Hartmann 2002, 2003, 2004,
 * 2005
 */
template <int dim>
class PolynomialSpace
{
public:
  /**
   * Access to the dimension of this object, for checking and automatic
   * setting of dimension in other classes.
   */
  static const unsigned int dimension = dim;

  /**
   * Constructor. <tt>pols</tt> is a vector of pointers to one-dimensional
   * polynomials and will be copied into a private member variable. The static
   * type of the template argument <tt>pols</tt> needs to be convertible to
   * Polynomials::Polynomial@<double@>, i.e. should usually be a derived class
   * of Polynomials::Polynomial@<double@>.
   */
  template <class Pol>
  PolynomialSpace(const std::vector<Pol> &pols);

  /**
   * Prints the list of the indices to <tt>out</tt>.
   */
  template <class StreamType>
  void
  output_indices(StreamType &out) const;

  /**
   * Set the ordering of the polynomials. Requires
   * <tt>renumber.size()==n()</tt>. Stores a copy of <tt>renumber</tt>.
   */
  void
  set_numbering(const std::vector<unsigned int> &renumber);

  /**
   * Compute the value and the first and second derivatives of each
   * polynomial at <tt>unit_point</tt>.
   *
   * The size of the vectors must either be equal 0 or equal n(). In the first
   * case, the function will not compute these values, i.e. you indicate what
   * you want to have computed by resizing those vectors which you want
   * filled.
   *
   * If you need values or derivatives of all polynomials then use this
   * function, rather than using any of the compute_value(), compute_grad() or
   * compute_grad_grad() functions, see below, in a loop over all polynomials.
   */
  void
  compute(const Point<dim> &           unit_point,
          std::vector<double> &        values,
          std::vector<Tensor<1, dim>> &grads,
          std::vector<Tensor<2, dim>> &grad_grads,
          std::vector<Tensor<3, dim>> &third_derivatives,
          std::vector<Tensor<4, dim>> &fourth_derivatives) const;

  /**
   * Compute the value of the <tt>i</tt>th polynomial at unit point
   * <tt>p</tt>.
   *
   * Consider using compute() instead.
   */
  double
  compute_value(const unsigned int i, const Point<dim> &p) const;

  /**
   * Compute the <tt>order</tt>th derivative of the <tt>i</tt>th polynomial
   * at unit point <tt>p</tt>.
   *
   * Consider using compute() instead.
   *
   * @tparam order The order of the derivative.
   */
  template <int order>
  Tensor<order, dim>
  compute_derivative(const unsigned int i, const Point<dim> &p) const;

  /**
   * Compute the gradient of the <tt>i</tt>th polynomial at unit point
   * <tt>p</tt>.
   *
   * Consider using compute() instead.
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const;

  /**
   * Compute the second derivative (grad_grad) of the <tt>i</tt>th polynomial
   * at unit point <tt>p</tt>.
   *
   * Consider using compute() instead.
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const;

  /**
   * Return the number of polynomials spanning the space represented by this
   * class. Here, if <tt>N</tt> is the number of one-dimensional polynomials
   * given, then the result of this function is <i>N</i> in 1d,
   * <i>N(N+1)/2</i> in 2d, and <i>N(N+1)(N+2)/6</i> in 3d.
   */
  unsigned int
  n() const;

  /**
   * Degree of the space. This is by definition the number of polynomials
   * given to the constructor, NOT the maximal degree of a polynomial in this
   * vector. The latter value is never checked and therefore left to the
   * application.
   */
  unsigned int
  degree() const;

  /**
   * Static function used in the constructor to compute the number of
   * polynomials.
   *
   * @warning The argument `n` is not the maximal degree, but the number of
   * onedimensional polynomials, thus the degree plus one.
   */
  static unsigned int
  compute_n_pols(const unsigned int n);

protected:
  /**
   * Compute numbers in x, y and z direction. Given an index <tt>n</tt> in the
   * d-dimensional polynomial space, return the indices i,j,k such that
   * <i>p<sub>n</sub>(x,y,z) =
   * p<sub>i</sub>(x)p<sub>j</sub>(y)p<sub>k</sub>(z)</i>.
   *
   * In 1d and 2d, obviously only i and i,j are returned.
   */
  std::array<unsigned int, dim>
  compute_index(const unsigned int n) const;

private:
  /**
   * Copy of the vector <tt>pols</tt> of polynomials given to the constructor.
   */
  const std::vector<Polynomials::Polynomial<double>> polynomials;

  /**
   * Store the precomputed value which the <tt>n()</tt> function returns.
   */
  const unsigned int n_pols;

  /**
   * Index map for reordering the polynomials.
   */
  std::vector<unsigned int> index_map;

  /**
   * Index map for reordering the polynomials.
   */
  std::vector<unsigned int> index_map_inverse;
};


/* -------------- declaration of explicit specializations --- */

template <>
std::array<unsigned int, 1>
PolynomialSpace<1>::compute_index(const unsigned int n) const;
template <>
std::array<unsigned int, 2>
PolynomialSpace<2>::compute_index(const unsigned int n) const;
template <>
std::array<unsigned int, 3>
PolynomialSpace<3>::compute_index(const unsigned int n) const;



/* -------------- inline and template functions ------------- */

template <int dim>
template <class Pol>
PolynomialSpace<dim>::PolynomialSpace(const std::vector<Pol> &pols)
  : polynomials(pols.begin(), pols.end())
  , n_pols(compute_n_pols(polynomials.size()))
  , index_map(n_pols)
  , index_map_inverse(n_pols)
{
  // per default set this index map
  // to identity. This map can be
  // changed by the user through the
  // set_numbering function
  for (unsigned int i = 0; i < n_pols; ++i)
    {
      index_map[i]         = i;
      index_map_inverse[i] = i;
    }
}


template <int dim>
inline unsigned int
PolynomialSpace<dim>::n() const
{
  return n_pols;
}



template <int dim>
inline unsigned int
PolynomialSpace<dim>::degree() const
{
  return polynomials.size();
}


template <int dim>
template <class StreamType>
void
PolynomialSpace<dim>::output_indices(StreamType &out) const
{
  for (unsigned int i = 0; i < n_pols; ++i)
    {
      const std::array<unsigned int, dim> ix = compute_index(i);
      out << i << "\t";
      for (unsigned int d = 0; d < dim; ++d)
        out << ix[d] << " ";
      out << std::endl;
    }
}

template <int dim>
template <int order>
Tensor<order, dim>
PolynomialSpace<dim>::compute_derivative(const unsigned int i,
                                         const Point<dim> & p) const
{
  const std::array<unsigned int, dim> indices = compute_index(i);

  double v[dim][order + 1];
  {
    std::vector<double> tmp(order + 1);
    for (unsigned int d = 0; d < dim; ++d)
      {
        polynomials[indices[d]].value(p(d), tmp);
        for (unsigned int j = 0; j < order + 1; ++j)
          v[d][j] = tmp[j];
      }
  }

  Tensor<order, dim> derivative;
  switch (order)
    {
      case 1:
        {
          Tensor<1, dim> &derivative_1 =
            *reinterpret_cast<Tensor<1, dim> *>(&derivative);
          for (unsigned int d = 0; d < dim; ++d)
            {
              derivative_1[d] = 1.;
              for (unsigned int x = 0; x < dim; ++x)
                {
                  unsigned int x_order = 0;
                  if (d == x)
                    ++x_order;

                  derivative_1[d] *= v[x][x_order];
                }
            }

          return derivative;
        }
      case 2:
        {
          Tensor<2, dim> &derivative_2 =
            *reinterpret_cast<Tensor<2, dim> *>(&derivative);
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              {
                derivative_2[d1][d2] = 1.;
                for (unsigned int x = 0; x < dim; ++x)
                  {
                    unsigned int x_order = 0;
                    if (d1 == x)
                      ++x_order;
                    if (d2 == x)
                      ++x_order;

                    derivative_2[d1][d2] *= v[x][x_order];
                  }
              }

          return derivative;
        }
      case 3:
        {
          Tensor<3, dim> &derivative_3 =
            *reinterpret_cast<Tensor<3, dim> *>(&derivative);
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              for (unsigned int d3 = 0; d3 < dim; ++d3)
                {
                  derivative_3[d1][d2][d3] = 1.;
                  for (unsigned int x = 0; x < dim; ++x)
                    {
                      unsigned int x_order = 0;
                      if (d1 == x)
                        ++x_order;
                      if (d2 == x)
                        ++x_order;
                      if (d3 == x)
                        ++x_order;

                      derivative_3[d1][d2][d3] *= v[x][x_order];
                    }
                }

          return derivative;
        }
      case 4:
        {
          Tensor<4, dim> &derivative_4 =
            *reinterpret_cast<Tensor<4, dim> *>(&derivative);
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              for (unsigned int d3 = 0; d3 < dim; ++d3)
                for (unsigned int d4 = 0; d4 < dim; ++d4)
                  {
                    derivative_4[d1][d2][d3][d4] = 1.;
                    for (unsigned int x = 0; x < dim; ++x)
                      {
                        unsigned int x_order = 0;
                        if (d1 == x)
                          ++x_order;
                        if (d2 == x)
                          ++x_order;
                        if (d3 == x)
                          ++x_order;
                        if (d4 == x)
                          ++x_order;

                        derivative_4[d1][d2][d3][d4] *= v[x][x_order];
                      }
                  }

          return derivative;
        }
      default:
        {
          Assert(false, ExcNotImplemented());
          return derivative;
        }
    }
}


DEAL_II_NAMESPACE_CLOSE

#endif
