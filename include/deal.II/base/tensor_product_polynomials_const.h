// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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

#ifndef __deal2__tensor_product_polynomials_const_h
#define __deal2__tensor_product_polynomials_const_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/utilities.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * @addtogroup Polynomials
 * @{
 */

/**
 * Tensor product of given polynomials and a locally constant function. This
 * class inherits most of its functionality from TensorProductPolynomials. It
 * works similarly to that class but adds a constant function for the last
 * index.
 *
 * @author Timo Heister, 2012
 */
template <int dim>
class TensorProductPolynomialsConst : public TensorProductPolynomials<dim>
{
public:
  /**
   * Access to the dimension of
   * this object, for checking and
   * automatic setting of dimension
   * in other classes.
   */
  static const unsigned int dimension = dim;

  /**
   * Constructor. <tt>pols</tt> is a vector of objects that should be derived
   * or otherwise convertible to one-dimensional polynomial objects. It will
   * be copied element by element into a private variable.
   */
  template <class Pol>
  TensorProductPolynomialsConst (const std::vector<Pol> &pols);

  /**
   * Computes the value and the first and second derivatives of each tensor
   * product polynomial at <tt>unit_point</tt>.
   *
   * The size of the vectors must either be equal 0 or equal n(). In the first
   * case, the function will not compute these values.
   *
   * If you need values or derivatives of all tensor product polynomials then
   * use this function, rather than using any of the compute_value(),
   * compute_grad() or compute_grad_grad() functions, see below, in a loop
   * over all tensor product polynomials.
   */
  void compute (const Point<dim>            &unit_point,
                std::vector<double>         &values,
                std::vector<Tensor<1,dim> > &grads,
                std::vector<Tensor<2,dim> > &grad_grads) const;

  /**
   * Computes the value of the <tt>i</tt>th tensor product polynomial at
   * <tt>unit_point</tt>. Here <tt>i</tt> is given in tensor product
   * numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each point value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the compute() function with
   * <tt>values.size()==</tt>n() to get the point values of all tensor
   * polynomials all at once and in a much more efficient way.
   */
  double compute_value (const unsigned int i,
                        const Point<dim> &p) const;

  /**
   * Computes the grad of the <tt>i</tt>th tensor product polynomial at
   * <tt>unit_point</tt>. Here <tt>i</tt> is given in tensor product
   * numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the compute() function, see above, with
   * <tt>grads.size()==</tt>n() to get the point value of all tensor
   * polynomials all at once and in a much more efficient way.
   */
  Tensor<1,dim> compute_grad (const unsigned int i,
                              const Point<dim> &p) const;

  /**
   * Computes the second derivative (grad_grad) of the <tt>i</tt>th tensor
   * product polynomial at <tt>unit_point</tt>. Here <tt>i</tt> is given in
   * tensor product numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the compute() function, see above, with
   * <tt>grad_grads.size()==</tt>n() to get the point value of all tensor
   * polynomials all at once and in a much more efficient way.
   */
  Tensor<2,dim> compute_grad_grad (const unsigned int i,
                                   const Point<dim> &p) const;

  /**
   * Returns the number of tensor product polynomials plus the constant
   * function. For <i>n</i> 1d polynomials this is <i>n<sup>dim</sup>+1</i>.
   */
  unsigned int n () const;
};

/** @} */


/* ---------------- template and inline functions ---------- */

#ifndef DOXYGEN

template <int dim>
template <class Pol>
inline
TensorProductPolynomialsConst<dim>::
TensorProductPolynomialsConst(const std::vector<Pol> &pols)
  :
  TensorProductPolynomials<dim>(pols)
{
  // append index for renumbering
  this->index_map.push_back(this->n_tensor_pols);
  this->index_map_inverse.push_back(this->n_tensor_pols);
}



template <int dim>
inline
unsigned int
TensorProductPolynomialsConst<dim>::n() const
{
  return this->n_tensor_pols+1;
}



template <>
inline
unsigned int
TensorProductPolynomialsConst<0>::n() const
{
  return numbers::invalid_unsigned_int;
}


#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
