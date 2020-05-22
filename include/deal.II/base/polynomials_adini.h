// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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


#ifndef dealii_polynomials_adini_h
#define dealii_polynomials_adini_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>

DEAL_II_NAMESPACE_OPEN

/**
 * The cubic polynomial space for the Adini element
 *
 * This space consists of the cubic space <i>P<sub>3</sub></i> augmented by
 * the functions <i>xy<sup>3</sup></i> and <i>x<sup>3</sup>y</i>.
 *
 * The basis of the space is chosen to match the node functionals of the Adini
 * element.
 *
 * @todo This polynomial space is implemented in 2D only and does not compute
 * derivatives of order 3 or higher.
 *
 * @ingroup Polynomials
 * @author Bärbel Janssen
 * @date 2007
 */
template <int dim>
class PolynomialsAdini : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * Constructor for the polynomials of the described space
   */
  PolynomialsAdini();

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
  evaluate(const Point<dim> &           unit_point,
           std::vector<double> &        values,
           std::vector<Tensor<1, dim>> &grads,
           std::vector<Tensor<2, dim>> &grad_grads,
           std::vector<Tensor<3, dim>> &third_derivatives,
           std::vector<Tensor<4, dim>> &fourth_derivatives) const override;

  /**
   * Compute the value of the <tt>i</tt>th polynomial at <tt>unit_point</tt>.
   *
   * Consider using evaluate() instead.
   */
  double
  compute_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_1st_derivative()
   */
  virtual Tensor<1, dim>
  compute_1st_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_2nd_derivative()
   */
  virtual Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_3rd_derivative()
   */
  virtual Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_4th_derivative()
   */
  virtual Tensor<4, dim>
  compute_4th_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * Compute the gradient of the <tt>i</tt>th polynomial at
   * <tt>unit_point</tt>.
   *
   * Consider using evaluate() instead.
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Compute the second derivative (grad_grad) of the <tt>i</tt>th polynomial
   * at <tt>unit_point</tt>.
   *
   * Consider using evaluate() instead.
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Return the name of the space, which is <tt>PolynomialsAdini</tt>.
   */
  std::string
  name() const override;

  /**
   * @copydoc ScalarPolynomialsBase::clone()
   */
  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;

private:
  /**
   * Store the coefficients of the polynomials in the order
   * $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   */
  Table<2, double> coef;

  /**
   * Store the coefficients of the x-derivative of the polynomials in the
   * order $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   */
  Table<2, double> dx;

  /**
   * Store the coefficients of the y-derivative of the polynomials in the
   * order $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   */
  Table<2, double> dy;

  /**
   * Store the coefficients of the second x-derivative of the polynomials in
   * the order $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   */
  Table<2, double> dxx;

  /**
   * Store the coefficients of the second y-derivative of the polynomials in
   * the order $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   */
  Table<2, double> dyy;

  /**
   * Store the coefficients of the second mixed derivative of the polynomials
   * in the order $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   */
  Table<2, double> dxy;
};



template <int dim>
inline Tensor<1, dim>
PolynomialsAdini<dim>::compute_1st_derivative(const unsigned int /*i*/,
                                              const Point<dim> & /*p*/) const
{
  Assert(false, ExcNotImplemented());
  return {};
}



template <int dim>
inline Tensor<2, dim>
PolynomialsAdini<dim>::compute_2nd_derivative(const unsigned int /*i*/,
                                              const Point<dim> & /*p*/) const
{
  Assert(false, ExcNotImplemented());
  return {};
}



template <int dim>
inline Tensor<3, dim>
PolynomialsAdini<dim>::compute_3rd_derivative(const unsigned int /*i*/,
                                              const Point<dim> & /*p*/) const
{
  Assert(false, ExcNotImplemented());
  return {};
}



template <int dim>
inline Tensor<4, dim>
PolynomialsAdini<dim>::compute_4th_derivative(const unsigned int /*i*/,
                                              const Point<dim> & /*p*/) const
{
  Assert(false, ExcNotImplemented());
  return {};
}



template <int dim>
inline std::string
PolynomialsAdini<dim>::name() const
{
  return "PolynomialsAdini";
}



DEAL_II_NAMESPACE_CLOSE

#endif
