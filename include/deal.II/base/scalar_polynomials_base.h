// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_scalar_polynomials_base_h
#define dealii_scalar_polynomials_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <memory>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * This class provides a framework for the finite element polynomial
 * classes for use with finite element classes that are derived from
 * FE_Poly. An object of this type (or rather of a type derived from
 * this class) is stored as a member variable in each object of type
 * FE_Poly.
 *
 * <h3>Deriving classes</h3>
 *
 * Any derived class must provide the most basic properties for shape
 * functions evaluated on the reference cell. This includes, but is not
 * limited to, implementing the evaluate(), name(), and
 * clone() member functions. These functions are necessary to store the
 * most basic information of how the polynomials in the derived class evaluate
 * at a given point on the reference cell. More information on each function can
 * be found in the corresponding function's documentation.
 *
 * Some classes that derive from this class include
 * <ul>
 *   <li> <tt>PolynomialsAdini</tt>
 *   <li> <tt>PolynomialsRannacherTurek</tt>
 *   <li> <tt>PolynomialsP</tt>
 *   <li> <tt>PolynomialSpace</tt>
 *   <li> <tt>TensorProductPolynomials</tt>
 *   <li> <tt>TensorProductPolynomialsConst</tt>
 *   <li> <tt>TensorProductPolynomialsBubbles</tt>
 * </ul>
 *
 * @ingroup Polynomials
 */
template <int dim>
class ScalarPolynomialsBase
{
public:
  /**
   * Constructor. This takes the degree of the space, @p deg from the finite element
   * class, and @p n, the number of polynomials for the space.
   */
  ScalarPolynomialsBase(const unsigned int deg,
                        const unsigned int n_polynomials);

  /**
   * Move constructor.
   */
  ScalarPolynomialsBase(ScalarPolynomialsBase<dim> &&) = default; // NOLINT

  /**
   * Copy constructor.
   */
  ScalarPolynomialsBase(const ScalarPolynomialsBase<dim> &) = default;

  /**
   * Virtual destructor. Makes sure that pointers to this class are deleted
   * properly.
   */
  virtual ~ScalarPolynomialsBase() = default;

  /**
   * Compute the value and the derivatives of the polynomials at
   * @p unit_point.
   *
   * The size of the vectors must either be zero or equal <tt>n()</tt>.  In
   * the first case, the function will not compute these values.
   *
   * If you need values or derivatives of all polynomials then use this
   * function, rather than using any of the <tt>compute_value</tt>,
   * <tt>compute_grad</tt> or <tt>compute_grad_grad</tt> functions, see below,
   * in a loop over all tensor product polynomials.
   */
  virtual void
  evaluate(const Point<dim>            &unit_point,
           std::vector<double>         &values,
           std::vector<Tensor<1, dim>> &grads,
           std::vector<Tensor<2, dim>> &grad_grads,
           std::vector<Tensor<3, dim>> &third_derivatives,
           std::vector<Tensor<4, dim>> &fourth_derivatives) const = 0;

  /**
   * Compute the value of the <tt>i</tt>th polynomial at unit point
   * <tt>p</tt>.
   *
   * Consider using evaluate() instead.
   */
  virtual double
  compute_value(const unsigned int i, const Point<dim> &p) const = 0;

  /**
   * Compute the <tt>order</tt>th derivative of the <tt>i</tt>th polynomial
   * at unit point <tt>p</tt>.
   *
   * Consider using evaluate() instead.
   *
   * @tparam order The order of the derivative.
   */
  template <int order>
  Tensor<order, dim>
  compute_derivative(const unsigned int i, const Point<dim> &p) const;

  /**
   * Compute the first derivative of the <tt>i</tt>th polynomial
   * at unit point <tt>p</tt>.
   *
   * Consider using evaluate() instead.
   */
  virtual Tensor<1, dim>
  compute_1st_derivative(const unsigned int i, const Point<dim> &p) const = 0;

  /**
   * Compute the second derivative of the <tt>i</tt>th polynomial
   * at unit point <tt>p</tt>.
   *
   * Consider using evaluate() instead.
   */
  virtual Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i, const Point<dim> &p) const = 0;

  /**
   * Compute the third derivative of the <tt>i</tt>th polynomial
   * at unit point <tt>p</tt>.
   *
   * Consider using evaluate() instead.
   */
  virtual Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i, const Point<dim> &p) const = 0;

  /**
   * Compute the fourth derivative of the <tt>i</tt>th polynomial
   * at unit point <tt>p</tt>.
   *
   * Consider using evaluate() instead.
   */
  virtual Tensor<4, dim>
  compute_4th_derivative(const unsigned int i, const Point<dim> &p) const = 0;

  /**
   * Compute the gradient of the <tt>i</tt>th polynomial at unit point
   * <tt>p</tt>.
   *
   * Consider using evaluate() instead.
   */
  virtual Tensor<1, dim>
  compute_grad(const unsigned int /*i*/, const Point<dim> & /*p*/) const = 0;

  /**
   * Compute the second derivative (grad_grad) of the <tt>i</tt>th polynomial
   * at unit point <tt>p</tt>.
   *
   * Consider using evaluate() instead.
   */
  virtual Tensor<2, dim>
  compute_grad_grad(const unsigned int /*i*/,
                    const Point<dim> & /*p*/) const = 0;

  /**
   * Return the number of polynomials.
   */
  unsigned int
  n() const;

  /**
   * Return the highest polynomial degree of polynomials represented by this
   * class. A derived class may override this if its value is different from
   * @p my_degree.
   */
  virtual unsigned int
  degree() const;

  /**
   * A sort of virtual copy constructor, this function returns a copy of
   * the polynomial space object. Derived classes need to override the function
   * here in this base class and return an object of the same type as the
   * derived class.
   *
   * Some places in the library, for example the constructors of FE_Poly,
   * need to make copies of polynomial spaces without knowing their exact type.
   * They do so through this function.
   */
  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const = 0;

  /**
   * Return the name of the space.
   */
  virtual std::string
  name() const = 0;

  /**
   * Return an estimate (in bytes) for the memory consumption of this object.
   */
  virtual std::size_t
  memory_consumption() const;

private:
  /**
   * The highest polynomial degree of this functions represented by this object.
   */
  const unsigned int polynomial_degree;

  /**
   * The number of polynomials represented by this object.
   */
  const unsigned int n_pols;
};



template <int dim>
inline unsigned int
ScalarPolynomialsBase<dim>::n() const
{
  return n_pols;
}



template <int dim>
inline unsigned int
ScalarPolynomialsBase<dim>::degree() const
{
  return polynomial_degree;
}



template <int dim>
template <int order>
inline Tensor<order, dim>
ScalarPolynomialsBase<dim>::compute_derivative(const unsigned int i,
                                               const Point<dim>  &p) const
{
  if constexpr (order == 1)
    {
      return compute_1st_derivative(i, p);
    }
  else if constexpr (order == 2)
    {
      return compute_2nd_derivative(i, p);
    }
  else if constexpr (order == 3)
    {
      return compute_3rd_derivative(i, p);
    }
  else if constexpr (order == 4)
    {
      return compute_4th_derivative(i, p);
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
}

DEAL_II_NAMESPACE_CLOSE

#endif
