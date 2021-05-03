// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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

#ifndef dealii_tensor_polynomials_base_h
#define dealii_tensor_polynomials_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * This class provides a framework for the finite element polynomial
 * classes for use with finite element classes that are derived from
 * FE_PolyTensor. An object of this type (or rather of a type derived
 * from this class) is stored as a member variable in each object of
 * type FE_PolyTensor.
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
 *   <li> <tt>PolynomialsABF</tt>
 *   <li> <tt>PolynomialsBDM</tt>
 *   <li> <tt>PolynomialsBernardiRaugel</tt>
 *   <li> <tt>PolynomialsNedelec</tt>
 *   <li> <tt>PolynomialsRaviartThomas</tt>
 *   <li> <tt>PolynomialsRT_Bubbles</tt>
 * </ul>
 *
 * @ingroup Polynomials
 */
template <int dim>
class TensorPolynomialsBase
{
public:
  /**
   * Constructor. This takes the degree of the space, @p deg from the finite element
   * class, and @p n, the number of polynomials for the space.
   */
  TensorPolynomialsBase(const unsigned int deg,
                        const unsigned int n_polynomials);

  /**
   * Move constructor.
   */
  TensorPolynomialsBase(TensorPolynomialsBase<dim> &&) = default; // NOLINT

  /**
   * Copy constructor.
   */
  TensorPolynomialsBase(const TensorPolynomialsBase<dim> &) = default;

  /**
   * Virtual destructor. Makes sure that pointers to this class are deleted
   * properly.
   */
  virtual ~TensorPolynomialsBase() = default;

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
  evaluate(const Point<dim> &           unit_point,
           std::vector<Tensor<1, dim>> &values,
           std::vector<Tensor<2, dim>> &grads,
           std::vector<Tensor<3, dim>> &grad_grads,
           std::vector<Tensor<4, dim>> &third_derivatives,
           std::vector<Tensor<5, dim>> &fourth_derivatives) const = 0;

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
  unsigned int
  degree() const;

  /**
   * A sort of virtual copy constructor, this function returns a copy of
   * the polynomial space object. Derived classes need to override the function
   * here in this base class and return an object of the same type as the
   * derived class.
   *
   * Some places in the library, for example the constructors of FE_PolyTensor,
   * need to make copies of polynomial spaces without knowing their exact type.
   * They do so through this function.
   */
  virtual std::unique_ptr<TensorPolynomialsBase<dim>>
  clone() const = 0;

  /**
   * Return the name of the space.
   */
  virtual std::string
  name() const = 0;

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
TensorPolynomialsBase<dim>::n() const
{
  return n_pols;
}



template <int dim>
inline unsigned int
TensorPolynomialsBase<dim>::degree() const
{
  return polynomial_degree;
}



DEAL_II_NAMESPACE_CLOSE

#endif
