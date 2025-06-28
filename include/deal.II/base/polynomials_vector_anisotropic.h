// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_polynomials_vector_anisotropic_h
#define dealii_polynomials_vector_anisotropic_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_polynomials_base.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <mutex>
#include <vector>


DEAL_II_NAMESPACE_OPEN

/**
 * This class implements vector-valued anisotropic polynomials and can be
 * used further for the generation of vector-valued polynomials with
 * degrees being different in the `dim - 1` tangential directions and the
 * normal direction, respectively. The known
 * examples include Raviart-Thomas and Nédélec polynomials.
 * The polynomial space can be renumbered as required by specific elements.
 * In fact, Raviart-Thomas and Nédélec elements will choose a different ordering
 * to ensure continuity in normal or tangential direction, respectively
 * @ingroup Polynomials
 */
template <int dim>
class PolynomialsVectorAnisotropic : public TensorPolynomialsBase<dim>
{
public:
  /**
   * Constructor. Creates all basis functions for anisotropic vector-valued
   * polynomials of given degrees in normal and tangential directions,
   * respectively.
   */
  PolynomialsVectorAnisotropic(
    const unsigned int               degree_normal,
    const unsigned int               degree_tangential,
    const std::vector<unsigned int> &polynomial_ordering);

  /**
   * Copy constructor.
   */
  PolynomialsVectorAnisotropic(const PolynomialsVectorAnisotropic &other) =
    default;

  /**
   *  Compute the value and certain derivatives of each anisotropic polynomial
   * at
   * @p unit_point.
   *
   * The size of the vectors must either be zero or equal <tt>n()</tt>. In
   * the first case, the function will not compute these values.
   */
  void
  evaluate(const Point<dim>            &unit_point,
           std::vector<Tensor<1, dim>> &values,
           std::vector<Tensor<2, dim>> &grads,
           std::vector<Tensor<3, dim>> &grad_grads,
           std::vector<Tensor<4, dim>> &third_derivatives,
           std::vector<Tensor<5, dim>> &fourth_derivatives) const override;

  /**
   * Return the name of the space, which is
   * <tt>PolynomialsVectorAnisotropic</tt>.
   */
  std::string
  name() const override;

  /**
   * Return the number of polynomials in the space without requiring to
   * build an object of PolynomialsVectorAnisotropic. This is required by the
   * FiniteElement classes.
   */
  static unsigned int
  n_polynomials(const unsigned int normal_degree,
                const unsigned int tangential_degree);

  /**
   * Return degree of polynomials in tangential direction.
   */
  unsigned int
  get_tangential_degree() const;

  /**
   * Return degree of polynomials in normal direction.
   */
  unsigned int
  get_normal_degree() const;

  /**
   * @copydoc TensorPolynomialsBase::clone()
   */
  virtual std::unique_ptr<TensorPolynomialsBase<dim>>
  clone() const override;

  /**
   * Compute the generalized support points in the ordering used by the
   * polynomial shape functions. Note that these points are not support points
   * in the classical sense as the Lagrange polynomials of the different
   * components have different points, which need to be combined in terms of
   * Piola transformations.
   */
  std::vector<Point<dim>>
  get_polynomial_support_points() const;

private:
  /**
   * The given degree in the normal direction.
   */
  const unsigned int normal_degree;

  /**
   * The given degree in the tangential direction.
   */
  const unsigned int tangential_degree;

  /**
   * An object representing the polynomial space for a single component.
   * We can re-use it for all components by rotating the coordinates of the
   * evaluation point.
   */
  const AnisotropicPolynomials<dim> polynomial_space;

  /**
   * Renumbering from lexicographic order that is used in the underlying
   * scalar anisotropic tensor-product polynomial space to the numbering
   * used by the finite element.
   */
  std::vector<unsigned int> lexicographic_to_hierarchic;

  /**
   * Renumbering from hierarchic to lexicographic order. Inverse of
   * lexicographic_to_hierarchic.
   */
  std::vector<unsigned int> hierarchic_to_lexicographic;

  /**
   * Renumbering from shifted polynomial spaces to lexicographic one.
   */
  std::array<std::vector<unsigned int>, dim> renumber_aniso;
};


DEAL_II_NAMESPACE_CLOSE

#endif
