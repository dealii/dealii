// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


#ifndef dealii_simplex_polynomials_h
#define dealii_simplex_polynomials_h

#include <deal.II/base/config.h>

#include <deal.II/base/scalar_polynomials_base.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace for functions and classes that provide support for simplex
 * reference cell entities, i.e., triangles and tetrahedrons.
 *
 *  @ingroup simplex
 */
namespace Simplex
{
  /**
   * Polynomials defined on dim-dimensional simplex entities. This class is
   * basis of Simplex::FE_P.
   *
   *  @ingroup simplex
   */
  template <int dim>
  class ScalarPolynomial : public ScalarPolynomialsBase<dim>
  {
    static_assert(dim == 2 || dim == 3, "Dimension not supported!");

  public:
    /**
     * Make the dimension available to the outside.
     */
    static const unsigned int dimension = dim;

    /*
     * Constructor taking the polynomial @p degree as input.
     *
     * @note Currently, only linear (degree=1) and quadratic polynomials
     *   (degree=2) are implemented.
     */
    ScalarPolynomial(const unsigned int degree);

    /**
     * @copydoc ScalarPolynomialsBase::evaluate()
     *
     * @note Currently, only the vectors @p values and @grads are filled.
     */
    void
    evaluate(const Point<dim> &           unit_point,
             std::vector<double> &        values,
             std::vector<Tensor<1, dim>> &grads,
             std::vector<Tensor<2, dim>> &grad_grads,
             std::vector<Tensor<3, dim>> &third_derivatives,
             std::vector<Tensor<4, dim>> &fourth_derivatives) const override;

    /**
     * @copydoc ScalarPolynomialsBase::compute_value()
     */
    double
    compute_value(const unsigned int i, const Point<dim> &p) const override;

    /**
     * @copydoc ScalarPolynomialsBase::compute_derivative()
     *
     * @note Currently, only implemented for first derivative.
     */
    template <int order>
    Tensor<order, dim>
    compute_derivative(const unsigned int i, const Point<dim> &p) const;

    /**
     * @copydoc ScalarPolynomialsBase::compute_1st_derivative()
     */
    Tensor<1, dim>
    compute_1st_derivative(const unsigned int i,
                           const Point<dim> & p) const override;

    /**
     * @copydoc ScalarPolynomialsBase::compute_2nd_derivative()
     *
     * @note Not implemented yet.
     */
    Tensor<2, dim>
    compute_2nd_derivative(const unsigned int i,
                           const Point<dim> & p) const override;

    /**
     * @copydoc ScalarPolynomialsBase::compute_3rd_derivative()
     *
     * @note Not implemented yet.
     */
    Tensor<3, dim>
    compute_3rd_derivative(const unsigned int i,
                           const Point<dim> & p) const override;

    /**
     * @copydoc ScalarPolynomialsBase::compute_4th_derivative()
     *
     * @note Not implemented yet.
     */
    Tensor<4, dim>
    compute_4th_derivative(const unsigned int i,
                           const Point<dim> & p) const override;

    /**
     * @copydoc ScalarPolynomialsBase::compute_grad()
     *
     * @note Not implemented yet.
     */
    Tensor<1, dim>
    compute_grad(const unsigned int i, const Point<dim> &p) const override;

    /**
     * @copydoc ScalarPolynomialsBase::compute_grad_grad()
     *
     * @note Not implemented yet.
     */
    Tensor<2, dim>
    compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

    /**
     * @copydoc ScalarPolynomialsBase::name()
     */
    std::string
    name() const override;

    /**
     * @copydoc ScalarPolynomialsBase::clone()
     */
    virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
    clone() const override;
  };

  // template functions
  template <int dim>
  template <int order>
  Tensor<order, dim>
  ScalarPolynomial<dim>::compute_derivative(const unsigned int i,
                                            const Point<dim> & p) const
  {
    Tensor<order, dim> der;

    AssertDimension(order, 1);
    const auto grad = compute_grad(i, p);

    for (unsigned int i = 0; i < dim; i++)
      der[i] = grad[i];

    return der;
  }
} // namespace Simplex

DEAL_II_NAMESPACE_CLOSE

#endif
