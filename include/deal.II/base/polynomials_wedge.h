// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#ifndef dealii_base_polynomials_wedge_h
#define dealii_base_polynomials_wedge_h

#include <deal.II/base/config.h>

#include <deal.II/base/ndarray.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/scalar_polynomials_vandermonde_base.h>
#include <deal.II/base/tensor.h>

#include <deal.II/grid/reference_cell.h>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  /**
   * Return the support points of a wedge in a way that can also be compiled
   * for dim==1 and dim==2.
   */
  template <int dim>
  std::vector<Point<dim>>
  get_wedge_support_points(const unsigned int degree)
  {
    if constexpr (dim == 3)
      {
        if (degree == 1)
          return {ReferenceCells::Wedge.vertex(0),
                  ReferenceCells::Wedge.vertex(1),
                  ReferenceCells::Wedge.vertex(2),
                  ReferenceCells::Wedge.vertex(3),
                  ReferenceCells::Wedge.vertex(4),
                  ReferenceCells::Wedge.vertex(5)};
        else if (degree == 2)
          return {// vertices
                  ReferenceCells::Wedge.vertex(0),
                  ReferenceCells::Wedge.vertex(1),
                  ReferenceCells::Wedge.vertex(2),
                  ReferenceCells::Wedge.vertex(3),
                  ReferenceCells::Wedge.vertex(4),
                  ReferenceCells::Wedge.vertex(5),
                  // edges
                  0.5 * (ReferenceCells::Wedge.vertex(0) +
                         ReferenceCells::Wedge.vertex(1)),
                  0.5 * (ReferenceCells::Wedge.vertex(1) +
                         ReferenceCells::Wedge.vertex(2)),
                  0.5 * (ReferenceCells::Wedge.vertex(2) +
                         ReferenceCells::Wedge.vertex(0)),
                  0.5 * (ReferenceCells::Wedge.vertex(3) +
                         ReferenceCells::Wedge.vertex(4)),
                  0.5 * (ReferenceCells::Wedge.vertex(4) +
                         ReferenceCells::Wedge.vertex(5)),
                  0.5 * (ReferenceCells::Wedge.vertex(5) +
                         ReferenceCells::Wedge.vertex(3)),
                  0.5 * (ReferenceCells::Wedge.vertex(3) +
                         ReferenceCells::Wedge.vertex(0)),
                  0.5 * (ReferenceCells::Wedge.vertex(1) +
                         ReferenceCells::Wedge.vertex(4)),
                  0.5 * (ReferenceCells::Wedge.vertex(2) +
                         ReferenceCells::Wedge.vertex(5)),
                  // quad midpoints
                  0.25 * (ReferenceCells::Wedge.vertex(0) +
                          ReferenceCells::Wedge.vertex(1) +
                          ReferenceCells::Wedge.vertex(3) +
                          ReferenceCells::Wedge.vertex(4)),
                  0.25 * (ReferenceCells::Wedge.vertex(1) +
                          ReferenceCells::Wedge.vertex(2) +
                          ReferenceCells::Wedge.vertex(4) +
                          ReferenceCells::Wedge.vertex(5)),
                  0.25 * (ReferenceCells::Wedge.vertex(2) +
                          ReferenceCells::Wedge.vertex(0) +
                          ReferenceCells::Wedge.vertex(5) +
                          ReferenceCells::Wedge.vertex(3))};
        else
          DEAL_II_NOT_IMPLEMENTED();
      }

    DEAL_II_ASSERT_UNREACHABLE();
    return {};
  }
} // namespace internal


/**
 * Polynomials defined on wedge entities. This class can be a basis of FE_WedgeP
 * and FE_WedgeDGP. Jacobi polynomials are used to construct a modal basis. The
 * polynomials are based on the implementation of triangles (see
 * ScalarLagrangePolynomialSimplex) using a tensor product structure. With the
 * modal basis a Vandermonde matrix is calculated which leads to a nodal basis.
 * For computing the values of the nodal basis the inverse of the Vandermonde
 * matrix is multiplied with the modal basis vector evaluated at the evaluation
 * point.
 */
template <int dim>
class ScalarLagrangePolynomialWedge
  : public ScalarPolynomialsVandermondeBase<dim>
{
public:
  /**
   * Make the dimension available to the outside.
   */
  static constexpr unsigned int dimension = dim;

  /*
   * Constructor taking the polynomial @p degree and the support points as
   * arguments.
   */
  ScalarLagrangePolynomialWedge(const unsigned int             degree,
                                const std::vector<Point<dim>> &support_points);

  /*
   * Legacy constructor taking the polynomial @p degree.
   * Note: only works for degrees one and two.
   */
  ScalarLagrangePolynomialWedge(const unsigned int degree);

  std::string
  name() const override;

  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;

private:
  /**
   * Evaluate the orthogonal basis at point @p p. The indices @p i, @p j
   * and @p k correspond to the polynomial degrees of the Jacobi polynomials.
   */
  double
  evaluate_orthogonal_basis_function_by_degree(
    const unsigned int i,
    const unsigned int j,
    const unsigned int k,
    const Point<dim>  &p) const override;

  /**
   * Evaluate the orthogonal basis function @p i at point @p p.
   * This function determines the corresponding indices for the Jacobi
   * polynomials and calls the function taking all indices as arguments.
   */
  double
  evaluate_orthogonal_basis_function(const unsigned int i,
                                     const Point<dim>  &p) const override;

  /**
   * Evaluate the derivative of the orthogonal basis at point @p p.
   * The indices @p i, @p j and @p k correspond to the polynomial degrees of
   * the Jacobi polynomials.
   */
  Tensor<1, dim>
  evaluate_orthogonal_basis_derivative_by_degree(
    const unsigned int i,
    const unsigned int j,
    const unsigned int k,
    const Point<dim>  &p) const override;

  /**
   * Evaluate the derivative of the orthogonal basis function @p i at point
   * @p p. This function determines the corresponding indices for the Jacobi
   * polynomials and calls the function taking all indices as arguments.
   */
  Tensor<1, dim>
  evaluate_orthogonal_basis_derivative(const unsigned int i,
                                       const Point<dim>  &p) const override;
};

DEAL_II_NAMESPACE_CLOSE

#endif
