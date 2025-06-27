// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/polynomials_vector_anisotropic.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>

#include <iomanip>
#include <iostream>
#include <memory>


DEAL_II_NAMESPACE_OPEN


namespace
{
  /**
   * Creating support points for normal and tangential directions.
   * Since the polynomial degrees differ for the 2 directions, number of
   * support points also differ, respectively.
   */
  std::vector<std::vector<Point<1>>>
  create_anisotropic_support_points(const unsigned int dim,
                                    const unsigned int normal_degree,
                                    const unsigned int tangential_degree)
  {
    Assert(dim > 0 && dim <= 3, ExcImpossibleInDim(dim));

    std::vector<std::vector<Point<1>>> points_normal_tangential(dim);
    points_normal_tangential[0] =
      normal_degree == 0 ? QMidpoint<1>().get_points() :
                           QGaussLobatto<1>(normal_degree + 1).get_points();
    for (unsigned int d = 1; d < dim; ++d)
      {
        points_normal_tangential[d] =
          tangential_degree == 0 ?
            QMidpoint<1>().get_points() :
            QGaussLobatto<1>(tangential_degree + 1).get_points();
      }
    return points_normal_tangential;
  }



  // Create nodal anisotropic polynomials as the tensor product of Lagrange
  // polynomials on Gauss-Lobatto points of the given degrees in the normal and
  // tangential directions, respectively (we could also choose Lagrange
  // polynomials on Gauss points but those are slightly more expensive to handle
  // in classes).
  std::vector<std::vector<Polynomials::Polynomial<double>>>
  create_aniso_polynomials(const unsigned int dim,
                           const unsigned int normal_degree,
                           const unsigned int tangential_degree)
  {
    std::vector<std::vector<Point<1>>> points_aniso;
    points_aniso =
      create_anisotropic_support_points(dim, normal_degree, tangential_degree);
    std::vector<std::vector<Polynomials::Polynomial<double>>> pols(dim);
    pols[0] = Polynomials::generate_complete_Lagrange_basis(points_aniso[0]);
    for (unsigned int d = 1; d < dim; ++d)
      pols[d] = Polynomials::generate_complete_Lagrange_basis(points_aniso[d]);
    return pols;
  }
} // namespace



template <int dim>
PolynomialsVectorAnisotropic<dim>::PolynomialsVectorAnisotropic(
  const unsigned int               normal_degree,
  const unsigned int               tangential_degree,
  const std::vector<unsigned int> &polynomial_ordering)
  : TensorPolynomialsBase<dim>(std::min(normal_degree, tangential_degree),
                               n_polynomials(normal_degree, tangential_degree))
  , normal_degree(normal_degree)
  , tangential_degree(tangential_degree)
  , polynomial_space(
      create_aniso_polynomials(dim, normal_degree, tangential_degree))
  , lexicographic_to_hierarchic(polynomial_ordering)
  , hierarchic_to_lexicographic(
      Utilities::invert_permutation(lexicographic_to_hierarchic))
{
  // create renumbering of the unknowns from the lexicographic order to the
  // actual order required by the finite element class with unknowns on
  // faces placed first
  const unsigned int n_pols = polynomial_space.n();

  // since we only store an anisotropic polynomial for the first component,
  // we set up a second numbering to switch out the actual coordinate
  // direction
  renumber_aniso[0].resize(n_pols);
  for (unsigned int i = 0; i < n_pols; ++i)
    renumber_aniso[0][i] = i;
  if (dim == 2)
    {
      // switch x and y component (i and j loops)
      renumber_aniso[1].resize(n_pols);
      for (unsigned int j = 0; j < normal_degree + 1; ++j)
        for (unsigned int i = 0; i < tangential_degree + 1; ++i)
          renumber_aniso[1][j * (tangential_degree + 1) + i] =
            j + i * (normal_degree + 1);
    }
  if (dim == 3)
    {
      // switch x, y, and z component (i, j, k) -> (j, k, i)
      renumber_aniso[1].resize(n_pols);
      for (unsigned int k = 0; k < tangential_degree + 1; ++k)
        for (unsigned int j = 0; j < normal_degree + 1; ++j)
          for (unsigned int i = 0; i < tangential_degree + 1; ++i)
            renumber_aniso[1][(k * (normal_degree + 1) + j) *
                                (tangential_degree + 1) +
                              i] =
              j + (normal_degree + 1) * (k + i * (tangential_degree + 1));

      // switch x, y, and z component (i, j, k) -> (k, i, j)
      renumber_aniso[2].resize(n_pols);
      for (unsigned int k = 0; k < normal_degree + 1; ++k)
        for (unsigned int j = 0; j < tangential_degree + 1; ++j)
          for (unsigned int i = 0; i < tangential_degree + 1; ++i)
            renumber_aniso[2][(k * (tangential_degree + 1) + j) *
                                (tangential_degree + 1) +
                              i] =
              k + (normal_degree + 1) * (i + j * (tangential_degree + 1));
    }
}



template <int dim>
void
PolynomialsVectorAnisotropic<dim>::evaluate(
  const Point<dim>            &unit_point,
  std::vector<Tensor<1, dim>> &values,
  std::vector<Tensor<2, dim>> &grads,
  std::vector<Tensor<3, dim>> &grad_grads,
  std::vector<Tensor<4, dim>> &third_derivatives,
  std::vector<Tensor<5, dim>> &fourth_derivatives) const
{
  Assert(values.size() == this->n() || values.empty(),
         ExcDimensionMismatch(values.size(), this->n()));
  Assert(grads.size() == this->n() || grads.empty(),
         ExcDimensionMismatch(grads.size(), this->n()));
  Assert(grad_grads.size() == this->n() || grad_grads.empty(),
         ExcDimensionMismatch(grad_grads.size(), this->n()));
  Assert(third_derivatives.size() == this->n() || third_derivatives.empty(),
         ExcDimensionMismatch(third_derivatives.size(), this->n()));
  Assert(fourth_derivatives.size() == this->n() || fourth_derivatives.empty(),
         ExcDimensionMismatch(fourth_derivatives.size(), this->n()));

  std::vector<double>         p_values;
  std::vector<Tensor<1, dim>> p_grads;
  std::vector<Tensor<2, dim>> p_grad_grads;
  std::vector<Tensor<3, dim>> p_third_derivatives;
  std::vector<Tensor<4, dim>> p_fourth_derivatives;

  const unsigned int n_sub = polynomial_space.n();
  p_values.resize((values.empty()) ? 0 : n_sub);
  p_grads.resize((grads.empty()) ? 0 : n_sub);
  p_grad_grads.resize((grad_grads.empty()) ? 0 : n_sub);
  p_third_derivatives.resize((third_derivatives.empty()) ? 0 : n_sub);
  p_fourth_derivatives.resize((fourth_derivatives.empty()) ? 0 : n_sub);

  for (unsigned int d = 0; d < dim; ++d)
    {
      // First we copy the point. The polynomial space for component d
      // consists of polynomials of degree k in x_d and degree k+1 in the
      // other variables. in order to simplify this, we use the same
      // AnisotropicPolynomial space and simply rotate the coordinates
      // through all directions.
      Point<dim> p;
      for (unsigned int c = 0; c < dim; ++c)
        p[c] = unit_point[(c + d) % dim];

      polynomial_space.evaluate(p,
                                p_values,
                                p_grads,
                                p_grad_grads,
                                p_third_derivatives,
                                p_fourth_derivatives);

      for (unsigned int i = 0; i < p_values.size(); ++i)
        values[lexicographic_to_hierarchic[i + d * n_sub]][d] =
          p_values[renumber_aniso[d][i]];

      for (unsigned int i = 0; i < p_grads.size(); ++i)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          grads[lexicographic_to_hierarchic[i + d * n_sub]][d][(d1 + d) % dim] =
            p_grads[renumber_aniso[d][i]][d1];

      for (unsigned int i = 0; i < p_grad_grads.size(); ++i)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            grad_grads[lexicographic_to_hierarchic[i + d * n_sub]][d]
                      [(d1 + d) % dim][(d2 + d) % dim] =
                        p_grad_grads[renumber_aniso[d][i]][d1][d2];

      for (unsigned int i = 0; i < p_third_derivatives.size(); ++i)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            for (unsigned int d3 = 0; d3 < dim; ++d3)
              third_derivatives[lexicographic_to_hierarchic[i + d * n_sub]][d]
                               [(d1 + d) % dim][(d2 + d) % dim]
                               [(d3 + d) % dim] =
                                 p_third_derivatives[renumber_aniso[d][i]][d1]
                                                    [d2][d3];

      for (unsigned int i = 0; i < p_fourth_derivatives.size(); ++i)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            for (unsigned int d3 = 0; d3 < dim; ++d3)
              for (unsigned int d4 = 0; d4 < dim; ++d4)
                fourth_derivatives[lexicographic_to_hierarchic[i + d * n_sub]]
                                  [d][(d1 + d) % dim][(d2 + d) % dim]
                                  [(d3 + d) % dim][(d4 + d) % dim] =
                                    p_fourth_derivatives[renumber_aniso[d][i]]
                                                        [d1][d2][d3][d4];
    }
}



template <int dim>
std::string
PolynomialsVectorAnisotropic<dim>::name() const
{
  return "VectorAnisotropic";
}



template <int dim>
unsigned int
PolynomialsVectorAnisotropic<dim>::n_polynomials(
  const unsigned int normal_degree,
  const unsigned int tangential_degree)
{
  return dim * (normal_degree + 1) *
         Utilities::pow(tangential_degree + 1, dim - 1);
}



template <int dim>
unsigned int
PolynomialsVectorAnisotropic<dim>::get_tangential_degree() const
{
  return tangential_degree;
}



template <int dim>
unsigned int
PolynomialsVectorAnisotropic<dim>::get_normal_degree() const
{
  return normal_degree;
}



template <int dim>
std::unique_ptr<TensorPolynomialsBase<dim>>
PolynomialsVectorAnisotropic<dim>::clone() const
{
  return std::make_unique<PolynomialsVectorAnisotropic<dim>>(*this);
}



template <int dim>
std::vector<Point<dim>>
PolynomialsVectorAnisotropic<dim>::get_polynomial_support_points() const
{
  Assert(dim > 0 && dim <= 3, ExcImpossibleInDim(dim));
  const std::vector<std::vector<Point<1>>> points_aniso =
    create_anisotropic_support_points(dim, normal_degree, tangential_degree);
  const unsigned int      n_sub = polynomial_space.n();
  std::vector<Point<dim>> points(dim * n_sub);
  points.resize(n_polynomials(normal_degree, tangential_degree));
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int i = 0; i < n_sub; ++i)
      {
        unsigned int                  renumbered_index = renumber_aniso[d][i];
        std::array<unsigned int, dim> indices_points_anisotropic;
        indices_points_anisotropic[0] = renumbered_index % (normal_degree + 1);
        if (dim > 1)
          {
            renumbered_index /= (normal_degree + 1);
            indices_points_anisotropic[1] =
              renumbered_index % (tangential_degree + 1);
          }
        if (dim > 2)
          indices_points_anisotropic[2] =
            renumbered_index / (tangential_degree + 1);
        for (unsigned int c = 0; c < dim; ++c)
          {
            points[lexicographic_to_hierarchic[i + d * n_sub]][(c + d) % dim] =
              points_aniso[c][indices_points_anisotropic[c]][0];
          }
      }
  return points;
}



template class PolynomialsVectorAnisotropic<1>;
template class PolynomialsVectorAnisotropic<2>;
template class PolynomialsVectorAnisotropic<3>;


DEAL_II_NAMESPACE_CLOSE
