// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2022 by the deal.II authors
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


#include <deal.II/base/polynomials_raviart_thomas.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>

#include <iomanip>
#include <iostream>
#include <memory>


DEAL_II_NAMESPACE_OPEN


namespace
{
  // Create nodal Raviart-Thomas polynomials as the tensor product of Lagrange
  // polynomials on Gauss-Lobatto points of the given degrees in the normal and
  // tangential directions, respectively (we could also choose Lagrange
  // polynomials on Gauss points but those are slightly more expensive to handle
  // in classes).
  std::vector<std::vector<Polynomials::Polynomial<double>>>
  create_rt_polynomials(const unsigned int dim,
                        const unsigned int normal_degree,
                        const unsigned int tangential_degree)
  {
    std::vector<std::vector<Polynomials::Polynomial<double>>> pols(dim);
    if (normal_degree > 0)
      pols[0] = Polynomials::generate_complete_Lagrange_basis(
        QGaussLobatto<1>(normal_degree + 1).get_points());
    else
      pols[0] = Polynomials::generate_complete_Lagrange_basis(
        QMidpoint<1>().get_points());
    if (tangential_degree > 0)
      for (unsigned int d = 1; d < dim; ++d)
        pols[d] = Polynomials::generate_complete_Lagrange_basis(
          QGaussLobatto<1>(tangential_degree + 1).get_points());
    else
      for (unsigned int d = 1; d < dim; ++d)
        pols[d] = Polynomials::generate_complete_Lagrange_basis(
          QMidpoint<1>().get_points());

    return pols;
  }
} // namespace



template <int dim>
PolynomialsRaviartThomas<dim>::PolynomialsRaviartThomas(
  const unsigned int normal_degree,
  const unsigned int tangential_degree)
  : TensorPolynomialsBase<dim>(std::min(normal_degree, tangential_degree),
                               n_polynomials(normal_degree, tangential_degree))
  , normal_degree(normal_degree)
  , tangential_degree(tangential_degree)
  , polynomial_space(
      create_rt_polynomials(dim, normal_degree, tangential_degree))
{
  // create renumbering of the unknowns from the lexicographic order to the
  // actual order required by the finite element class with unknowns on
  // faces placed first
  const unsigned int n_pols = polynomial_space.n();
  lexicographic_to_hierarchic =
    get_lexicographic_numbering(normal_degree, tangential_degree);

  hierarchic_to_lexicographic =
    Utilities::invert_permutation(lexicographic_to_hierarchic);

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
PolynomialsRaviartThomas<dim>::PolynomialsRaviartThomas(const unsigned int k)
  : PolynomialsRaviartThomas(k + 1, k)
{}



template <int dim>
void
PolynomialsRaviartThomas<dim>::evaluate(
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
        p(c) = unit_point((c + d) % dim);

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
PolynomialsRaviartThomas<dim>::name() const
{
  return "RaviartThomas";
}



template <int dim>
unsigned int
PolynomialsRaviartThomas<dim>::n_polynomials(const unsigned int degree)
{
  return n_polynomials(degree + 1, degree);
}



template <int dim>
unsigned int
PolynomialsRaviartThomas<dim>::n_polynomials(
  const unsigned int normal_degree,
  const unsigned int tangential_degree)
{
  return dim * (normal_degree + 1) *
         Utilities::pow(tangential_degree + 1, dim - 1);
}



template <int dim>
std::vector<unsigned int>
PolynomialsRaviartThomas<dim>::get_lexicographic_numbering(
  const unsigned int normal_degree,
  const unsigned int tangential_degree)
{
  const unsigned int n_dofs_face =
    Utilities::pow(tangential_degree + 1, dim - 1);
  std::vector<unsigned int> lexicographic_numbering;
  // component 1
  for (unsigned int j = 0; j < n_dofs_face; ++j)
    {
      lexicographic_numbering.push_back(j);
      if (normal_degree > 1)
        for (unsigned int i = n_dofs_face * 2 * dim;
             i < n_dofs_face * 2 * dim + normal_degree - 1;
             ++i)
          lexicographic_numbering.push_back(i + j * (normal_degree - 1));
      lexicographic_numbering.push_back(n_dofs_face + j);
    }

  // component 2
  unsigned int layers = (dim == 3) ? tangential_degree + 1 : 1;
  for (unsigned int k = 0; k < layers; ++k)
    {
      unsigned int k_add = k * (tangential_degree + 1);
      for (unsigned int j = n_dofs_face * 2;
           j < n_dofs_face * 2 + tangential_degree + 1;
           ++j)
        lexicographic_numbering.push_back(j + k_add);

      if (normal_degree > 1)
        for (unsigned int i = n_dofs_face * (2 * dim + (normal_degree - 1));
             i < n_dofs_face * (2 * dim + (normal_degree - 1)) +
                   (normal_degree - 1) * (tangential_degree + 1);
             ++i)
          {
            lexicographic_numbering.push_back(i + k_add * tangential_degree);
          }
      for (unsigned int j = n_dofs_face * 3;
           j < n_dofs_face * 3 + tangential_degree + 1;
           ++j)
        lexicographic_numbering.push_back(j + k_add);
    }

  // component 3
  if (dim == 3)
    {
      for (unsigned int i = 4 * n_dofs_face; i < 5 * n_dofs_face; ++i)
        lexicographic_numbering.push_back(i);
      if (normal_degree > 1)
        for (unsigned int i =
               6 * n_dofs_face + n_dofs_face * 2 * (normal_degree - 1);
             i < 6 * n_dofs_face + n_dofs_face * 3 * (normal_degree - 1);
             ++i)
          lexicographic_numbering.push_back(i);
      for (unsigned int i = 5 * n_dofs_face; i < 6 * n_dofs_face; ++i)
        lexicographic_numbering.push_back(i);
    }

  return lexicographic_numbering;
}



template <int dim>
std::unique_ptr<TensorPolynomialsBase<dim>>
PolynomialsRaviartThomas<dim>::clone() const
{
  return std::make_unique<PolynomialsRaviartThomas<dim>>(*this);
}



template <int dim>
std::vector<Point<dim>>
PolynomialsRaviartThomas<dim>::get_polynomial_support_points() const
{
  Assert(dim > 0 && dim <= 3, ExcImpossibleInDim(dim));
  const Quadrature<1> tangential(
    tangential_degree == 0 ?
      static_cast<Quadrature<1>>(QMidpoint<1>()) :
      static_cast<Quadrature<1>>(QGaussLobatto<1>(tangential_degree + 1)));
  const Quadrature<1> normal(
    normal_degree == 0 ?
      static_cast<Quadrature<1>>(QMidpoint<1>()) :
      static_cast<Quadrature<1>>(QGaussLobatto<1>(normal_degree + 1)));
  const QAnisotropic<dim> quad =
    (dim == 1 ? QAnisotropic<dim>(normal) :
                (dim == 2 ? QAnisotropic<dim>(normal, tangential) :
                            QAnisotropic<dim>(normal, tangential, tangential)));

  const unsigned int      n_sub = polynomial_space.n();
  std::vector<Point<dim>> points(dim * n_sub);
  points.resize(n_polynomials(normal_degree, tangential_degree));
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int i = 0; i < n_sub; ++i)
      for (unsigned int c = 0; c < dim; ++c)
        points[lexicographic_to_hierarchic[i + d * n_sub]][(c + d) % dim] =
          quad.point(renumber_aniso[d][i])[c];
  return points;
}



template class PolynomialsRaviartThomas<1>;
template class PolynomialsRaviartThomas<2>;
template class PolynomialsRaviartThomas<3>;


DEAL_II_NAMESPACE_CLOSE
