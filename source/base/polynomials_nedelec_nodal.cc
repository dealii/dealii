// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_nedelec_nodal.h>
#include <deal.II/base/quadrature_lib.h>

#include <iomanip>
#include <iostream>
#include <memory>

DEAL_II_NAMESPACE_OPEN


template <int dim>
PolynomialsNedelecNodal<dim>::PolynomialsNedelecNodal(const unsigned int degree)
  : TensorPolynomialsBase<dim>(degree, n_polynomials(degree))
  , polynomial_space(create_polynomials(degree))
  , deg(degree)
{
  renumber_lexicographic_to_hierarchic = get_lexicographic_numbering(degree);

  renumber_hierarchic_to_lexicographic =
    Utilities::invert_permutation(renumber_lexicographic_to_hierarchic);

  const unsigned int n_pols = polynomial_space.n();

  renumber_aniso[0].resize(n_pols);
  for (unsigned int i = 0; i < n_pols; ++i)
    renumber_aniso[0][i] = i;

  if (dim == 2)
    {
      // switch x and y component (i and j loops)
      renumber_aniso[1].resize(n_pols);
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 2; ++i)
          renumber_aniso[1][j * (degree + 2) + i] = j + i * (degree + 1);
    }
  if (dim == 3)
    {
      // switch x, y, and z component (i, j, k) -> (j, k, i)
      renumber_aniso[1].resize(n_pols);
      for (unsigned int k = 0; k < degree + 2; ++k)
        for (unsigned int j = 0; j < degree + 1; ++j)
          for (unsigned int i = 0; i < degree + 2; ++i)
            renumber_aniso[1][(k * (degree + 1) + j) * (degree + 2) + i] =
              j + (degree + 1) * (i * (degree + 2) + k);

      // switch x, y, and z component (i, j, k) -> (k, i, j)
      renumber_aniso[2].resize(n_pols);
      for (unsigned int k = 0; k < degree + 1; ++k)
        for (unsigned int j = 0; j < degree + 2; ++j)
          for (unsigned int i = 0; i < degree + 2; ++i)
            {
              // std::cout<<"index: ("<<k<<","<<j<<","<<i<<"): ";
              // std::cout<<(k * (degree + 2) + j) * (degree + 2) + i<<"\n";
              renumber_aniso[2][(k * (degree + 2) + j) * (degree + 2) + i] =
                (j * (degree + 2) + i) * (degree + 1) + k;
            }
    }
}

template <int dim>
std::vector<std::vector<Polynomials::Polynomial<double>>>
PolynomialsNedelecNodal<dim>::create_polynomials(const unsigned int degree)
{
  std::vector<std::vector<Polynomials::Polynomial<double>>> pols(dim);
  pols[0] = Polynomials::generate_complete_Lagrange_basis(
    QGauss<1>(degree + 1).get_points());
  for (unsigned int d = 1; d < dim; ++d)
    pols[d] = Polynomials::generate_complete_Lagrange_basis(
      QGaussLobatto<1>(degree + 2).get_points());

  return pols;
}


// Compute the values, gradients
// and double gradients of the
// polynomial at the given point.
template <int dim>
void
PolynomialsNedelecNodal<dim>::evaluate(
  const Point<dim>            &unit_point,
  std::vector<Tensor<1, dim>> &values,
  std::vector<Tensor<2, dim>> &grads,
  std::vector<Tensor<3, dim>> &grad_grads,
  std::vector<Tensor<4, dim>> &third_derivatives,
  std::vector<Tensor<5, dim>> &fourth_derivatives) const
{
  Assert(values.size() == this->n() || values.size() == 0,
         ExcDimensionMismatch(values.size(), this->n()));
  Assert(grads.size() == this->n() || grads.size() == 0,
         ExcDimensionMismatch(grads.size(), this->n()));
  Assert(grad_grads.size() == this->n() || grad_grads.size() == 0,
         ExcDimensionMismatch(grad_grads.size(), this->n()));
  Assert(third_derivatives.size() == this->n() || third_derivatives.size() == 0,
         ExcDimensionMismatch(third_derivatives.size(), this->n()));
  Assert(fourth_derivatives.size() == this->n() ||
           fourth_derivatives.size() == 0,
         ExcDimensionMismatch(fourth_derivatives.size(), this->n()));

  std::vector<double>         p_values;
  std::vector<Tensor<1, dim>> p_grads;
  std::vector<Tensor<2, dim>> p_grad_grads;
  std::vector<Tensor<3, dim>> p_third_derivatives;
  std::vector<Tensor<4, dim>> p_fourth_derivatives;

  const unsigned int n_sub = polynomial_space.n();
  // std::cout<<"n_sub\t"<<n_sub<<std::endl;
  p_values.resize((values.size() == 0) ? 0 : n_sub);
  p_grads.resize((grads.size() == 0) ? 0 : n_sub);
  p_grad_grads.resize((grad_grads.size() == 0) ? 0 : n_sub);
  p_third_derivatives.resize((third_derivatives.size() == 0) ? 0 : n_sub);
  p_fourth_derivatives.resize((fourth_derivatives.size() == 0) ? 0 : n_sub);

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
        values[renumber_lexicographic_to_hierarchic[i + d * n_sub]][d] =
          p_values[renumber_aniso[d][i]];

      for (unsigned int i = 0; i < p_grads.size(); ++i)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          grads[renumber_lexicographic_to_hierarchic[i + d * n_sub]][d]
               [(d1 + d) % dim] = p_grads[renumber_aniso[d][i]][d1];

      for (unsigned int i = 0; i < p_grad_grads.size(); ++i)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            grad_grads[renumber_lexicographic_to_hierarchic[i + d * n_sub]][d]
                      [(d1 + d) % dim][(d2 + d) % dim] =
                        p_grad_grads[renumber_aniso[d][i]][d1][d2];

      for (unsigned int i = 0; i < p_third_derivatives.size(); ++i)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            for (unsigned int d3 = 0; d3 < dim; ++d3)
              third_derivatives
                [renumber_lexicographic_to_hierarchic[i + d * n_sub]][d]
                [(d1 + d) % dim][(d2 + d) % dim][(d3 + d) % dim] =
                  p_third_derivatives[renumber_aniso[d][i]][d1][d2][d3];

      for (unsigned int i = 0; i < p_fourth_derivatives.size(); ++i)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            for (unsigned int d3 = 0; d3 < dim; ++d3)
              for (unsigned int d4 = 0; d4 < dim; ++d4)
                fourth_derivatives
                  [renumber_lexicographic_to_hierarchic[i + d * n_sub]][d]
                  [(d1 + d) % dim][(d2 + d) % dim][(d3 + d) % dim]
                  [(d4 + d) % dim] =
                    p_fourth_derivatives[renumber_aniso[d][i]][d1][d2][d3][d4];
    }
}


template <int dim>
unsigned int
PolynomialsNedelecNodal<dim>::n_polynomials(const unsigned int k)
{
  switch (dim)
    {
      case 1:
        return k + 1;

      case 2:
        return 2 * (k + 1) * (k + 2);

      case 3:
        return 3 * (k + 1) * (k + 2) * (k + 2);

      default:
        {
          DEAL_II_NOT_IMPLEMENTED();
          return 0;
        }
    }
}

template <int dim>
std::vector<unsigned int>
PolynomialsNedelecNodal<dim>::get_lexicographic_numbering(
  const unsigned int degree)
{
  const unsigned int n_pols = (dim == 2) ?
                                (degree + 1) * (degree + 2) :
                                (degree + 1) * (degree + 2) * (degree + 2);

  std::vector<unsigned int> renumber_hierarchic_to_lexicographic;

  // edges 0-3
  for (unsigned int i = 0; i < degree + 1; ++i)
    renumber_hierarchic_to_lexicographic.push_back(n_pols + i * (degree + 2));
  for (unsigned int i = 0; i < degree + 1; ++i)
    renumber_hierarchic_to_lexicographic.push_back(n_pols + degree + 1 +
                                                   i * (degree + 2));
  for (unsigned int i = 0; i < degree + 1; ++i)
    renumber_hierarchic_to_lexicographic.push_back(i);
  for (unsigned int i = 0; i < degree + 1; ++i)
    renumber_hierarchic_to_lexicographic.push_back((degree + 1) * (degree + 1) +
                                                   i);
  if (dim == 2)
    {
      // quads for 2D
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(j * (degree + 1) + i);
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(n_pols +
                                                         j * (degree + 2) + i);
    }
  else if (dim == 3)
    {
      // edges 4-7
      const unsigned int offset_z = (degree + 1) * (degree + 1) * (degree + 2);
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          n_pols + i * (degree + 2) + offset_z);
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          n_pols + degree + 1 + i * (degree + 2) + offset_z);
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(i + offset_z);
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          (degree + 1) * (degree + 1) + i + offset_z);

      // edges 8-11
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          2 * n_pols + i * (degree + 2) * (degree + 2));
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          2 * n_pols + degree + 1 + i * (degree + 2) * (degree + 2));
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          2 * n_pols + (degree + 1) * (degree + 2) +
          i * (degree + 2) * (degree + 2));
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          2 * n_pols + (degree + 2) * (degree + 2) - 1 +
          i * (degree + 2) * (degree + 2));

      // quad 0
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            n_pols + i * (degree + 2) + j * (degree + 2) * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            2 * n_pols + i * (degree + 2) + j * (degree + 2) * (degree + 2));
      // quad 1
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            n_pols + degree + 1 + i * (degree + 2) +
            j * (degree + 2) * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            2 * n_pols + degree + 1 + i * (degree + 2) +
            j * (degree + 2) * (degree + 2));
      // quad 2
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(i + j * (degree + 2) *
                                                               (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            2 * n_pols + i + j * (degree + 2) * (degree + 2));
      // quad 3
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            (degree + 1) * (degree + 1) + i + j * (degree + 2) * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            2 * n_pols + (degree + 1) * (degree + 2) + i +
            j * (degree + 2) * (degree + 2));
      // quad 4
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(i + j * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(n_pols + i +
                                                         j * (degree + 2));
      // quad 5
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            (degree + 1) * (degree + 1) * (degree + 2) + i + j * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            n_pols + (degree + 1) * (degree + 1) * (degree + 2) + i +
            j * (degree + 2));

      // hexes
      for (unsigned int k = 1; k < degree + 1; ++k)
        for (unsigned int j = 1; j < degree + 1; ++j)
          for (unsigned int i = 0; i < degree + 1; ++i)
            renumber_hierarchic_to_lexicographic.push_back(
              k * (degree + 1) * (degree + 2) + j * (degree + 1) + i);
      for (unsigned int k = 1; k < degree + 1; ++k)
        for (unsigned int j = 0; j < degree + 1; ++j)
          for (unsigned int i = 1; i < degree + 1; ++i)
            renumber_hierarchic_to_lexicographic.push_back(
              n_pols + k * (degree + 2) * (degree + 1) + j * (degree + 2) + i);
      for (unsigned int k = 0; k < degree + 1; ++k)
        for (unsigned int j = 1; j < degree + 1; ++j)
          for (unsigned int i = 1; i < degree + 1; ++i)
            renumber_hierarchic_to_lexicographic.push_back(
              2 * n_pols + k * (degree + 2) * (degree + 2) + j * (degree + 2) +
              i);
    }

  return Utilities::invert_permutation(renumber_hierarchic_to_lexicographic);
}

template <int dim>
std::vector<Point<dim>>
PolynomialsNedelecNodal<dim>::get_polynomial_support_points() const
{
  Assert(dim > 0 && dim <= 3, ExcImpossibleInDim(dim));
  const Quadrature<1> tangential(
    static_cast<Quadrature<1>>(QGaussLobatto<1>(deg + 2)));
  const Quadrature<1> normal(static_cast<Quadrature<1>>(QGauss<1>(deg + 1)));
  const QAnisotropic<dim> quad =
    (dim == 1 ? QAnisotropic<dim>(normal) :
                (dim == 2 ? QAnisotropic<dim>(normal, tangential) :
                            QAnisotropic<dim>(normal, tangential, tangential)));

  const unsigned int      n_sub = polynomial_space.n();
  std::vector<Point<dim>> points(dim * n_sub);
  points.resize(n_polynomials(deg));
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int i = 0; i < n_sub; ++i)
      for (unsigned int c = 0; c < dim; ++c)
        points[renumber_lexicographic_to_hierarchic[i + d * n_sub]]
              [(c + d) % dim] = quad.point(renumber_aniso[d][i])[c];
  return points;
}


template <int dim>
std::unique_ptr<TensorPolynomialsBase<dim>>
PolynomialsNedelecNodal<dim>::clone() const
{
  return std::make_unique<PolynomialsNedelecNodal<dim>>(*this);
}

template <int dim>
inline std::string
PolynomialsNedelecNodal<dim>::name() const
{
  return "NedelecNodal";
}



template class PolynomialsNedelecNodal<1>;
template class PolynomialsNedelecNodal<2>;
template class PolynomialsNedelecNodal<3>;


DEAL_II_NAMESPACE_CLOSE
