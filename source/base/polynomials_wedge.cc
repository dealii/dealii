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


#include <deal.II/base/ndarray.h>
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/polynomials_wedge.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  unsigned int
  compute_n_polynomials_wedge(const unsigned int dim, const unsigned int degree)
  {
    if (dim == 3)
      {
        if (degree == 1)
          return 6;
        if (degree == 2)
          return 18;
      }

    Assert(false, ExcNotImplemented());

    return 0;
  }
} // namespace



template <int dim>
ScalarLagrangePolynomialWedge<dim>::ScalarLagrangePolynomialWedge(
  const unsigned int degree)
  : ScalarPolynomialsBase<dim>(degree, compute_n_polynomials_wedge(dim, degree))
  , poly_tri(BarycentricPolynomials<2>::get_fe_p_basis(degree))
  , poly_line(BarycentricPolynomials<1>::get_fe_p_basis(degree))
{}


namespace
{
  /**
   * Decompose the shape-function index of a linear wedge into an index
   * to access the right shape function within the triangle and and within
   * the line.
   */
  static const constexpr ndarray<unsigned int, 6, 2> wedge_table_1{
    {{{0, 0}}, {{1, 0}}, {{2, 0}}, {{0, 1}}, {{1, 1}}, {{2, 1}}}};

  /**
   * Decompose the shape-function index of a quadratic wedge into an index
   * to access the right shape function within the triangle and and within
   * the line.
   */
  static const constexpr ndarray<unsigned int, 18, 2> wedge_table_2{{{{0, 0}},
                                                                     {{1, 0}},
                                                                     {{2, 0}},
                                                                     {{0, 1}},
                                                                     {{1, 1}},
                                                                     {{2, 1}},
                                                                     {{3, 0}},
                                                                     {{4, 0}},
                                                                     {{5, 0}},
                                                                     {{3, 1}},
                                                                     {{4, 1}},
                                                                     {{5, 1}},
                                                                     {{0, 2}},
                                                                     {{1, 2}},
                                                                     {{2, 2}},
                                                                     {{3, 2}},
                                                                     {{4, 2}},
                                                                     {{5, 2}}}};
} // namespace



template <int dim>
double
ScalarLagrangePolynomialWedge<dim>::compute_value(const unsigned int i,
                                                  const Point<dim> & p) const
{
  const auto pair = this->degree() == 1 ? wedge_table_1[i] : wedge_table_2[i];

  const Point<2> p_tri(p[0], p[1]);
  const auto     v_tri = poly_tri.compute_value(pair[0], p_tri);

  const Point<1> p_line(p[2]);
  const auto     v_line = poly_line.compute_value(pair[1], p_line);

  return v_tri * v_line;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialWedge<dim>::compute_grad(const unsigned int i,
                                                 const Point<dim> & p) const
{
  const auto pair = this->degree() == 1 ? wedge_table_1[i] : wedge_table_2[i];

  const Point<2> p_tri(p[0], p[1]);
  const auto     v_tri = poly_tri.compute_value(pair[0], p_tri);
  const auto     g_tri = poly_tri.compute_grad(pair[0], p_tri);

  const Point<1> p_line(p[2]);
  const auto     v_line = poly_line.compute_value(pair[1], p_line);
  const auto     g_line = poly_line.compute_grad(pair[1], p_line);

  Tensor<1, dim> grad;
  grad[0] = g_tri[0] * v_line;
  grad[1] = g_tri[1] * v_line;
  grad[2] = v_tri * g_line[0];

  return grad;
}



template <int dim>
Tensor<2, dim>
ScalarLagrangePolynomialWedge<dim>::compute_grad_grad(const unsigned int i,
                                                      const Point<dim> &p) const
{
  (void)i;
  (void)p;

  Assert(false, ExcNotImplemented());
  return Tensor<2, dim>();
}



template <int dim>
void
ScalarLagrangePolynomialWedge<dim>::evaluate(
  const Point<dim> &           unit_point,
  std::vector<double> &        values,
  std::vector<Tensor<1, dim>> &grads,
  std::vector<Tensor<2, dim>> &grad_grads,
  std::vector<Tensor<3, dim>> &third_derivatives,
  std::vector<Tensor<4, dim>> &fourth_derivatives) const
{
  (void)grads;
  (void)grad_grads;
  (void)third_derivatives;
  (void)fourth_derivatives;

  if (values.size() == this->n())
    for (unsigned int i = 0; i < this->n(); i++)
      values[i] = compute_value(i, unit_point);

  if (grads.size() == this->n())
    for (unsigned int i = 0; i < this->n(); i++)
      grads[i] = compute_grad(i, unit_point);
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialWedge<dim>::compute_1st_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_grad(i, p);
}



template <int dim>
Tensor<2, dim>
ScalarLagrangePolynomialWedge<dim>::compute_2nd_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  (void)i;
  (void)p;

  Assert(false, ExcNotImplemented());

  return {};
}



template <int dim>
Tensor<3, dim>
ScalarLagrangePolynomialWedge<dim>::compute_3rd_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  (void)i;
  (void)p;

  Assert(false, ExcNotImplemented());

  return {};
}



template <int dim>
Tensor<4, dim>
ScalarLagrangePolynomialWedge<dim>::compute_4th_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  (void)i;
  (void)p;

  Assert(false, ExcNotImplemented());

  return {};
}



template <int dim>
std::string
ScalarLagrangePolynomialWedge<dim>::name() const
{
  return "ScalarLagrangePolynomialWedge";
}



template <int dim>
std::unique_ptr<ScalarPolynomialsBase<dim>>
ScalarLagrangePolynomialWedge<dim>::clone() const
{
  return std::make_unique<ScalarLagrangePolynomialWedge<dim>>(*this);
}



template class ScalarLagrangePolynomialWedge<1>;
template class ScalarLagrangePolynomialWedge<2>;
template class ScalarLagrangePolynomialWedge<3>;

DEAL_II_NAMESPACE_CLOSE
