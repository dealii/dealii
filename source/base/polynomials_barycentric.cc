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


#include <deal.II/base/polynomials_barycentric.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * Get the highest degree of the barycentric polynomial (in Cartesian
   * coordinates).
   */
  template <int dim>
  unsigned int
  get_degree(const std::vector<BarycentricPolynomial<dim>> &polys)
  {
    // Since the first variable in a simplex polynomial is, e.g., in 2D,
    //
    // t0 = 1 - x - y
    //
    // (that is, it depends on the Cartesian variables), we have to compute
    // its degree separately. An example: t0*t1*t2 has degree 1 in the affine
    // polynomial basis but is degree 2 in the Cartesian polynomial basis.
    std::size_t max_degree = 0;
    for (const auto &poly : polys)
      {
        const TableIndices<dim + 1> degrees = poly.degrees();

        const auto  degree_0 = degrees[0];
        std::size_t degree_d = 0;
        for (unsigned int d = 1; d < dim + 1; ++d)
          degree_d = std::max(degree_d, degrees[d]);

        max_degree = std::max(max_degree, degree_d + degree_0);
      }

    return max_degree;
  }
} // namespace internal


template <int dim>
BarycentricPolynomials<dim>
BarycentricPolynomials<dim>::get_fe_p_basis(const unsigned int degree)
{
  std::vector<BarycentricPolynomial<dim>> polys;

  auto M = [](const unsigned int d) {
    return BarycentricPolynomial<dim, double>::monomial(d);
  };
  switch (degree)
    {
      case 0:
        polys.push_back(0 * M(0) + 1);
        break;
      case 1:
        {
          for (unsigned int d = 0; d < dim + 1; ++d)
            polys.push_back(M(d));
          break;
        }
      case 2:
        {
          for (unsigned int d = 0; d < dim + 1; ++d)
            polys.push_back(M(d) * (2 * M(d) - 1));
          polys.push_back(4 * M(1) * M(0));
          if (dim >= 2)
            {
              polys.push_back(4 * M(1) * M(2));
              polys.push_back(4 * M(2) * M(0));
            }
          if (dim == 3)
            {
              polys.push_back(4 * M(3) * M(0));
              polys.push_back(4 * M(1) * M(3));
              polys.push_back(4 * M(2) * M(3));
            }
          break;
        }
      default:
        Assert(false, ExcNotImplemented());
    }

  return BarycentricPolynomials<dim>(polys);
}



template <int dim>
BarycentricPolynomials<dim>::BarycentricPolynomials(
  const std::vector<BarycentricPolynomial<dim>> &polynomials)
  : ScalarPolynomialsBase<dim>(internal::get_degree(polynomials),
                               polynomials.size())
{
  polys = polynomials;

  poly_grads.reinit({polynomials.size(), dim});
  poly_hessians.reinit({polynomials.size(), dim, dim});
  poly_third_derivatives.reinit({polynomials.size(), dim, dim, dim});
  poly_fourth_derivatives.reinit({polynomials.size(), dim, dim, dim, dim});

  for (std::size_t i = 0; i < polynomials.size(); ++i)
    {
      // gradients
      for (unsigned int d = 0; d < dim; ++d)
        poly_grads[i][d] = polynomials[i].derivative(d);

      // hessians
      for (unsigned int d0 = 0; d0 < dim; ++d0)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          poly_hessians[i][d0][d1] = poly_grads[i][d0].derivative(d1);

      // third derivatives
      for (unsigned int d0 = 0; d0 < dim; ++d0)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            poly_third_derivatives[i][d0][d1][d2] =
              poly_hessians[i][d0][d1].derivative(d2);

      // fourth derivatives
      for (unsigned int d0 = 0; d0 < dim; ++d0)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            for (unsigned int d3 = 0; d3 < dim; ++d3)
              poly_fourth_derivatives[i][d0][d1][d2][d3] =
                poly_third_derivatives[i][d0][d1][d2].derivative(d3);
    }
}



template <int dim>
void
BarycentricPolynomials<dim>::evaluate(
  const Point<dim> &           unit_point,
  std::vector<double> &        values,
  std::vector<Tensor<1, dim>> &grads,
  std::vector<Tensor<2, dim>> &grad_grads,
  std::vector<Tensor<3, dim>> &third_derivatives,
  std::vector<Tensor<4, dim>> &fourth_derivatives) const
{
  Assert(values.size() == this->n() || values.size() == 0,
         ExcDimensionMismatch2(values.size(), this->n(), 0));
  Assert(grads.size() == this->n() || grads.size() == 0,
         ExcDimensionMismatch2(grads.size(), this->n(), 0));
  Assert(grad_grads.size() == this->n() || grad_grads.size() == 0,
         ExcDimensionMismatch2(grad_grads.size(), this->n(), 0));
  Assert(third_derivatives.size() == this->n() || third_derivatives.size() == 0,
         ExcDimensionMismatch2(third_derivatives.size(), this->n(), 0));
  Assert(fourth_derivatives.size() == this->n() ||
           fourth_derivatives.size() == 0,
         ExcDimensionMismatch2(fourth_derivatives.size(), this->n(), 0));

  for (std::size_t i = 0; i < polys.size(); ++i)
    {
      if (values.size() == this->n())
        values[i] = polys[i].value(unit_point);

      // gradients
      if (grads.size() == this->n())
        for (unsigned int d = 0; d < dim; ++d)
          grads[i][d] = poly_grads[i][d].value(unit_point);

      // hessians
      if (grad_grads.size() == this->n())
        for (unsigned int d0 = 0; d0 < dim; ++d0)
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            grad_grads[i][d0][d1] = poly_hessians[i][d0][d1].value(unit_point);

      // third derivatives
      if (third_derivatives.size() == this->n())
        for (unsigned int d0 = 0; d0 < dim; ++d0)
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              third_derivatives[i][d0][d1][d2] =
                poly_third_derivatives[i][d0][d1][d2].value(unit_point);

      // fourth derivatives
      if (fourth_derivatives.size() == this->n())
        for (unsigned int d0 = 0; d0 < dim; ++d0)
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              for (unsigned int d3 = 0; d3 < dim; ++d3)
                fourth_derivatives[i][d0][d1][d2][d3] =
                  poly_fourth_derivatives[i][d0][d1][d2][d3].value(unit_point);
    }
}



template <int dim>
double
BarycentricPolynomials<dim>::compute_value(const unsigned int i,
                                           const Point<dim> & p) const
{
  AssertIndexRange(i, this->n());
  return polys[i].value(p);
}



template <int dim>
Tensor<1, dim>
BarycentricPolynomials<dim>::compute_1st_derivative(const unsigned int i,
                                                    const Point<dim> & p) const
{
  Tensor<1, dim> result;
  for (unsigned int d = 0; d < dim; ++d)
    result[d] = poly_grads[i][d].value(p);
  return result;
}



template <int dim>
Tensor<2, dim>
BarycentricPolynomials<dim>::compute_2nd_derivative(const unsigned int i,
                                                    const Point<dim> & p) const
{
  Tensor<2, dim> result;
  for (unsigned int d0 = 0; d0 < dim; ++d0)
    for (unsigned int d1 = 0; d1 < dim; ++d1)
      result[d0][d1] = poly_hessians[i][d0][d1].value(p);

  return result;
}



template <int dim>
Tensor<3, dim>
BarycentricPolynomials<dim>::compute_3rd_derivative(const unsigned int i,
                                                    const Point<dim> & p) const
{
  Tensor<3, dim> result;
  for (unsigned int d0 = 0; d0 < dim; ++d0)
    for (unsigned int d1 = 0; d1 < dim; ++d1)
      for (unsigned int d2 = 0; d2 < dim; ++d2)
        result[d0][d1][d2] = poly_third_derivatives[i][d0][d1][d2].value(p);

  return result;
}



template <int dim>
Tensor<4, dim>
BarycentricPolynomials<dim>::compute_4th_derivative(const unsigned int i,
                                                    const Point<dim> & p) const
{
  Tensor<4, dim> result;
  for (unsigned int d0 = 0; d0 < dim; ++d0)
    for (unsigned int d1 = 0; d1 < dim; ++d1)
      for (unsigned int d2 = 0; d2 < dim; ++d2)
        for (unsigned int d3 = 0; d3 < dim; ++d3)
          result[d0][d1][d2][d3] =
            poly_fourth_derivatives[i][d0][d1][d2][d3].value(p);

  return result;
}



template <int dim>
Tensor<1, dim>
BarycentricPolynomials<dim>::compute_grad(const unsigned int i,
                                          const Point<dim> & p) const
{
  return compute_1st_derivative(i, p);
}



template <int dim>
Tensor<2, dim>
BarycentricPolynomials<dim>::compute_grad_grad(const unsigned int i,
                                               const Point<dim> & p) const
{
  return compute_2nd_derivative(i, p);
}



template <int dim>
std::unique_ptr<ScalarPolynomialsBase<dim>>
BarycentricPolynomials<dim>::clone() const
{
  return std::make_unique<BarycentricPolynomials<dim>>(*this);
}



template <int dim>
std::string
BarycentricPolynomials<dim>::name() const
{
  return "BarycentricPolynomials<" + std::to_string(dim) + ">";
}



template <int dim>
std::size_t
BarycentricPolynomials<dim>::memory_consumption() const
{
  std::size_t poly_memory = 0;
  for (const auto &poly : polys)
    poly_memory += poly.memory_consumption();
  return ScalarPolynomialsBase<dim>::memory_consumption() + poly_memory +
         poly_grads.memory_consumption() + poly_hessians.memory_consumption() +
         poly_third_derivatives.memory_consumption() +
         poly_fourth_derivatives.memory_consumption();
}

template class BarycentricPolynomials<1>;
template class BarycentricPolynomials<2>;
template class BarycentricPolynomials<3>;

DEAL_II_NAMESPACE_CLOSE
