// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2019 by the deal.II authors
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


#include <deal.II/base/polynomials_abf.h>
#include <deal.II/base/quadrature_lib.h>

#include <iomanip>
#include <iostream>
#include <memory>


DEAL_II_NAMESPACE_OPEN



namespace
{
  template <int dim>
  std::vector<std::vector<Polynomials::Polynomial<double>>>
  get_abf_polynomials(const unsigned int k)
  {
    std::vector<std::vector<Polynomials::Polynomial<double>>> pols(dim);
    pols[0] = Polynomials::LagrangeEquidistant::generate_complete_basis(k + 2);

    if (k == 0)
      for (unsigned int d = 1; d < dim; ++d)
        pols[d] = Polynomials::Legendre::generate_complete_basis(0);
    else
      for (unsigned int d = 1; d < dim; ++d)
        pols[d] = Polynomials::LagrangeEquidistant::generate_complete_basis(k);

    return pols;
  }
} // namespace

template <int dim>
PolynomialsABF<dim>::PolynomialsABF(const unsigned int k)
  : TensorPolynomialsBase<dim>(k, n_polynomials(k))
  , polynomial_space(get_abf_polynomials<dim>(k))
{
  // check that the dimensions match. we only store one of the 'dim'
  // anisotropic polynomials that make up the vector-valued space, so
  // multiply by 'dim'
  Assert(dim * polynomial_space.n() == n_polynomials(k), ExcInternalError());
}



template <int dim>
void
PolynomialsABF<dim>::evaluate(
  const Point<dim> &           unit_point,
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

  const unsigned int n_sub = polynomial_space.n();
  // guard access to the scratch
  // arrays in the following block
  // using a mutex to make sure they
  // are not used by multiple threads
  // at once
  std::lock_guard<std::mutex> lock(mutex);

  p_values.resize((values.size() == 0) ? 0 : n_sub);
  p_grads.resize((grads.size() == 0) ? 0 : n_sub);
  p_grad_grads.resize((grad_grads.size() == 0) ? 0 : n_sub);
  p_third_derivatives.resize((third_derivatives.size() == 0) ? 0 : n_sub);
  p_fourth_derivatives.resize((fourth_derivatives.size() == 0) ? 0 : n_sub);

  for (unsigned int d = 0; d < dim; ++d)
    {
      // First we copy the point. The
      // polynomial space for
      // component d consists of
      // polynomials of degree k+1 in
      // x_d and degree k in the
      // other variables. in order to
      // simplify this, we use the
      // same AnisotropicPolynomial
      // space and simply rotate the
      // coordinates through all
      // directions.
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
        values[i + d * n_sub][d] = p_values[i];

      for (unsigned int i = 0; i < p_grads.size(); ++i)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          grads[i + d * n_sub][d][(d1 + d) % dim] = p_grads[i][d1];

      for (unsigned int i = 0; i < p_grad_grads.size(); ++i)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            grad_grads[i + d * n_sub][d][(d1 + d) % dim][(d2 + d) % dim] =
              p_grad_grads[i][d1][d2];

      for (unsigned int i = 0; i < p_third_derivatives.size(); ++i)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            for (unsigned int d3 = 0; d3 < dim; ++d3)
              third_derivatives[i + d * n_sub][d][(d1 + d) % dim]
                               [(d2 + d) % dim][(d3 + d) % dim] =
                                 p_third_derivatives[i][d1][d2][d3];

      for (unsigned int i = 0; i < p_fourth_derivatives.size(); ++i)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            for (unsigned int d3 = 0; d3 < dim; ++d3)
              for (unsigned int d4 = 0; d4 < dim; ++d4)
                fourth_derivatives[i + d * n_sub][d][(d1 + d) % dim]
                                  [(d2 + d) % dim][(d3 + d) % dim]
                                  [(d4 + d) % dim] =
                                    p_fourth_derivatives[i][d1][d2][d3][d4];
    }
}


template <int dim>
unsigned int
PolynomialsABF<dim>::n_polynomials(const unsigned int k)
{
  switch (dim)
    {
      case 1:
        // in 1d, we simply have Q_{k+2}, which has dimension k+3
        return k + 3;

      case 2:
        // the polynomial space is Q_{k+2,k} \times Q_{k,k+2}, which has
        // 2(k+3)(k+1) DoFs
        return 2 * (k + 3) * (k + 1);

      case 3:
        // the polynomial space is Q_{k+2,k,k} \times Q_{k,k+2,k} \times
        // Q_{k,k,k+2}, which has 3(k+3)(k+1)(k+1) DoFs
        return 3 * (k + 3) * (k + 1) * (k + 1);

      default:
        Assert(false, ExcNotImplemented());
    }

  return 0;
}


template <int dim>
std::unique_ptr<TensorPolynomialsBase<dim>>
PolynomialsABF<dim>::clone() const
{
  return std::make_unique<PolynomialsABF<dim>>(*this);
}


template class PolynomialsABF<1>;
template class PolynomialsABF<2>;
template class PolynomialsABF<3>;


DEAL_II_NAMESPACE_CLOSE
