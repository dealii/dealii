// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2024 by the deal.II authors
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
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/quadrature_lib.h>

#include <iomanip>
#include <iostream>
#include <memory>

DEAL_II_NAMESPACE_OPEN


template <int dim>
PolynomialsBDM<dim>::PolynomialsBDM(const unsigned int k)
  : TensorPolynomialsBase<dim>(k + 1, n_polynomials(k))
  , polynomial_space(Polynomials::Legendre::generate_complete_basis(k))
  , monomials((dim == 2) ? (1) : (k + 2))
  , p_values(polynomial_space.n())
  , p_grads(polynomial_space.n())
  , p_grad_grads(polynomial_space.n())
{
  switch (dim)
    {
      case 2:
        monomials[0] = Polynomials::Monomial<double>(k + 1);
        break;
      case 3:
        for (unsigned int i = 0; i < monomials.size(); ++i)
          monomials[i] = Polynomials::Monomial<double>(i);
        break;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim>
void
PolynomialsBDM<dim>::evaluate(
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

  // third and fourth derivatives not implemented
  (void)third_derivatives;
  Assert(third_derivatives.empty(), ExcNotImplemented());
  (void)fourth_derivatives;
  Assert(fourth_derivatives.empty(), ExcNotImplemented());

  const unsigned int n_sub = polynomial_space.n();

  // guard access to the scratch arrays in the following block using a
  // mutex to make sure they are not used by multiple threads at once
  {
    std::lock_guard<std::mutex> lock(mutex);

    p_values.resize((values.empty()) ? 0 : n_sub);
    p_grads.resize((grads.empty()) ? 0 : n_sub);
    p_grad_grads.resize((grad_grads.empty()) ? 0 : n_sub);

    // Compute values of complete space and insert into tensors.  Result
    // will have first all polynomials in the x-component, then y and z.
    polynomial_space.evaluate(unit_point,
                              p_values,
                              p_grads,
                              p_grad_grads,
                              p_third_derivatives,
                              p_fourth_derivatives);

    std::fill(values.begin(), values.end(), Tensor<1, dim>());
    for (unsigned int i = 0; i < p_values.size(); ++i)
      for (unsigned int j = 0; j < dim; ++j)
        values[i + j * n_sub][j] = p_values[i];

    std::fill(grads.begin(), grads.end(), Tensor<2, dim>());
    for (unsigned int i = 0; i < p_grads.size(); ++i)
      for (unsigned int j = 0; j < dim; ++j)
        grads[i + j * n_sub][j] = p_grads[i];

    std::fill(grad_grads.begin(), grad_grads.end(), Tensor<3, dim>());
    for (unsigned int i = 0; i < p_grad_grads.size(); ++i)
      for (unsigned int j = 0; j < dim; ++j)
        grad_grads[i + j * n_sub][j] = p_grad_grads[i];
  }

  // This is the first polynomial not covered by the P_k subspace
  unsigned int start = dim * n_sub;

  // Store values of auxiliary polynomials and their three derivatives
  std::vector<std::vector<double>> monovali(dim, std::vector<double>(4));
  std::vector<std::vector<double>> monovalk(dim, std::vector<double>(4));

  if (dim == 1)
    {
      // Despite the fact that we are instantiating this class for 1, 2 and
      // 3 space dimensions we only support dimension 2 and 3.
      Assert(false,
             dealii::ExcMessage("PolynomialsBDF::evaluate is only "
                                "available for dim == 2, or dim == 3"));
    }
  else if (dim == 2)
    {
      for (unsigned int d = 0; d < dim; ++d)
        monomials[0].value(unit_point[d], monovali[d]);
      if (values.size() != 0)
        {
          values[start][0]     = monovali[0][0];
          values[start][1]     = -unit_point[1] * monovali[0][1];
          values[start + 1][0] = unit_point[0] * monovali[1][1];
          values[start + 1][1] = -monovali[1][0];
        }
      if (grads.size() != 0)
        {
          grads[start][0][0]     = monovali[0][1];
          grads[start][0][1]     = 0.;
          grads[start][1][0]     = -unit_point[1] * monovali[0][2];
          grads[start][1][1]     = -monovali[0][1];
          grads[start + 1][0][0] = monovali[1][1];
          grads[start + 1][0][1] = unit_point[0] * monovali[1][2];
          grads[start + 1][1][0] = 0.;
          grads[start + 1][1][1] = -monovali[1][1];
        }
      if (grad_grads.size() != 0)
        {
          grad_grads[start][0][0][0]     = monovali[0][2];
          grad_grads[start][0][0][1]     = 0.;
          grad_grads[start][0][1][0]     = 0.;
          grad_grads[start][0][1][1]     = 0.;
          grad_grads[start][1][0][0]     = -unit_point[1] * monovali[0][3];
          grad_grads[start][1][0][1]     = -monovali[0][2];
          grad_grads[start][1][1][0]     = -monovali[0][2];
          grad_grads[start][1][1][1]     = 0.;
          grad_grads[start + 1][0][0][0] = 0;
          grad_grads[start + 1][0][0][1] = monovali[1][2];
          grad_grads[start + 1][0][1][0] = monovali[1][2];
          grad_grads[start + 1][0][1][1] = unit_point[0] * monovali[1][3];
          grad_grads[start + 1][1][0][0] = 0.;
          grad_grads[start + 1][1][0][1] = 0.;
          grad_grads[start + 1][1][1][0] = 0.;
          grad_grads[start + 1][1][1][1] = -monovali[1][2];
        }
    }
  else if (dim == 3)
    {
      // The number of curls in each component. Note that the table in
      // BrezziFortin91 has a typo, but the text has the right basis

      // Note that the next basis function is always obtained from the
      // previous by cyclic rotation of the coordinates
      const unsigned int n_curls = monomials.size() - 1;
      for (unsigned int i = 0; i < n_curls; ++i, start += dim)
        {
          for (unsigned int d = 0; d < dim; ++d)
            {
              // p(t) = t^(i+1)
              monomials[i + 1].value(unit_point[d], monovali[d]);
              // q(t) = t^(k-i)
              monomials[this->degree() - 1 - i].value(unit_point[d],
                                                      monovalk[d]);
            }

          if (values.size() != 0)
            {
              // x p'(y) q(z)
              values[start][0] =
                unit_point[0] * monovali[1][1] * monovalk[2][0];
              // - p(y) q(z)
              values[start][1] = -monovali[1][0] * monovalk[2][0];
              values[start][2] = 0.;

              // y p'(z) q(x)
              values[start + 1][1] =
                unit_point[1] * monovali[2][1] * monovalk[0][0];
              // - p(z) q(x)
              values[start + 1][2] = -monovali[2][0] * monovalk[0][0];
              values[start + 1][0] = 0.;

              // z p'(x) q(y)
              values[start + 2][2] =
                unit_point[2] * monovali[0][1] * monovalk[1][0];
              // -p(x) q(y)
              values[start + 2][0] = -monovali[0][0] * monovalk[1][0];
              values[start + 2][1] = 0.;
            }

          if (grads.size() != 0)
            {
              grads[start][0][0] = monovali[1][1] * monovalk[2][0];
              grads[start][0][1] =
                unit_point[0] * monovali[1][2] * monovalk[2][0];
              grads[start][0][2] =
                unit_point[0] * monovali[1][1] * monovalk[2][1];
              grads[start][1][0] = 0.;
              grads[start][1][1] = -monovali[1][1] * monovalk[2][0];
              grads[start][1][2] = -monovali[1][0] * monovalk[2][1];
              grads[start][2][0] = 0.;
              grads[start][2][1] = 0.;
              grads[start][2][2] = 0.;

              grads[start + 1][1][1] = monovali[2][1] * monovalk[0][0];
              grads[start + 1][1][2] =
                unit_point[1] * monovali[2][2] * monovalk[0][0];
              grads[start + 1][1][0] =
                unit_point[1] * monovali[2][1] * monovalk[0][1];
              grads[start + 1][2][1] = 0.;
              grads[start + 1][2][2] = -monovali[2][1] * monovalk[0][0];
              grads[start + 1][2][0] = -monovali[2][0] * monovalk[0][1];
              grads[start + 1][0][1] = 0.;
              grads[start + 1][0][2] = 0.;
              grads[start + 1][0][0] = 0.;

              grads[start + 2][2][2] = monovali[0][1] * monovalk[1][0];
              grads[start + 2][2][0] =
                unit_point[2] * monovali[0][2] * monovalk[1][0];
              grads[start + 2][2][1] =
                unit_point[2] * monovali[0][1] * monovalk[1][1];
              grads[start + 2][0][2] = 0.;
              grads[start + 2][0][0] = -monovali[0][1] * monovalk[1][0];
              grads[start + 2][0][1] = -monovali[0][0] * monovalk[1][1];
              grads[start + 2][1][2] = 0.;
              grads[start + 2][1][0] = 0.;
              grads[start + 2][1][1] = 0.;
            }

          if (grad_grads.size() != 0)
            {
              grad_grads[start][0][0][0] = 0.;
              grad_grads[start][0][0][1] = monovali[1][2] * monovalk[2][0];
              grad_grads[start][0][0][2] = monovali[1][1] * monovalk[2][1];
              grad_grads[start][0][1][0] = monovali[1][2] * monovalk[2][0];
              grad_grads[start][0][1][1] =
                unit_point[0] * monovali[1][3] * monovalk[2][0];
              grad_grads[start][0][1][2] =
                unit_point[0] * monovali[1][2] * monovalk[2][1];
              grad_grads[start][0][2][0] = monovali[1][1] * monovalk[2][1];
              grad_grads[start][0][2][1] =
                unit_point[0] * monovali[1][2] * monovalk[2][1];
              grad_grads[start][0][2][2] =
                unit_point[0] * monovali[1][1] * monovalk[2][2];
              grad_grads[start][1][0][0] = 0.;
              grad_grads[start][1][0][1] = 0.;
              grad_grads[start][1][0][2] = 0.;
              grad_grads[start][1][1][0] = 0.;
              grad_grads[start][1][1][1] = -monovali[1][2] * monovalk[2][0];
              grad_grads[start][1][1][2] = -monovali[1][1] * monovalk[2][1];
              grad_grads[start][1][2][0] = 0.;
              grad_grads[start][1][2][1] = -monovali[1][1] * monovalk[2][1];
              grad_grads[start][1][2][2] = -monovali[1][0] * monovalk[2][2];
              grad_grads[start][2][0][0] = 0.;
              grad_grads[start][2][0][1] = 0.;
              grad_grads[start][2][0][2] = 0.;
              grad_grads[start][2][1][0] = 0.;
              grad_grads[start][2][1][1] = 0.;
              grad_grads[start][2][1][2] = 0.;
              grad_grads[start][2][2][0] = 0.;
              grad_grads[start][2][2][1] = 0.;
              grad_grads[start][2][2][2] = 0.;

              grad_grads[start + 1][0][0][0] = 0.;
              grad_grads[start + 1][0][0][1] = 0.;
              grad_grads[start + 1][0][0][2] = 0.;
              grad_grads[start + 1][0][1][0] = 0.;
              grad_grads[start + 1][0][1][1] = 0.;
              grad_grads[start + 1][0][1][2] = 0.;
              grad_grads[start + 1][0][2][0] = 0.;
              grad_grads[start + 1][0][2][1] = 0.;
              grad_grads[start + 1][0][2][2] = 0.;
              grad_grads[start + 1][1][0][0] =
                unit_point[1] * monovali[2][1] * monovalk[0][2];
              grad_grads[start + 1][1][0][1] = monovali[2][1] * monovalk[0][1];
              grad_grads[start + 1][1][0][2] =
                unit_point[1] * monovali[2][2] * monovalk[0][1];
              grad_grads[start + 1][1][1][0] = monovalk[0][1] * monovali[2][1];
              grad_grads[start + 1][1][1][1] = 0.;
              grad_grads[start + 1][1][1][2] = monovalk[0][0] * monovali[2][2];
              grad_grads[start + 1][1][2][0] =
                unit_point[1] * monovalk[0][1] * monovali[2][2];
              grad_grads[start + 1][1][2][1] = monovalk[0][0] * monovali[2][2];
              grad_grads[start + 1][1][2][2] =
                unit_point[1] * monovalk[0][0] * monovali[2][3];
              grad_grads[start + 1][2][0][0] = -monovalk[0][2] * monovali[2][0];
              grad_grads[start + 1][2][0][1] = 0.;
              grad_grads[start + 1][2][0][2] = -monovalk[0][1] * monovali[2][1];
              grad_grads[start + 1][2][1][0] = 0.;
              grad_grads[start + 1][2][1][1] = 0.;
              grad_grads[start + 1][2][1][2] = 0.;
              grad_grads[start + 1][2][2][0] = -monovalk[0][1] * monovali[2][1];
              grad_grads[start + 1][2][2][1] = 0.;
              grad_grads[start + 1][2][2][2] = -monovalk[0][0] * monovali[2][2];

              grad_grads[start + 2][0][0][0] = -monovali[0][2] * monovalk[1][0];
              grad_grads[start + 2][0][0][1] = -monovali[0][1] * monovalk[1][1];
              grad_grads[start + 2][0][0][2] = 0.;
              grad_grads[start + 2][0][1][0] = -monovali[0][1] * monovalk[1][1];
              grad_grads[start + 2][0][1][1] = -monovali[0][0] * monovalk[1][2];
              grad_grads[start + 2][0][1][2] = 0.;
              grad_grads[start + 2][0][2][0] = 0.;
              grad_grads[start + 2][0][2][1] = 0.;
              grad_grads[start + 2][0][2][2] = 0.;
              grad_grads[start + 2][1][0][0] = 0.;
              grad_grads[start + 2][1][0][1] = 0.;
              grad_grads[start + 2][1][0][2] = 0.;
              grad_grads[start + 2][1][1][0] = 0.;
              grad_grads[start + 2][1][1][1] = 0.;
              grad_grads[start + 2][1][1][2] = 0.;
              grad_grads[start + 2][1][2][0] = 0.;
              grad_grads[start + 2][1][2][1] = 0.;
              grad_grads[start + 2][1][2][2] = 0.;
              grad_grads[start + 2][2][0][0] =
                unit_point[2] * monovali[0][3] * monovalk[1][0];
              grad_grads[start + 2][2][0][1] =
                unit_point[2] * monovali[0][2] * monovalk[1][1];
              grad_grads[start + 2][2][0][2] = monovali[0][2] * monovalk[1][0];
              grad_grads[start + 2][2][1][0] =
                unit_point[2] * monovali[0][2] * monovalk[1][1];
              grad_grads[start + 2][2][1][1] =
                unit_point[2] * monovali[0][1] * monovalk[1][2];
              grad_grads[start + 2][2][1][2] = monovali[0][1] * monovalk[1][1];
              grad_grads[start + 2][2][2][0] = monovali[0][2] * monovalk[1][0];
              grad_grads[start + 2][2][2][1] = monovali[0][1] * monovalk[1][1];
              grad_grads[start + 2][2][2][2] = 0.;
            }
        }
      Assert(start == this->n(), ExcInternalError());
    }
}


template <int dim>
unsigned int
PolynomialsBDM<dim>::n_polynomials(const unsigned int k)
{
  if (dim == 1)
    return k + 1;
  if (dim == 2)
    return (k + 1) * (k + 2) + 2;
  if (dim == 3)
    return ((k + 1) * (k + 2) * (k + 3)) / 2 + 3 * (k + 1);
  DEAL_II_NOT_IMPLEMENTED();
  return 0;
}


template <int dim>
std::unique_ptr<TensorPolynomialsBase<dim>>
PolynomialsBDM<dim>::clone() const
{
  return std::make_unique<PolynomialsBDM<dim>>(*this);
}


template class PolynomialsBDM<1>;
template class PolynomialsBDM<2>;
template class PolynomialsBDM<3>;


DEAL_II_NAMESPACE_CLOSE
