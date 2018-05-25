// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/base/exceptions.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor_product_polynomials_bubbles.h>

DEAL_II_NAMESPACE_OPEN



/* ------------------- TensorProductPolynomialsBubbles -------------- */


template <int dim>
double
TensorProductPolynomialsBubbles<dim>::compute_value(const unsigned int i,
                                                    const Point<dim> & p) const
{
  const unsigned int q_degree      = this->polynomials.size() - 1;
  const unsigned int max_q_indices = this->n_tensor_pols;
  const unsigned int n_bubbles     = ((q_degree <= 1) ? 1 : dim);
  (void)n_bubbles;
  Assert(i < max_q_indices + n_bubbles, ExcInternalError());

  // treat the regular basis functions
  if (i < max_q_indices)
    return this->TensorProductPolynomials<dim>::compute_value(i, p);

  const unsigned int comp = i - this->n_tensor_pols;

  // compute \prod_{i=1}^d 4*(1-x_i^2)(p)
  double value = 1.;
  for (unsigned int j = 0; j < dim; ++j)
    value *= 4 * p(j) * (1 - p(j));
  // and multiply with (2x_i-1)^{r-1}
  for (unsigned int i = 0; i < q_degree - 1; ++i)
    value *= 2 * p(comp) - 1;
  return value;
}



template <>
double
TensorProductPolynomialsBubbles<0>::compute_value(const unsigned int,
                                                  const Point<0> &) const
{
  Assert(false, ExcNotImplemented());
  return 0.;
}


template <int dim>
Tensor<1, dim>
TensorProductPolynomialsBubbles<dim>::compute_grad(const unsigned int i,
                                                   const Point<dim> & p) const
{
  const unsigned int q_degree      = this->polynomials.size() - 1;
  const unsigned int max_q_indices = this->n_tensor_pols;
  const unsigned int n_bubbles     = ((q_degree <= 1) ? 1 : dim);
  (void)n_bubbles;
  Assert(i < max_q_indices + n_bubbles, ExcInternalError());

  // treat the regular basis functions
  if (i < max_q_indices)
    return this->TensorProductPolynomials<dim>::compute_grad(i, p);

  const unsigned int comp = i - this->n_tensor_pols;
  Tensor<1, dim>     grad;

  for (unsigned int d = 0; d < dim; ++d)
    {
      grad[d] = 1.;
      // compute grad(4*\prod_{i=1}^d (x_i(1-x_i)))(p)
      for (unsigned j = 0; j < dim; ++j)
        grad[d] *= (d == j ? 4 * (1 - 2 * p(j)) : 4 * p(j) * (1 - p(j)));
      // and multiply with (2*x_i-1)^{r-1}
      for (unsigned int i = 0; i < q_degree - 1; ++i)
        grad[d] *= 2 * p(comp) - 1;
    }

  if (q_degree >= 2)
    {
      // add \prod_{i=1}^d 4*(x_i(1-x_i))(p)
      double value = 1.;
      for (unsigned int j = 0; j < dim; ++j)
        value *= 4 * p(j) * (1 - p(j));
      // and multiply with grad(2*x_i-1)^{r-1}
      double tmp = value * 2 * (q_degree - 1);
      for (unsigned int i = 0; i < q_degree - 2; ++i)
        tmp *= 2 * p(comp) - 1;
      grad[comp] += tmp;
    }

  return grad;
}



template <int dim>
Tensor<2, dim>
TensorProductPolynomialsBubbles<dim>::compute_grad_grad(
  const unsigned int i,
  const Point<dim> & p) const
{
  const unsigned int q_degree      = this->polynomials.size() - 1;
  const unsigned int max_q_indices = this->n_tensor_pols;
  const unsigned int n_bubbles     = ((q_degree <= 1) ? 1 : dim);
  (void)n_bubbles;
  Assert(i < max_q_indices + n_bubbles, ExcInternalError());

  // treat the regular basis functions
  if (i < max_q_indices)
    return this->TensorProductPolynomials<dim>::compute_grad_grad(i, p);

  const unsigned int comp = i - this->n_tensor_pols;

  double v[dim + 1][3];
  {
    for (unsigned int c = 0; c < dim; ++c)
      {
        v[c][0] = 4 * p(c) * (1 - p(c));
        v[c][1] = 4 * (1 - 2 * p(c));
        v[c][2] = -8;
      }

    double tmp = 1.;
    for (unsigned int i = 0; i < q_degree - 1; ++i)
      tmp *= 2 * p(comp) - 1;
    v[dim][0] = tmp;

    if (q_degree >= 2)
      {
        double tmp = 2 * (q_degree - 1);
        for (unsigned int i = 0; i < q_degree - 2; ++i)
          tmp *= 2 * p(comp) - 1;
        v[dim][1] = tmp;
      }
    else
      v[dim][1] = 0.;

    if (q_degree >= 3)
      {
        double tmp = 4 * (q_degree - 2) * (q_degree - 1);
        for (unsigned int i = 0; i < q_degree - 3; ++i)
          tmp *= 2 * p(comp) - 1;
        v[dim][2] = tmp;
      }
    else
      v[dim][2] = 0.;
  }

  // calculate (\partial_j \partial_k \psi) * monomial
  Tensor<2, dim> grad_grad_1;
  for (unsigned int d1 = 0; d1 < dim; ++d1)
    for (unsigned int d2 = 0; d2 < dim; ++d2)
      {
        grad_grad_1[d1][d2] = v[dim][0];
        for (unsigned int x = 0; x < dim; ++x)
          {
            unsigned int derivative = 0;
            if (d1 == x || d2 == x)
              {
                if (d1 == d2)
                  derivative = 2;
                else
                  derivative = 1;
              }
            grad_grad_1[d1][d2] *= v[x][derivative];
          }
      }

  // calculate (\partial_j  \psi) *(\partial_k monomial)
  // and (\partial_k  \psi) *(\partial_j monomial)
  Tensor<2, dim> grad_grad_2;
  Tensor<2, dim> grad_grad_3;
  for (unsigned int d = 0; d < dim; ++d)
    {
      grad_grad_2[d][comp] = v[dim][1];
      grad_grad_3[comp][d] = v[dim][1];
      for (unsigned int x = 0; x < dim; ++x)
        {
          grad_grad_2[d][comp] *= v[x][d == x];
          grad_grad_3[comp][d] *= v[x][d == x];
        }
    }

  // calculate \psi *(\partial j \partial_k monomial) and sum
  Tensor<2, dim> grad_grad;
  double         psi_value = 1.;
  for (unsigned int x = 0; x < dim; ++x)
    psi_value *= v[x][0];

  for (unsigned int d1 = 0; d1 < dim; ++d1)
    for (unsigned int d2 = 0; d2 < dim; ++d2)
      grad_grad[d1][d2] =
        grad_grad_1[d1][d2] + grad_grad_2[d1][d2] + grad_grad_3[d1][d2];
  grad_grad[comp][comp] += psi_value * v[dim][2];

  return grad_grad;
}

template <int dim>
void
TensorProductPolynomialsBubbles<dim>::compute(
  const Point<dim> &           p,
  std::vector<double> &        values,
  std::vector<Tensor<1, dim>> &grads,
  std::vector<Tensor<2, dim>> &grad_grads,
  std::vector<Tensor<3, dim>> &third_derivatives,
  std::vector<Tensor<4, dim>> &fourth_derivatives) const
{
  const unsigned int q_degree      = this->polynomials.size() - 1;
  const unsigned int max_q_indices = this->n_tensor_pols;
  (void)max_q_indices;
  const unsigned int n_bubbles = ((q_degree <= 1) ? 1 : dim);
  Assert(values.size() == max_q_indices + n_bubbles || values.size() == 0,
         ExcDimensionMismatch2(values.size(), max_q_indices + n_bubbles, 0));
  Assert(grads.size() == max_q_indices + n_bubbles || grads.size() == 0,
         ExcDimensionMismatch2(grads.size(), max_q_indices + n_bubbles, 0));
  Assert(
    grad_grads.size() == max_q_indices + n_bubbles || grad_grads.size() == 0,
    ExcDimensionMismatch2(grad_grads.size(), max_q_indices + n_bubbles, 0));
  Assert(third_derivatives.size() == max_q_indices + n_bubbles ||
           third_derivatives.size() == 0,
         ExcDimensionMismatch2(
           third_derivatives.size(), max_q_indices + n_bubbles, 0));
  Assert(fourth_derivatives.size() == max_q_indices + n_bubbles ||
           fourth_derivatives.size() == 0,
         ExcDimensionMismatch2(
           fourth_derivatives.size(), max_q_indices + n_bubbles, 0));

  bool do_values = false, do_grads = false, do_grad_grads = false;
  bool do_3rd_derivatives = false, do_4th_derivatives = false;
  if (values.empty() == false)
    {
      values.resize(this->n_tensor_pols);
      do_values = true;
    }
  if (grads.empty() == false)
    {
      grads.resize(this->n_tensor_pols);
      do_grads = true;
    }
  if (grad_grads.empty() == false)
    {
      grad_grads.resize(this->n_tensor_pols);
      do_grad_grads = true;
    }
  if (third_derivatives.empty() == false)
    {
      third_derivatives.resize(this->n_tensor_pols);
      do_3rd_derivatives = true;
    }
  if (fourth_derivatives.empty() == false)
    {
      fourth_derivatives.resize(this->n_tensor_pols);
      do_4th_derivatives = true;
    }

  this->TensorProductPolynomials<dim>::compute(
    p, values, grads, grad_grads, third_derivatives, fourth_derivatives);

  for (unsigned int i = this->n_tensor_pols;
       i < this->n_tensor_pols + n_bubbles;
       ++i)
    {
      if (do_values)
        values.push_back(compute_value(i, p));
      if (do_grads)
        grads.push_back(compute_grad(i, p));
      if (do_grad_grads)
        grad_grads.push_back(compute_grad_grad(i, p));
      if (do_3rd_derivatives)
        third_derivatives.push_back(compute_derivative<3>(i, p));
      if (do_4th_derivatives)
        fourth_derivatives.push_back(compute_derivative<4>(i, p));
    }
}


/* ------------------- explicit instantiations -------------- */
template class TensorProductPolynomialsBubbles<1>;
template class TensorProductPolynomialsBubbles<2>;
template class TensorProductPolynomialsBubbles<3>;

DEAL_II_NAMESPACE_CLOSE
