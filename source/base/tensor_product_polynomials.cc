// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2019 by the deal.II authors
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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/base/std_cxx14/memory.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <boost/container/small_vector.hpp>

#include <array>

DEAL_II_NAMESPACE_OPEN



/* ------------------- TensorProductPolynomials -------------- */


namespace internal
{
  namespace
  {
    template <int dim>
    inline void
    compute_tensor_index(const unsigned int,
                         const unsigned int,
                         const unsigned int,
                         unsigned int (&)[dim])
    {
      Assert(false, ExcNotImplemented());
    }

    inline void
    compute_tensor_index(const unsigned int n,
                         const unsigned int,
                         const unsigned int,
                         unsigned int (&indices)[1])
    {
      indices[0] = n;
    }

    inline void
    compute_tensor_index(const unsigned int n,
                         const unsigned int n_pols_0,
                         const unsigned int,
                         unsigned int (&indices)[2])
    {
      indices[0] = n % n_pols_0;
      indices[1] = n / n_pols_0;
    }

    inline void
    compute_tensor_index(const unsigned int n,
                         const unsigned int n_pols_0,
                         const unsigned int n_pols_1,
                         unsigned int (&indices)[3])
    {
      indices[0] = n % n_pols_0;
      indices[1] = (n / n_pols_0) % n_pols_1;
      indices[2] = n / (n_pols_0 * n_pols_1);
    }
  } // namespace
} // namespace internal



template <int dim, typename PolynomialType>
inline void
TensorProductPolynomials<dim, PolynomialType>::compute_index(
  const unsigned int i,
  unsigned int (&indices)[(dim > 0 ? dim : 1)]) const
{
  Assert(i < Utilities::fixed_power<dim>(polynomials.size()),
         ExcInternalError());
  internal::compute_tensor_index(index_map[i],
                                 polynomials.size(),
                                 polynomials.size(),
                                 indices);
}



template <int dim, typename PolynomialType>
void
TensorProductPolynomials<dim, PolynomialType>::output_indices(
  std::ostream &out) const
{
  unsigned int ix[(dim > 0) ? dim : 1];
  for (unsigned int i = 0; i < this->n(); ++i)
    {
      compute_index(i, ix);
      out << i << "\t";
      for (unsigned int d = 0; d < dim; ++d)
        out << ix[d] << " ";
      out << std::endl;
    }
}



template <int dim, typename PolynomialType>
void
TensorProductPolynomials<dim, PolynomialType>::set_numbering(
  const std::vector<unsigned int> &renumber)
{
  Assert(renumber.size() == index_map.size(),
         ExcDimensionMismatch(renumber.size(), index_map.size()));

  index_map = renumber;
  for (unsigned int i = 0; i < index_map.size(); ++i)
    index_map_inverse[index_map[i]] = i;
}



template <>
double
TensorProductPolynomials<0, Polynomials::Polynomial<double>>::compute_value(
  const unsigned int,
  const Point<0> &) const
{
  Assert(false, ExcNotImplemented());
  return 0;
}



template <int dim, typename PolynomialType>
double
TensorProductPolynomials<dim, PolynomialType>::compute_value(
  const unsigned int i,
  const Point<dim> & p) const
{
  Assert(dim > 0, ExcNotImplemented());

  unsigned int indices[dim];
  compute_index(i, indices);

  double value = 1.;
  for (unsigned int d = 0; d < dim; ++d)
    value *= polynomials[indices[d]].value(p(d));

  return value;
}



template <int dim, typename PolynomialType>
Tensor<1, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_grad(
  const unsigned int i,
  const Point<dim> & p) const
{
  unsigned int indices[dim];
  compute_index(i, indices);

  // compute values and
  // uni-directional derivatives at
  // the given point in each
  // co-ordinate direction
  double v[dim][2];
  {
    std::vector<double> tmp(2);
    for (unsigned int d = 0; d < dim; ++d)
      {
        polynomials[indices[d]].value(p(d), tmp);
        v[d][0] = tmp[0];
        v[d][1] = tmp[1];
      }
  }

  Tensor<1, dim> grad;
  for (unsigned int d = 0; d < dim; ++d)
    {
      grad[d] = 1.;
      for (unsigned int x = 0; x < dim; ++x)
        grad[d] *= v[x][d == x];
    }

  return grad;
}



template <>
Tensor<1, 0>
TensorProductPolynomials<0, Polynomials::Polynomial<double>>::compute_grad(
  const unsigned int,
  const Point<0> &) const
{
  return Tensor<1, 0>();
}



template <int dim, typename PolynomialType>
Tensor<2, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_grad_grad(
  const unsigned int i,
  const Point<dim> & p) const
{
  unsigned int indices[(dim > 0) ? dim : 1];
  compute_index(i, indices);

  double v[dim][3];
  {
    std::vector<double> tmp(3);
    for (unsigned int d = 0; d < dim; ++d)
      {
        polynomials[indices[d]].value(p(d), tmp);
        v[d][0] = tmp[0];
        v[d][1] = tmp[1];
        v[d][2] = tmp[2];
      }
  }

  Tensor<2, dim> grad_grad;
  for (unsigned int d1 = 0; d1 < dim; ++d1)
    for (unsigned int d2 = 0; d2 < dim; ++d2)
      {
        grad_grad[d1][d2] = 1.;
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
            grad_grad[d1][d2] *= v[x][derivative];
          }
      }

  return grad_grad;
}



template <>
Tensor<2, 0>
TensorProductPolynomials<0, Polynomials::Polynomial<double>>::compute_grad_grad(
  const unsigned int,
  const Point<0> &) const
{
  return Tensor<2, 0>();
}



template <int dim, typename PolynomialType>
void
TensorProductPolynomials<dim, PolynomialType>::evaluate(
  const Point<dim> &           p,
  std::vector<double> &        values,
  std::vector<Tensor<1, dim>> &grads,
  std::vector<Tensor<2, dim>> &grad_grads,
  std::vector<Tensor<3, dim>> &third_derivatives,
  std::vector<Tensor<4, dim>> &fourth_derivatives) const
{
  Assert(dim <= 3, ExcNotImplemented());
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

  const bool update_values          = (values.size() == this->n()),
             update_grads           = (grads.size() == this->n()),
             update_grad_grads      = (grad_grads.size() == this->n()),
             update_3rd_derivatives = (third_derivatives.size() == this->n()),
             update_4th_derivatives = (fourth_derivatives.size() == this->n());

  // check how many values/derivatives we have to compute
  unsigned int n_values_and_derivatives = 0;
  if (update_values)
    n_values_and_derivatives = 1;
  if (update_grads)
    n_values_and_derivatives = 2;
  if (update_grad_grads)
    n_values_and_derivatives = 3;
  if (update_3rd_derivatives)
    n_values_and_derivatives = 4;
  if (update_4th_derivatives)
    n_values_and_derivatives = 5;

  // Compute the values (and derivatives, if necessary) of all 1D polynomials
  // at this evaluation point. We need to compute dim*n_polynomials
  // evaluations, involving an evaluation of each polynomial for each
  // coordinate direction. Once we have those values, we perform the
  // multiplications for the tensor product in the arbitrary dimension.
  const unsigned int n_polynomials = polynomials.size();
  boost::container::small_vector<std::array<std::array<double, 5>, dim>, 20>
    values_1d(n_polynomials);
  if (n_values_and_derivatives == 1)
    for (unsigned int i = 0; i < n_polynomials; ++i)
      for (unsigned int d = 0; d < dim; ++d)
        values_1d[i][d][0] = polynomials[i].value(p(d));
  else
    for (unsigned int i = 0; i < n_polynomials; ++i)
      for (unsigned d = 0; d < dim; ++d)
        polynomials[i].value(p(d),
                             n_values_and_derivatives,
                             values_1d[i][d].data());

  unsigned int indices[3];
  unsigned int ind = 0;
  for (indices[2] = 0; indices[2] < (dim > 2 ? n_polynomials : 1); ++indices[2])
    for (indices[1] = 0; indices[1] < (dim > 1 ? n_polynomials : 1);
         ++indices[1])
      if (n_values_and_derivatives == 1)
        for (indices[0] = 0; indices[0] < n_polynomials; ++indices[0], ++ind)
          {
            double value = values_1d[indices[0]][0][0];
            for (unsigned int d = 1; d < dim; ++d)
              value *= values_1d[indices[d]][d][0];
            values[index_map_inverse[ind]] = value;
          }
      else
        for (indices[0] = 0; indices[0] < n_polynomials; ++indices[0], ++ind)
          {
            unsigned int i = index_map_inverse[ind];

            if (update_values)
              {
                double value = values_1d[indices[0]][0][0];
                for (unsigned int x = 1; x < dim; ++x)
                  value *= values_1d[indices[x]][x][0];
                values[i] = value;
              }

            if (update_grads)
              for (unsigned int d = 0; d < dim; ++d)
                {
                  double grad = 1.;
                  for (unsigned int x = 0; x < dim; ++x)
                    grad *= values_1d[indices[x]][x][(d == x) ? 1 : 0];
                  grads[i][d] = grad;
                }

            if (update_grad_grads)
              for (unsigned int d1 = 0; d1 < dim; ++d1)
                for (unsigned int d2 = 0; d2 < dim; ++d2)
                  {
                    double der2 = 1.;
                    for (unsigned int x = 0; x < dim; ++x)
                      {
                        unsigned int derivative = 0;
                        if (d1 == x)
                          ++derivative;
                        if (d2 == x)
                          ++derivative;

                        der2 *= values_1d[indices[x]][x][derivative];
                      }
                    grad_grads[i][d1][d2] = der2;
                  }

            if (update_3rd_derivatives)
              for (unsigned int d1 = 0; d1 < dim; ++d1)
                for (unsigned int d2 = 0; d2 < dim; ++d2)
                  for (unsigned int d3 = 0; d3 < dim; ++d3)
                    {
                      double der3 = 1.;
                      for (unsigned int x = 0; x < dim; ++x)
                        {
                          unsigned int derivative = 0;
                          if (d1 == x)
                            ++derivative;
                          if (d2 == x)
                            ++derivative;
                          if (d3 == x)
                            ++derivative;

                          der3 *= values_1d[indices[x]][x][derivative];
                        }
                      third_derivatives[i][d1][d2][d3] = der3;
                    }

            if (update_4th_derivatives)
              for (unsigned int d1 = 0; d1 < dim; ++d1)
                for (unsigned int d2 = 0; d2 < dim; ++d2)
                  for (unsigned int d3 = 0; d3 < dim; ++d3)
                    for (unsigned int d4 = 0; d4 < dim; ++d4)
                      {
                        double der4 = 1.;
                        for (unsigned int x = 0; x < dim; ++x)
                          {
                            unsigned int derivative = 0;
                            if (d1 == x)
                              ++derivative;
                            if (d2 == x)
                              ++derivative;
                            if (d3 == x)
                              ++derivative;
                            if (d4 == x)
                              ++derivative;

                            der4 *= values_1d[indices[x]][x][derivative];
                          }
                        fourth_derivatives[i][d1][d2][d3][d4] = der4;
                      }
          }
}



template <int dim, typename PolynomialType>
std::unique_ptr<ScalarPolynomialsBase<dim>>
TensorProductPolynomials<dim, PolynomialType>::clone() const
{
  return std_cxx14::make_unique<TensorProductPolynomials<dim, PolynomialType>>(
    *this);
}



/* ------------------- AnisotropicPolynomials -------------- */


template <int dim>
AnisotropicPolynomials<dim>::AnisotropicPolynomials(
  const std::vector<std::vector<Polynomials::Polynomial<double>>> &pols)
  : ScalarPolynomialsBase<dim>(1, get_n_tensor_pols(pols))
  , polynomials(pols)
{
  Assert(pols.size() == dim, ExcDimensionMismatch(pols.size(), dim));
  for (unsigned int d = 0; d < dim; ++d)
    Assert(pols[d].size() > 0,
           ExcMessage("The number of polynomials must be larger than zero "
                      "for all coordinate directions."));
}



template <int dim>
void
AnisotropicPolynomials<dim>::compute_index(
  const unsigned int i,
  unsigned int (&indices)[(dim > 0 ? dim : 1)]) const
{
#ifdef DEBUG
  unsigned int n_poly = 1;
  for (unsigned int d = 0; d < dim; ++d)
    n_poly *= polynomials[d].size();
  Assert(i < n_poly, ExcInternalError());
#endif

  if (dim == 0)
    {
    }
  else if (dim == 1)
    internal::compute_tensor_index(i,
                                   polynomials[0].size(),
                                   0 /*not used*/,
                                   indices);
  else
    internal::compute_tensor_index(i,
                                   polynomials[0].size(),
                                   polynomials[1].size(),
                                   indices);
}



template <int dim>
double
AnisotropicPolynomials<dim>::compute_value(const unsigned int i,
                                           const Point<dim> & p) const
{
  unsigned int indices[(dim > 0) ? dim : 1];
  compute_index(i, indices);

  double value = 1.;
  for (unsigned int d = 0; d < dim; ++d)
    value *= polynomials[d][indices[d]].value(p(d));

  return value;
}


template <int dim>
Tensor<1, dim>
AnisotropicPolynomials<dim>::compute_grad(const unsigned int i,
                                          const Point<dim> & p) const
{
  unsigned int indices[(dim > 0) ? dim : 1];
  compute_index(i, indices);

  // compute values and
  // uni-directional derivatives at
  // the given point in each
  // co-ordinate direction
  std::vector<std::vector<double>> v(dim, std::vector<double>(2));
  for (unsigned int d = 0; d < dim; ++d)
    polynomials[d][indices[d]].value(p(d), v[d]);

  Tensor<1, dim> grad;
  for (unsigned int d = 0; d < dim; ++d)
    {
      grad[d] = 1.;
      for (unsigned int x = 0; x < dim; ++x)
        grad[d] *= v[x][d == x];
    }

  return grad;
}


template <int dim>
Tensor<2, dim>
AnisotropicPolynomials<dim>::compute_grad_grad(const unsigned int i,
                                               const Point<dim> & p) const
{
  unsigned int indices[(dim > 0) ? dim : 1];
  compute_index(i, indices);

  std::vector<std::vector<double>> v(dim, std::vector<double>(3));
  for (unsigned int d = 0; d < dim; ++d)
    polynomials[d][indices[d]].value(p(d), v[d]);

  Tensor<2, dim> grad_grad;
  for (unsigned int d1 = 0; d1 < dim; ++d1)
    for (unsigned int d2 = 0; d2 < dim; ++d2)
      {
        grad_grad[d1][d2] = 1.;
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
            grad_grad[d1][d2] *= v[x][derivative];
          }
      }

  return grad_grad;
}



template <int dim>
void
AnisotropicPolynomials<dim>::evaluate(
  const Point<dim> &           p,
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

  const bool update_values          = (values.size() == this->n()),
             update_grads           = (grads.size() == this->n()),
             update_grad_grads      = (grad_grads.size() == this->n()),
             update_3rd_derivatives = (third_derivatives.size() == this->n()),
             update_4th_derivatives = (fourth_derivatives.size() == this->n());

  // check how many
  // values/derivatives we have to
  // compute
  unsigned int n_values_and_derivatives = 0;
  if (update_values)
    n_values_and_derivatives = 1;
  if (update_grads)
    n_values_and_derivatives = 2;
  if (update_grad_grads)
    n_values_and_derivatives = 3;
  if (update_3rd_derivatives)
    n_values_and_derivatives = 4;
  if (update_4th_derivatives)
    n_values_and_derivatives = 5;

  // compute the values (and
  // derivatives, if necessary) of
  // all polynomials at this
  // evaluation point
  std::vector<std::vector<std::vector<double>>> v(dim);
  for (unsigned int d = 0; d < dim; ++d)
    {
      v[d].resize(polynomials[d].size());
      for (unsigned int i = 0; i < polynomials[d].size(); ++i)
        {
          v[d][i].resize(n_values_and_derivatives, 0.);
          polynomials[d][i].value(p(d), v[d][i]);
        }
    }

  for (unsigned int i = 0; i < this->n(); ++i)
    {
      // first get the
      // one-dimensional indices of
      // this particular tensor
      // product polynomial
      unsigned int indices[(dim > 0) ? dim : 1];
      compute_index(i, indices);

      if (update_values)
        {
          values[i] = 1;
          for (unsigned int x = 0; x < dim; ++x)
            values[i] *= v[x][indices[x]][0];
        }

      if (update_grads)
        for (unsigned int d = 0; d < dim; ++d)
          {
            grads[i][d] = 1.;
            for (unsigned int x = 0; x < dim; ++x)
              grads[i][d] *= v[x][indices[x]][d == x ? 1 : 0];
          }

      if (update_grad_grads)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            {
              grad_grads[i][d1][d2] = 1.;
              for (unsigned int x = 0; x < dim; ++x)
                {
                  unsigned int derivative = 0;
                  if (d1 == x)
                    ++derivative;
                  if (d2 == x)
                    ++derivative;

                  grad_grads[i][d1][d2] *= v[x][indices[x]][derivative];
                }
            }

      if (update_3rd_derivatives)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            for (unsigned int d3 = 0; d3 < dim; ++d3)
              {
                third_derivatives[i][d1][d2][d3] = 1.;
                for (unsigned int x = 0; x < dim; ++x)
                  {
                    unsigned int derivative = 0;
                    if (d1 == x)
                      ++derivative;
                    if (d2 == x)
                      ++derivative;
                    if (d3 == x)
                      ++derivative;

                    third_derivatives[i][d1][d2][d3] *=
                      v[x][indices[x]][derivative];
                  }
              }

      if (update_4th_derivatives)
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          for (unsigned int d2 = 0; d2 < dim; ++d2)
            for (unsigned int d3 = 0; d3 < dim; ++d3)
              for (unsigned int d4 = 0; d4 < dim; ++d4)
                {
                  fourth_derivatives[i][d1][d2][d3][d4] = 1.;
                  for (unsigned int x = 0; x < dim; ++x)
                    {
                      unsigned int derivative = 0;
                      if (d1 == x)
                        ++derivative;
                      if (d2 == x)
                        ++derivative;
                      if (d3 == x)
                        ++derivative;
                      if (d4 == x)
                        ++derivative;

                      fourth_derivatives[i][d1][d2][d3][d4] *=
                        v[x][indices[x]][derivative];
                    }
                }
    }
}


template <int dim>
unsigned int
AnisotropicPolynomials<dim>::get_n_tensor_pols(
  const std::vector<std::vector<Polynomials::Polynomial<double>>> &pols)
{
  unsigned int y = 1;
  for (unsigned int d = 0; d < dim; ++d)
    y *= pols[d].size();
  return y;
}


template <int dim>
std::unique_ptr<ScalarPolynomialsBase<dim>>
AnisotropicPolynomials<dim>::clone() const
{
  return std_cxx14::make_unique<AnisotropicPolynomials<dim>>(*this);
}



/* ------------------- explicit instantiations -------------- */
template class TensorProductPolynomials<0, Polynomials::Polynomial<double>>;
template class TensorProductPolynomials<1, Polynomials::Polynomial<double>>;
template class TensorProductPolynomials<2, Polynomials::Polynomial<double>>;
template class TensorProductPolynomials<3, Polynomials::Polynomial<double>>;

template class TensorProductPolynomials<
  1,
  Polynomials::PiecewisePolynomial<double>>;
template class TensorProductPolynomials<
  2,
  Polynomials::PiecewisePolynomial<double>>;
template class TensorProductPolynomials<
  3,
  Polynomials::PiecewisePolynomial<double>>;

template class AnisotropicPolynomials<0>;
template class AnisotropicPolynomials<1>;
template class AnisotropicPolynomials<2>;
template class AnisotropicPolynomials<3>;

DEAL_II_NAMESPACE_CLOSE
