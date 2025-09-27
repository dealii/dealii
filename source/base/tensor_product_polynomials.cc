// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <boost/container/small_vector.hpp>

#include <array>
#include <memory>

DEAL_II_NAMESPACE_OPEN



/* ------------------- TensorProductPolynomials -------------- */


namespace internal
{
  namespace
  {
    template <std::size_t dim>
    inline void
    compute_tensor_index(const unsigned int,
                         const unsigned int,
                         const unsigned int,
                         std::array<unsigned int, dim> &)
    {
      DEAL_II_NOT_IMPLEMENTED();
    }

    inline void
    compute_tensor_index(const unsigned int n,
                         const unsigned int,
                         const unsigned int,
                         std::array<unsigned int, 1> &indices)
    {
      indices[0] = n;
    }

    inline void
    compute_tensor_index(const unsigned int n,
                         const unsigned int n_pols_0,
                         const unsigned int,
                         std::array<unsigned int, 2> &indices)
    {
      indices[0] = n % n_pols_0;
      indices[1] = n / n_pols_0;
    }

    inline void
    compute_tensor_index(const unsigned int           n,
                         const unsigned int           n_pols_0,
                         const unsigned int           n_pols_1,
                         std::array<unsigned int, 3> &indices)
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
  const unsigned int             i,
  std::array<unsigned int, dim> &indices) const
{
  if constexpr (dim == 0)
    {
      (void)i;
      (void)indices;
      DEAL_II_NOT_IMPLEMENTED();
    }
  else
    {
      Assert(i < Utilities::fixed_power<dim>(polynomials.size()),
             ExcInternalError());
      internal::compute_tensor_index(index_map[i],
                                     polynomials.size(),
                                     polynomials.size(),
                                     indices);
    }
}



template <int dim, typename PolynomialType>
void
TensorProductPolynomials<dim, PolynomialType>::output_indices(
  std::ostream &out) const
{
  if constexpr (dim == 0)
    {
      (void)out;
      DEAL_II_NOT_IMPLEMENTED();
    }
  else
    {
      std::array<unsigned int, dim> ix;
      for (unsigned int i = 0; i < this->n(); ++i)
        {
          compute_index(i, ix);
          out << i << "\t";
          for (unsigned int d = 0; d < dim; ++d)
            out << ix[d] << " ";
          out << std::endl;
        }
    }
}



template <int dim>
inline const std::vector<unsigned int> &
AnisotropicPolynomials<dim>::get_numbering() const
{
  return index_map;
}



template <int dim>
inline const std::vector<unsigned int> &
AnisotropicPolynomials<dim>::get_numbering_inverse() const
{
  return index_map_inverse;
}



template <int dim>
void
AnisotropicPolynomials<dim>::set_numbering(
  const std::vector<unsigned int> &renumber)
{
  Assert(renumber.size() == index_map.size(),
         ExcDimensionMismatch(renumber.size(), index_map.size()));

  index_map = renumber;
  for (unsigned int i = 0; i < index_map.size(); ++i)
    index_map_inverse[index_map[i]] = i;
}



template <>
void
AnisotropicPolynomials<0>::set_numbering(const std::vector<unsigned int> &)
{
  AssertThrow(false, ExcNotImplemented("This function does not work in 0-d!"));
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
void
TensorProductPolynomials<0, Polynomials::Polynomial<double>>::set_numbering(
  const std::vector<unsigned int> &)
{
  AssertThrow(false, ExcNotImplemented("This function does not work in 0-d!"));
}



template <int dim, typename PolynomialType>
double
TensorProductPolynomials<dim, PolynomialType>::compute_value(
  const unsigned int i,
  const Point<dim>  &p) const
{
  if constexpr (dim == 0)
    {
      (void)i;
      (void)p;
      DEAL_II_NOT_IMPLEMENTED();
      return 0;
    }
  else
    {
      std::array<unsigned int, dim> indices;
      compute_index(i, indices);

      double value = 1.;
      for (unsigned int d = 0; d < dim; ++d)
        value *= polynomials[indices[d]].value(p[d]);

      return value;
    }
}



template <int dim, typename PolynomialType>
Tensor<1, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_grad(
  const unsigned int i,
  const Point<dim>  &p) const
{
  if constexpr (dim == 0)
    {
      (void)i;
      (void)p;
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
  else
    {
      std::array<unsigned int, dim> indices;
      compute_index(i, indices);

      // compute values and
      // uni-directional derivatives at
      // the given point in each
      // coordinate direction
      ndarray<double, dim, 2> v;
      {
        std::vector<double> tmp(2);
        for (unsigned int d = 0; d < dim; ++d)
          {
            polynomials[indices[d]].value(p[d], tmp);
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
}



template <int dim, typename PolynomialType>
Tensor<2, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_grad_grad(
  const unsigned int i,
  const Point<dim>  &p) const
{
  if constexpr (dim == 0)
    {
      (void)i;
      (void)p;
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
  else
    {
      std::array<unsigned int, dim> indices;
      compute_index(i, indices);

      ndarray<double, dim, 3> v;
      {
        std::vector<double> tmp(3);
        for (unsigned int d = 0; d < dim; ++d)
          {
            polynomials[indices[d]].value(p[d], tmp);
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
}



namespace internal
{
  namespace TensorProductPolynomials
  {
    // This function computes the tensor product of some tabulated
    // one-dimensional polynomials (also the anisotropic case is supported)
    // with tensor product indices of all dimensions except the first one
    // tabulated in the 'indices' array; the first dimension is manually
    // iterated through because these are possibly performance-critical loops,
    // so we want to avoid indirect addressing.
    template <int dim, std::size_t dim1>
    void
    evaluate_tensor_product(
      const unsigned int n_derivatives,
      const boost::container::small_vector<dealii::ndarray<double, 5, dim>, 10>
                        &values_1d,
      const unsigned int size_x,
      const boost::container::small_vector<std::array<unsigned int, dim1>, 64>
                                      &indices,
      const std::vector<unsigned int> &index_map,
      std::vector<double>             &values,
      std::vector<Tensor<1, dim>>     &grads,
      std::vector<Tensor<2, dim>>     &grad_grads,
      std::vector<Tensor<3, dim>>     &third_derivatives,
      std::vector<Tensor<4, dim>>     &fourth_derivatives)
    {
      const bool update_values = (values.size() == indices.size() * size_x);
      const bool update_grads  = (grads.size() == indices.size() * size_x);
      const bool update_grad_grads =
        (grad_grads.size() == indices.size() * size_x);
      const bool update_3rd_derivatives =
        (third_derivatives.size() == indices.size() * size_x);
      const bool update_4th_derivatives =
        (fourth_derivatives.size() == indices.size() * size_x);

      // For values, 1st and 2nd derivatives use a more lengthy code that
      // minimizes the number of arithmetic operations and memory accesses
      if (n_derivatives == 0)
        for (unsigned int i = 0, i1 = 0; i1 < indices.size(); ++i1)
          {
            double value_outer = 1.;
            if constexpr (dim > 1)
              for (unsigned int d = 1; d < dim; ++d)
                value_outer *= values_1d[indices[i1][d - 1]][0][d];
            if (index_map.empty())
              for (unsigned int ix = 0; ix < size_x; ++ix, ++i)
                values[i] = value_outer * values_1d[ix][0][0];
            else
              for (unsigned int ix = 0; ix < size_x; ++ix, ++i)
                values[index_map[i]] = value_outer * values_1d[ix][0][0];
          }
      else
        for (unsigned int iy = 0, i1 = 0; i1 < indices.size(); ++i1)
          {
            // prepare parts of products in y (and z) directions
            std::array<double, dim + (dim * (dim - 1)) / 2> value_outer;
            value_outer[0] = 1.;
            if constexpr (dim > 1)
              {
                for (unsigned int x = 1; x < dim; ++x)
                  value_outer[0] *= values_1d[indices[i1][x - 1]][0][x];
                for (unsigned int d = 1; d < dim; ++d)
                  {
                    value_outer[d] = values_1d[indices[i1][d - 1]][1][d];
                    for (unsigned int x = 1; x < dim; ++x)
                      if (x != d)
                        value_outer[d] *= values_1d[indices[i1][x - 1]][0][x];
                  }
                for (unsigned int d1 = 1, count = dim; d1 < dim; ++d1)
                  for (unsigned int d2 = d1; d2 < dim; ++d2, ++count)
                    {
                      value_outer[count] = 1.;
                      for (unsigned int x = 1; x < dim; ++x)
                        {
                          unsigned int derivative = 0;
                          if (d1 == x)
                            ++derivative;
                          if (d2 == x)
                            ++derivative;

                          value_outer[count] *=
                            values_1d[indices[i1][x - 1]][derivative][x];
                        }
                    }
              }

            // now run the loop over x and multiply by the values/derivatives
            // in x direction
            for (unsigned int ix = 0, i = iy; ix < size_x; ++ix, ++i)
              {
                std::array<double, 3> val_x{{values_1d[ix][0][0],
                                             values_1d[ix][1][0],
                                             values_1d[ix][2][0]}};
                const unsigned int    index =
                  (index_map.empty() ? i : index_map[i]);

                if (update_values)
                  values[index] = value_outer[0] * val_x[0];

                if (update_grads)
                  {
                    grads[index][0] = value_outer[0] * val_x[1];
                    if constexpr (dim > 1)
                      for (unsigned int d = 1; d < dim; ++d)
                        grads[index][d] = value_outer[d] * val_x[0];
                  }

                if (update_grad_grads)
                  {
                    grad_grads[index][0][0] = value_outer[0] * val_x[2];
                    if constexpr (dim > 1)
                      {
                        for (unsigned int d = 1; d < dim; ++d)
                          grad_grads[index][0][d] = grad_grads[index][d][0] =
                            value_outer[d] * val_x[1];
                        for (unsigned int d1 = 1, count = dim; d1 < dim; ++d1)
                          for (unsigned int d2 = d1; d2 < dim; ++d2, ++count)
                            grad_grads[index][d1][d2] =
                              grad_grads[index][d2][d1] =
                                value_outer[count] * val_x[0];
                      }
                  }
              }

            // Use slower code for 3rd and 4th derivatives
            if (update_3rd_derivatives)
              for (unsigned int ix = 0, i = iy; ix < size_x; ++ix, ++i)
                {
                  const unsigned int index =
                    (index_map.empty() ? i : index_map[i]);
                  std::array<unsigned int, dim> my_indices;
                  my_indices[0] = ix;
                  if constexpr (dim > 1)
                    for (unsigned int d = 1; d < dim; ++d)
                      my_indices[d] = indices[i1][d - 1];
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

                              der3 *= values_1d[my_indices[x]][derivative][x];
                            }
                          third_derivatives[index][d1][d2][d3] = der3;
                        }
                }

            if (update_4th_derivatives)
              for (unsigned int ix = 0, i = iy; ix < size_x; ++ix, ++i)
                {
                  const unsigned int index =
                    (index_map.empty() ? i : index_map[i]);
                  std::array<unsigned int, dim> my_indices;
                  my_indices[0] = ix;
                  if constexpr (dim > 1)
                    for (unsigned int d = 1; d < dim; ++d)
                      my_indices[d] = indices[i1][d - 1];
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

                                der4 *= values_1d[my_indices[x]][derivative][x];
                              }
                            fourth_derivatives[index][d1][d2][d3][d4] = der4;
                          }
                }

            iy += size_x;
          }
    }
  } // namespace TensorProductPolynomials
} // namespace internal



template <int dim, typename PolynomialType>
void
TensorProductPolynomials<dim, PolynomialType>::evaluate(
  const Point<dim>            &p,
  std::vector<double>         &values,
  std::vector<Tensor<1, dim>> &grads,
  std::vector<Tensor<2, dim>> &grad_grads,
  std::vector<Tensor<3, dim>> &third_derivatives,
  std::vector<Tensor<4, dim>> &fourth_derivatives) const
{
  if constexpr (dim == 0)
    {
      (void)p;
      (void)values;
      (void)grads;
      (void)grad_grads;
      (void)third_derivatives;
      (void)fourth_derivatives;
      DEAL_II_NOT_IMPLEMENTED();
    }
  else
    {
      Assert(dim <= 3, ExcNotImplemented());
      Assert(values.size() == this->n() || values.empty(),
             ExcDimensionMismatch2(values.size(), this->n(), 0));
      Assert(grads.size() == this->n() || grads.empty(),
             ExcDimensionMismatch2(grads.size(), this->n(), 0));
      Assert(grad_grads.size() == this->n() || grad_grads.empty(),
             ExcDimensionMismatch2(grad_grads.size(), this->n(), 0));
      Assert(third_derivatives.size() == this->n() || third_derivatives.empty(),
             ExcDimensionMismatch2(third_derivatives.size(), this->n(), 0));
      Assert(fourth_derivatives.size() == this->n() ||
               fourth_derivatives.empty(),
             ExcDimensionMismatch2(fourth_derivatives.size(), this->n(), 0));

      // check how many values/derivatives we have to compute
      unsigned int n_derivatives = 0;
      if (values.size() == this->n())
        n_derivatives = 0;
      if (grads.size() == this->n())
        n_derivatives = 1;
      if (grad_grads.size() == this->n())
        n_derivatives = 2;
      if (third_derivatives.size() == this->n())
        n_derivatives = 3;
      if (fourth_derivatives.size() == this->n())
        n_derivatives = 4;

      // Compute the values (and derivatives, if necessary) of all 1d
      // polynomials at this evaluation point. We can use the more optimized
      // values_of_array function to compute 'dim' polynomials at once
      const unsigned int n_polynomials = polynomials.size();
      boost::container::small_vector<ndarray<double, 5, dim>, 10> values_1d(
        n_polynomials);
      if constexpr (std::is_same_v<PolynomialType,
                                   dealii::Polynomials::Polynomial<double>>)
        {
          std::array<double, dim> point_array;
          for (unsigned int d = 0; d < dim; ++d)
            point_array[d] = p[d];
          for (unsigned int i = 0; i < n_polynomials; ++i)
            polynomials[i].values_of_array(point_array,
                                           n_derivatives,
                                           values_1d[i].data());
        }
      else
        for (unsigned int i = 0; i < n_polynomials; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            {
              std::array<double, 5> derivatives;
              polynomials[i].value(p[d], n_derivatives, derivatives.data());
              for (unsigned int j = 0; j <= n_derivatives; ++j)
                values_1d[i][j][d] = derivatives[j];
            }

      // Unroll the tensor product indices of all but the first dimension in
      // arbitrary dimension
      constexpr unsigned int dim1 = dim > 1 ? dim - 1 : 1;
      boost::container::small_vector<std::array<unsigned int, dim1>, 64>
        indices(1);
      if constexpr (dim > 1)
        for (unsigned int d = 1; d < dim; ++d)
          {
            const unsigned int size = indices.size();
            for (unsigned int i = 1; i < n_polynomials; ++i)
              for (unsigned int j = 0; j < size; ++j)
                {
                  std::array<unsigned int, dim1> next_index = indices[j];
                  next_index[d - 1]                         = i;
                  indices.push_back(next_index);
                }
          }
      AssertDimension(indices.size(), Utilities::pow(n_polynomials, dim - 1));

      internal::TensorProductPolynomials::evaluate_tensor_product<dim>(
        n_derivatives,
        values_1d,
        n_polynomials,
        indices,
        index_map_inverse,
        values,
        grads,
        grad_grads,
        third_derivatives,
        fourth_derivatives);
    }
}



template <int dim, typename PolynomialType>
std::unique_ptr<ScalarPolynomialsBase<dim>>
TensorProductPolynomials<dim, PolynomialType>::clone() const
{
  return std::make_unique<TensorProductPolynomials<dim, PolynomialType>>(*this);
}



template <int dim, typename PolynomialType>
std::size_t
TensorProductPolynomials<dim, PolynomialType>::memory_consumption() const
{
  return (MemoryConsumption::memory_consumption(polynomials) +
          MemoryConsumption::memory_consumption(index_map) +
          MemoryConsumption::memory_consumption(index_map_inverse));
}



template <int dim, typename PolynomialType>
std::vector<PolynomialType>
TensorProductPolynomials<dim, PolynomialType>::get_underlying_polynomials()
  const
{
  return polynomials;
}



/* ------------------- AnisotropicPolynomials -------------- */


template <int dim>
AnisotropicPolynomials<dim>::AnisotropicPolynomials(
  const std::vector<std::vector<Polynomials::Polynomial<double>>> &pols)
  : ScalarPolynomialsBase<dim>(1, get_n_tensor_pols(pols))
  , polynomials(pols)
  , index_map(this->n())
  , index_map_inverse(this->n())
{
  Assert(pols.size() == dim, ExcDimensionMismatch(pols.size(), dim));
  for (const auto &pols_d : pols)
    {
      (void)pols_d;
      Assert(pols_d.size() > 0,
             ExcMessage("The number of polynomials must be larger than zero "
                        "for all coordinate directions."));
    }

  // per default set this index map to identity. This map can be changed by
  // the user through the set_numbering() function
  for (unsigned int i = 0; i < this->n(); ++i)
    {
      index_map[i]         = i;
      index_map_inverse[i] = i;
    }
}



template <int dim>
void
AnisotropicPolynomials<dim>::compute_index(
  const unsigned int             i,
  std::array<unsigned int, dim> &indices) const
{
  if constexpr (dim == 0)
    {
      (void)i;
      (void)indices;
      DEAL_II_NOT_IMPLEMENTED();
    }
  else
    {
      if constexpr (running_in_debug_mode())
        {
          unsigned int n_poly = 1;
          for (unsigned int d = 0; d < dim; ++d)
            n_poly *= polynomials[d].size();
          Assert(i < n_poly, ExcInternalError());
        }

      if (dim == 0)
        {
        }
      else if (dim == 1)
        internal::compute_tensor_index(index_map[i],
                                       polynomials[0].size(),
                                       0 /*not used*/,
                                       indices);
      else
        internal::compute_tensor_index(index_map[i],
                                       polynomials[0].size(),
                                       polynomials[1].size(),
                                       indices);
    }
}



template <int dim>
double
AnisotropicPolynomials<dim>::compute_value(const unsigned int i,
                                           const Point<dim>  &p) const
{
  if constexpr (dim == 0)
    {
      (void)i;
      (void)p;
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
  else
    {
      std::array<unsigned int, dim> indices;
      compute_index(i, indices);

      double value = 1.;
      for (unsigned int d = 0; d < dim; ++d)
        value *= polynomials[d][indices[d]].value(p[d]);

      return value;
    }
}



template <int dim>
Tensor<1, dim>
AnisotropicPolynomials<dim>::compute_grad(const unsigned int i,
                                          const Point<dim>  &p) const
{
  if constexpr (dim == 0)
    {
      (void)i;
      (void)p;
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
  else
    {
      std::array<unsigned int, dim> indices;
      compute_index(i, indices);

      // compute values and
      // uni-directional derivatives at
      // the given point in each
      // coordinate direction
      ndarray<double, dim, 2> v;
      for (unsigned int d = 0; d < dim; ++d)
        polynomials[d][indices[d]].value(p[d], 1, v[d].data());

      Tensor<1, dim> grad;
      for (unsigned int d = 0; d < dim; ++d)
        {
          grad[d] = 1.;
          for (unsigned int x = 0; x < dim; ++x)
            grad[d] *= v[x][d == x];
        }

      return grad;
    }
}



template <int dim>
Tensor<2, dim>
AnisotropicPolynomials<dim>::compute_grad_grad(const unsigned int i,
                                               const Point<dim>  &p) const
{
  if constexpr (dim == 0)
    {
      (void)i;
      (void)p;
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
  else
    {
      std::array<unsigned int, dim> indices;
      compute_index(i, indices);

      ndarray<double, dim, 3> v;
      for (unsigned int d = 0; d < dim; ++d)
        polynomials[d][indices[d]].value(p[d], 2, v[d].data());

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
}



template <int dim>
void
AnisotropicPolynomials<dim>::evaluate(
  const Point<dim>            &p,
  std::vector<double>         &values,
  std::vector<Tensor<1, dim>> &grads,
  std::vector<Tensor<2, dim>> &grad_grads,
  std::vector<Tensor<3, dim>> &third_derivatives,
  std::vector<Tensor<4, dim>> &fourth_derivatives) const
{
  if constexpr (dim == 0)
    {
      (void)p;
      (void)values;
      (void)grads;
      (void)grad_grads;
      (void)third_derivatives;
      (void)fourth_derivatives;
      DEAL_II_NOT_IMPLEMENTED();
    }
  else
    {
      Assert(values.size() == this->n() || values.empty(),
             ExcDimensionMismatch2(values.size(), this->n(), 0));
      Assert(grads.size() == this->n() || grads.empty(),
             ExcDimensionMismatch2(grads.size(), this->n(), 0));
      Assert(grad_grads.size() == this->n() || grad_grads.empty(),
             ExcDimensionMismatch2(grad_grads.size(), this->n(), 0));
      Assert(third_derivatives.size() == this->n() || third_derivatives.empty(),
             ExcDimensionMismatch2(third_derivatives.size(), this->n(), 0));
      Assert(fourth_derivatives.size() == this->n() ||
               fourth_derivatives.empty(),
             ExcDimensionMismatch2(fourth_derivatives.size(), this->n(), 0));

      // check how many values/derivatives we have to compute
      unsigned int n_derivatives = 0;
      if (values.size() == this->n())
        n_derivatives = 0;
      if (grads.size() == this->n())
        n_derivatives = 1;
      if (grad_grads.size() == this->n())
        n_derivatives = 2;
      if (third_derivatives.size() == this->n())
        n_derivatives = 3;
      if (fourth_derivatives.size() == this->n())
        n_derivatives = 4;

      // compute the values (and derivatives, if necessary) of all polynomials
      // at this evaluation point
      std::size_t max_n_polynomials = 0;
      for (unsigned int d = 0; d < dim; ++d)
        max_n_polynomials = std::max(max_n_polynomials, polynomials[d].size());

      // 5 is enough to store values and derivatives in all supported cases
      boost::container::small_vector<ndarray<double, 5, dim>, 10> values_1d(
        max_n_polynomials);
      if (n_derivatives == 0)
        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int i = 0; i < polynomials[d].size(); ++i)
            values_1d[i][0][d] = polynomials[d][i].value(p[d]);
      else
        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int i = 0; i < polynomials[d].size(); ++i)
            {
              // The isotropic tensor product function wants us to use a
              // different innermost index, so we cannot pass the values_1d
              // array into the function directly
              std::array<double, 5> derivatives;
              polynomials[d][i].value(p[d], n_derivatives, derivatives.data());
              for (unsigned int j = 0; j <= n_derivatives; ++j)
                values_1d[i][j][d] = derivatives[j];
            }

      // Unroll the tensor product indices in arbitrary dimension
      constexpr unsigned int dim1 = dim > 1 ? dim - 1 : 1;
      boost::container::small_vector<std::array<unsigned int, dim1>, 64>
        indices(1);
      for (unsigned int d = 1; d < dim; ++d)
        {
          const unsigned int size = indices.size();
          for (unsigned int i = 1; i < polynomials[d].size(); ++i)
            for (unsigned int j = 0; j < size; ++j)
              {
                std::array<unsigned int, dim1> next_index = indices[j];
                next_index[d - 1]                         = i;
                indices.push_back(next_index);
              }
        }

      internal::TensorProductPolynomials::evaluate_tensor_product<dim>(
        n_derivatives,
        values_1d,
        polynomials[0].size(),
        indices,
        index_map_inverse,
        values,
        grads,
        grad_grads,
        third_derivatives,
        fourth_derivatives);
    }
}



template <int dim>
unsigned int
AnisotropicPolynomials<dim>::get_n_tensor_pols(
  const std::vector<std::vector<Polynomials::Polynomial<double>>> &pols)
{
  if constexpr (dim == 0)
    {
      (void)pols;
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
  else
    {
      unsigned int y = 1;
      for (unsigned int d = 0; d < dim; ++d)
        y *= pols[d].size();
      return y;
    }
}



template <int dim>
std::unique_ptr<ScalarPolynomialsBase<dim>>
AnisotropicPolynomials<dim>::clone() const
{
  return std::make_unique<AnisotropicPolynomials<dim>>(*this);
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
