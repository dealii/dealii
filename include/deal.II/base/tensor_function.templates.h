// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tensor_function_templates_h
#define dealii_tensor_function_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/vector.h>

#include <cmath>
#include <vector>

DEAL_II_NAMESPACE_OPEN


template <int rank, int dim, typename Number>
TensorFunction<rank, dim, Number>::TensorFunction(
  const typename TensorFunction<rank, dim, Number>::time_type initial_time)
  : FunctionTime<typename TensorFunction<rank, dim, Number>::time_type>(
      initial_time)
{}



template <int rank, int dim, typename Number>
typename TensorFunction<rank, dim, Number>::value_type
TensorFunction<rank, dim, Number>::value(const Point<dim> &) const
{
  Assert(false, ExcPureFunctionCalled());
  return Tensor<rank, dim, Number>();
}


template <int rank, int dim, typename Number>
void
TensorFunction<rank, dim, Number>::value_list(
  const std::vector<Point<dim>> &points,
  std::vector<value_type>       &values) const
{
  Assert(values.size() == points.size(),
         ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i = 0; i < points.size(); ++i)
    values[i] = this->value(points[i]);
}


template <int rank, int dim, typename Number>
typename TensorFunction<rank, dim, Number>::gradient_type
TensorFunction<rank, dim, Number>::gradient(const Point<dim> &) const
{
  Assert(false, ExcPureFunctionCalled());
  return Tensor<rank + 1, dim, Number>();
}


template <int rank, int dim, typename Number>
void
TensorFunction<rank, dim, Number>::gradient_list(
  const std::vector<Point<dim>> &points,
  std::vector<gradient_type>    &gradients) const
{
  Assert(gradients.size() == points.size(),
         ExcDimensionMismatch(gradients.size(), points.size()));

  for (unsigned int i = 0; i < points.size(); ++i)
    gradients[i] = gradient(points[i]);
}



template <int rank, int dim, typename Number>
ConstantTensorFunction<rank, dim, Number>::ConstantTensorFunction(
  const Tensor<rank, dim, Number> &value,
  const typename ConstantTensorFunction<rank, dim, Number>::time_type
    initial_time)
  : TensorFunction<rank, dim, Number>(initial_time)
  , _value(value)
{}



template <int rank, int dim, typename Number>
typename TensorFunction<rank, dim, Number>::value_type
ConstantTensorFunction<rank, dim, Number>::value(
  const Point<dim> & /*point*/) const
{
  return _value;
}


template <int rank, int dim, typename Number>
void
ConstantTensorFunction<rank, dim, Number>::value_list(
  const std::vector<Point<dim>>                                       &points,
  std::vector<typename TensorFunction<rank, dim, Number>::value_type> &values)
  const
{
  (void)points;
  Assert(values.size() == points.size(),
         ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = _value;
}


template <int rank, int dim, typename Number>
typename TensorFunction<rank, dim, Number>::gradient_type
ConstantTensorFunction<rank, dim, Number>::gradient(const Point<dim> &) const
{
  // Return a zero (=default initialized) tensor
  return {};
}


template <int rank, int dim, typename Number>
void
ConstantTensorFunction<rank, dim, Number>::gradient_list(
  const std::vector<Point<dim>> &points,
  std::vector<typename TensorFunction<rank, dim, Number>::gradient_type>
    &gradients) const
{
  (void)points;
  Assert(gradients.size() == points.size(),
         ExcDimensionMismatch(gradients.size(), points.size()));

  // Return an array of zero tensors.
  std::fill(gradients.begin(),
            gradients.end(),
            typename TensorFunction<rank, dim, Number>::gradient_type());
}



template <int rank, int dim, typename Number>
ZeroTensorFunction<rank, dim, Number>::ZeroTensorFunction(
  const typename ZeroTensorFunction<rank, dim, Number>::time_type initial_time)
  : ConstantTensorFunction<rank, dim, Number>(
      dealii::Tensor<rank, dim, Number>(),
      initial_time)
{}


DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_tensor_function_templates_h */
