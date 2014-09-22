// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#include <deal.II/base/tensor_function.h>
#include <vector>
#include <deal.II/base/tensor.h>
#include <cmath>
#include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN



//////////////////////////////////////////////////////////////////////
// TensorFunction
//////////////////////////////////////////////////////////////////////

template <int rank, int dim, typename Number>
TensorFunction<rank, dim, Number>::TensorFunction (const Number initial_time)
  :
  FunctionTime<Number> (initial_time)
{}



template <int rank, int dim, typename Number>
TensorFunction<rank, dim, Number>::~TensorFunction ()
{}


template <int rank, int dim, typename Number>
typename TensorFunction<rank, dim, Number>::value_type
TensorFunction<rank, dim, Number>::value (const Point<dim,Number> &) const
{
  Assert (false, ExcPureFunctionCalled());
  return Tensor<rank,dim>();
}


template <int rank, int dim, typename Number>
void
TensorFunction<rank, dim, Number>::value_list (const std::vector<Point<dim,Number> > &points,
                                               std::vector<value_type>        &values) const
{
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    values[i]  = this->value (points[i]);
}


template <int rank, int dim, typename Number>
typename TensorFunction<rank, dim, Number>::gradient_type
TensorFunction<rank, dim, Number>::gradient (const Point<dim,Number> &) const
{
  Assert (false, ExcPureFunctionCalled());
  return Tensor<rank+1,dim>();
}


template <int rank, int dim, typename Number>
void
TensorFunction<rank, dim, Number>::gradient_list (const std::vector<Point<dim,Number> >   &points,
                                                  std::vector<gradient_type> &gradients) const
{
  Assert (gradients.size() == points.size(),
          ExcDimensionMismatch(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    gradients[i] = gradient(points[i]);
}



//////////////////////////////////////////////////////////////////////
// ConstantTensorFunction
//////////////////////////////////////////////////////////////////////


template <int rank, int dim, typename Number>
ConstantTensorFunction<rank, dim, Number>::ConstantTensorFunction (const Tensor<rank, dim, Number> &value,
    const Number initial_time)
  :
  TensorFunction<rank, dim, Number> (initial_time),
  _value(value)
{}


template <int rank, int dim, typename Number>
ConstantTensorFunction<rank, dim, Number>::~ConstantTensorFunction ()
{}


template <int rank, int dim, typename Number>
typename TensorFunction<rank, dim, Number>::value_type
ConstantTensorFunction<rank, dim, Number>::value (const Point<dim,Number> &/*point*/) const
{
  return _value;
}


template <int rank, int dim, typename Number>
void
ConstantTensorFunction<rank, dim, Number>::value_list (const std::vector<Point<dim,Number> > &points,
                                                       std::vector<typename TensorFunction<rank, dim, Number>::value_type> &values) const
{
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0; i<values.size(); ++i)
    values[i]  = _value;
}


template <int rank, int dim, typename Number>
typename TensorFunction<rank, dim, Number>::gradient_type
ConstantTensorFunction<rank, dim, Number>::gradient (const Point<dim,Number> &) const
{
  static const Tensor<rank+1,dim> zero;

  return zero;
}


template <int rank, int dim, typename Number>
void
ConstantTensorFunction<rank, dim, Number>::gradient_list (const std::vector<Point<dim,Number> >   &points,
                                                          std::vector<typename TensorFunction<rank, dim, Number>::gradient_type> &gradients) const
{
  Assert (gradients.size() == points.size(),
          ExcDimensionMismatch(gradients.size(), points.size()));

  static const Tensor<rank+1,dim> zero;

  for (unsigned int i=0; i<gradients.size(); ++i)
    gradients[i] = zero;
}



//////////////////////////////////////////////////////////////////////
// ZeroTensorFunction
//////////////////////////////////////////////////////////////////////


template <int rank, int dim, typename Number>
ZeroTensorFunction<rank, dim, Number>::ZeroTensorFunction (const Number initial_time)
  :
  ConstantTensorFunction<rank, dim, Number> (dealii::Tensor<rank, dim, Number>(), initial_time)
{}



//////////////////////////////////////////////////////////////////////
// Explicit instantiations:
//////////////////////////////////////////////////////////////////////

template class TensorFunction<1,1>;
template class TensorFunction<2,1>;
template class TensorFunction<3,1>;
template class TensorFunction<4,1>;
template class TensorFunction<1,2>;
template class TensorFunction<2,2>;
template class TensorFunction<3,2>;
template class TensorFunction<4,2>;
template class TensorFunction<1,3>;
template class TensorFunction<2,3>;
template class TensorFunction<3,3>;
template class TensorFunction<4,3>;

template class ConstantTensorFunction<1,1>;
template class ConstantTensorFunction<2,1>;
template class ConstantTensorFunction<3,1>;
template class ConstantTensorFunction<4,1>;
template class ConstantTensorFunction<1,2>;
template class ConstantTensorFunction<2,2>;
template class ConstantTensorFunction<3,2>;
template class ConstantTensorFunction<4,2>;
template class ConstantTensorFunction<1,3>;
template class ConstantTensorFunction<2,3>;
template class ConstantTensorFunction<3,3>;
template class ConstantTensorFunction<4,3>;

template class ZeroTensorFunction<1,1>;
template class ZeroTensorFunction<2,1>;
template class ZeroTensorFunction<3,1>;
template class ZeroTensorFunction<4,1>;
template class ZeroTensorFunction<1,2>;
template class ZeroTensorFunction<2,2>;
template class ZeroTensorFunction<3,2>;
template class ZeroTensorFunction<4,2>;
template class ZeroTensorFunction<1,3>;
template class ZeroTensorFunction<2,3>;
template class ZeroTensorFunction<3,3>;
template class ZeroTensorFunction<4,3>;


DEAL_II_NAMESPACE_CLOSE
