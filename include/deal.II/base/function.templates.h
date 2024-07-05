// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_function_templates_h
#define dealii_function_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/vector.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


template <int dim, typename RangeNumberType, typename PointNumberType>
const unsigned int Function<dim, RangeNumberType, PointNumberType>::dimension;


template <int dim, typename RangeNumberType, typename PointNumberType>
Function<dim, RangeNumberType, PointNumberType>::Function(
  const unsigned int n_components,
  const typename Function<dim, RangeNumberType, PointNumberType>::time_type
    initial_time)
  : FunctionTime<
      typename Function<dim, RangeNumberType, PointNumberType>::time_type>(
      initial_time)
  , n_components(n_components)
{
  // avoid the construction of function objects that don't return any
  // values. This doesn't make much sense in the first place, but will lead
  // to odd errors later on (happened to me in fact :-)
  Assert(n_components > 0, ExcZero());
}



template <int dim, typename RangeNumberType, typename PointNumberType>
Function<dim, RangeNumberType, PointNumberType> &
Function<dim, RangeNumberType, PointNumberType>::operator=(const Function &f)
{
  (void)f;
  AssertDimension(n_components, f.n_components);
  return *this;
}


template <int dim, typename RangeNumberType, typename PointNumberType>
RangeNumberType
Function<dim, RangeNumberType, PointNumberType>::value(
  const Point<dim, PointNumberType> &,
  const unsigned int) const
{
  Assert(false, ExcPureFunctionCalled());
  return 0;
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::vector_value(
  const Point<dim, PointNumberType> &p,
  Vector<RangeNumberType>           &v) const
{
  AssertDimension(v.size(), this->n_components);
  for (unsigned int i = 0; i < this->n_components; ++i)
    v(i) = value(p, i);
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::value_list(
  const std::vector<Point<dim, PointNumberType>> &points,
  std::vector<RangeNumberType>                   &values,
  const unsigned int                              component) const
{
  // check whether component is in the valid range is up to the derived
  // class
  Assert(values.size() == points.size(),
         ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i = 0; i < points.size(); ++i)
    values[i] = this->value(points[i], component);
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::vector_value_list(
  const std::vector<Point<dim, PointNumberType>> &points,
  std::vector<Vector<RangeNumberType>>           &values) const
{
  // check whether component is in the valid range is up to the derived
  // class
  Assert(values.size() == points.size(),
         ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i = 0; i < points.size(); ++i)
    this->vector_value(points[i], values[i]);
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::vector_values(
  const std::vector<Point<dim, PointNumberType>> &points,
  std::vector<std::vector<RangeNumberType>>      &values) const
{
  const unsigned int n = this->n_components;
  AssertDimension(values.size(), n);
  for (unsigned int i = 0; i < n; ++i)
    value_list(points, values[i], i);
}


template <int dim, typename RangeNumberType, typename PointNumberType>
Tensor<1, dim, RangeNumberType>
Function<dim, RangeNumberType, PointNumberType>::gradient(
  const Point<dim, PointNumberType> &,
  const unsigned int) const
{
  Assert(false, ExcPureFunctionCalled());
  return Tensor<1, dim, RangeNumberType>();
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::vector_gradient(
  const Point<dim, PointNumberType>            &p,
  std::vector<Tensor<1, dim, RangeNumberType>> &v) const
{
  AssertDimension(v.size(), this->n_components);
  for (unsigned int i = 0; i < this->n_components; ++i)
    v[i] = gradient(p, i);
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::gradient_list(
  const std::vector<Point<dim, PointNumberType>> &points,
  std::vector<Tensor<1, dim, RangeNumberType>>   &gradients,
  const unsigned int                              component) const
{
  Assert(gradients.size() == points.size(),
         ExcDimensionMismatch(gradients.size(), points.size()));

  for (unsigned int i = 0; i < points.size(); ++i)
    gradients[i] = gradient(points[i], component);
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::vector_gradient_list(
  const std::vector<Point<dim, PointNumberType>>            &points,
  std::vector<std::vector<Tensor<1, dim, RangeNumberType>>> &gradients) const
{
  Assert(gradients.size() == points.size(),
         ExcDimensionMismatch(gradients.size(), points.size()));

  for (unsigned int i = 0; i < points.size(); ++i)
    {
      Assert(gradients[i].size() == n_components,
             ExcDimensionMismatch(gradients[i].size(), n_components));
      vector_gradient(points[i], gradients[i]);
    }
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::vector_gradients(
  const std::vector<Point<dim, PointNumberType>>            &points,
  std::vector<std::vector<Tensor<1, dim, RangeNumberType>>> &values) const
{
  const unsigned int n = this->n_components;
  AssertDimension(values.size(), n);
  for (unsigned int i = 0; i < n; ++i)
    gradient_list(points, values[i], i);
}



template <int dim, typename RangeNumberType, typename PointNumberType>
RangeNumberType
Function<dim, RangeNumberType, PointNumberType>::laplacian(
  const Point<dim, PointNumberType> &,
  const unsigned int) const
{
  Assert(false, ExcPureFunctionCalled());
  return 0;
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::vector_laplacian(
  const Point<dim, PointNumberType> &,
  Vector<RangeNumberType> &) const
{
  Assert(false, ExcPureFunctionCalled());
}



template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::laplacian_list(
  const std::vector<Point<dim, PointNumberType>> &points,
  std::vector<RangeNumberType>                   &laplacians,
  const unsigned int                              component) const
{
  // check whether component is in the valid range is up to the derived
  // class
  Assert(laplacians.size() == points.size(),
         ExcDimensionMismatch(laplacians.size(), points.size()));

  for (unsigned int i = 0; i < points.size(); ++i)
    laplacians[i] = this->laplacian(points[i], component);
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::vector_laplacian_list(
  const std::vector<Point<dim, PointNumberType>> &points,
  std::vector<Vector<RangeNumberType>>           &laplacians) const
{
  // check whether component is in the valid range is up to the derived
  // class
  Assert(laplacians.size() == points.size(),
         ExcDimensionMismatch(laplacians.size(), points.size()));

  for (unsigned int i = 0; i < points.size(); ++i)
    this->vector_laplacian(points[i], laplacians[i]);
}


template <int dim, typename RangeNumberType, typename PointNumberType>
SymmetricTensor<2, dim, RangeNumberType>
Function<dim, RangeNumberType, PointNumberType>::hessian(
  const Point<dim, PointNumberType> &,
  const unsigned int) const
{
  Assert(false, ExcPureFunctionCalled());
  return SymmetricTensor<2, dim, RangeNumberType>();
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::vector_hessian(
  const Point<dim, PointNumberType>                     &p,
  std::vector<SymmetricTensor<2, dim, RangeNumberType>> &v) const
{
  AssertDimension(v.size(), this->n_components);
  for (unsigned int i = 0; i < this->n_components; ++i)
    v[i] = hessian(p, i);
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::hessian_list(
  const std::vector<Point<dim, PointNumberType>>        &points,
  std::vector<SymmetricTensor<2, dim, RangeNumberType>> &hessians,
  const unsigned int                                     component) const
{
  Assert(hessians.size() == points.size(),
         ExcDimensionMismatch(hessians.size(), points.size()));

  for (unsigned int i = 0; i < points.size(); ++i)
    hessians[i] = hessian(points[i], component);
}


template <int dim, typename RangeNumberType, typename PointNumberType>
void
Function<dim, RangeNumberType, PointNumberType>::vector_hessian_list(
  const std::vector<Point<dim, PointNumberType>>                     &points,
  std::vector<std::vector<SymmetricTensor<2, dim, RangeNumberType>>> &hessians)
  const
{
  Assert(hessians.size() == points.size(),
         ExcDimensionMismatch(hessians.size(), points.size()));

  for (unsigned int i = 0; i < points.size(); ++i)
    {
      Assert(hessians[i].size() == n_components,
             ExcDimensionMismatch(hessians[i].size(), n_components));
      vector_hessian(points[i], hessians[i]);
    }
}



template <int dim, typename RangeNumberType, typename PointNumberType>
std::size_t
Function<dim, RangeNumberType, PointNumberType>::memory_consumption() const
{
  // only simple data elements, so use sizeof operator
  return sizeof(*this);
}



//---------------------------------------------------------------------------

namespace Functions
{
  template <int dim, typename RangeNumberType, typename PointNumberType>
  ConstantFunction<dim, RangeNumberType, PointNumberType>::ConstantFunction(
    const RangeNumberType value,
    const unsigned int    n_components)
    : Function<dim, RangeNumberType, PointNumberType>(n_components)
    , function_value_vector(n_components, value)
  {}



  template <int dim, typename RangeNumberType, typename PointNumberType>
  ConstantFunction<dim, RangeNumberType, PointNumberType>::ConstantFunction(
    const std::vector<RangeNumberType> &values)
    : Function<dim, RangeNumberType, PointNumberType>(values.size())
    , function_value_vector(values)
  {}



  template <int dim, typename RangeNumberType, typename PointNumberType>
  ConstantFunction<dim, RangeNumberType, PointNumberType>::ConstantFunction(
    const Vector<RangeNumberType> &values)
    : Function<dim, RangeNumberType, PointNumberType>(values.size())
    , function_value_vector(values.size())
  {
    Assert(values.size() == function_value_vector.size(),
           ExcDimensionMismatch(values.size(), function_value_vector.size()));
    std::copy(values.begin(), values.end(), function_value_vector.begin());
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  ConstantFunction<dim, RangeNumberType, PointNumberType>::ConstantFunction(
    const RangeNumberType *begin_ptr,
    const unsigned int     n_components)
    : Function<dim, RangeNumberType, PointNumberType>(n_components)
    , function_value_vector(n_components)
  {
    Assert(begin_ptr != nullptr, ExcMessage("Null pointer encountered!"));
    std::copy(begin_ptr,
              begin_ptr + n_components,
              function_value_vector.begin());
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  RangeNumberType
  ConstantFunction<dim, RangeNumberType, PointNumberType>::value(
    const Point<dim, PointNumberType> &,
    const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    return function_value_vector[component];
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  void
  ConstantFunction<dim, RangeNumberType, PointNumberType>::vector_value(
    const Point<dim, PointNumberType> &,
    Vector<RangeNumberType> &return_value) const
  {
    Assert(return_value.size() == this->n_components,
           ExcDimensionMismatch(return_value.size(), this->n_components));

    std::copy(function_value_vector.begin(),
              function_value_vector.end(),
              return_value.begin());
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  void
  ConstantFunction<dim, RangeNumberType, PointNumberType>::value_list(
    const std::vector<Point<dim, PointNumberType>> &points,
    std::vector<RangeNumberType>                   &return_values,
    const unsigned int                              component) const
  {
    // To avoid warning of unused parameter
    (void)points;
    AssertIndexRange(component, this->n_components);
    Assert(return_values.size() == points.size(),
           ExcDimensionMismatch(return_values.size(), points.size()));

    std::fill(return_values.begin(),
              return_values.end(),
              function_value_vector[component]);
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  void
  ConstantFunction<dim, RangeNumberType, PointNumberType>::vector_value_list(
    const std::vector<Point<dim, PointNumberType>> &points,
    std::vector<Vector<RangeNumberType>>           &return_values) const
  {
    Assert(return_values.size() == points.size(),
           ExcDimensionMismatch(return_values.size(), points.size()));

    for (unsigned int i = 0; i < points.size(); ++i)
      {
        Assert(return_values[i].size() == this->n_components,
               ExcDimensionMismatch(return_values[i].size(),
                                    this->n_components));
        std::copy(function_value_vector.begin(),
                  function_value_vector.end(),
                  return_values[i].begin());
      };
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  std::size_t
  ConstantFunction<dim, RangeNumberType, PointNumberType>::memory_consumption()
    const
  {
    // Here we assume RangeNumberType is a simple type.
    return (sizeof(*this) + this->n_components * sizeof(RangeNumberType));
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  Tensor<1, dim, RangeNumberType>
  ConstantFunction<dim, RangeNumberType, PointNumberType>::gradient(
    const Point<dim, PointNumberType> &,
    const unsigned int) const
  {
    return Tensor<1, dim, RangeNumberType>();
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  void
  ConstantFunction<dim, RangeNumberType, PointNumberType>::vector_gradient(
    const Point<dim, PointNumberType> &,
    std::vector<Tensor<1, dim, RangeNumberType>> &gradients) const
  {
    Assert(gradients.size() == this->n_components,
           ExcDimensionMismatch(gradients.size(), this->n_components));

    for (unsigned int c = 0; c < this->n_components; ++c)
      gradients[c].clear();
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  void
  ConstantFunction<dim, RangeNumberType, PointNumberType>::gradient_list(
    const std::vector<Point<dim, PointNumberType>> &points,
    std::vector<Tensor<1, dim, RangeNumberType>>   &gradients,
    const unsigned int /*component*/) const
  {
    Assert(gradients.size() == points.size(),
           ExcDimensionMismatch(gradients.size(), points.size()));

    for (unsigned int i = 0; i < points.size(); ++i)
      gradients[i].clear();
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  void
  ConstantFunction<dim, RangeNumberType, PointNumberType>::vector_gradient_list(
    const std::vector<Point<dim, PointNumberType>>            &points,
    std::vector<std::vector<Tensor<1, dim, RangeNumberType>>> &gradients) const
  {
    Assert(gradients.size() == points.size(),
           ExcDimensionMismatch(gradients.size(), points.size()));
    for (unsigned int i = 0; i < points.size(); ++i)
      {
        Assert(gradients[i].size() == this->n_components,
               ExcDimensionMismatch(gradients[i].size(), this->n_components));
        for (unsigned int c = 0; c < this->n_components; ++c)
          gradients[i][c].clear();
      };
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  SymmetricTensor<2, dim, RangeNumberType>
  ConstantFunction<dim, RangeNumberType, PointNumberType>::hessian(
    const Point<dim, PointNumberType> &,
    const unsigned int) const
  {
    return SymmetricTensor<2, dim, RangeNumberType>();
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  RangeNumberType
  ConstantFunction<dim, RangeNumberType, PointNumberType>::laplacian(
    const Point<dim, PointNumberType> &,
    const unsigned int) const
  {
    return 0;
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  ZeroFunction<dim, RangeNumberType, PointNumberType>::ZeroFunction(
    const unsigned int n_components)
    : ConstantFunction<dim, RangeNumberType, PointNumberType>(RangeNumberType(),
                                                              n_components)
  {}



  template <int dim, typename RangeNumberType, typename PointNumberType>
  IdentityFunction<dim, RangeNumberType, PointNumberType>::IdentityFunction()
    : Function<dim, RangeNumberType, PointNumberType>(dim)
  {}



  template <int dim, typename RangeNumberType, typename PointNumberType>
  RangeNumberType
  IdentityFunction<dim, RangeNumberType, PointNumberType>::value(
    const Point<dim, PointNumberType> &p,
    const unsigned int                 component) const
  {
    AssertIndexRange(component, this->n_components);
    return p[component];
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  Tensor<1, dim, RangeNumberType>
  IdentityFunction<dim, RangeNumberType, PointNumberType>::gradient(
    const Point<dim, PointNumberType> &,
    const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    Tensor<1, dim, RangeNumberType> result;
    result[component] = RangeNumberType(1);
    return result;
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  SymmetricTensor<2, dim, RangeNumberType>
  IdentityFunction<dim, RangeNumberType, PointNumberType>::hessian(
    const Point<dim, PointNumberType> &,
    const unsigned int) const
  {
    return SymmetricTensor<2, dim, RangeNumberType>();
  }



  template <int dim, typename RangeNumberType, typename PointNumberType>
  RangeNumberType
  IdentityFunction<dim, RangeNumberType, PointNumberType>::laplacian(
    const Point<dim, PointNumberType> &,
    const unsigned int) const
  {
    return 0;
  }
} // namespace Functions

//---------------------------------------------------------------------------

template <int dim, typename RangeNumberType, typename PointNumberType>
ComponentSelectFunction<dim, RangeNumberType, PointNumberType>::
  ComponentSelectFunction(const unsigned int    selected,
                          const RangeNumberType value,
                          const unsigned int    n_components)
  : Functions::ConstantFunction<dim, RangeNumberType, PointNumberType>(
      value,
      n_components)
  , selected_components(std::make_pair(selected, selected + 1))
{}



template <int dim, typename RangeNumberType, typename PointNumberType>
ComponentSelectFunction<dim, RangeNumberType, PointNumberType>::
  ComponentSelectFunction(const unsigned int selected,
                          const unsigned int n_components)
  : Functions::ConstantFunction<dim, RangeNumberType, PointNumberType>(
      1.,
      n_components)
  , selected_components(std::make_pair(selected, selected + 1))
{
  AssertIndexRange(selected, n_components);
}



template <int dim, typename RangeNumberType, typename PointNumberType>
ComponentSelectFunction<dim, RangeNumberType, PointNumberType>::
  ComponentSelectFunction(const std::pair<unsigned int, unsigned int> &selected,
                          const unsigned int n_components)
  : Functions::ConstantFunction<dim, RangeNumberType, PointNumberType>(
      1.,
      n_components)
  , selected_components(selected)
{
  Assert(selected_components.first < selected_components.second,
         ExcMessage("The upper bound of the interval must be larger than "
                    "the lower bound"));
  Assert(selected_components.second <= n_components,
         ExcMessage("The upper bound of the interval must be less than "
                    "or equal to the total number of vector components"));
}



template <int dim, typename RangeNumberType, typename PointNumberType>
void
ComponentSelectFunction<dim, RangeNumberType, PointNumberType>::
  substitute_function_value_with(
    const Functions::ConstantFunction<dim, RangeNumberType, PointNumberType> &f)
{
  Point<dim, PointNumberType> p;
  for (unsigned int i = 0; i < this->function_value_vector.size(); ++i)
    this->function_value_vector[i] = f.value(p, i);
}



template <int dim, typename RangeNumberType, typename PointNumberType>
void
ComponentSelectFunction<dim, RangeNumberType, PointNumberType>::vector_value(
  const Point<dim, PointNumberType> &,
  Vector<RangeNumberType> &return_value) const
{
  Assert(return_value.size() == this->n_components,
         ExcDimensionMismatch(return_value.size(), this->n_components));

  return_value = 0;
  std::copy(this->function_value_vector.begin() + selected_components.first,
            this->function_value_vector.begin() + selected_components.second,
            return_value.begin() + selected_components.first);
}



template <int dim, typename RangeNumberType, typename PointNumberType>
void
ComponentSelectFunction<dim, RangeNumberType, PointNumberType>::
  vector_value_list(const std::vector<Point<dim, PointNumberType>> &points,
                    std::vector<Vector<RangeNumberType>> &values) const
{
  Assert(values.size() == points.size(),
         ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i = 0; i < points.size(); ++i)
    ComponentSelectFunction<dim, RangeNumberType, PointNumberType>::
      vector_value(points[i], values[i]);
}



template <int dim, typename RangeNumberType, typename PointNumberType>
std::size_t
ComponentSelectFunction<dim, RangeNumberType, PointNumberType>::
  memory_consumption() const
{
  // No new complex data structure is introduced here, just evaluate how much
  // more memory is used *inside* the class via sizeof() and add that value to
  // parent class's memory_consumption()
  return (
    sizeof(*this) -
    sizeof(Functions::ConstantFunction<dim, RangeNumberType, PointNumberType>) +
    Functions::ConstantFunction<dim, RangeNumberType, PointNumberType>::
      memory_consumption());
}

//---------------------------------------------------------------------------

template <int dim, typename RangeNumberType, typename PointNumberType>
ScalarFunctionFromFunctionObject<dim, RangeNumberType, PointNumberType>::
  ScalarFunctionFromFunctionObject(
    const std::function<RangeNumberType(const Point<dim, PointNumberType> &)>
      &function_object)
  : Function<dim, RangeNumberType, PointNumberType>(1)
  , function_object(function_object)
{}



template <int dim, typename RangeNumberType, typename PointNumberType>
RangeNumberType
ScalarFunctionFromFunctionObject<dim, RangeNumberType, PointNumberType>::value(
  const Point<dim, PointNumberType> &p,
  const unsigned int                 component) const
{
  (void)component;
  Assert(component == 0,
         ExcMessage("This object represents only scalar functions"));
  return function_object(p);
}



template <int dim, typename RangeNumberType, typename PointNumberType>
VectorFunctionFromScalarFunctionObject<dim, RangeNumberType, PointNumberType>::
  VectorFunctionFromScalarFunctionObject(
    const std::function<RangeNumberType(const Point<dim, PointNumberType> &)>
                      &function_object,
    const unsigned int selected_component,
    const unsigned int n_components)
  : Function<dim, RangeNumberType, PointNumberType>(n_components)
  , function_object(function_object)
  , selected_component(selected_component)
{
  AssertIndexRange(selected_component, this->n_components);
}



template <int dim, typename RangeNumberType, typename PointNumberType>
RangeNumberType
VectorFunctionFromScalarFunctionObject<dim, RangeNumberType, PointNumberType>::
  value(const Point<dim, PointNumberType> &p,
        const unsigned int                 component) const
{
  AssertIndexRange(component, this->n_components);

  if (component == selected_component)
    return function_object(p);
  else
    return 0;
}



template <int dim, typename RangeNumberType, typename PointNumberType>
void
VectorFunctionFromScalarFunctionObject<dim, RangeNumberType, PointNumberType>::
  vector_value(const Point<dim, PointNumberType> &p,
               Vector<RangeNumberType>           &values) const
{
  AssertDimension(values.size(), this->n_components);

  // set everything to zero, and then the right component to its correct
  // value
  values                     = 0;
  values(selected_component) = function_object(p);
}



/**
 * The constructor for <tt>VectorFunctionFromTensorFunction</tt> which
 * initiates the return vector to be size <tt>n_components</tt>.
 */
template <int dim, typename RangeNumberType, typename PointNumberType>
VectorFunctionFromTensorFunction<dim, RangeNumberType, PointNumberType>::
  VectorFunctionFromTensorFunction(
    const TensorFunction<1, dim, RangeNumberType, PointNumberType>
                      &tensor_function,
    const unsigned int selected_component,
    const unsigned int n_components)
  : Function<dim, RangeNumberType, PointNumberType>(n_components)
  , tensor_function(tensor_function)
  , selected_component(selected_component)
{
  // Verify that the Tensor<1,dim,RangeNumberType> will fit in the given length
  // selected_components and not hang over the end of the vector.
  AssertIndexRange(selected_component + dim - 1, this->n_components);
}



template <int dim, typename RangeNumberType, typename PointNumberType>
inline RangeNumberType
VectorFunctionFromTensorFunction<dim, RangeNumberType, PointNumberType>::value(
  const Point<dim, PointNumberType> &p,
  const unsigned int                 component) const
{
  AssertIndexRange(component, this->n_components);

  // if the requested component is out of the range selected, then we can
  // return early
  if ((component < selected_component) ||
      (component >= selected_component + dim))
    return 0;

  // otherwise retrieve the values from the <tt>tensor_function</tt> to be
  // placed at the <tt>selected_component</tt> to
  // <tt>selected_component + dim - 1</tt> elements of the <tt>Vector</tt>
  // values and pick the correct one
  const Tensor<1, dim, RangeNumberType> tensor_value = tensor_function.value(p);

  return tensor_value[component - selected_component];
}


template <int dim, typename RangeNumberType, typename PointNumberType>
inline void
VectorFunctionFromTensorFunction<dim, RangeNumberType, PointNumberType>::
  vector_value(const Point<dim, PointNumberType> &p,
               Vector<RangeNumberType>           &values) const
{
  Assert(values.size() == this->n_components,
         ExcDimensionMismatch(values.size(), this->n_components));

  // Retrieve the values from the <tt>tensor_function</tt> to be placed at
  // the <tt>selected_component</tt> to
  // <tt>selected_component + dim - 1</tt> elements of the <tt>Vector</tt>
  // values.
  const Tensor<1, dim, RangeNumberType> tensor_value = tensor_function.value(p);

  // First we make all elements of values = 0
  values = 0;

  // Second we adjust the desired components to take on the values in
  // <tt>tensor_value</tt>.
  for (unsigned int i = 0; i < dim; ++i)
    values(i + selected_component) = tensor_value[i];
}


/**
 * Member function <tt>vector_value_list </tt> is the interface for giving a
 * list of points (<code>vector<Point<dim, PointNumberType> ></code>) of which
 * to evaluate using the <tt>vector_value</tt> member function.  Again, this
 * function is written so as to not replicate the function definition but passes
 * each point on to <tt>vector_value</tt> to be evaluated.
 */
template <int dim, typename RangeNumberType, typename PointNumberType>
void
VectorFunctionFromTensorFunction<dim, RangeNumberType, PointNumberType>::
  vector_value_list(const std::vector<Point<dim, PointNumberType>> &points,
                    std::vector<Vector<RangeNumberType>> &value_list) const
{
  Assert(value_list.size() == points.size(),
         ExcDimensionMismatch(value_list.size(), points.size()));

  const unsigned int n_points = points.size();

  for (unsigned int p = 0; p < n_points; ++p)
    VectorFunctionFromTensorFunction<dim, RangeNumberType, PointNumberType>::
      vector_value(points[p], value_list[p]);
}

template <int dim, typename RangeNumberType, typename PointNumberType>
inline Tensor<1, dim, RangeNumberType>
VectorFunctionFromTensorFunction<dim, RangeNumberType, PointNumberType>::
  gradient(const Point<dim, PointNumberType> &p,
           const unsigned int                 component) const
{
  AssertIndexRange(component, this->n_components);

  // if the requested component is out of the range selected, then we can
  // return early
  if ((component < selected_component) ||
      (component >= selected_component + dim))
    return Tensor<1, dim>();

  // // otherwise retrieve the values from the <tt>tensor_function</tt> to be
  // // placed at the <tt>selected_component</tt> to
  // // <tt>selected_component + dim - 1</tt> elements of the <tt>Vector</tt>
  // // values and pick the correct one
  const Tensor<2, dim, RangeNumberType> tensor_gradient =
    tensor_function.gradient(p);

  return tensor_gradient[component - selected_component];
}

template <int dim, typename RangeNumberType, typename PointNumberType>
void
VectorFunctionFromTensorFunction<dim, RangeNumberType, PointNumberType>::
  vector_gradient(const Point<dim, PointNumberType>            &p,
                  std::vector<Tensor<1, dim, RangeNumberType>> &gradients) const
{
  AssertDimension(gradients.size(), this->n_components);


  for (unsigned int i = 0; i < this->n_components; ++i)
    gradients[i] = gradient(p, i);
}

template <int dim, typename RangeNumberType, typename PointNumberType>
void
VectorFunctionFromTensorFunction<dim, RangeNumberType, PointNumberType>::
  gradient_list(const std::vector<Point<dim, PointNumberType>> &points,
                std::vector<Tensor<1, dim, RangeNumberType>>   &gradients,
                const unsigned int                              component) const
{
  Assert(gradients.size() == points.size(),
         ExcDimensionMismatch(gradients.size(), points.size()));


  for (unsigned int i = 0; i < points.size(); ++i)
    gradients[i] = gradient(points[i], component);
}



template <int dim, typename RangeNumberType, typename PointNumberType>
void
VectorFunctionFromTensorFunction<dim, RangeNumberType, PointNumberType>::
  vector_gradient_list(
    const std::vector<Point<dim, PointNumberType>>            &points,
    std::vector<std::vector<Tensor<1, dim, RangeNumberType>>> &gradients) const
{
  Assert(gradients.size() == points.size(),
         ExcDimensionMismatch(gradients.size(), points.size()));


  for (unsigned int i = 0; i < points.size(); ++i)
    {
      Assert(gradients[i].size() == this->n_components,
             ExcDimensionMismatch(gradients[i].size(), this->n_components));
      vector_gradient(points[i], gradients[i]);
    }
}

template <int dim, typename RangeNumberType, typename PointNumberType>
void
VectorFunctionFromTensorFunction<dim, RangeNumberType, PointNumberType>::
  vector_gradients(
    const std::vector<Point<dim, PointNumberType>>            &points,
    std::vector<std::vector<Tensor<1, dim, RangeNumberType>>> &gradients) const
{
  const unsigned int n = this->n_components;
  AssertDimension(gradients.size(), n);
  for (unsigned int i = 0; i < n; ++i)
    gradient_list(points, gradients[i], i);
}


template <int dim, typename RangeNumberType, typename PointNumberType>
FunctionFromFunctionObjects<dim, RangeNumberType, PointNumberType>::
  FunctionFromFunctionObjects(const unsigned int n_components,
                              const double       initial_time)
  : Function<dim, RangeNumberType, PointNumberType>(n_components, initial_time)
  , function_values(n_components)
  , function_gradients(n_components)
{}



template <int dim, typename RangeNumberType, typename PointNumberType>
FunctionFromFunctionObjects<dim, RangeNumberType, PointNumberType>::
  FunctionFromFunctionObjects(const std::vector<std::function<RangeNumberType(
                                const Point<dim, PointNumberType> &)>> &values,
                              const double initial_time)
  : Function<dim, RangeNumberType, PointNumberType>(values.size(), initial_time)
  , function_values(values)
  , function_gradients(values.size())
{}



template <int dim, typename RangeNumberType, typename PointNumberType>
FunctionFromFunctionObjects<dim, RangeNumberType, PointNumberType>::
  FunctionFromFunctionObjects(
    const std::vector<
      std::function<RangeNumberType(const Point<dim, PointNumberType> &)>>
                                             &values,
    const std::vector<std::function<Tensor<1, dim, RangeNumberType>(
      const Point<dim, PointNumberType> &)>> &gradients,
    const double                              initial_time)
  : Function<dim, RangeNumberType, PointNumberType>(values.size(), initial_time)
  , function_values(values)
  , function_gradients(gradients)
{}



template <int dim, typename RangeNumberType, typename PointNumberType>
RangeNumberType
FunctionFromFunctionObjects<dim, RangeNumberType, PointNumberType>::value(
  const Point<dim, PointNumberType> &p,
  const unsigned int                 component) const
{
  AssertIndexRange(component, this->n_components);
  Assert(function_values[component],
         ExcMessage("Accessing value() in FunctionFromFunctionObjects requires "
                    "setting the std::function objects for the value"));
  return function_values[component](p);
}



template <int dim, typename RangeNumberType, typename PointNumberType>
Tensor<1, dim, RangeNumberType>
FunctionFromFunctionObjects<dim, RangeNumberType, PointNumberType>::gradient(
  const Point<dim, PointNumberType> &p,
  const unsigned int                 component) const
{
  AssertIndexRange(component, this->n_components);
  Assert(function_gradients[component],
         ExcMessage(
           "Accessing gradient() in FunctionFromFunctionObjects "
           "requires setting the std::function objects for the gradient"));
  return function_gradients[component](p);
}



template <int dim, typename RangeNumberType, typename PointNumberType>
void
FunctionFromFunctionObjects<dim, RangeNumberType, PointNumberType>::
  set_function_values(const std::vector<std::function<RangeNumberType(
                        const Point<dim, PointNumberType> &)>> &values)
{
  AssertDimension(this->n_components, values.size());
  function_values = values;
}



template <int dim, typename RangeNumberType, typename PointNumberType>
void
FunctionFromFunctionObjects<dim, RangeNumberType, PointNumberType>::
  set_function_gradients(
    const std::vector<std::function<Tensor<1, dim, RangeNumberType>(
      const Point<dim, PointNumberType> &)>> &gradients)
{
  AssertDimension(this->n_components, gradients.size());
  function_gradients = gradients;
}

DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_function_templates_h */
