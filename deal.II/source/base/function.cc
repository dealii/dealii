// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <vector>

DEAL_II_NAMESPACE_OPEN


template <int dim>
const unsigned int Function<dim>::dimension;


template <int dim>
Function<dim>::Function (const unsigned int n_components,
                         const double       initial_time)
  :
  FunctionTime(initial_time),
  n_components(n_components)
{
  // avoid the construction of function
  // objects that don't return any
  // values. This doesn't make much sense in
  // the first place, but will lead to odd
  // errors later on (happened to me in fact
  // :-)
  Assert (n_components > 0,
          ExcZero());
}


template <int dim>
Function<dim>::~Function ()
{}



template <int dim>
Function<dim> &Function<dim>::operator= (const Function &f)
{
  AssertDimension (n_components, f.n_components);
  return *this;
}


template <int dim>
double Function<dim>::value (const Point<dim> &,
                             const unsigned int) const
{
  Assert (false, ExcPureFunctionCalled());
  return 0;
}


template <int dim>
void Function<dim>::vector_value (const Point<dim> &p,
                                  Vector<double> &v) const
{
  AssertDimension(v.size(), this->n_components);
  for (unsigned int i=0; i<this->n_components; ++i)
    v(i) = value(p, i);
}


template <int dim>
void Function<dim>::value_list (const std::vector<Point<dim> > &points,
                                std::vector<double>                     &values,
                                const unsigned int                       component) const
{
  // check whether component is in
  // the valid range is up to the
  // derived class
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    values[i]  = this->value (points[i], component);
}


template <int dim>
void Function<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                       std::vector<Vector<double> >   &values) const
{
  // check whether component is in
  // the valid range is up to the
  // derived class
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    this->vector_value (points[i], values[i]);
}


template <int dim>
void Function<dim>::vector_values (
  const std::vector<Point<dim> > &points,
  std::vector<std::vector<double> > &values) const
{
  const unsigned int n = this->n_components;
  AssertDimension (values.size(), n);
  for (unsigned int i=0; i<n; ++i)
    value_list(points, values[i], i);
}


template <int dim>
Tensor<1,dim> Function<dim>::gradient (const Point<dim> &,
                                       const unsigned int) const
{
  Assert (false, ExcPureFunctionCalled());
  return Point<dim>();
}


template <int dim>
void Function<dim>::vector_gradient (
  const Point<dim> &p,
  std::vector<Tensor<1,dim> > &v) const
{
  AssertDimension(v.size(), this->n_components);
  for (unsigned int i=0; i<this->n_components; ++i)
    v[i] = gradient(p, i);
}


template <int dim>
void Function<dim>::gradient_list (const std::vector<Point<dim> > &points,
                                   std::vector<Tensor<1,dim> >    &gradients,
                                   const unsigned int              component) const
{
  Assert (gradients.size() == points.size(),
          ExcDimensionMismatch(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    gradients[i] = gradient(points[i], component);
}


template <int dim>
void Function<dim>::vector_gradient_list (const std::vector<Point<dim> >            &points,
                                          std::vector<std::vector<Tensor<1,dim> > > &gradients) const
{
  Assert (gradients.size() == points.size(),
          ExcDimensionMismatch(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    {
      Assert (gradients[i].size() == n_components,
              ExcDimensionMismatch(gradients[i].size(), n_components));
      vector_gradient (points[i], gradients[i]);
    }
}


template <int dim>
void Function<dim>::vector_gradients (
  const std::vector<Point<dim> > &points,
  std::vector<std::vector<Tensor<1,dim> > > &values) const
{
  const unsigned int n = this->n_components;
  AssertDimension (values.size(), n);
  for (unsigned int i=0; i<n; ++i)
    gradient_list(points, values[i], i);
}



template <int dim>
double Function<dim>::laplacian (const Point<dim> &,
                                 const unsigned int) const
{
  Assert (false, ExcPureFunctionCalled());
  return 0;
}


template <int dim>
void Function<dim>::vector_laplacian (const Point<dim> &,
                                      Vector<double> &) const
{
  Assert (false, ExcPureFunctionCalled());
}



template <int dim>
void Function<dim>::laplacian_list (const std::vector<Point<dim> > &points,
                                    std::vector<double>            &laplacians,
                                    const unsigned int              component) const
{
  // check whether component is in
  // the valid range is up to the
  // derived class
  Assert (laplacians.size() == points.size(),
          ExcDimensionMismatch(laplacians.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    laplacians[i]  = this->laplacian (points[i], component);
}


template <int dim>
void Function<dim>::vector_laplacian_list (const std::vector<Point<dim> > &points,
                                           std::vector<Vector<double> >   &laplacians) const
{
  // check whether component is in
  // the valid range is up to the
  // derived class
  Assert (laplacians.size() == points.size(),
          ExcDimensionMismatch(laplacians.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    this->vector_laplacian (points[i], laplacians[i]);
}



template <int dim>
std::size_t
Function<dim>::memory_consumption () const
{
  // only simple data elements, so
  // use sizeof operator
  return sizeof (*this);
}


//---------------------------------------------------------------------------

template <int dim>
ZeroFunction<dim>::ZeroFunction (const unsigned int n_components)
  :
  Function<dim> (n_components)
{}


template <int dim>
ZeroFunction<dim>::~ZeroFunction ()
{}


template <int dim>
double ZeroFunction<dim>::value (const Point<dim> &,
                                 const unsigned int) const
{
  return 0.;
}


template <int dim>
void ZeroFunction<dim>::vector_value (const Point<dim> &,
                                      Vector<double>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
          ExcDimensionMismatch (return_value.size(), this->n_components));

  std::fill (return_value.begin(), return_value.end(), 0.0);
}


template <int dim>
void ZeroFunction<dim>::value_list (const std::vector<Point<dim> > &points,
                                    std::vector<double>            &values,
                                    const unsigned int         /*component*/) const
{
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  std::fill (values.begin(), values.end(), 0.);
}


template <int dim>
void ZeroFunction<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                           std::vector<Vector<double> >   &values) const
{
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    {
      Assert (values[i].size() == this->n_components,
              ExcDimensionMismatch(values[i].size(), this->n_components));
      std::fill (values[i].begin(), values[i].end(), 0.);
    };
}


template <int dim>
Tensor<1,dim> ZeroFunction<dim>::gradient (const Point<dim> &,
                                           const unsigned int) const
{
  return Tensor<1,dim>();
}


template <int dim>
void ZeroFunction<dim>::vector_gradient (const Point<dim> &,
                                         std::vector<Tensor<1,dim> > &gradients) const
{
  Assert (gradients.size() == this->n_components,
          ExcDimensionMismatch(gradients.size(), this->n_components));

  for (unsigned int c=0; c<this->n_components; ++c)
    gradients[c].clear ();
}


template <int dim>
void ZeroFunction<dim>::gradient_list (const std::vector<Point<dim> > &points,
                                       std::vector<Tensor<1,dim> >    &gradients,
                                       const unsigned int              /*component*/) const
{
  Assert (gradients.size() == points.size(),
          ExcDimensionMismatch(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    gradients[i].clear ();
}


template <int dim>
void ZeroFunction<dim>::vector_gradient_list (const std::vector<Point<dim> >            &points,
                                              std::vector<std::vector<Tensor<1,dim> > > &gradients) const
{
  Assert (gradients.size() == points.size(),
          ExcDimensionMismatch(gradients.size(), points.size()));
  for (unsigned int i=0; i<points.size(); ++i)
    {
      Assert (gradients[i].size() == this->n_components,
              ExcDimensionMismatch(gradients[i].size(), this->n_components));
      for (unsigned int c=0; c<this->n_components; ++c)
        gradients[i][c].clear ();
    };
}

//---------------------------------------------------------------------------


template <int dim>
ConstantFunction<dim>::ConstantFunction (const double value,
                                         const unsigned int n_components)
  :
  ZeroFunction<dim> (n_components),
  function_value    (value)
{}



template <int dim>
ConstantFunction<dim>::~ConstantFunction () {}



template <int dim>
double ConstantFunction<dim>::value (const Point<dim> &,
                                     const unsigned int component) const
{
  Assert (component < this->n_components,
          ExcIndexRange (component, 0, this->n_components));
  return function_value;
}



template <int dim>
void ConstantFunction<dim>::vector_value (const Point<dim> &,
                                          Vector<double>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
          ExcDimensionMismatch (return_value.size(), this->n_components));

  std::fill (return_value.begin(), return_value.end(), function_value);
}



template <int dim>
void ConstantFunction<dim>::value_list (const std::vector<Point<dim> > &points,
                                        std::vector<double>            &values,
                                        const unsigned int              component) const
{
  Assert (component < this->n_components,
          ExcIndexRange (component, 0, this->n_components));
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  std::fill (values.begin(), values.end(), function_value);
}



template <int dim>
void ConstantFunction<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                               std::vector<Vector<double> >   &values) const
{
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    {
      Assert (values[i].size() == this->n_components,
              ExcDimensionMismatch(values[i].size(), this->n_components));
      std::fill (values[i].begin(), values[i].end(), function_value);
    };
}



template <int dim>
std::size_t
ConstantFunction<dim>::memory_consumption () const
{
  // only simple data elements, so
  // use sizeof operator
  return sizeof (*this);
}

//---------------------------------------------------------------------------

template <int dim>
ComponentSelectFunction<dim>::
ComponentSelectFunction (const unsigned int selected,
                         const double value,
                         const unsigned int n_components)
  :
  ConstantFunction<dim> (value, n_components),
  selected_components(std::make_pair(selected,selected+1))
{}



template <int dim>
ComponentSelectFunction<dim>::
ComponentSelectFunction (const unsigned int selected,
                         const unsigned int n_components)
  :
  ConstantFunction<dim> (1., n_components),
  selected_components(std::make_pair(selected,selected+1))
{
  Assert (selected < n_components,
          ExcIndexRange (selected, 0, n_components));
}



template <int dim>
ComponentSelectFunction<dim>::
ComponentSelectFunction (const std::pair<unsigned int,unsigned int> &selected,
                         const unsigned int n_components)
  :
  ConstantFunction<dim> (1., n_components),
  selected_components(selected)
{
  Assert (selected_components.first < selected_components.second,
          ExcMessage ("The upper bound of the interval must be larger than "
                      "the lower bound"));
  Assert (selected_components.second <= n_components,
          ExcMessage ("The upper bound of the interval must be less than "
                      "or equal to the total number of vector components"));
}



template <int dim>
void ComponentSelectFunction<dim>::vector_value (const Point<dim> &,
                                                 Vector<double>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
          ExcDimensionMismatch (return_value.size(), this->n_components));

  return_value = 0;
  std::fill (return_value.begin()+selected_components.first,
             return_value.begin()+selected_components.second,
             this->function_value);
}



template <int dim>
void ComponentSelectFunction<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                                      std::vector<Vector<double> >   &values) const
{
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    ComponentSelectFunction<dim>::vector_value (points[i],
                                                values[i]);
}



template <int dim>
std::size_t
ComponentSelectFunction<dim>::memory_consumption () const
{
  // only simple data elements, so
  // use sizeof operator
  return sizeof (*this);
}


template <int dim>
ScalarFunctionFromFunctionObject<dim>::
ScalarFunctionFromFunctionObject (const std_cxx1x::function<double (const Point<dim> &)> &function_object)
  :
  Function<dim>(1),
  function_object (function_object)
{}



template <int dim>
double
ScalarFunctionFromFunctionObject<dim>::value (const Point<dim> &p,
                                              const unsigned int component) const
{
  Assert (component == 0,
          ExcMessage ("This object represents only scalar functions"));
  return function_object (p);
}



template <int dim>
VectorFunctionFromScalarFunctionObject<dim>::
VectorFunctionFromScalarFunctionObject (const std_cxx1x::function<double (const Point<dim> &)> &function_object,
                                        const unsigned int selected_component,
                                        const unsigned int n_components)
  :
  Function<dim>(n_components),
  function_object (function_object),
  selected_component (selected_component)
{
  Assert (selected_component < this->n_components,
          ExcIndexRange (selected_component, 0, this->n_components));
}



template <int dim>
double
VectorFunctionFromScalarFunctionObject<dim>::value (const Point<dim> &p,
                                                    const unsigned int component) const
{
  Assert (component < this->n_components,
          ExcIndexRange (component, 0, this->n_components));

  if (component == selected_component)
    return function_object (p);
  else
    return 0;
}



template <int dim>
void
VectorFunctionFromScalarFunctionObject<dim>::
vector_value (const Point<dim>   &p,
              Vector<double>     &values) const
{
  AssertDimension(values.size(), this->n_components);

  // set everything to zero, and then the right component to its correct value
  values = 0;
  values(selected_component) = function_object (p);
}



/**
 * The constructor for <tt>VectorFunctionFromTensorFunction</tt> which initiates the return vector
 * to be size <tt>n_components</tt>.
 */
template <int dim>
VectorFunctionFromTensorFunction<dim>::VectorFunctionFromTensorFunction (const TensorFunction<1,dim> &tensor_function,
    const unsigned int selected_component,
    const unsigned int n_components)
  :
  Function<dim> (n_components),
  tensor_function (tensor_function),
  selected_component (selected_component)
{

  // Verify that the Tensor<1,dim> will fit in the given length selected_components
  // and not hang over the end of the vector.
  Assert (selected_component + dim - 1 < this->n_components,
          ExcIndexRange (selected_component, 0, this->n_components));
}


template <int dim>
VectorFunctionFromTensorFunction<dim>::~VectorFunctionFromTensorFunction ()
{}


template <int dim>
inline
double VectorFunctionFromTensorFunction<dim>::value (const Point<dim> &p,
                                                     const unsigned int component) const
{
  Assert (component<this->n_components,
          ExcIndexRange(component, 0, this->n_components));

  // if the requested component is out of the range selected, then we
  // can return early
  if ((component < selected_component)
      ||
      (component >= selected_component+dim))
    return 0;

  // otherwise retrieve the values from the <tt>tensor_function</tt> to be
  // placed at the <tt>selected_component</tt> to <tt>selected_component + dim - 1</tt>
  // elements of the <tt>Vector</tt> values and pick the correct one
  const Tensor<1,dim> tensor_value = tensor_function.value (p);

  return tensor_value[component-selected_component];
}


template <int dim>
inline
void VectorFunctionFromTensorFunction<dim>::vector_value (const Point<dim> &p,
                                                          Vector<double>   &values) const
{
  Assert(values.size() == this->n_components,
         ExcDimensionMismatch(values.size(),this->n_components));

  // Retrieve the values from the <tt>tensor_function</tt> to be
  // placed at the <tt>selected_component</tt> to <tt>selected_component + dim - 1</tt>
  // elements of the <tt>Vector</tt> values.
  const Tensor<1,dim> tensor_value = tensor_function.value (p);

  // First we make all elements of values = 0
  values = 0;

  // Second we adjust the desired components to take on the values in <tt>tensor_value</tt>.
  for (unsigned int i=0; i<dim; ++i)
    values(i+selected_component) = tensor_value[i];
}


/** Member function <tt>vector_value_list </tt> is the interface for giving a list of
 * points (<code>vector<Point<dim> ></code>) of which to evaluate using the <tt>vector_value</tt>
 * member function.  Again, this function is written so as to not replicate the function definition
 * but passes each point on to <tt>vector_value</tt> to be evaluated.
 */
template <int dim>
void VectorFunctionFromTensorFunction<dim>::vector_value_list (const std::vector<Point<dim> > &points,
    std::vector<Vector<double> > &value_list) const
{
  Assert (value_list.size() == points.size(),
          ExcDimensionMismatch (value_list.size(), points.size()));

  const unsigned int n_points = points.size();

  for (unsigned int p=0; p<n_points; ++p)
    VectorFunctionFromTensorFunction<dim>::vector_value (points[p], value_list[p]);
}




// explicit instantiations

template class Function<1>;
template class ZeroFunction<1>;
template class ConstantFunction<1>;
template class ComponentSelectFunction<1>;
template class ScalarFunctionFromFunctionObject<1>;
template class VectorFunctionFromScalarFunctionObject<1>;
template class VectorFunctionFromTensorFunction<1>;

template class Function<2>;
template class ZeroFunction<2>;
template class ConstantFunction<2>;
template class ComponentSelectFunction<2>;
template class ScalarFunctionFromFunctionObject<2>;
template class VectorFunctionFromScalarFunctionObject<2>;
template class VectorFunctionFromTensorFunction<2>;

template class Function<3>;
template class ZeroFunction<3>;
template class ConstantFunction<3>;
template class ComponentSelectFunction<3>;
template class ScalarFunctionFromFunctionObject<3>;
template class VectorFunctionFromScalarFunctionObject<3>;
template class VectorFunctionFromTensorFunction<3>;

DEAL_II_NAMESPACE_CLOSE
