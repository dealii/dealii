// ---------------------------------------------------------------------
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


template <int dim, typename Number>
const unsigned int Function<dim, Number>::dimension;


template <int dim, typename Number>
Function<dim, Number>::Function (const unsigned int n_components,
                                 const Number       initial_time)
  :
  FunctionTime<Number>(initial_time),
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


template <int dim, typename Number>
Function<dim, Number>::~Function ()
{}



template <int dim, typename Number>
Function<dim, Number> &Function<dim, Number>::operator= (const Function &f)
{
  AssertDimension (n_components, f.n_components);
  return *this;
}


template <int dim, typename Number>
Number Function<dim, Number>::value (const Point<dim, Number> &,
                                     const unsigned int) const
{
  Assert (false, ExcPureFunctionCalled());
  return 0;
}


template <int dim, typename Number>
void Function<dim, Number>::vector_value (const Point<dim, Number> &p,
                                          Vector<Number> &v) const
{
  AssertDimension(v.size(), this->n_components);
  for (unsigned int i=0; i<this->n_components; ++i)
    v(i) = value(p, i);
}


template <int dim, typename Number>
void Function<dim, Number>::value_list (const std::vector<Point<dim, Number> > &points,
                                        std::vector<Number>                     &values,
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


template <int dim, typename Number>
void Function<dim, Number>::vector_value_list (const std::vector<Point<dim, Number> > &points,
                                               std::vector<Vector<Number> >   &values) const
{
  // check whether component is in
  // the valid range is up to the
  // derived class
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    this->vector_value (points[i], values[i]);
}


template <int dim, typename Number>
void Function<dim, Number>::vector_values (
  const std::vector<Point<dim, Number> > &points,
  std::vector<std::vector<Number> > &values) const
{
  const unsigned int n = this->n_components;
  AssertDimension (values.size(), n);
  for (unsigned int i=0; i<n; ++i)
    value_list(points, values[i], i);
}


template <int dim, typename Number>
Tensor<1,dim,Number> Function<dim, Number>::gradient (const Point<dim, Number> &,
                                                      const unsigned int) const
{
  Assert (false, ExcPureFunctionCalled());
  return Point<dim, Number>();
}


template <int dim, typename Number>
void Function<dim, Number>::vector_gradient (
  const Point<dim, Number> &p,
  std::vector<Tensor<1,dim,Number> > &v) const
{
  AssertDimension(v.size(), this->n_components);
  for (unsigned int i=0; i<this->n_components; ++i)
    v[i] = gradient(p, i);
}


template <int dim, typename Number>
void Function<dim, Number>::gradient_list (const std::vector<Point<dim, Number> > &points,
                                           std::vector<Tensor<1,dim,Number> >    &gradients,
                                           const unsigned int              component) const
{
  Assert (gradients.size() == points.size(),
          ExcDimensionMismatch(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    gradients[i] = gradient(points[i], component);
}


template <int dim, typename Number>
void Function<dim, Number>::vector_gradient_list (const std::vector<Point<dim, Number> >            &points,
                                                  std::vector<std::vector<Tensor<1,dim,Number> > > &gradients) const
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


template <int dim, typename Number>
void Function<dim, Number>::vector_gradients (
  const std::vector<Point<dim, Number> > &points,
  std::vector<std::vector<Tensor<1,dim,Number> > > &values) const
{
  const unsigned int n = this->n_components;
  AssertDimension (values.size(), n);
  for (unsigned int i=0; i<n; ++i)
    gradient_list(points, values[i], i);
}



template <int dim, typename Number>
Number Function<dim, Number>::laplacian (const Point<dim, Number> &,
                                         const unsigned int) const
{
  Assert (false, ExcPureFunctionCalled());
  return 0;
}


template <int dim, typename Number>
void Function<dim, Number>::vector_laplacian (const Point<dim, Number> &,
                                              Vector<Number> &) const
{
  Assert (false, ExcPureFunctionCalled());
}



template <int dim, typename Number>
void Function<dim, Number>::laplacian_list (const std::vector<Point<dim, Number> > &points,
                                            std::vector<Number>            &laplacians,
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


template <int dim, typename Number>
void Function<dim, Number>::vector_laplacian_list (const std::vector<Point<dim, Number> > &points,
                                                   std::vector<Vector<Number> >   &laplacians) const
{
  // check whether component is in
  // the valid range is up to the
  // derived class
  Assert (laplacians.size() == points.size(),
          ExcDimensionMismatch(laplacians.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    this->vector_laplacian (points[i], laplacians[i]);
}



template <int dim, typename Number>
std::size_t
Function<dim, Number>::memory_consumption () const
{
  // only simple data elements, so
  // use sizeof operator
  return sizeof (*this);
}


//---------------------------------------------------------------------------

template <int dim, typename Number>
ZeroFunction<dim, Number>::ZeroFunction (const unsigned int n_components)
  :
  Function<dim, Number> (n_components)
{}


template <int dim, typename Number>
ZeroFunction<dim, Number>::~ZeroFunction ()
{}


template <int dim, typename Number>
Number ZeroFunction<dim, Number>::value (const Point<dim, Number> &,
                                         const unsigned int) const
{
  return 0.;
}


template <int dim, typename Number>
void ZeroFunction<dim, Number>::vector_value (const Point<dim, Number> &,
                                              Vector<Number>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
          ExcDimensionMismatch (return_value.size(), this->n_components));

  std::fill (return_value.begin(), return_value.end(), 0.0);
}


template <int dim, typename Number>
void ZeroFunction<dim, Number>::value_list (const std::vector<Point<dim, Number> > &points,
                                            std::vector<Number>            &values,
                                            const unsigned int         /*component*/) const
{
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  std::fill (values.begin(), values.end(), 0.);
}


template <int dim, typename Number>
void ZeroFunction<dim, Number>::vector_value_list (const std::vector<Point<dim, Number> > &points,
                                                   std::vector<Vector<Number> >   &values) const
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


template <int dim, typename Number>
Tensor<1,dim,Number> ZeroFunction<dim, Number>::gradient (const Point<dim, Number> &,
                                                          const unsigned int) const
{
  return Tensor<1,dim,Number>();
}


template <int dim, typename Number>
void ZeroFunction<dim, Number>::vector_gradient (const Point<dim, Number> &,
                                                 std::vector<Tensor<1,dim,Number> > &gradients) const
{
  Assert (gradients.size() == this->n_components,
          ExcDimensionMismatch(gradients.size(), this->n_components));

  for (unsigned int c=0; c<this->n_components; ++c)
    gradients[c].clear ();
}


template <int dim, typename Number>
void ZeroFunction<dim, Number>::gradient_list (const std::vector<Point<dim, Number> > &points,
                                               std::vector<Tensor<1,dim,Number> >    &gradients,
                                               const unsigned int              /*component*/) const
{
  Assert (gradients.size() == points.size(),
          ExcDimensionMismatch(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    gradients[i].clear ();
}


template <int dim, typename Number>
void ZeroFunction<dim, Number>::vector_gradient_list (const std::vector<Point<dim, Number> >            &points,
                                                      std::vector<std::vector<Tensor<1,dim,Number> > > &gradients) const
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


template <int dim, typename Number>
ConstantFunction<dim, Number>::ConstantFunction (const Number value,
                                                 const unsigned int n_components)
  :
  ZeroFunction<dim, Number> (n_components),
  function_value    (value)
{}



template <int dim, typename Number>
ConstantFunction<dim, Number>::~ConstantFunction () {}



template <int dim, typename Number>
Number ConstantFunction<dim, Number>::value (const Point<dim, Number> &,
                                             const unsigned int component) const
{
  Assert (component < this->n_components,
          ExcIndexRange (component, 0, this->n_components));
  return function_value;
}



template <int dim, typename Number>
void ConstantFunction<dim, Number>::vector_value (const Point<dim, Number> &,
                                                  Vector<Number>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
          ExcDimensionMismatch (return_value.size(), this->n_components));

  std::fill (return_value.begin(), return_value.end(), function_value);
}



template <int dim, typename Number>
void ConstantFunction<dim, Number>::value_list (const std::vector<Point<dim, Number> > &points,
                                                std::vector<Number>            &values,
                                                const unsigned int              component) const
{
  Assert (component < this->n_components,
          ExcIndexRange (component, 0, this->n_components));
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  std::fill (values.begin(), values.end(), function_value);
}



template <int dim, typename Number>
void ConstantFunction<dim, Number>::vector_value_list (const std::vector<Point<dim, Number> > &points,
                                                       std::vector<Vector<Number> >   &values) const
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



template <int dim, typename Number>
std::size_t
ConstantFunction<dim, Number>::memory_consumption () const
{
  // only simple data elements, so
  // use sizeof operator
  return sizeof (*this);
}

//---------------------------------------------------------------------------

template <int dim, typename Number>
ComponentSelectFunction<dim, Number>::
ComponentSelectFunction (const unsigned int selected,
                         const Number value,
                         const unsigned int n_components)
  :
  ConstantFunction<dim, Number> (value, n_components),
  selected_components(std::make_pair(selected,selected+1))
{}



template <int dim, typename Number>
ComponentSelectFunction<dim, Number>::
ComponentSelectFunction (const unsigned int selected,
                         const unsigned int n_components)
  :
  ConstantFunction<dim, Number> (1., n_components),
  selected_components(std::make_pair(selected,selected+1))
{
  Assert (selected < n_components,
          ExcIndexRange (selected, 0, n_components));
}



template <int dim, typename Number>
ComponentSelectFunction<dim, Number>::
ComponentSelectFunction (const std::pair<unsigned int,unsigned int> &selected,
                         const unsigned int n_components)
  :
  ConstantFunction<dim, Number> (1., n_components),
  selected_components(selected)
{
  Assert (selected_components.first < selected_components.second,
          ExcMessage ("The upper bound of the interval must be larger than "
                      "the lower bound"));
  Assert (selected_components.second <= n_components,
          ExcMessage ("The upper bound of the interval must be less than "
                      "or equal to the total number of vector components"));
}



template <int dim, typename Number>
void ComponentSelectFunction<dim, Number>::vector_value (const Point<dim, Number> &,
                                                         Vector<Number>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
          ExcDimensionMismatch (return_value.size(), this->n_components));

  return_value = 0;
  std::fill (return_value.begin()+selected_components.first,
             return_value.begin()+selected_components.second,
             this->function_value);
}



template <int dim, typename Number>
void ComponentSelectFunction<dim, Number>::vector_value_list (const std::vector<Point<dim, Number> > &points,
    std::vector<Vector<Number> >   &values) const
{
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    ComponentSelectFunction<dim, Number>::vector_value (points[i],
                                                        values[i]);
}



template <int dim, typename Number>
std::size_t
ComponentSelectFunction<dim, Number>::memory_consumption () const
{
  // only simple data elements, so
  // use sizeof operator
  return sizeof (*this);
}


template <int dim, typename Number>
ScalarFunctionFromFunctionObject<dim, Number>::
ScalarFunctionFromFunctionObject (const std_cxx11::function<Number (const Point<dim, Number> &)> &function_object)
  :
  Function<dim, Number>(1),
  function_object (function_object)
{}



template <int dim, typename Number>
Number
ScalarFunctionFromFunctionObject<dim, Number>::value (const Point<dim, Number> &p,
                                                      const unsigned int component) const
{
  Assert (component == 0,
          ExcMessage ("This object represents only scalar functions"));
  return function_object (p);
}



template <int dim, typename Number>
VectorFunctionFromScalarFunctionObject<dim, Number>::
VectorFunctionFromScalarFunctionObject (const std_cxx11::function<Number (const Point<dim, Number> &)> &function_object,
                                        const unsigned int selected_component,
                                        const unsigned int n_components)
  :
  Function<dim, Number>(n_components),
  function_object (function_object),
  selected_component (selected_component)
{
  Assert (selected_component < this->n_components,
          ExcIndexRange (selected_component, 0, this->n_components));
}



template <int dim, typename Number>
Number
VectorFunctionFromScalarFunctionObject<dim, Number>::value (const Point<dim, Number> &p,
                                                            const unsigned int component) const
{
  Assert (component < this->n_components,
          ExcIndexRange (component, 0, this->n_components));

  if (component == selected_component)
    return function_object (p);
  else
    return 0;
}



template <int dim, typename Number>
void
VectorFunctionFromScalarFunctionObject<dim, Number>::
vector_value (const Point<dim, Number>   &p,
              Vector<Number>     &values) const
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
template <int dim, typename Number>
VectorFunctionFromTensorFunction<dim, Number>::VectorFunctionFromTensorFunction (const TensorFunction<1,dim,Number> &tensor_function,
    const unsigned int selected_component,
    const unsigned int n_components)
  :
  Function<dim, Number> (n_components),
  tensor_function (tensor_function),
  selected_component (selected_component)
{

  // Verify that the Tensor<1,dim,Number> will fit in the given length selected_components
  // and not hang over the end of the vector.
  Assert (selected_component + dim - 1 < this->n_components,
          ExcIndexRange (selected_component, 0, this->n_components));
}


template <int dim, typename Number>
VectorFunctionFromTensorFunction<dim, Number>::~VectorFunctionFromTensorFunction ()
{}


template <int dim, typename Number>
inline
Number VectorFunctionFromTensorFunction<dim, Number>::value (const Point<dim, Number> &p,
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
  const Tensor<1,dim,Number> tensor_value = tensor_function.value (p);

  return tensor_value[component-selected_component];
}


template <int dim, typename Number>
inline
void VectorFunctionFromTensorFunction<dim, Number>::vector_value (const Point<dim, Number> &p,
    Vector<Number>   &values) const
{
  Assert(values.size() == this->n_components,
         ExcDimensionMismatch(values.size(),this->n_components));

  // Retrieve the values from the <tt>tensor_function</tt> to be
  // placed at the <tt>selected_component</tt> to <tt>selected_component + dim - 1</tt>
  // elements of the <tt>Vector</tt> values.
  const Tensor<1,dim,Number> tensor_value = tensor_function.value (p);

  // First we make all elements of values = 0
  values = 0;

  // Second we adjust the desired components to take on the values in <tt>tensor_value</tt>.
  for (unsigned int i=0; i<dim; ++i)
    values(i+selected_component) = tensor_value[i];
}


/** Member function <tt>vector_value_list </tt> is the interface for giving a list of
 * points (<code>vector<Point<dim, Number> ></code>) of which to evaluate using the <tt>vector_value</tt>
 * member function.  Again, this function is written so as to not replicate the function definition
 * but passes each point on to <tt>vector_value</tt> to be evaluated.
 */
template <int dim, typename Number>
void VectorFunctionFromTensorFunction<dim, Number>::vector_value_list (const std::vector<Point<dim, Number> > &points,
    std::vector<Vector<Number> > &value_list) const
{
  Assert (value_list.size() == points.size(),
          ExcDimensionMismatch (value_list.size(), points.size()));

  const unsigned int n_points = points.size();

  for (unsigned int p=0; p<n_points; ++p)
    VectorFunctionFromTensorFunction<dim, Number>::vector_value (points[p], value_list[p]);
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
