/*      $Id$                 */


#include <base/function.h>
#include <base/point.h>
#include <lac/vector.h>
#include <vector>


template <int dim>
Function<dim>::Function (const unsigned int n_components,
			 const double       initial_time) :
		FunctionTime(initial_time),
		n_components(n_components)
{};



template <int dim>
Function<dim>::~Function ()
{};



template <int dim>
double Function<dim>::value (const Point<dim> &,
			     const unsigned int) const
{
  Assert (false, ExcPureFunctionCalled());
  return 0;
};



template <int dim>
void Function<dim>::vector_value (const Point<dim> &,
				  Vector<double>   &) const
{
  Assert (false, ExcPureFunctionCalled());
};



template <int dim>
void Function<dim>::value_list (const vector<Point<dim> > &points,
				vector<double>            &values,
				const unsigned int         component) const
{
				   // check whether component is in
				   // the valid range is up to the
				   // derived class
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    values[i]  = this->value (points[i], component);
};



template <int dim>
void Function<dim>::vector_value_list (const vector<Point<dim> > &points,
				       vector<Vector<double> >   &values) const
{
				   // check whether component is in
				   // the valid range is up to the
				   // derived class
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    this->vector_value (points[i], values[i]);
};




template <int dim>
Tensor<1,dim> Function<dim>::gradient (const Point<dim> &,
				       const unsigned int) const
{
  Assert (false, ExcPureFunctionCalled());
  return Point<dim>();
};



template <int dim>
void Function<dim>::vector_gradient (const Point<dim>       &,
				     vector<Tensor<1,dim> > &) const
{
  Assert (false, ExcPureFunctionCalled());
};



template <int dim>
void Function<dim>::gradient_list (const vector<Point<dim> > &points,
				   vector<Tensor<1,dim> >    &gradients,
				   const unsigned int         component) const
{
  Assert (gradients.size() == points.size(),
	  ExcVectorHasWrongSize(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    gradients[i] = gradient(points[i], component);
};



template <int dim>
void Function<dim>::vector_gradient_list (const vector<Point<dim> >       &points,
					  vector<vector<Tensor<1,dim> > > &gradients) const
{
  Assert (gradients.size() == points.size(),
	  ExcVectorHasWrongSize(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    {
      Assert (gradients[i].size() == n_components,
	      ExcVectorHasWrongSize(gradients[i].size(), n_components));
      for (unsigned int component=0; component<n_components; ++component)
	gradients[i][component] = gradient(points[i], component);
    };
};





template <int dim>
ZeroFunction<dim>::ZeroFunction (const unsigned int n_components) :
		Function<dim> (n_components)
{};



template <int dim>
ZeroFunction<dim>::~ZeroFunction () {};



template <int dim>
double ZeroFunction<dim>::value (const Point<dim> &,
				       const unsigned int) const
{
  return 0.;
};



template <int dim>
void ZeroFunction<dim>::vector_value (const Point<dim> &,
				      Vector<double>   &return_value) const
{
  Assert (return_value.size() == n_components,
	  ExcVectorHasWrongSize (return_value.size(), n_components));

  fill_n (return_value.begin(), n_components, 0.0);
};



template <int dim>
void ZeroFunction<dim>::value_list (const vector<Point<dim> > &points,
				    vector<double>            &values,
				    const unsigned int         /*component*/) const {
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  fill_n (values.begin(), points.size(), 0.);
};



template <int dim>
void ZeroFunction<dim>::vector_value_list (const vector<Point<dim> > &points,
					   vector<Vector<double> >   &values) const
{
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    {
      Assert (values[i].size() == n_components,
	      ExcVectorHasWrongSize(values[i].size(), n_components));
      fill_n (values[i].begin(), n_components, 0.);
    };
};



template <int dim>
Tensor<1,dim> ZeroFunction<dim>::gradient (const Point<dim> &,
					   const unsigned int) const
{
  return Tensor<1,dim>();
};



template <int dim>
void ZeroFunction<dim>::vector_gradient (const Point<dim>       &,
					 vector<Tensor<1,dim> > &gradients) const
{
  Assert (gradients.size() == n_components,
	  ExcVectorHasWrongSize(gradients.size(), n_components));

  for (unsigned int c=0; c<n_components; ++c)
    gradients[c].clear ();
};



template <int dim>
void ZeroFunction<dim>::gradient_list (const vector<Point<dim> > &points,
				       vector<Tensor<1,dim> >    &gradients,
				       const unsigned int         /*component*/) const
{
  Assert (gradients.size() == points.size(),
	  ExcVectorHasWrongSize(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    gradients[i].clear ();
};



template <int dim>
void ZeroFunction<dim>::vector_gradient_list (const vector<Point<dim> >       &points,
					      vector<vector<Tensor<1,dim> > > &gradients) const
{
  Assert (gradients.size() == points.size(),
	  ExcVectorHasWrongSize(gradients.size(), points.size()));
  for (unsigned int i=0; i<points.size(); ++i)
    {
      Assert (gradients[i].size() == n_components,
	      ExcVectorHasWrongSize(gradients[i].size(), n_components));
      for (unsigned int c=0; c<n_components; ++c)
	gradients[i][c].clear ();
    };
};





template <int dim>
ConstantFunction<dim>::ConstantFunction (const double value,
					 const unsigned int n_components) :
		ZeroFunction<dim> (n_components),
		function_value    (value) {};


template <int dim>
ConstantFunction<dim>::~ConstantFunction () {};



template <int dim>
double ConstantFunction<dim>::value (const Point<dim> &,
					   const unsigned int) const
{
  return function_value;
};



template <int dim>
void ConstantFunction<dim>::vector_value (const Point<dim> &,
					  Vector<double>   &return_value) const
{
  Assert (return_value.size() == n_components,
	  ExcVectorHasWrongSize (return_value.size(), n_components));

  fill_n (return_value.begin(), n_components, function_value);
};



template <int dim>
void ConstantFunction<dim>::value_list (const vector<Point<dim> > &points,
					vector<double>            &values,
					const unsigned int         /*component*/) const {
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  fill_n (values.begin(), points.size(), function_value);
};



template <int dim>
void ConstantFunction<dim>::vector_value_list (const vector<Point<dim> > &points,
					       vector<Vector<double> >   &values) const
{
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    {
      Assert (values[i].size() == n_components,
	      ExcVectorHasWrongSize(values[i].size(), n_components));
      fill_n (values[i].begin(), n_components, function_value);
    };
};




// explicit instantiations

template class Function<1>;
template class ZeroFunction<1>;
template class ConstantFunction<1>;
template class Function<2>;
template class ZeroFunction<2>;
template class ConstantFunction<2>;
template class Function<3>;
template class ZeroFunction<3>;
template class ConstantFunction<3>;
