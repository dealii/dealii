/*      $Id$                 */


#include <basic/function.h>
#include <vector.h>



template <int dim>
Function<dim>::~Function () {};


template <int dim>
double Function<dim>::operator () (const Point<dim> &) const {
  Assert (false, ExcPureFunctionCalled());
  return 0;
};



template <int dim>
void Function<dim>::value_list (const vector<Point<dim> > &points,
				vector<double>            &values) const {
  Assert (values.size() == 0,
	  ExcVectorNotEmpty());

  values.reserve (points.size());
  for (unsigned int i=0; i<points.size(); ++i)
    values.push_back (this->operator() (points[i]));
};




template <int dim>
Point<dim> Function<dim>::gradient (const Point<dim> &) const {
  Assert (false, ExcPureFunctionCalled());
  return Point<dim>();
};



template <int dim>
void Function<dim>::gradient_list (const vector<Point<dim> > &points,
				   vector<Point<dim> >       &gradients) const {
  Assert (gradients.size() == 0,
	  ExcVectorNotEmpty());

  gradients.reserve (points.size());
  for (unsigned int i=0; i<points.size(); ++i)
    gradients.push_back (gradient(points[i]));
};







template <int dim>
ZeroFunction<dim>::~ZeroFunction () {};



template <int dim>
double ZeroFunction<dim>::operator () (const Point<dim> &) const {
  return 0.;
};



template <int dim>
void ZeroFunction<dim>::value_list (const vector<Point<dim> > &points,
				    vector<double>            &values) const {
  Assert (values.size() == 0,
	  ExcVectorNotEmpty());

  values.reserve (points.size());
  values.insert (values.begin(), points.size(), 0.);
};



template <int dim>
Point<dim> ZeroFunction<dim>::gradient (const Point<dim> &) const {
  return Point<dim>();
};



template <int dim>
void ZeroFunction<dim>::gradient_list (const vector<Point<dim> > &points,
				       vector<Point<dim> >       &values) const {
  Assert (values.size() == 0,
	  ExcVectorNotEmpty());

  values.reserve (points.size());
  values.insert (values.begin(), points.size(), Point<dim>());
};




template <int dim>
ConstantFunction<dim>::ConstantFunction (const double value) :
		function_value(value) {};


template <int dim>
ConstantFunction<dim>::~ConstantFunction () {};



template <int dim>
double ConstantFunction<dim>::operator () (const Point<dim> &) const {
  return function_value;
};



template <int dim>
void ConstantFunction<dim>::value_list (const vector<Point<dim> > &points,
				    vector<double>            &values) const {
  Assert (values.size() == 0,
	  ExcVectorNotEmpty());

  values.reserve (points.size());
  values.insert (values.begin(), points.size(), function_value);
};



// explicit instantiations
template class Function<1>;
template class Function<2>;

template class ZeroFunction<1>;
template class ZeroFunction<2>;

template class ConstantFunction<1>;
template class ConstantFunction<2>;
