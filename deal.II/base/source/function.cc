/*      $Id$                 */


#include <basic/function.h>
#include <vector>


template <int dim>
Function<dim>::Function (const double initial_time) :
		time(initial_time) {};



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
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    values[i]  = this->operator() (points[i]);
};




template <int dim>
Point<dim> Function<dim>::gradient (const Point<dim> &) const {
  Assert (false, ExcPureFunctionCalled());
  return Point<dim>();
};



template <int dim>
void Function<dim>::gradient_list (const vector<Point<dim> > &points,
				   vector<Point<dim> >       &gradients) const {
  Assert (gradients.size() == points.size(),
	  ExcVectorHasWrongSize(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    gradients[i] = gradient(points[i]);
};



template <int dim>
double Function<dim>::get_time () const {
  return time;
};



template <int dim>
void Function<dim>::set_time (const double new_time) {
  time = new_time;
};



template <int dim>
void Function<dim>::advance_time (const double delta_t) {
  set_time (time+delta_t);
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
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  fill_n (values.begin(), points.size(), 0.);
};



template <int dim>
Point<dim> ZeroFunction<dim>::gradient (const Point<dim> &) const {
  return Point<dim>();
};



template <int dim>
void ZeroFunction<dim>::gradient_list (const vector<Point<dim> > &points,
				       vector<Point<dim> >       &gradients) const {
  Assert (gradients.size() == points.size(),
	  ExcVectorHasWrongSize(gradients.size(), points.size()));

  fill_n (gradients.begin(), points.size(), Point<dim>());
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
  Assert (values.size() == points.size(),
	  ExcVectorHasWrongSize(values.size(), points.size()));

  fill_n (values.begin(), points.size(), function_value);
};



// explicit instantiations

template class Function<deal_II_dimension>;
template class ZeroFunction<deal_II_dimension>;
template class ConstantFunction<deal_II_dimension>;
