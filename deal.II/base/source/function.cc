/*      $Id$                 */


#include <basic/function.h>
#include <vector.h>


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




// explicit instantiations
template class Function<1>;
template class Function<2>;

template class ZeroFunction<1>;
template class ZeroFunction<2>;
