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





// explicit instantiations
template class Function<1>;
template class Function<2>;
