/* $Id$ */

#include <fe/quadrature.h>


template <int dim>
Quadrature<dim>::Quadrature (const unsigned int n_q) :
		n_quadrature_points(n_q),
		quadrature_points (n_q, Point<dim>()),
		weights (n_q, 0) {};



template <int dim>
Quadrature<dim>::~Quadrature () {};



template <int dim>
const Point<dim> & Quadrature<dim>::quad_point (const unsigned int i) const {
  Assert (i<n_quadrature_points, ExcInvalidIndex(i, n_quadrature_points));
  return quadrature_points[i];
};



template <int dim>
double Quadrature<dim>::weight (const unsigned int i) const {
  Assert (i<n_quadrature_points, ExcInvalidIndex(i, n_quadrature_points));
  return weights[i];
};





// explicite instantiations
template class Quadrature<1>;
template class Quadrature<2>;
