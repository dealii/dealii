//----------------------------  quadrature_lib.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  quadrature_lib.cc  ---------------------------


#include <base/quadrature_lib.h>
#include <cmath>


// please note: for a given dimension, we need the quadrature formulae
// for all lower dimensions as well. That is why in this file the check
// is for deal_II_dimension >= any_number and not for ==

template <>
QGauss2<1>::QGauss2 () :
		Quadrature<1> (2)
{
				   // points on [-1,1]
  static const double xpts_normal[] = { -sqrt(1./3.), sqrt(1./3.) };
				   // weights on [-1,1]
  static const double wts_normal[]  = { 1., 1. };

				   // points and weights on [0,1]
  static const double xpts[] = { (xpts_normal[0]+1)/2.,
				 (xpts_normal[1]+1)/2.  };
  static const double wts[]  = { wts_normal[0]/2.,
				 wts_normal[1]/2.  };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};


template <>
QGauss3<1>::QGauss3 () :
		Quadrature<1> (3)
{
				   // points on [-1,1]
  static const double xpts_normal[] = { -sqrt(3./5.),
					0.,
					sqrt(3./5.) };
				   // weights on [-1,1]
  static const double wts_normal[]  = { 5./9.,
					8./9.,
					5./9. };

				   // points and weights on [0,1]
  static const double xpts[] = { (xpts_normal[0]+1)/2.,
				 (xpts_normal[1]+1)/2.,
				 (xpts_normal[2]+1)/2. };
  static const double wts[]  = { wts_normal[0]/2.,
				 wts_normal[1]/2.,
				 wts_normal[2]/2. };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};


template <>
QGauss4<1>::QGauss4 () :
		Quadrature<1> (4)
{
				   // points on [-1,1]
  static const double xpts_normal[] = { -sqrt(1./7.*(3-4*sqrt(0.3))),
					-sqrt(1./7.*(3+4*sqrt(0.3))),
					+sqrt(1./7.*(3-4*sqrt(0.3))),
					+sqrt(1./7.*(3+4*sqrt(0.3)))  };
				   // weights on [-1,1]
  static const double wts_normal[]  = { 1./2. + 1./12.*sqrt(10./3.),
					1./2. - 1./12.*sqrt(10./3.),
					1./2. + 1./12.*sqrt(10./3.),
					1./2. - 1./12.*sqrt(10./3.)  };

				   // points and weights on [0,1]
  static const double xpts[] = { (xpts_normal[0]+1)/2.,
				 (xpts_normal[1]+1)/2.,
				 (xpts_normal[2]+1)/2.,
				 (xpts_normal[3]+1)/2. };
  static const double wts[]  = { wts_normal[0]/2.,
				 wts_normal[1]/2.,
				 wts_normal[2]/2.,
				 wts_normal[3]/2. };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};


template <>
QGauss5<1>::QGauss5 () :
		Quadrature<1> (5)
{
				   // points on [-1,1]
  static const double xpts_normal[] = { -sqrt(1./9.*(5.-2*sqrt(10./7.))),
					-sqrt(1./9.*(5.+2*sqrt(10./7.))),
					0,
					+sqrt(1./9.*(5.-2*sqrt(10./7.))),
					+sqrt(1./9.*(5.+2*sqrt(10./7.)))  };
				   // weights on [-1,1]
  static const double wts_normal[]  = { 0.3*(-0.7+5.*sqrt(0.7))/(-2.+5.*sqrt(0.7)),
					0.3*(+0.7+5.*sqrt(0.7))/(+2.+5.*sqrt(0.7)),
					128./225.,
					0.3*(-0.7+5.*sqrt(0.7))/(-2.+5.*sqrt(0.7)),
					0.3*(+0.7+5.*sqrt(0.7))/(+2.+5.*sqrt(0.7)) };

				   // points and weights on [0,1]
  static const double xpts[] = { (xpts_normal[0]+1)/2.,
				 (xpts_normal[1]+1)/2.,
				 (xpts_normal[2]+1)/2.,
				 (xpts_normal[3]+1)/2.,
				 (xpts_normal[4]+1)/2. };
  static const double wts[]  = { wts_normal[0]/2.,
				 wts_normal[1]/2.,
				 wts_normal[2]/2.,
				 wts_normal[3]/2.,
				 wts_normal[4]/2. };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};


template <>
QGauss6<1>::QGauss6 () :
		Quadrature<1> (6)
{
				   // points on [-1,1]
  static const double xpts_normal[] = { -0.932469514203152,
					-0.661209386466265,
					-0.238619186083197,
					+0.238619186083197,
					+0.661209386466265,
					+0.932469514203152  };
				   // weights on [-1,1]
  static const double wts_normal[]  = { 0.171324492379170,
					0.360761573048139,
					0.467913934572691,
					0.467913934572691,
					0.360761573048139,
					0.171324492379170  };

				   // points and weights on [0,1]
  static const double xpts[] = { (xpts_normal[0]+1)/2.,
				 (xpts_normal[1]+1)/2.,
				 (xpts_normal[2]+1)/2.,
				 (xpts_normal[3]+1)/2.,
				 (xpts_normal[4]+1)/2.,
				 (xpts_normal[5]+1)/2. };
  static const double wts[]  = { wts_normal[0]/2.,
				 wts_normal[1]/2.,
				 wts_normal[2]/2.,
				 wts_normal[3]/2.,
				 wts_normal[4]/2.,
				 wts_normal[5]/2. };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};


template <>
QGauss7<1>::QGauss7 () :
		Quadrature<1> (7)
{
				   // points on [-1,1]
  static const double xpts_normal[] = { -0.949107912342759,
					-0.741531185599394,
					-0.405845151377397,
					0,
					+0.405845151377397,
					+0.741531185599394,
					+0.949107912342759 };
				   // weights on [-1,1]
  static const double wts_normal[]  = { 0.129484966168870,
					0.279705391489277,
					0.381830050505119,
					0.417959183673469,
					0.381830050505119,
					0.279705391489277,
					0.129484966168870  };

				   // points and weights on [0,1]
  static const double xpts[] = { (xpts_normal[0]+1)/2.,
				 (xpts_normal[1]+1)/2.,
				 (xpts_normal[2]+1)/2.,
				 (xpts_normal[3]+1)/2.,
				 (xpts_normal[4]+1)/2.,
				 (xpts_normal[5]+1)/2.,
				 (xpts_normal[6]+1)/2. };
  static const double wts[]  = { wts_normal[0]/2.,
				 wts_normal[1]/2.,
				 wts_normal[2]/2.,
				 wts_normal[3]/2.,
				 wts_normal[4]/2.,
				 wts_normal[5]/2.,
				 wts_normal[6]/2. };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};


template <>
QGauss8<1>::QGauss8 () :
		Quadrature<1> (8)
{
				   // points on [-1,1]
  static const double xpts_normal[] = { -0.960289856497536,
					-0.796666477413627,
					-0.525532409916329,
					-0.183434642495650,
					+0.183434642495650,
					+0.525532409916329,
					+0.796666477413627,
					+0.960289856497536 };
				   // weights on [-1,1]
  static const double wts_normal[]  = { 0.101228536200376,
					0.222381034453374,
					0.313706645877887,
					0.362683783378362,
					0.362683783378362,
					0.313706645877887,
					0.222381034453374,
					0.101228536200376  };

				   // points and weights on [0,1]
  static const double xpts[] = { (xpts_normal[0]+1)/2.,
				 (xpts_normal[1]+1)/2.,
				 (xpts_normal[2]+1)/2.,
				 (xpts_normal[3]+1)/2.,
				 (xpts_normal[4]+1)/2.,
				 (xpts_normal[5]+1)/2.,
				 (xpts_normal[6]+1)/2.,
				 (xpts_normal[7]+1)/2. };
  static const double wts[]  = { wts_normal[0]/2.,
				 wts_normal[1]/2.,
				 wts_normal[2]/2.,
				 wts_normal[3]/2.,
				 wts_normal[4]/2.,
				 wts_normal[5]/2.,
				 wts_normal[6]/2.,
				 wts_normal[7]/2. };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};


template <>
QMidpoint<1>::QMidpoint () :
		Quadrature<1>(1)
{
  quadrature_points[0] = Point<1>(0.5);
  weights[0] = 1.0;
};


template <>
QSimpson<1>::QSimpson () :
		Quadrature<1> (3)
{
  static const double xpts[] = { 0.0, 0.5, 1.0 };
  static const double wts[]  = { 1./6., 2./3., 1./6. };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};


template <>
QTrapez<1>::QTrapez () :
		Quadrature<1> (2)
{
  static const double xpts[] = { 0.0, 1.0 };
  static const double wts[]  = { 0.5, 0.5 };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};


// construct the quadrature formulae in higher dimensions by
// tensor product of lower dimensions
template <int dim>
QGauss2<dim>::QGauss2 () :  Quadrature<dim> (QGauss2<dim-1>(), QGauss2<1>())  {};

template <int dim>
QGauss3<dim>::QGauss3 () :
		Quadrature<dim> (QGauss3<dim-1>(), QGauss3<1>()){};

template <int dim>
QGauss4<dim>::QGauss4 () :
		Quadrature<dim> (QGauss4<dim-1>(), QGauss4<1>()){};

template <int dim>
QGauss5<dim>::QGauss5 () :
		Quadrature<dim> (QGauss5<dim-1>(), QGauss5<1>()){};

template <int dim>
QGauss6<dim>::QGauss6 () :
		Quadrature<dim> (QGauss6<dim-1>(), QGauss6<1>()){};

template <int dim>
QGauss7<dim>::QGauss7 () :
		Quadrature<dim> (QGauss7<dim-1>(), QGauss7<1>()){};

template <int dim>
QGauss8<dim>::QGauss8 () :
		Quadrature<dim> (QGauss8<dim-1>(), QGauss8<1>()){};


template <int dim>
QMidpoint<dim>::QMidpoint () :
		Quadrature<dim> (QMidpoint<dim-1>(), QMidpoint<1>()){};

template <int dim>
QSimpson<dim>::QSimpson () :
		Quadrature<dim> (QSimpson<dim-1>(), QSimpson<1>()){};

template <int dim>
QTrapez<dim>::QTrapez () :
		Quadrature<dim> (QTrapez<dim-1>(), QTrapez<1>()){};


// explicite specialization
// note that 1d formulae are specialized by implementation above
template class QGauss2<2>;
template class QGauss3<2>;
template class QGauss4<2>;
template class QGauss5<2>;
template class QGauss6<2>;
template class QGauss7<2>;
template class QGauss8<2>;
template class QMidpoint<2>;
template class QSimpson<2>;
template class QTrapez<2>;

template class QGauss2<3>;
template class QGauss3<3>;
template class QGauss4<3>;
template class QGauss5<3>;
template class QGauss6<3>;
template class QGauss7<3>;
template class QGauss8<3>;
template class QMidpoint<3>;
template class QSimpson<3>;
template class QTrapez<3>;
