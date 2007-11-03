//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/quadrature_lib.h>
#include <cmath>

#ifdef HAVE_STD_NUMERIC_LIMITS
#  include <limits>
#endif

DEAL_II_NAMESPACE_OPEN


// please note: for a given dimension, we need the quadrature formulae
// for all lower dimensions as well. That is why in this file the check
// is for deal_II_dimension >= any_number and not for ==

namespace
{
  template <typename number>
  number abs (const number a)
  {
    return ((a>0) ? a : -a);
  }
}



template <>
QGauss<0>::QGauss (const unsigned int)
                :
                Quadrature<0> (0)
{
                                   // this function has to be provided to
                                   // avoid certain linker failures, but it
                                   // should never be called
  Assert (false, ExcInternalError());
}



template <>
QGauss<1>::QGauss (const unsigned int n)
                :
                Quadrature<1> (n)
{
  if (n == 0)
    return;
  
  const unsigned int m = (n+1)/2;

                                   // tolerance for the Newton
                                   // iteration below. we need to make
                                   // it adaptive since on some
                                   // machines (for example PowerPC)
                                   // long double is the same as
                                   // double -- in that case we can
                                   // only get to a certain multiple
                                   // of the accuracy of double there,
                                   // while on other machines we'd
                                   // like to go further down
				   //
				   // the situation is complicated by
				   // the fact that even if long
				   // double exists and is described
				   // by std::numeric_limits, we may
				   // not actually get the additional
				   // precission. One case where this
				   // happens is on x86, where one can
				   // set hardware flags that disable
				   // long double precision even for
				   // long double variables. these
				   // flags are not usually set, but
				   // for example matlab sets them and
				   // this then breaks deal.II code
				   // that is run as a subroutine to
				   // matlab...
                                   //
                                   // a similar situation exists, btw,
                                   // when running programs under
                                   // valgrind up to and including at
                                   // least version 3.1: valgrind's
                                   // emulator only supports 64 bit
                                   // arithmetic, even for 80 bit long
                                   // doubles.
#ifdef HAVE_STD_NUMERIC_LIMITS
  const long double
    long_double_eps = static_cast<long double>(std::numeric_limits<long double>::epsilon()),
    double_eps      = static_cast<long double>(std::numeric_limits<double>::epsilon());
#else
  const long double
    long_double_eps = 1.09-19L,
    double_eps      = 2.23-16;
#endif

				   // now check whether long double is
				   // more accurate than double, and
				   // set tolerances accordingly
  const long double tolerance
    = (static_cast<long double>(1.0) + long_double_eps != static_cast<long double>(1.0)
       ?
       std::max (double_eps / 100, long_double_eps * 5)
       :
       double_eps * 5
       );

  
  for (unsigned int i=1; i<=m; ++i)
    {
      long double z = std::cos(numbers::PI * (i-.25)/(n+.5));

      long double pp;
      long double p1, p2, p3;

                                       // Newton iteration
      do
	{
					   // compute L_n (z)
	  p1 = 1.;
	  p2 = 0.;
	  for (unsigned int j=0;j<n;++j)
	    {
	      p3 = p2;
	      p2 = p1;
	      p1 = ((2.*j+1.)*z*p2-j*p3)/(j+1);
	    }
	  pp = n*(z*p1-p2)/(z*z-1);
	  z = z-p1/pp;
	}
      while (abs(p1/pp) > tolerance);

      double x = .5*z;
      this->quadrature_points[i-1] = Point<1>(.5-x);
      this->quadrature_points[n-i] = Point<1>(.5+x);
      
      double w = 1./((1.-z*z)*pp*pp);
      this->weights[i-1] = w;
      this->weights[n-i] = w;
    }
}


template <>
QGaussLobatto<1>::QGaussLobatto (unsigned int n)
		:
		Quadrature<1> (n)
{
  Assert (n >= 2, ExcNotImplemented());
  
  std::vector<long double> points  = compute_quadrature_points(n, 1, 1);
  std::vector<long double> w       = compute_quadrature_weights(points, 0, 0);

				   // scale points to the interval
				   // [0.0, 1.0]:
  for (unsigned int i=0; i<points.size(); ++i)
    {
      this->quadrature_points[i] = Point<1>(0.5 + 0.5*static_cast<double>(points[i]));
      this->weights[i]           = 0.5*w[i];
    }
}



template <>
std::vector<long double> QGaussLobatto<1>::
compute_quadrature_points(const unsigned int q,
			  const int alpha,
			  const int beta) const
{
  const unsigned int m = q-2; // no. of inner points
  std::vector<long double> x(m);
  
				   // compute quadrature points with
				   // a Newton algorithm.

				   // Set tolerance. See class QGauss
				   // for detailed explanation.
#ifdef HAVE_STD_NUMERIC_LIMITS
  const long double
    long_double_eps = static_cast<long double>(std::numeric_limits<long double>::epsilon()),
    double_eps      = static_cast<long double>(std::numeric_limits<double>::epsilon());
#else
  const long double
    long_double_eps = 1.09-19L,
    double_eps      = 2.23-16;
#endif

				   // check whether long double is
				   // more accurate than double, and
				   // set tolerances accordingly
  const long double epsilon
    = (static_cast<long double>(1.0) + long_double_eps != static_cast<long double>(1.0)
       ?
       std::max (double_eps / 100, long_double_eps * 5)
       :
       double_eps * 5
       );
  
				   // we take the zeros of the Chebyshev
				   // polynomial (alpha=beta=-0.5) as
				   // initial values:
  for (unsigned int i=0; i<m; ++i)
    x[i] = - std::cos( (long double) (2*i+1)/(2*m) * M_PI );
  
  long double r, s, J_x, f, delta;
  for (unsigned int k=0; k<m; ++k)
    {
      r = x[k];
      if (k>0)
	r = (r + x[k-1])/2;

      do 
	{
	  s = 1.;
	  for (unsigned int i=0; i<k; ++i)
	    s /= r - x[i];

	  J_x   =  0.5*(alpha + beta + m + 1)*JacobiP(r, alpha+1, beta+1, m-1);
	  f     = JacobiP(r, alpha, beta, m);
	  delta = f/(f*s- J_x);
	  r += delta;
	}
      while (std::fabs(delta) >= epsilon);
      
      x[k] = r;
    } // for

				   // add boundary points:
  x.insert(x.begin(), -1.L);
  x.push_back(+1.L);

  return x;
}



template <>
std::vector<long double> QGaussLobatto<1>::
compute_quadrature_weights(const std::vector<long double> &x,
			   const int alpha,
			   const int beta) const
{
  const unsigned int q = x.size();
  std::vector<long double> w(q);
  long double s = 0.L;
  
  const long double factor = std::pow(2., alpha+beta+1) *
			     gamma(alpha+q) *
			     gamma(beta+q) /
			     ((q-1)*gamma(q)*gamma(alpha+beta+q+1));
  for (unsigned int i=0; i<q; ++i)
    {
      s = JacobiP(x[i], alpha, beta, q-1);
      w[i] = factor/(s*s);
    }
  w[0]   *= (beta + 1);
  w[q-1] *= (alpha + 1);

  return w;
}



template <>
long double QGaussLobatto<1>::JacobiP(const long double x,
				      const int alpha,
				      const int beta,
				      const unsigned int n) const
{
				   // the Jacobi polynomial is evaluated
				   // using a recursion formula.
  std::vector<long double> p(n+1);
  int v, a1, a2, a3, a4;
  
				   // initial values P_0(x), P_1(x):
  p[0] = 1.0L;
  if (n==0) return p[0];
  p[1] = ((alpha+beta+2)*x + (alpha-beta))/2;
  if (n==1) return p[1];
  
  for (unsigned int i=1; i<=(n-1); ++i)
    {
      v  = 2*i + alpha + beta;
      a1 = 2*(i+1)*(i + alpha + beta + 1)*v;
      a2 = (v + 1)*(alpha*alpha - beta*beta);
      a3 = v*(v + 1)*(v + 2);
      a4 = 2*(i+alpha)*(i+beta)*(v + 2);

      p[i+1] = static_cast<long double>( (a2 + a3*x)*p[i] - a4*p[i-1])/a1;
    } // for
  return p[n];
}



template <>
unsigned int QGaussLobatto<1>::gamma(const unsigned int n) const
{
  unsigned int result = n - 1;
  for (int i=n-2; i>1; --i)
    result *= i;
  return result;
}



template <>
QGauss2<1>::QGauss2 ()
                :
		Quadrature<1> (2)
{
				   // points on [-1,1]
  static const double xpts_normal[] = { -std::sqrt(1./3.), std::sqrt(1./3.) };
				   // weights on [-1,1]
  static const double wts_normal[]  = { 1., 1. };

				   // points and weights on [0,1]
  static const double xpts[] = { (xpts_normal[0]+1)/2.,
				 (xpts_normal[1]+1)/2.  };
  static const double wts[]  = { wts_normal[0]/2.,
				 wts_normal[1]/2.  };

  for (unsigned int i=0; i<this->n_quadrature_points; ++i) 
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i] = wts[i];
    };
}



template <>
QGauss3<1>::QGauss3 ()
                :
		Quadrature<1> (3)
{
				   // points on [-1,1]
  static const double xpts_normal[] = { -std::sqrt(3./5.),
					0.,
					std::sqrt(3./5.) };
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

  for (unsigned int i=0; i<this->n_quadrature_points; ++i) 
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i] = wts[i];
    };
}



template <>
QGauss4<1>::QGauss4 () 
                :
		Quadrature<1> (4)
{
				   // points on [-1,1]
  static const double xpts_normal[] = { -std::sqrt(1./7.*(3+4*std::sqrt(0.3))),
					-std::sqrt(1./7.*(3-4*std::sqrt(0.3))),
					+std::sqrt(1./7.*(3-4*std::sqrt(0.3))),
					+std::sqrt(1./7.*(3+4*std::sqrt(0.3)))  };
				   // weights on [-1,1]
  static const double wts_normal[]  = { 1./2. - 1./12.*std::sqrt(10./3.),
					1./2. + 1./12.*std::sqrt(10./3.),
					1./2. + 1./12.*std::sqrt(10./3.),
					1./2. - 1./12.*std::sqrt(10./3.)  };

				   // points and weights on [0,1]
  static const double xpts[] = { (xpts_normal[0]+1)/2.,
				 (xpts_normal[1]+1)/2.,
				 (xpts_normal[2]+1)/2.,
				 (xpts_normal[3]+1)/2. };
  static const double wts[]  = { wts_normal[0]/2.,
				 wts_normal[1]/2.,
				 wts_normal[2]/2.,
				 wts_normal[3]/2. };

  for (unsigned int i=0; i<this->n_quadrature_points; ++i) 
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i] = wts[i];
    };
}


template <>
QGauss5<1>::QGauss5 ()
                :
		Quadrature<1> (5)
{
				   // points on [-1,1]
  static const double xpts_normal[] = { -std::sqrt(1./9.*(5.+2*std::sqrt(10./7.))),
					-std::sqrt(1./9.*(5.-2*std::sqrt(10./7.))),
					0,
					+std::sqrt(1./9.*(5.-2*std::sqrt(10./7.))),
					+std::sqrt(1./9.*(5.+2*std::sqrt(10./7.)))  };
				   // weights on [-1,1]
  static const double wts_normal[]  = { 0.3*(+0.7+5.*std::sqrt(0.7))/(+2.+5.*std::sqrt(0.7)),
					0.3*(-0.7+5.*std::sqrt(0.7))/(-2.+5.*std::sqrt(0.7)),
					128./225.,
					0.3*(-0.7+5.*std::sqrt(0.7))/(-2.+5.*std::sqrt(0.7)),
					0.3*(+0.7+5.*std::sqrt(0.7))/(+2.+5.*std::sqrt(0.7)) };

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

  for (unsigned int i=0; i<this->n_quadrature_points; ++i) 
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i] = wts[i];
    };
}



template <>
QGauss6<1>::QGauss6 ()
                :
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

  for (unsigned int i=0; i<this->n_quadrature_points; ++i) 
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i] = wts[i];
    };
}



template <>
QGauss7<1>::QGauss7 ()
                :
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

  for (unsigned int i=0; i<this->n_quadrature_points; ++i) 
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i] = wts[i];
    };
}



template <>
QMidpoint<1>::QMidpoint ()
                :
		Quadrature<1>(1)
{
  this->quadrature_points[0] = Point<1>(0.5);
  this->weights[0] = 1.0;
}



template <>
QTrapez<1>::QTrapez ()
                :
		Quadrature<1> (2)
{
  static const double xpts[] = { 0.0, 1.0 };
  static const double wts[]  = { 0.5, 0.5 };

  for (unsigned int i=0; i<this->n_quadrature_points; ++i) 
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i] = wts[i];
    };
}



template <>
QSimpson<1>::QSimpson ()
                :
		Quadrature<1> (3)
{
  static const double xpts[] = { 0.0, 0.5, 1.0 };
  static const double wts[]  = { 1./6., 2./3., 1./6. };

  for (unsigned int i=0; i<this->n_quadrature_points; ++i) 
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i] = wts[i];
    };
}



template <>
QMilne<1>::QMilne ()
                :
		Quadrature<1> (5)
{
  static const double xpts[] = { 0.0, .25, .5, .75, 1.0 };
  static const double wts[]  = { 7./90., 32./90., 12./90., 32./90., 7./90. };

  for (unsigned int i=0; i<this->n_quadrature_points; ++i) 
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i] = wts[i];
    };
}



template <>
QWeddle<1>::QWeddle ()
                :
		Quadrature<1> (7)
{
  static const double xpts[] = { 0.0, 1./6., 1./3., .5, 2./3., 5./6., 1.0 };
  static const double wts[]  = { 41./840., 216./840., 27./840., 272./840.,
				 27./840., 216./840., 41./840.
  };

  for (unsigned int i=0; i<this->n_quadrature_points; ++i) 
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i] = wts[i];
    };
}



// construct the quadrature formulae in higher dimensions by
// tensor product of lower dimensions

template <int dim>
QGauss<dim>::QGauss (const unsigned int n)
  :  Quadrature<dim> (QGauss<dim-1>(n), QGauss<1>(n))
{}



template <int dim>
QGaussLobatto<dim>::QGaussLobatto (const unsigned int n)
  :  Quadrature<dim> (QGaussLobatto<dim-1>(n), QGaussLobatto<1>(n))
{}



template <int dim>
QGauss2<dim>::QGauss2 ()
                :
                Quadrature<dim> (QGauss2<dim-1>(), QGauss2<1>())
{}



template <int dim>
QGauss3<dim>::QGauss3 ()
                :
		Quadrature<dim> (QGauss3<dim-1>(), QGauss3<1>())
{}



template <int dim>
QGauss4<dim>::QGauss4 ()
                :
		Quadrature<dim> (QGauss4<dim-1>(), QGauss4<1>())
{}



template <int dim>
QGauss5<dim>::QGauss5 ()
                :
		Quadrature<dim> (QGauss5<dim-1>(), QGauss5<1>())
{}



template <int dim>
QGauss6<dim>::QGauss6 ()
                :
		Quadrature<dim> (QGauss6<dim-1>(), QGauss6<1>())
{}



template <int dim>
QGauss7<dim>::QGauss7 ()
                :
		Quadrature<dim> (QGauss7<dim-1>(), QGauss7<1>())
{}



template <int dim>
QMidpoint<dim>::QMidpoint ()
                :
		Quadrature<dim> (QMidpoint<dim-1>(), QMidpoint<1>())
{}



template <int dim>
QTrapez<dim>::QTrapez ()
                :
		Quadrature<dim> (QTrapez<dim-1>(), QTrapez<1>())
{}



template <int dim>
QSimpson<dim>::QSimpson ()
                :
		Quadrature<dim> (QSimpson<dim-1>(), QSimpson<1>())
{}



template <int dim>
QMilne<dim>::QMilne ()
                :
		Quadrature<dim> (QMilne<dim-1>(), QMilne<1>())
{}


template <int dim>
QWeddle<dim>::QWeddle ()
                :
		Quadrature<dim> (QWeddle<dim-1>(), QWeddle<1>())
{}


// explicit specialization
// note that 1d formulae are specialized by implementation above
template class QGauss<2>;
template class QGaussLobatto<2>;
template class QGauss2<2>;
template class QGauss3<2>;
template class QGauss4<2>;
template class QGauss5<2>;
template class QGauss6<2>;
template class QGauss7<2>;
template class QMidpoint<2>;
template class QTrapez<2>;
template class QSimpson<2>;
template class QMilne<2>;
template class QWeddle<2>;

template class QGauss<3>;
template class QGaussLobatto<3>;
template class QGauss2<3>;
template class QGauss3<3>;
template class QGauss4<3>;
template class QGauss5<3>;
template class QGauss6<3>;
template class QGauss7<3>;
template class QMidpoint<3>;
template class QTrapez<3>;
template class QSimpson<3>;
template class QMilne<3>;
template class QWeddle<3>;

DEAL_II_NAMESPACE_CLOSE
