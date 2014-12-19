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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/geometry_info.h>

#include <cmath>
#include <limits>
#include <algorithm>


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
  // there are n_q^dim == 1
  // points
  Quadrature<0> (1)
{
  // the single quadrature point gets unit
  // weight
  this->weights[0] = 1;
}



template <>
QGaussLobatto<0>::QGaussLobatto (const unsigned int)
  :
  // there are n_q^dim == 1
  // points
  Quadrature<0> (1)
{
  // the single quadrature point gets unit
  // weight
  this->weights[0] = 1;
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
  // precision. One case where this
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
  // least version 3.3: valgrind's
  // emulator only supports 64 bit
  // arithmetic, even for 80 bit long
  // doubles.
  const long double
  long_double_eps = static_cast<long double>(std::numeric_limits<long double>::epsilon()),
  double_eps      = static_cast<long double>(std::numeric_limits<double>::epsilon());

  // now check whether long double is more
  // accurate than double, and set
  // tolerances accordingly. generate a one
  // that really is generated at run-time
  // and is not optimized away by the
  // compiler. that makes sure that the
  // tolerance is set at run-time with the
  // current behavior, not at compile-time
  // (not doing so leads to trouble with
  // valgrind for example).
  volatile long double runtime_one = 1.0;
  const long double tolerance
    = (runtime_one + long_double_eps != runtime_one
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
          for (unsigned int j=0; j<n; ++j)
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
QGaussLobatto<1>::QGaussLobatto (const unsigned int n)
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
  const long double
  long_double_eps = static_cast<long double>(std::numeric_limits<long double>::epsilon()),
  double_eps      = static_cast<long double>(std::numeric_limits<double>::epsilon());

  // check whether long double is
  // more accurate than double, and
  // set tolerances accordingly
  volatile long double runtime_one = 1.0;
  const long double tolerance
    = (runtime_one + long_double_eps != runtime_one
       ?
       std::max (double_eps / 100, long_double_eps * 5)
       :
       double_eps * 5
      );

  // The following implementation
  // follows closely the one given in
  // the appendix of the book by
  // Karniadakis and Sherwin:
  // Spectral/hp element methods for
  // computational fluid dynamics
  // (Oxford University Press, 2005)

  // we take the zeros of the Chebyshev
  // polynomial (alpha=beta=-0.5) as
  // initial values:
  for (unsigned int i=0; i<m; ++i)
    x[i] = - std::cos( (long double) (2*i+1)/(2*m) * numbers::PI );

  long double r, s, J_x, f, delta;

  for (unsigned int k=0; k<m; ++k)
    {
      r = x[k];
      if (k>0)
        r = (r + x[k-1])/2;

      do
        {
          s = 0.;
          for (unsigned int i=0; i<k; ++i)
            s += 1./(r - x[i]);

          J_x   =  0.5*(alpha + beta + m + 1)*JacobiP(r, alpha+1, beta+1, m-1);
          f     = JacobiP(r, alpha, beta, m);
          delta = f/(f*s- J_x);
          r += delta;
        }
      while (std::fabs(delta) >= tolerance);

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
long double QGaussLobatto<1>::gamma(const unsigned int n) const
{
  long double result = n - 1;
  for (int i=n-2; i>1; --i)
    result *= i;
  return result;
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

  for (unsigned int i=0; i<this->size(); ++i)
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

  for (unsigned int i=0; i<this->size(); ++i)
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

  for (unsigned int i=0; i<this->size(); ++i)
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

  for (unsigned int i=0; i<this->size(); ++i)
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i] = wts[i];
    };
}


template <>
QGaussLog<1>::QGaussLog(const unsigned int n,
                        const bool revert)
  :
  Quadrature<1> (n)
{

  std::vector<double> p=set_quadrature_points(n);
  std::vector<double> w=set_quadrature_weights(n);

  for (unsigned int i=0; i<this->size(); ++i)
    {
      // Using the change of variables x=1-t, it's possible to show
      // that int f(x)ln|1-x| = int f(1-t) ln|t|, which implies that
      // we can use this quadrature formula also with weight ln|1-x|.
      this->quadrature_points[i] = revert ? Point<1>(1-p[n-1-i]) : Point<1>(p[i]);
      this->weights[i]           = revert ? w[n-1-i] : w[i];
    }

}

template <>
std::vector<double>
QGaussLog<1>::set_quadrature_points(const unsigned int n) const
{

  std::vector<double> points(n);

  switch (n)
    {
    case 1:
      points[0] = 0.3333333333333333;
      break;

    case 2:
      points[0] = 0.1120088061669761;
      points[1] = 0.6022769081187381;
      break;

    case 3:
      points[0] = 0.06389079308732544;
      points[1] = 0.3689970637156184;
      points[2] = 0.766880303938942;
      break;

    case 4:
      points[0] = 0.04144848019938324;
      points[1] = 0.2452749143206022;
      points[2] = 0.5561654535602751;
      points[3] = 0.848982394532986;
      break;

    case 5:
      points[0] = 0.02913447215197205;
      points[1] = 0.1739772133208974;
      points[2] =  0.4117025202849029;
      points[3] = 0.6773141745828183;
      points[4] = 0.89477136103101;
      break;

    case 6:
      points[0] = 0.02163400584411693;
      points[1] = 0.1295833911549506;
      points[2] = 0.3140204499147661;
      points[3] = 0.5386572173517997;
      points[4] = 0.7569153373774084;
      points[5] = 0.922668851372116;
      break;


    case 7:
      points[0] = 0.0167193554082585;
      points[1] = 0.100185677915675;
      points[2] = 0.2462942462079286;
      points[3] = 0.4334634932570557;
      points[4] = 0.6323509880476823;
      points[5] = 0.81111862674023;
      points[6] = 0.940848166743287;
      break;

    case 8:
      points[0] = 0.01332024416089244;
      points[1] = 0.07975042901389491;
      points[2] = 0.1978710293261864;
      points[3] =   0.354153994351925;
      points[4] =   0.5294585752348643;
      points[5] = 0.7018145299391673;
      points[6] = 0.849379320441094;
      points[7] = 0.953326450056343;
      break;

    case 9:
      points[0] = 0.01086933608417545;
      points[1] = 0.06498366633800794;
      points[2] = 0.1622293980238825;
      points[3] = 0.2937499039716641;
      points[4] = 0.4466318819056009;
      points[5] = 0.6054816627755208;
      points[6] = 0.7541101371585467;
      points[7] = 0.877265828834263;
      points[8] = 0.96225055941096;
      break;

    case 10:
      points[0] = 0.00904263096219963;
      points[1] = 0.05397126622250072;
      points[2] =  0.1353118246392511;
      points[3] = 0.2470524162871565;
      points[4] = 0.3802125396092744;
      points[5] = 0.5237923179723384;
      points[6] = 0.6657752055148032;
      points[7] = 0.7941904160147613;
      points[8] = 0.898161091216429;
      points[9] = 0.9688479887196;
      break;


    case 11:
      points[0] = 0.007643941174637681;
      points[1] = 0.04554182825657903;
      points[2] = 0.1145222974551244;
      points[3] = 0.2103785812270227;
      points[4] = 0.3266955532217897;
      points[5] = 0.4554532469286375;
      points[6] = 0.5876483563573721;
      points[7] = 0.7139638500230458;
      points[8] = 0.825453217777127;
      points[9] = 0.914193921640008;
      points[10] = 0.973860256264123;
      break;

    case 12:
      points[0] = 0.006548722279080035;
      points[1] = 0.03894680956045022;
      points[2] = 0.0981502631060046;
      points[3] = 0.1811385815906331;
      points[4] = 0.2832200676673157;
      points[5] = 0.398434435164983;
      points[6] = 0.5199526267791299;
      points[7] = 0.6405109167754819;
      points[8] = 0.7528650118926111;
      points[9] = 0.850240024421055;
      points[10] = 0.926749682988251;
      points[11] = 0.977756129778486;
      break;

    default:
      Assert(false, ExcNotImplemented());
      break;
    }

  return points;
}


template <>
std::vector<double>
QGaussLog<1>::set_quadrature_weights(const unsigned int n) const
{

  std::vector<double> weights(n);

  switch (n)
    {
    case 1:
      weights[0] = -1.0;
      break;
    case 2:
      weights[0] = -0.7185393190303845;
      weights[1] = -0.2814606809696154;
      break;

    case 3:
      weights[0] = -0.5134045522323634;
      weights[1] = -0.3919800412014877;
      weights[2] = -0.0946154065661483;
      break;

    case 4:
      weights[0] =-0.3834640681451353;
      weights[1] =-0.3868753177747627;
      weights[2] =-0.1904351269501432;
      weights[3] =-0.03922548712995894;
      break;

    case 5:
      weights[0] =-0.2978934717828955;
      weights[1] =-0.3497762265132236;
      weights[2] =-0.234488290044052;
      weights[3] =-0.0989304595166356;
      weights[4] =-0.01891155214319462;
      break;

    case 6:
      weights[0] = -0.2387636625785478;
      weights[1] = -0.3082865732739458;
      weights[2] = -0.2453174265632108;
      weights[3] = -0.1420087565664786;
      weights[4] = -0.05545462232488041;
      weights[5] = -0.01016895869293513;
      break;


    case 7:
      weights[0] = -0.1961693894252476;
      weights[1] = -0.2703026442472726;
      weights[2] = -0.239681873007687;
      weights[3] = -0.1657757748104267;
      weights[4] = -0.0889432271377365;
      weights[5] = -0.03319430435645653;
      weights[6] = -0.005932787015162054;
      break;

    case 8:
      weights[0] = -0.164416604728002;
      weights[1] = -0.2375256100233057;
      weights[2] = -0.2268419844319134;
      weights[3] = -0.1757540790060772;
      weights[4] = -0.1129240302467932;
      weights[5] = -0.05787221071771947;
      weights[6] = -0.02097907374214317;
      weights[7] = -0.003686407104036044;
      break;

    case 9:
      weights[0] = -0.1400684387481339;
      weights[1] = -0.2097722052010308;
      weights[2] = -0.211427149896601;
      weights[3] = -0.1771562339380667;
      weights[4] = -0.1277992280331758;
      weights[5] = -0.07847890261203835;
      weights[6] = -0.0390225049841783;
      weights[7] = -0.01386729555074604;
      weights[8] = -0.002408041036090773;
      break;

    case 10:
      weights[0] = -0.12095513195457;
      weights[1] = -0.1863635425640733;
      weights[2] = -0.1956608732777627;
      weights[3] = -0.1735771421828997;
      weights[4] = -0.135695672995467;
      weights[5] = -0.0936467585378491;
      weights[6] = -0.05578772735275126;
      weights[7] = -0.02715981089692378;
      weights[8] = -0.00951518260454442;
      weights[9] = -0.001638157633217673;
      break;


    case 11:
      weights[0] = -0.1056522560990997;
      weights[1] = -0.1665716806006314;
      weights[2] = -0.1805632182877528;
      weights[3] = -0.1672787367737502;
      weights[4] = -0.1386970574017174;
      weights[5] = -0.1038334333650771;
      weights[6] = -0.06953669788988512;
      weights[7] = -0.04054160079499477;
      weights[8] = -0.01943540249522013;
      weights[9] = -0.006737429326043388;
      weights[10] = -0.001152486965101561;
      break;

    case 12:
      weights[0] = -0.09319269144393;
      weights[1] = -0.1497518275763289;
      weights[2] = -0.166557454364573;
      weights[3] = -0.1596335594369941;
      weights[4] = -0.1384248318647479;
      weights[5] = -0.1100165706360573;
      weights[6] = -0.07996182177673273;
      weights[7] = -0.0524069547809709;
      weights[8] = -0.03007108900074863;
      weights[9] = -0.01424924540252916;
      weights[10] = -0.004899924710875609;
      weights[11] = -0.000834029009809656;
      break;

    default:
      Assert(false, ExcNotImplemented());
      break;
    }

  return weights;

}


template<>
QGaussLogR<1>::QGaussLogR(const unsigned int n,
                          const Point<1> origin,
                          const double alpha,
                          const bool factor_out_singularity) :
  Quadrature<1>( ( (origin[0] == 0) || (origin[0] == 1) ) ?
                 (alpha == 1 ? n : 2*n ) : 4*n ),
  fraction( ( (origin[0] == 0) || (origin[0] == 1.) ) ? 1. : origin[0] )
{
  // The three quadrature formulas that make this one up. There are
  // at most two when the origin is one of the extremes, and there is
  // only one if the origin is one of the extremes and alpha is
  // equal to one.
  //
  // If alpha is different from one, then we need a correction which
  // is performed with a standard Gauss quadrature rule on each
  // segment. This is not needed in the standard case where alpha is
  // equal to one and the origin is on one of the extremes. We
  // integrate with weight ln|(x-o)/alpha|. In the easy cases, we
  // only need n quadrature points. In the most difficult one, we
  // need 2*n points for the first segment, and 2*n points for the
  // second segment.
  QGaussLog<1> quad1(n, origin[0] != 0);
  QGaussLog<1> quad2(n);
  QGauss<1> quad(n);

  // Check that the origin is inside 0,1
  Assert( (fraction >= 0) && (fraction <= 1),
          ExcMessage("Origin is outside [0,1]."));

  // Non singular offset. This is the start of non singular quad
  // points.
  unsigned int ns_offset = (fraction == 1) ? n : 2*n;

  for (unsigned int i=0, j=ns_offset; i<n; ++i, ++j)
    {
      // The first i quadrature points are the same as quad1, and
      // are by default singular.
      this->quadrature_points[i] = quad1.point(i)*fraction;
      this->weights[i] = quad1.weight(i)*fraction;

      // We need to scale with -log|fraction*alpha|
      if ( (alpha != 1) || (fraction != 1) )
        {
          this->quadrature_points[j] = quad.point(i)*fraction;
          this->weights[j] = -std::log(alpha/fraction)*quad.weight(i)*fraction;
        }
      // In case we need the second quadrature as well, do it now.
      if (fraction != 1)
        {
          this->quadrature_points[i+n] = quad2.point(i)*(1-fraction)+Point<1>(fraction);
          this->weights[i+n] = quad2.weight(i)*(1-fraction);

          // We need to scale with -log|fraction*alpha|
          this->quadrature_points[j+n] = quad.point(i)*(1-fraction)+Point<1>(fraction);
          this->weights[j+n] = -std::log(alpha/(1-fraction))*quad.weight(i)*(1-fraction);
        }
    }
  if (factor_out_singularity == true)
    for (unsigned int i=0; i<size(); ++i)
      {
        Assert( this->quadrature_points[i] != origin,
                ExcMessage("The singularity cannot be on a Gauss point of the same order!") );
        this->weights[i] /= std::log(std::abs( (this->quadrature_points[i]-origin)[0] )/alpha );
      }
}


template<>
unsigned int QGaussOneOverR<2>::quad_size(const Point<2> singularity,
                                          const unsigned int n)
{
  double eps=1e-8;
  bool on_edge=false;
  bool on_vertex=false;
  for (unsigned int i=0; i<2; ++i)
    if ( ( std::abs(singularity[i]  ) < eps ) ||
         ( std::abs(singularity[i]-1) < eps ) )
      on_edge = true;
  if (on_edge && (std::abs( (singularity-Point<2>(.5, .5)).square()-.5)
                  < eps) )
    on_vertex = true;
  if (on_vertex) return (2*n*n);
  if (on_edge) return (4*n*n);
  return (8*n*n);
}

template<>
QGaussOneOverR<2>::QGaussOneOverR(const unsigned int n,
                                  const Point<2> singularity,
                                  const bool factor_out_singularity) :
  Quadrature<2>(quad_size(singularity, n))
{
  // We treat all the cases in the
  // same way. Split the element in 4
  // pieces, measure the area, if
  // it's relevant, add the
  // quadrature connected to that
  // singularity.
  std::vector<QGaussOneOverR<2> > quads;
  std::vector<Point<2> > origins;
  // Id of the corner with a
  // singularity
  quads.push_back(QGaussOneOverR(n, 3, factor_out_singularity));
  quads.push_back(QGaussOneOverR(n, 2, factor_out_singularity));
  quads.push_back(QGaussOneOverR(n, 1, factor_out_singularity));
  quads.push_back(QGaussOneOverR(n, 0, factor_out_singularity));

  origins.push_back(Point<2>(0.,0.));
  origins.push_back(Point<2>(singularity[0],0.));
  origins.push_back(Point<2>(0.,singularity[1]));
  origins.push_back(singularity);

  // Lexycographical ordering.

  double eps = 1e-8;
  unsigned int q_id = 0; // Current quad point index.
  double area = 0;
  Point<2> dist;

  for (unsigned int box=0; box<4; ++box)
    {
      dist = (singularity-GeometryInfo<2>::unit_cell_vertex(box));
      dist = Point<2>(std::abs(dist[0]), std::abs(dist[1]));
      area = dist[0]*dist[1];
      if (area > eps)
        for (unsigned int q=0; q<quads[box].size(); ++q, ++q_id)
          {
            const Point<2> &qp = quads[box].point(q);
            this->quadrature_points[q_id] =
              origins[box]+
              Point<2>(dist[0]*qp[0], dist[1]*qp[1]);
            this->weights[q_id] = quads[box].weight(q)*area;
          }
    }
}


template<>
QGaussOneOverR<2>::QGaussOneOverR(const unsigned int n,
                                  const unsigned int vertex_index,
                                  const bool factor_out_singularity) :
  Quadrature<2>(2*n *n)
{
  // This version of the constructor
  // works only for the 4
  // vertices. If you need a more
  // general one, you should use the
  // one with the Point<2> in the
  // constructor.
  Assert(vertex_index <4, ExcIndexRange(vertex_index, 0, 4));

  // Start with the gauss quadrature formula on the (u,v) reference
  // element.
  QGauss<2> gauss(n);

  Assert(gauss.size() == n*n, ExcInternalError());

  // For the moment we only implemented this for the vertices of a
  // quadrilateral. We are planning to do this also for the support
  // points of arbitrary FE_Q elements, to allow the use of this
  // class in boundary element programs with higher order mappings.
  Assert(vertex_index < 4, ExcIndexRange(vertex_index, 0, 4));

  // We create only the first one. All other pieces are rotation of
  // this one.
  // In this case the transformation is
  //
  // (x,y) = (u, u tan(pi/4 v))
  //
  // with Jacobian
  //
  // J = pi/4 R / cos(pi/4 v)
  //
  // And we get rid of R to take into account the singularity,
  // unless specified differently in the constructor.
  std::vector<Point<2> >      &ps = this->quadrature_points;
  std::vector<double>         &ws = this->weights;
  double pi4 = numbers::PI/4;

  for (unsigned int q=0; q<gauss.size(); ++q)
    {
      const Point<2> &gp = gauss.point(q);
      ps[q][0] = gp[0];
      ps[q][1] = gp[0] * std::tan(pi4*gp[1]);
      ws[q]    = gauss.weight(q)*pi4/std::cos(pi4 *gp[1]);
      if (factor_out_singularity)
        ws[q] *= (ps[q]-GeometryInfo<2>::unit_cell_vertex(0)).norm();
      // The other half of the quadrilateral is symmetric with
      // respect to xy plane.
      ws[gauss.size()+q]    = ws[q];
      ps[gauss.size()+q][0] = ps[q][1];
      ps[gauss.size()+q][1] = ps[q][0];
    }

  // Now we distribute these vertices in the correct manner
  double theta = 0;
  switch (vertex_index)
    {
    case 0:
      theta = 0;
      break;
    case 1:
      //
      theta = numbers::PI/2;
      break;
    case 2:
      theta = -numbers::PI/2;
      break;
    case 3:
      theta = numbers::PI;
      break;
    }

  double R00 =  std::cos(theta), R01 = -std::sin(theta);
  double R10 =  std::sin(theta), R11 =  std::cos(theta);

  if (vertex_index != 0)
    for (unsigned int q=0; q<size(); ++q)
      {
        double x = ps[q][0]-.5,  y = ps[q][1]-.5;

        ps[q][0] = R00*x + R01*y + .5;
        ps[q][1] = R10*x + R11*y + .5;
      }
}


template <int dim>
QSorted<dim>::QSorted(Quadrature<dim> quad) :
  Quadrature<dim>(quad.size())
{
  std::vector< std::pair<double, Point<dim> > > wp;
  for (unsigned int i=0; i<quad.size(); ++i)
    wp.push_back(std::pair<double, Point<dim> >(quad.weight(i),
                                                quad.point(i)));
  sort(wp.begin(), wp.end(), *this);
  for (unsigned int i=0; i<quad.size(); ++i)
    {
      this->weights[i] = wp[i].first;
      this->quadrature_points[i] = wp[i].second;
    }
}


template <int dim>
bool QSorted<dim>::operator()(const std::pair<double, Point<dim> > &a,
                              const std::pair<double, Point<dim> > &b)
{
  return (a.first < b.first);
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
template class QMidpoint<2>;
template class QTrapez<2>;
template class QSimpson<2>;
template class QMilne<2>;
template class QWeddle<2>;

template class QGauss<3>;
template class QGaussLobatto<3>;
template class QMidpoint<3>;
template class QTrapez<3>;
template class QSimpson<3>;
template class QMilne<3>;
template class QWeddle<3>;

template class QSorted<1>;
template class QSorted<2>;
template class QSorted<3>;

DEAL_II_NAMESPACE_CLOSE
