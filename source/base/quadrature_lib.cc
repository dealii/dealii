// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>


DEAL_II_NAMESPACE_OPEN


// please note: for a given dimension, we need the quadrature formulae
// for all lower dimensions as well. That is why in this file the check
// is for deal_II_dimension >= any_number and not for ==



template <>
QGauss<0>::QGauss(const unsigned int)
  : // there are n_q^dim == 1
    // points
  Quadrature<0>(1)
{
  // the single quadrature point gets unit
  // weight
  this->weights[0] = 1;
}



template <>
QGaussLobatto<0>::QGaussLobatto(const unsigned int)
  : // there are n_q^dim == 1
    // points
  Quadrature<0>(1)
{
  // the single quadrature point gets unit
  // weight
  this->weights[0] = 1;
}



template <>
QGauss<1>::QGauss(const unsigned int n)
  : Quadrature<1>(n)
{
  if (n == 0)
    return;

  std::vector<long double> points =
    Polynomials::jacobi_polynomial_roots<long double>(n, 0, 0);

  for (unsigned int i = 0; i < (points.size() + 1) / 2; ++i)
    {
      this->quadrature_points[i][0]         = points[i];
      this->quadrature_points[n - i - 1][0] = 1. - points[i];

      // derivative of Jacobi polynomial
      const long double pp =
        0.5 * (n + 1) *
        Polynomials::jacobi_polynomial_value(n - 1, 1, 1, points[i]);
      const long double x      = -1. + 2. * points[i];
      const double      w      = 1. / ((1. - x * x) * pp * pp);
      this->weights[i]         = w;
      this->weights[n - i - 1] = w;
    }
}

namespace internal
{
  namespace QGaussLobatto
  {
    /**
     * Evaluate the Gamma function $ \Gamma(n) = (n-1)! $.
     * @param n  point of evaluation (integer).
     */
    long double
    gamma(const unsigned int n)
    {
      long double result = n - 1;
      for (int i = n - 2; i > 1; --i)
        result *= i;
      return result;
    }



    /**
     * Compute Legendre-Gauss-Lobatto quadrature weights. The quadrature points
     * and weights are related to Jacobi polynomial specified by @p alpha, @p
     * beta. @p x denotes the quadrature points.
     *
     * @return Vector containing weights.
     */
    std::vector<long double>
    compute_quadrature_weights(const std::vector<long double> &x,
                               const int                       alpha,
                               const int                       beta)
    {
      const unsigned int       q = x.size();
      std::vector<long double> w(q);

      const long double factor =
        std::pow(2., alpha + beta + 1) * gamma(alpha + q) * gamma(beta + q) /
        ((q - 1) * gamma(q) * gamma(alpha + beta + q + 1));
      for (unsigned int i = 0; i < q; ++i)
        {
          const long double s =
            Polynomials::jacobi_polynomial_value(q - 1, alpha, beta, x[i]);
          w[i] = factor / (s * s);
        }
      w[0] *= (beta + 1);
      w[q - 1] *= (alpha + 1);

      return w;
    }
  } // namespace QGaussLobatto
} // namespace internal


#ifndef DOXYGEN
template <>
QGaussLobatto<1>::QGaussLobatto(const unsigned int n)
  : Quadrature<1>(n)
{
  Assert(n >= 2, ExcNotImplemented());

  std::vector<long double> points =
    Polynomials::jacobi_polynomial_roots<long double>(n - 2, 1, 1);
  points.insert(points.begin(), 0);
  points.push_back(1.);
  std::vector<long double> w =
    internal::QGaussLobatto::compute_quadrature_weights(points, 0, 0);

  // scale weights to the interval [0.0, 1.0]:
  for (unsigned int i = 0; i < points.size(); ++i)
    {
      this->quadrature_points[i][0] = points[i];
      this->weights[i]              = 0.5 * w[i];
    }
}
#endif


template <>
QMidpoint<1>::QMidpoint()
  : Quadrature<1>(1)
{
  this->quadrature_points[0] = Point<1>(0.5);
  this->weights[0]           = 1.0;
}



template <>
QTrapezoid<1>::QTrapezoid()
  : Quadrature<1>(2)
{
  static const double xpts[] = {0.0, 1.0};
  static const double wts[]  = {0.5, 0.5};

  for (unsigned int i = 0; i < this->size(); ++i)
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i]           = wts[i];
    }
}



template <>
QSimpson<1>::QSimpson()
  : Quadrature<1>(3)
{
  static const double xpts[] = {0.0, 0.5, 1.0};
  static const double wts[]  = {1. / 6., 2. / 3., 1. / 6.};

  for (unsigned int i = 0; i < this->size(); ++i)
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i]           = wts[i];
    }
}



template <>
QMilne<1>::QMilne()
  : Quadrature<1>(5)
{
  static const double xpts[] = {0.0, .25, .5, .75, 1.0};
  static const double wts[]  = {
    7. / 90., 32. / 90., 12. / 90., 32. / 90., 7. / 90.};

  for (unsigned int i = 0; i < this->size(); ++i)
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i]           = wts[i];
    }
}



template <>
QWeddle<1>::QWeddle()
  : Quadrature<1>(7)
{
  static const double xpts[] = {
    0.0, 1. / 6., 1. / 3., .5, 2. / 3., 5. / 6., 1.0};
  static const double wts[] = {41. / 840.,
                               216. / 840.,
                               27. / 840.,
                               272. / 840.,
                               27. / 840.,
                               216. / 840.,
                               41. / 840.};

  for (unsigned int i = 0; i < this->size(); ++i)
    {
      this->quadrature_points[i] = Point<1>(xpts[i]);
      this->weights[i]           = wts[i];
    }
}


template <>
QGaussLog<1>::QGaussLog(const unsigned int n, const bool revert)
  : Quadrature<1>(n)
{
  std::vector<double> p = get_quadrature_points(n);
  std::vector<double> w = get_quadrature_weights(n);

  for (unsigned int i = 0; i < this->size(); ++i)
    {
      // Using the change of variables x=1-t, it's possible to show
      // that int f(x)ln|1-x| = int f(1-t) ln|t|, which implies that
      // we can use this quadrature formula also with weight ln|1-x|.
      this->quadrature_points[i] =
        revert ? Point<1>(1 - p[n - 1 - i]) : Point<1>(p[i]);
      this->weights[i] = revert ? w[n - 1 - i] : w[i];
    }
}

template <>
std::vector<double>
QGaussLog<1>::get_quadrature_points(const unsigned int n)
{
  std::vector<double> q_points(n);

  switch (n)
    {
      case 1:
        q_points[0] = 0.3333333333333333;
        break;

      case 2:
        q_points[0] = 0.1120088061669761;
        q_points[1] = 0.6022769081187381;
        break;

      case 3:
        q_points[0] = 0.06389079308732544;
        q_points[1] = 0.3689970637156184;
        q_points[2] = 0.766880303938942;
        break;

      case 4:
        q_points[0] = 0.04144848019938324;
        q_points[1] = 0.2452749143206022;
        q_points[2] = 0.5561654535602751;
        q_points[3] = 0.848982394532986;
        break;

      case 5:
        q_points[0] = 0.02913447215197205;
        q_points[1] = 0.1739772133208974;
        q_points[2] = 0.4117025202849029;
        q_points[3] = 0.6773141745828183;
        q_points[4] = 0.89477136103101;
        break;

      case 6:
        q_points[0] = 0.02163400584411693;
        q_points[1] = 0.1295833911549506;
        q_points[2] = 0.3140204499147661;
        q_points[3] = 0.5386572173517997;
        q_points[4] = 0.7569153373774084;
        q_points[5] = 0.922668851372116;
        break;


      case 7:
        q_points[0] = 0.0167193554082585;
        q_points[1] = 0.100185677915675;
        q_points[2] = 0.2462942462079286;
        q_points[3] = 0.4334634932570557;
        q_points[4] = 0.6323509880476823;
        q_points[5] = 0.81111862674023;
        q_points[6] = 0.940848166743287;
        break;

      case 8:
        q_points[0] = 0.01332024416089244;
        q_points[1] = 0.07975042901389491;
        q_points[2] = 0.1978710293261864;
        q_points[3] = 0.354153994351925;
        q_points[4] = 0.5294585752348643;
        q_points[5] = 0.7018145299391673;
        q_points[6] = 0.849379320441094;
        q_points[7] = 0.953326450056343;
        break;

      case 9:
        q_points[0] = 0.01086933608417545;
        q_points[1] = 0.06498366633800794;
        q_points[2] = 0.1622293980238825;
        q_points[3] = 0.2937499039716641;
        q_points[4] = 0.4466318819056009;
        q_points[5] = 0.6054816627755208;
        q_points[6] = 0.7541101371585467;
        q_points[7] = 0.877265828834263;
        q_points[8] = 0.96225055941096;
        break;

      case 10:
        q_points[0] = 0.00904263096219963;
        q_points[1] = 0.05397126622250072;
        q_points[2] = 0.1353118246392511;
        q_points[3] = 0.2470524162871565;
        q_points[4] = 0.3802125396092744;
        q_points[5] = 0.5237923179723384;
        q_points[6] = 0.6657752055148032;
        q_points[7] = 0.7941904160147613;
        q_points[8] = 0.898161091216429;
        q_points[9] = 0.9688479887196;
        break;


      case 11:
        q_points[0]  = 0.007643941174637681;
        q_points[1]  = 0.04554182825657903;
        q_points[2]  = 0.1145222974551244;
        q_points[3]  = 0.2103785812270227;
        q_points[4]  = 0.3266955532217897;
        q_points[5]  = 0.4554532469286375;
        q_points[6]  = 0.5876483563573721;
        q_points[7]  = 0.7139638500230458;
        q_points[8]  = 0.825453217777127;
        q_points[9]  = 0.914193921640008;
        q_points[10] = 0.973860256264123;
        break;

      case 12:
        q_points[0]  = 0.006548722279080035;
        q_points[1]  = 0.03894680956045022;
        q_points[2]  = 0.0981502631060046;
        q_points[3]  = 0.1811385815906331;
        q_points[4]  = 0.2832200676673157;
        q_points[5]  = 0.398434435164983;
        q_points[6]  = 0.5199526267791299;
        q_points[7]  = 0.6405109167754819;
        q_points[8]  = 0.7528650118926111;
        q_points[9]  = 0.850240024421055;
        q_points[10] = 0.926749682988251;
        q_points[11] = 0.977756129778486;
        break;

      default:
        Assert(false, ExcNotImplemented());
        break;
    }

  return q_points;
}


template <>
std::vector<double>
QGaussLog<1>::get_quadrature_weights(const unsigned int n)
{
  std::vector<double> quadrature_weights(n);

  switch (n)
    {
      case 1:
        quadrature_weights[0] = -1.0;
        break;
      case 2:
        quadrature_weights[0] = -0.7185393190303845;
        quadrature_weights[1] = -0.2814606809696154;
        break;

      case 3:
        quadrature_weights[0] = -0.5134045522323634;
        quadrature_weights[1] = -0.3919800412014877;
        quadrature_weights[2] = -0.0946154065661483;
        break;

      case 4:
        quadrature_weights[0] = -0.3834640681451353;
        quadrature_weights[1] = -0.3868753177747627;
        quadrature_weights[2] = -0.1904351269501432;
        quadrature_weights[3] = -0.03922548712995894;
        break;

      case 5:
        quadrature_weights[0] = -0.2978934717828955;
        quadrature_weights[1] = -0.3497762265132236;
        quadrature_weights[2] = -0.234488290044052;
        quadrature_weights[3] = -0.0989304595166356;
        quadrature_weights[4] = -0.01891155214319462;
        break;

      case 6:
        quadrature_weights[0] = -0.2387636625785478;
        quadrature_weights[1] = -0.3082865732739458;
        quadrature_weights[2] = -0.2453174265632108;
        quadrature_weights[3] = -0.1420087565664786;
        quadrature_weights[4] = -0.05545462232488041;
        quadrature_weights[5] = -0.01016895869293513;
        break;


      case 7:
        quadrature_weights[0] = -0.1961693894252476;
        quadrature_weights[1] = -0.2703026442472726;
        quadrature_weights[2] = -0.239681873007687;
        quadrature_weights[3] = -0.1657757748104267;
        quadrature_weights[4] = -0.0889432271377365;
        quadrature_weights[5] = -0.03319430435645653;
        quadrature_weights[6] = -0.005932787015162054;
        break;

      case 8:
        quadrature_weights[0] = -0.164416604728002;
        quadrature_weights[1] = -0.2375256100233057;
        quadrature_weights[2] = -0.2268419844319134;
        quadrature_weights[3] = -0.1757540790060772;
        quadrature_weights[4] = -0.1129240302467932;
        quadrature_weights[5] = -0.05787221071771947;
        quadrature_weights[6] = -0.02097907374214317;
        quadrature_weights[7] = -0.003686407104036044;
        break;

      case 9:
        quadrature_weights[0] = -0.1400684387481339;
        quadrature_weights[1] = -0.2097722052010308;
        quadrature_weights[2] = -0.211427149896601;
        quadrature_weights[3] = -0.1771562339380667;
        quadrature_weights[4] = -0.1277992280331758;
        quadrature_weights[5] = -0.07847890261203835;
        quadrature_weights[6] = -0.0390225049841783;
        quadrature_weights[7] = -0.01386729555074604;
        quadrature_weights[8] = -0.002408041036090773;
        break;

      case 10:
        quadrature_weights[0] = -0.12095513195457;
        quadrature_weights[1] = -0.1863635425640733;
        quadrature_weights[2] = -0.1956608732777627;
        quadrature_weights[3] = -0.1735771421828997;
        quadrature_weights[4] = -0.135695672995467;
        quadrature_weights[5] = -0.0936467585378491;
        quadrature_weights[6] = -0.05578772735275126;
        quadrature_weights[7] = -0.02715981089692378;
        quadrature_weights[8] = -0.00951518260454442;
        quadrature_weights[9] = -0.001638157633217673;
        break;


      case 11:
        quadrature_weights[0]  = -0.1056522560990997;
        quadrature_weights[1]  = -0.1665716806006314;
        quadrature_weights[2]  = -0.1805632182877528;
        quadrature_weights[3]  = -0.1672787367737502;
        quadrature_weights[4]  = -0.1386970574017174;
        quadrature_weights[5]  = -0.1038334333650771;
        quadrature_weights[6]  = -0.06953669788988512;
        quadrature_weights[7]  = -0.04054160079499477;
        quadrature_weights[8]  = -0.01943540249522013;
        quadrature_weights[9]  = -0.006737429326043388;
        quadrature_weights[10] = -0.001152486965101561;
        break;

      case 12:
        quadrature_weights[0]  = -0.09319269144393;
        quadrature_weights[1]  = -0.1497518275763289;
        quadrature_weights[2]  = -0.166557454364573;
        quadrature_weights[3]  = -0.1596335594369941;
        quadrature_weights[4]  = -0.1384248318647479;
        quadrature_weights[5]  = -0.1100165706360573;
        quadrature_weights[6]  = -0.07996182177673273;
        quadrature_weights[7]  = -0.0524069547809709;
        quadrature_weights[8]  = -0.03007108900074863;
        quadrature_weights[9]  = -0.01424924540252916;
        quadrature_weights[10] = -0.004899924710875609;
        quadrature_weights[11] = -0.000834029009809656;
        break;

      default:
        Assert(false, ExcNotImplemented());
        break;
    }

  return quadrature_weights;
}


template <>
QGaussLogR<1>::QGaussLogR(const unsigned int n,
                          const Point<1>     origin,
                          const double       alpha,
                          const bool         factor_out_singularity)
  : Quadrature<1>(
      ((origin[0] == 0) || (origin[0] == 1)) ? (alpha == 1 ? n : 2 * n) : 4 * n)
  , fraction(((origin[0] == 0) || (origin[0] == 1.)) ? 1. : origin[0])
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
  QGauss<1>    quad(n);

  // Check that the origin is inside 0,1
  Assert((fraction >= 0) && (fraction <= 1),
         ExcMessage("Origin is outside [0,1]."));

  // Non singular offset. This is the start of non singular quad
  // points.
  unsigned int ns_offset = (fraction == 1) ? n : 2 * n;

  for (unsigned int i = 0, j = ns_offset; i < n; ++i, ++j)
    {
      // The first i quadrature points are the same as quad1, and
      // are by default singular.
      this->quadrature_points[i] = quad1.point(i) * fraction;
      this->weights[i]           = quad1.weight(i) * fraction;

      // We need to scale with -log|fraction*alpha|
      if ((alpha != 1) || (fraction != 1))
        {
          this->quadrature_points[j] = quad.point(i) * fraction;
          this->weights[j] =
            -std::log(alpha / fraction) * quad.weight(i) * fraction;
        }
      // In case we need the second quadrature as well, do it now.
      if (fraction != 1)
        {
          this->quadrature_points[i + n] =
            quad2.point(i) * (1 - fraction) + Point<1>(fraction);
          this->weights[i + n] = quad2.weight(i) * (1 - fraction);

          // We need to scale with -log|fraction*alpha|
          this->quadrature_points[j + n] =
            quad.point(i) * (1 - fraction) + Point<1>(fraction);
          this->weights[j + n] =
            -std::log(alpha / (1 - fraction)) * quad.weight(i) * (1 - fraction);
        }
    }
  if (factor_out_singularity == true)
    for (unsigned int i = 0; i < size(); ++i)
      {
        Assert(
          this->quadrature_points[i] != origin,
          ExcMessage(
            "The singularity cannot be on a Gauss point of the same order!"));
        double denominator =
          std::log(std::abs((this->quadrature_points[i] - origin)[0]) / alpha);
        Assert(denominator != 0.0,
               ExcMessage(
                 "The quadrature formula you are using does not allow to "
                 "factor out the singularity, which is zero at one point."));
        this->weights[i] /= denominator;
      }
}


template <>
unsigned int
QGaussOneOverR<2>::quad_size(const Point<2> singularity, const unsigned int n)
{
  const double eps = 1e-8;
  const bool   on_edge =
    std::any_of(singularity.begin_raw(),
                singularity.end_raw(),
                [eps](double coord) {
                  return std::abs(coord) < eps || std::abs(coord - 1.) < eps;
                });
  const bool on_vertex =
    on_edge &&
    std::abs((singularity - Point<2>(.5, .5)).norm_square() - .5) < eps;
  if (on_vertex)
    return 2 * n * n;
  else if (on_edge)
    return 4 * n * n;
  else
    return 8 * n * n;
}

template <>
QGaussOneOverR<2>::QGaussOneOverR(const unsigned int n,
                                  const Point<2>     singularity,
                                  const bool         factor_out_singularity)
  : Quadrature<2>(quad_size(singularity, n))
{
  // We treat all the cases in the
  // same way. Split the element in 4
  // pieces, measure the area, if
  // it's relevant, add the
  // quadrature connected to that
  // singularity.
  std::vector<QGaussOneOverR<2>> quads;
  std::vector<Point<2>>          origins;
  // Id of the corner with a
  // singularity
  quads.emplace_back(n, 3, factor_out_singularity);
  quads.emplace_back(n, 2, factor_out_singularity);
  quads.emplace_back(n, 1, factor_out_singularity);
  quads.emplace_back(n, 0, factor_out_singularity);

  origins.emplace_back(0., 0.);
  origins.emplace_back(singularity[0], 0.);
  origins.emplace_back(0., singularity[1]);
  origins.push_back(singularity);

  // Lexicographical ordering.

  double       eps  = 1e-8;
  unsigned int q_id = 0; // Current quad point index.
  Tensor<1, 2> dist;

  for (unsigned int box = 0; box < 4; ++box)
    {
      dist        = (singularity - GeometryInfo<2>::unit_cell_vertex(box));
      dist        = Point<2>(std::abs(dist[0]), std::abs(dist[1]));
      double area = dist[0] * dist[1];
      if (area > eps)
        for (unsigned int q = 0; q < quads[box].size(); ++q, ++q_id)
          {
            const Point<2> &qp = quads[box].point(q);
            this->quadrature_points[q_id] =
              origins[box] + Point<2>(dist[0] * qp[0], dist[1] * qp[1]);
            this->weights[q_id] = quads[box].weight(q) * area;
          }
    }
}


template <>
QGaussOneOverR<2>::QGaussOneOverR(const unsigned int n,
                                  const unsigned int vertex_index,
                                  const bool         factor_out_singularity)
  : Quadrature<2>(2 * n * n)
{
  // This version of the constructor
  // works only for the 4
  // vertices. If you need a more
  // general one, you should use the
  // one with the Point<2> in the
  // constructor.
  AssertIndexRange(vertex_index, 4);

  // Start with the gauss quadrature formula on the (u,v) reference
  // element.
  QGauss<2> gauss(n);

  Assert(gauss.size() == n * n, ExcInternalError());

  // For the moment we only implemented this for the vertices of a
  // quadrilateral. We are planning to do this also for the support
  // points of arbitrary FE_Q elements, to allow the use of this
  // class in boundary element programs with higher order mappings.
  AssertIndexRange(vertex_index, 4);

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
  std::vector<Point<2>> &ps  = this->quadrature_points;
  std::vector<double> &  ws  = this->weights;
  double                 pi4 = numbers::PI / 4;

  for (unsigned int q = 0; q < gauss.size(); ++q)
    {
      const Point<2> &gp = gauss.point(q);
      ps[q][0]           = gp[0];
      ps[q][1]           = gp[0] * std::tan(pi4 * gp[1]);
      ws[q]              = gauss.weight(q) * pi4 / std::cos(pi4 * gp[1]);
      if (factor_out_singularity)
        ws[q] *= (ps[q] - GeometryInfo<2>::unit_cell_vertex(0)).norm();
      // The other half of the quadrilateral is symmetric with
      // respect to xy plane.
      ws[gauss.size() + q]    = ws[q];
      ps[gauss.size() + q][0] = ps[q][1];
      ps[gauss.size() + q][1] = ps[q][0];
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
        theta = numbers::PI / 2;
        break;
      case 2:
        theta = -numbers::PI / 2;
        break;
      case 3:
        theta = numbers::PI;
        break;
    }

  double R00 = std::cos(theta), R01 = -std::sin(theta);
  double R10 = std::sin(theta), R11 = std::cos(theta);

  if (vertex_index != 0)
    for (unsigned int q = 0; q < size(); ++q)
      {
        double x = ps[q][0] - .5, y = ps[q][1] - .5;

        ps[q][0] = R00 * x + R01 * y + .5;
        ps[q][1] = R10 * x + R11 * y + .5;
      }
}


template <int dim>
QSorted<dim>::QSorted(const Quadrature<dim> &quad)
  : Quadrature<dim>(quad)
{
  std::vector<unsigned int> permutation(quad.size());
  for (unsigned int i = 0; i < quad.size(); ++i)
    permutation[i] = i;

  std::sort(permutation.begin(),
            permutation.end(),
            [this](const unsigned int x, const unsigned int y) {
              return this->compare_weights(x, y);
            });

  // At this point, the variable is_tensor_product_flag is set
  // to the respective value of the given Quadrature in the base
  // class copy constructor.
  // We only call a quadrature formula 'tensor product'
  // if the quadrature points are also sorted lexicographically.
  // In particular, any reordering destroys that property
  // and we might need to modify the variable accordingly.
  for (unsigned int i = 0; i < quad.size(); ++i)
    {
      this->weights[i]           = quad.weight(permutation[i]);
      this->quadrature_points[i] = quad.point(permutation[i]);
      if (permutation[i] != i)
        this->is_tensor_product_flag = false;
    }
}


template <int dim>
bool
QSorted<dim>::compare_weights(const unsigned int a, const unsigned int b) const
{
  return (this->weights[a] < this->weights[b]);
}


// construct the quadrature formulae in higher dimensions by
// tensor product of lower dimensions

template <int dim>
QGauss<dim>::QGauss(const unsigned int n)
  : Quadrature<dim>(QGauss<dim - 1>(n), QGauss<1>(n))
{}



template <int dim>
QGaussLobatto<dim>::QGaussLobatto(const unsigned int n)
  : Quadrature<dim>(QGaussLobatto<dim - 1>(n), QGaussLobatto<1>(n))
{}



template <int dim>
QMidpoint<dim>::QMidpoint()
  : Quadrature<dim>(QMidpoint<dim - 1>(), QMidpoint<1>())
{}



template <int dim>
QTrapezoid<dim>::QTrapezoid()
  : Quadrature<dim>(QTrapezoid<dim - 1>(), QTrapezoid<1>())
{}



template <int dim>
QSimpson<dim>::QSimpson()
  : Quadrature<dim>(QSimpson<dim - 1>(), QSimpson<1>())
{}



template <int dim>
QMilne<dim>::QMilne()
  : Quadrature<dim>(QMilne<dim - 1>(), QMilne<1>())
{}


template <int dim>
QWeddle<dim>::QWeddle()
  : Quadrature<dim>(QWeddle<dim - 1>(), QWeddle<1>())
{}

template <int dim>
QTelles<dim>::QTelles(const Quadrature<1> &base_quad,
                      const Point<dim> &   singularity)
  : // We need the explicit implementation if dim == 1. If dim > 1 we use the
    // former implementation and apply a tensorial product to obtain the higher
    // dimensions.
  Quadrature<dim>(
    dim == 2 ?
      QAnisotropic<dim>(QTelles<1>(base_quad, Point<1>(singularity[0])),
                        QTelles<1>(base_quad, Point<1>(singularity[1]))) :
      dim == 3 ?
      QAnisotropic<dim>(QTelles<1>(base_quad, Point<1>(singularity[0])),
                        QTelles<1>(base_quad, Point<1>(singularity[1])),
                        QTelles<1>(base_quad, Point<1>(singularity[2]))) :
      Quadrature<dim>())
{}

template <int dim>
QTelles<dim>::QTelles(const unsigned int n, const Point<dim> &singularity)
  : // In this case we map the standard Gauss Legendre formula using the given
    // singularity point coordinates.
  Quadrature<dim>(QTelles<dim>(QGauss<1>(n), singularity))
{}



template <>
QTelles<1>::QTelles(const Quadrature<1> &base_quad, const Point<1> &singularity)
  : // We explicitly implement the Telles' variable change if dim == 1.
  Quadrature<1>(base_quad)
{
  // We define all the constants to be used in the implementation of
  // Telles' rule
  const double eta_bar  = singularity[0] * 2. - 1.;
  const double eta_star = eta_bar * eta_bar - 1.;
  double       gamma_bar;

  std::vector<Point<1>> quadrature_points_dummy(quadrature_points.size());
  std::vector<double>   weights_dummy(weights.size());
  unsigned int          cont = 0;
  const double          tol  = 1e-10;
  for (unsigned int d = 0; d < quadrature_points.size(); ++d)
    {
      if (std::abs(quadrature_points[d][0] - singularity[0]) > tol)
        {
          quadrature_points_dummy[d - cont] = quadrature_points[d];
          weights_dummy[d - cont]           = weights[d];
        }
      else
        {
          // We need to remove the singularity point from the quadrature point
          // list. To do so we use the variable cont.
          cont = 1;
        }
    }
  if (cont == 1)
    {
      quadrature_points.resize(quadrature_points_dummy.size() - 1);
      weights.resize(weights_dummy.size() - 1);
      for (unsigned int d = 0; d < quadrature_points.size(); ++d)
        {
          quadrature_points[d] = quadrature_points_dummy[d];
          weights[d]           = weights_dummy[d];
        }
    }
  // We need to check if the singularity is at the boundary of the interval.
  if (std::abs(eta_star) <= tol)
    {
      gamma_bar =
        std::pow((eta_bar * eta_star + std::abs(eta_star)), 1.0 / 3.0) +
        std::pow((eta_bar * eta_star - std::abs(eta_star)), 1.0 / 3.0) +
        eta_bar;
    }
  else
    {
      gamma_bar = (eta_bar * eta_star + std::abs(eta_star)) /
                    std::abs(eta_bar * eta_star + std::abs(eta_star)) *
                    std::pow(std::abs(eta_bar * eta_star + std::abs(eta_star)),
                             1.0 / 3.0) +
                  (eta_bar * eta_star - std::abs(eta_star)) /
                    std::abs(eta_bar * eta_star - std::abs(eta_star)) *
                    std::pow(std::abs(eta_bar * eta_star - std::abs(eta_star)),
                             1.0 / 3.0) +
                  eta_bar;
    }
  for (unsigned int q = 0; q < quadrature_points.size(); ++q)
    {
      double gamma = quadrature_points[q][0] * 2 - 1;
      double eta   = (std::pow(gamma - gamma_bar, 3.0) +
                    gamma_bar * (gamma_bar * gamma_bar + 3)) /
                   (1 + 3 * gamma_bar * gamma_bar);

      double J = 3 * ((gamma - gamma_bar) * (gamma - gamma_bar)) /
                 (1 + 3 * gamma_bar * gamma_bar);

      quadrature_points[q][0] = (eta + 1) / 2.0;
      weights[q]              = J * weights[q];
    }
}

namespace internal
{
  namespace QGaussChebyshev
  {
    /**
     * Computes the points of the quadrature formula.
     */
    std::vector<double>
    get_quadrature_points(const unsigned int n)
    {
      std::vector<double> points(n);
      // n point quadrature: index from 0 to n-1
      for (unsigned short i = 0; i < n; ++i)
        // would be cos((2i+1)Pi/(2N+2))
        // put + Pi so we start from the smallest point
        // then map from [-1,1] to [0,1]
        points[i] =
          1. / 2. *
          (1. + std::cos(numbers::PI *
                         (1. + double(2 * i + 1) / double(2 * (n - 1) + 2))));

      return points;
    }



    /**
     * Computes the weights of the quadrature formula.
     */
    std::vector<double>
    get_quadrature_weights(const unsigned int n)
    {
      std::vector<double> weights(n);

      for (unsigned short i = 0; i < n; ++i)
        {
          // same weights as on [-1,1]
          weights[i] = numbers::PI / double(n);
        }

      return weights;
    }
  } // namespace QGaussChebyshev
} // namespace internal


template <>
QGaussChebyshev<1>::QGaussChebyshev(const unsigned int n)
  : Quadrature<1>(n)
{
  Assert(n > 0, ExcMessage("Need at least one point for the quadrature rule"));
  std::vector<double> p = internal::QGaussChebyshev::get_quadrature_points(n);
  std::vector<double> w = internal::QGaussChebyshev::get_quadrature_weights(n);

  for (unsigned int i = 0; i < this->size(); ++i)
    {
      this->quadrature_points[i] = Point<1>(p[i]);
      this->weights[i]           = w[i];
    }
}


template <int dim>
QGaussChebyshev<dim>::QGaussChebyshev(const unsigned int n)
  : Quadrature<dim>(QGaussChebyshev<1>(n))
{}


namespace internal
{
  namespace QGaussRadauChebyshev
  {
    // Computes the points of the quadrature formula.
    std::vector<double>
    get_quadrature_points(const unsigned int                          n,
                          ::dealii::QGaussRadauChebyshev<1>::EndPoint ep)
    {
      std::vector<double> points(n);
      // n point quadrature: index from 0 to n-1
      for (unsigned short i = 0; i < n; ++i)
        // would be -cos(2i Pi/(2N+1))
        // put + Pi so we start from the smallest point
        // then map from [-1,1] to [0,1]
        switch (ep)
          {
            case ::dealii::QGaussRadauChebyshev<1>::left:
              {
                points[i] =
                  1. / 2. *
                  (1. -
                   std::cos(numbers::PI *
                            (1 + 2 * double(i) / (2 * double(n - 1) + 1.))));
                break;
              }

            case ::dealii::QGaussRadauChebyshev<1>::right:
              {
                points[i] =
                  1. / 2. *
                  (1. - std::cos(numbers::PI * (2 * double(n - 1 - i) /
                                                (2 * double(n - 1) + 1.))));
                break;
              }

            default:
              Assert(
                false,
                ExcMessage(
                  "This constructor can only be called with either "
                  "QGaussRadauChebyshev::left or QGaussRadauChebyshev::right as "
                  "second argument."));
          }

      return points;
    }



    // Computes the weights of the quadrature formula.
    std::vector<double>
    get_quadrature_weights(const unsigned int                          n,
                           ::dealii::QGaussRadauChebyshev<1>::EndPoint ep)
    {
      std::vector<double> weights(n);

      for (unsigned short i = 0; i < n; ++i)
        {
          // same weights as on [-1,1]
          weights[i] = 2. * numbers::PI / double(2 * (n - 1) + 1.);
          if (ep == ::dealii::QGaussRadauChebyshev<1>::left && i == 0)
            weights[i] /= 2.;
          else if (ep == ::dealii::QGaussRadauChebyshev<1>::right &&
                   i == (n - 1))
            weights[i] /= 2.;
        }

      return weights;
    }
  } // namespace QGaussRadauChebyshev
} // namespace internal


template <>
QGaussRadauChebyshev<1>::QGaussRadauChebyshev(const unsigned int n, EndPoint ep)
  : Quadrature<1>(n)
  , ep(ep)
{
  Assert(n > 0, ExcMessage("Need at least one point for quadrature rules"));
  std::vector<double> p =
    internal::QGaussRadauChebyshev::get_quadrature_points(n, ep);
  std::vector<double> w =
    internal::QGaussRadauChebyshev::get_quadrature_weights(n, ep);

  for (unsigned int i = 0; i < this->size(); ++i)
    {
      this->quadrature_points[i] = Point<1>(p[i]);
      this->weights[i]           = w[i];
    }
}


template <int dim>
QGaussRadauChebyshev<dim>::QGaussRadauChebyshev(const unsigned int n,
                                                EndPoint           ep)
  : Quadrature<dim>(QGaussRadauChebyshev<1>(
      n,
      static_cast<QGaussRadauChebyshev<1>::EndPoint>(ep)))
  , ep(ep)
{}



namespace internal
{
  namespace QGaussLobattoChebyshev
  {
    // Computes the points of the quadrature formula.
    std::vector<double>
    get_quadrature_points(const unsigned int n)
    {
      std::vector<double> points(n);
      // n point quadrature: index from 0 to n-1
      for (unsigned short i = 0; i < n; ++i)
        // would be cos(i Pi/N)
        // put + Pi so we start from the smallest point
        // then map from [-1,1] to [0,1]
        points[i] =
          1. / 2. *
          (1. + std::cos(numbers::PI * (1 + double(i) / double(n - 1))));

      return points;
    }

    // Computes the weights of the quadrature formula.
    std::vector<double>
    get_quadrature_weights(const unsigned int n)
    {
      std::vector<double> weights(n);

      for (unsigned short i = 0; i < n; ++i)
        {
          // same weights as on [-1,1]
          weights[i] = numbers::PI / double((n - 1));
          if (i == 0 || i == (n - 1))
            weights[i] /= 2.;
        }

      return weights;
    }
  } // namespace QGaussLobattoChebyshev
} // namespace internal



template <>
QGaussLobattoChebyshev<1>::QGaussLobattoChebyshev(const unsigned int n)
  : Quadrature<1>(n)
{
  Assert(n > 1,
         ExcMessage(
           "Need at least two points for Gauss-Lobatto quadrature rule"));
  std::vector<double> p =
    internal::QGaussLobattoChebyshev::get_quadrature_points(n);
  std::vector<double> w =
    internal::QGaussLobattoChebyshev::get_quadrature_weights(n);

  for (unsigned int i = 0; i < this->size(); ++i)
    {
      this->quadrature_points[i] = Point<1>(p[i]);
      this->weights[i]           = w[i];
    }
}


template <int dim>
QGaussLobattoChebyshev<dim>::QGaussLobattoChebyshev(const unsigned int n)
  : Quadrature<dim>(QGaussLobattoChebyshev<1>(n))
{}



template <int dim>
QSimplex<dim>::QSimplex(const Quadrature<dim> &quad)
{
  std::vector<Point<dim>> qpoints;
  std::vector<double>     weights;

  for (unsigned int i = 0; i < quad.size(); ++i)
    {
      double r = 0;
      /* Use "int d" instead of the more natural "unsigned int d" to work
       * around a wrong diagnostic in gcc-10.3.0 that warns about that the
       * comparison "d < dim" is always false in case of "dim == 0".
       * MM 2021 */
      for (int d = 0; d < dim; ++d)
        r += quad.point(i)[d];
      if (r <= 1 + 1e-10)
        {
          this->quadrature_points.push_back(quad.point(i));
          this->weights.push_back(quad.weight(i));
        }
    }
}



template <int dim>
Quadrature<dim>
QSimplex<dim>::compute_affine_transformation(
  const std::array<Point<dim>, dim + 1> &vertices) const
{
  Tensor<2, dim> B;
  for (unsigned int d = 0; d < dim; ++d)
    B[d] = vertices[d + 1] - vertices[0];

  B              = transpose(B);
  const double J = std::abs(determinant(B));

  // if the determinant is zero, we return an empty quadrature
  if (J < 1e-12)
    return Quadrature<dim>();

  std::vector<Point<dim>> qp(this->size());
  std::vector<double>     w(this->size());

  for (unsigned int i = 0; i < this->size(); ++i)
    {
      qp[i] = Point<dim>(vertices[0] + B * this->point(i));
      w[i]  = J * this->weight(i);
    }

  return Quadrature<dim>(qp, w);
}



QTrianglePolar::QTrianglePolar(const Quadrature<1> &radial_quadrature,
                               const Quadrature<1> &angular_quadrature)
  : QSimplex<2>(Quadrature<2>())
{
  const QAnisotropic<2> base(radial_quadrature, angular_quadrature);
  this->quadrature_points.resize(base.size());
  this->weights.resize(base.size());
  for (unsigned int i = 0; i < base.size(); ++i)
    {
      const auto q = base.point(i);
      const auto w = base.weight(i);

      const auto xhat = q[0];
      const auto yhat = q[1];

      const double t  = numbers::PI_2 * yhat;
      const double pi = numbers::PI;
      const double st = std::sin(t);
      const double ct = std::cos(t);
      const double r  = xhat / (st + ct);

      const double J = pi * xhat / (2 * (std::sin(pi * yhat) + 1));

      this->quadrature_points[i] = Point<2>(r * ct, r * st);
      this->weights[i]           = w * J;
    }
}



QTrianglePolar::QTrianglePolar(const unsigned int n)
  : QTrianglePolar(QGauss<1>(n), QGauss<1>(n))
{}



QDuffy::QDuffy(const Quadrature<1> &radial_quadrature,
               const Quadrature<1> &angular_quadrature,
               const double         beta)
  : QSimplex<2>(Quadrature<2>())
{
  const QAnisotropic<2> base(radial_quadrature, angular_quadrature);
  this->quadrature_points.resize(base.size());
  this->weights.resize(base.size());
  for (unsigned int i = 0; i < base.size(); ++i)
    {
      const auto q = base.point(i);
      const auto w = base.weight(i);

      const auto xhat = q[0];
      const auto yhat = q[1];

      const double x = std::pow(xhat, beta) * (1 - yhat);
      const double y = std::pow(xhat, beta) * yhat;

      const double J = beta * std::pow(xhat, 2. * beta - 1.);

      this->quadrature_points[i] = Point<2>(x, y);
      this->weights[i]           = w * J;
    }
}



QDuffy::QDuffy(const unsigned int n, const double beta)
  : QDuffy(QGauss<1>(n), QGauss<1>(n), beta)
{}



template <int dim>
QSplit<dim>::QSplit(const QSimplex<dim> &base, const Point<dim> &split_point)
{
  AssertThrow(GeometryInfo<dim>::is_inside_unit_cell(split_point, 1e-12),
              ExcMessage(
                "The split point should be inside the unit reference cell."));

  std::array<Point<dim>, dim + 1> vertices;
  vertices[0] = split_point;

  // Make a simplex from the split_point and the first dim vertices of each
  // face. In dimension three, we need to split the face in two triangles, so
  // we use once the first dim vertices of each face, and the second time the
  // the dim vertices of each face starting from 1.
  for (auto f : GeometryInfo<dim>::face_indices())
    for (unsigned int start = 0; start < (dim > 2 ? 2 : 1); ++start)
      {
        for (unsigned int i = 0; i < dim; ++i)
          vertices[i + 1] = GeometryInfo<dim>::unit_cell_vertex(
            GeometryInfo<dim>::face_to_cell_vertices(f, i + start));
        const auto quad = base.compute_affine_transformation(vertices);
        if (quad.size())
          {
            this->quadrature_points.insert(this->quadrature_points.end(),
                                           quad.get_points().begin(),
                                           quad.get_points().end());
            this->weights.insert(this->weights.end(),
                                 quad.get_weights().begin(),
                                 quad.get_weights().end());
          }
      }
}



template <int dim>
QGaussSimplex<dim>::QGaussSimplex(const unsigned int n_points_1D)
  : QSimplex<dim>(Quadrature<dim>())
{
  // fill quadrature points and quadrature weights
  if (dim == 0 || dim == 1)
    {
      const dealii::QGauss<dim> quad(n_points_1D);

      this->quadrature_points = quad.get_points();
      this->weights           = quad.get_weights();
    }
  else if (dim == 2)
    {
      if (n_points_1D == 1)
        {
          const double p = 1.0 / 3.0;
          this->quadrature_points.emplace_back(p, p);
          this->weights.emplace_back(0.5);
        }
      else if (n_points_1D == 2)
        {
          const double Q23 = 2.0 / 3.0;
          const double Q16 = 1.0 / 6.0;

          this->quadrature_points.emplace_back(Q23, Q16);
          this->quadrature_points.emplace_back(Q16, Q23);
          this->quadrature_points.emplace_back(Q16, Q16);
          this->weights.emplace_back(Q16);
          this->weights.emplace_back(Q16);
          this->weights.emplace_back(Q16);
        }
      else if (n_points_1D == 3)
        {
          const double q12 = 0.5;

          // clang-format off
            this->quadrature_points.emplace_back(0.3333333333330, 0.3333333333330);
            this->quadrature_points.emplace_back(0.7974269853530, 0.1012865073230);
            this->quadrature_points.emplace_back(0.1012865073230, 0.7974269853530);
            this->quadrature_points.emplace_back(0.1012865073230, 0.1012865073230);
            this->quadrature_points.emplace_back(0.0597158717898, 0.4701420641050);
            this->quadrature_points.emplace_back(0.4701420641050, 0.0597158717898);
            this->quadrature_points.emplace_back(0.4701420641050, 0.4701420641050);
          // clang-format on

          this->weights.emplace_back(q12 * 0.225);
          this->weights.emplace_back(q12 * 0.125939180545);
          this->weights.emplace_back(q12 * 0.125939180545);
          this->weights.emplace_back(q12 * 0.125939180545);
          this->weights.emplace_back(q12 * 0.132394152789);
          this->weights.emplace_back(q12 * 0.132394152789);
          this->weights.emplace_back(q12 * 0.132394152789);
        }
      else if (n_points_1D == 4)
        {
          Quadrature<dim>::operator=(
            QWitherdenVincentSimplex<dim>(n_points_1D));
        }
    }
  else if (dim == 3)
    {
      if (n_points_1D == 1)
        {
          const double Q14 = 1.0 / 4.0;
          const double Q16 = 1.0 / 6.0;

          this->quadrature_points.emplace_back(Q14, Q14, Q14);
          this->weights.emplace_back(Q16);
        }
      else if (n_points_1D == 2)
        {
          const double Q124 = 1.0 / 6.0 / 4.0;

          const double palpha = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
          const double pbeta  = (5.0 - sqrt(5.0)) / 20.0;
          this->quadrature_points.emplace_back(pbeta, pbeta, pbeta);
          this->quadrature_points.emplace_back(palpha, pbeta, pbeta);
          this->quadrature_points.emplace_back(pbeta, palpha, pbeta);
          this->quadrature_points.emplace_back(pbeta, pbeta, palpha);
          this->weights.emplace_back(Q124);
          this->weights.emplace_back(Q124);
          this->weights.emplace_back(Q124);
          this->weights.emplace_back(Q124);
        }
      else if (n_points_1D == 3)
        {
          const double Q16 = 1.0 / 6.0;

          // clang-format off
            this->quadrature_points.emplace_back(0.5684305841968444, 0.1438564719343852, 0.1438564719343852);
            this->quadrature_points.emplace_back(0.1438564719343852, 0.1438564719343852, 0.1438564719343852);
            this->quadrature_points.emplace_back(0.1438564719343852, 0.1438564719343852, 0.5684305841968444);
            this->quadrature_points.emplace_back(0.1438564719343852, 0.5684305841968444, 0.1438564719343852);
            this->quadrature_points.emplace_back(0.0000000000000000, 0.5000000000000000, 0.5000000000000000);
            this->quadrature_points.emplace_back(0.5000000000000000, 0.0000000000000000, 0.5000000000000000);
            this->quadrature_points.emplace_back(0.5000000000000000, 0.5000000000000000, 0.0000000000000000);
            this->quadrature_points.emplace_back(0.5000000000000000, 0.0000000000000000, 0.0000000000000000);
            this->quadrature_points.emplace_back(0.0000000000000000, 0.5000000000000000, 0.0000000000000000);
            this->quadrature_points.emplace_back(0.0000000000000000, 0.0000000000000000, 0.5000000000000000);
          // clang-format on

          this->weights.emplace_back(0.2177650698804054 * Q16);
          this->weights.emplace_back(0.2177650698804054 * Q16);
          this->weights.emplace_back(0.2177650698804054 * Q16);
          this->weights.emplace_back(0.2177650698804054 * Q16);
          this->weights.emplace_back(0.0214899534130631 * Q16);
          this->weights.emplace_back(0.0214899534130631 * Q16);
          this->weights.emplace_back(0.0214899534130631 * Q16);
          this->weights.emplace_back(0.0214899534130631 * Q16);
          this->weights.emplace_back(0.0214899534130631 * Q16);
          this->weights.emplace_back(0.0214899534130631 * Q16);
        }
      else if (n_points_1D == 4)
        {
          Quadrature<dim>::operator=(
            QWitherdenVincentSimplex<dim>(n_points_1D));
        }
    }

  AssertDimension(this->quadrature_points.size(), this->weights.size());
  Assert(this->quadrature_points.size() > 0,
         ExcNotImplemented(
           "QGaussSimplex is currently only implemented for "
           "n_points_1D = 1, 2, 3, and 4 while you are asking for "
           "n_points_1D = " +
           Utilities::to_string(n_points_1D)));
}

namespace
{
  template <std::size_t b_dim>
  std::vector<std::array<double, b_dim>>
  all_permutations(const std::array<double, b_dim> &b_point)
  {
    std::vector<std::array<double, b_dim>> output;

    // We want all possible permutations of the barycentric coordinates.
    // The easiest way to get all of them is to sort the input first and
    // then use next_permutation to cycle through them all.
    std::array<double, b_dim> temp = b_point;
    std::sort(temp.begin(), temp.end());
    do
      {
        output.push_back(temp);
      }
    while (std::next_permutation(temp.begin(), temp.end()));

    return output;
  }
} // namespace



template <int dim>
QWitherdenVincentSimplex<dim>::QWitherdenVincentSimplex(
  const unsigned int n_points_1D)
  : QSimplex<dim>(Quadrature<dim>())
{
  Assert(1 <= dim && dim <= 3, ExcNotImplemented());
  // Just use Gauss in 1D: this is a high-order open rule so this is a
  // reasonable equivalent for generic programming.
  if (dim == 1)
    {
      Quadrature<dim>::operator=(dealii::QGauss<dim>(n_points_1D));
      return;
    }

  std::array<double, dim + 1> centroid;
  std::fill(centroid.begin(), centroid.end(), 1.0 / (dim + 1.0));
  std::vector<std::vector<std::array<double, dim + 1>>> b_point_permutations;
  std::vector<double>                                   b_weights;

  // We can simplify the implementation of these quadrature rules
  // by quite a bit by exploiting symmetry - we do essentially the
  // same thing for each barycentric coordinate, so we can express
  // our quadrature rule as permutations of barycentric points
  // instead of writing things out explicitly.

  // Apply a Barycentric permutation where one point is different.
  auto process_point_1 = [&](const double a, const double w) {
    const double                b = 1.0 - dim * a;
    std::array<double, dim + 1> b_point;
    std::fill(b_point.begin(), b_point.begin() + dim, a);
    b_point[dim] = b;

    b_weights.push_back(w);
    b_point_permutations.push_back(all_permutations(b_point));
  };

  // Apply a Barycentric permutation where two points (in 3D) are different.
  auto process_point_2 = [&](const double a, const double w) {
    Assert(dim == 3, ExcInternalError());
    const double                b = (1.0 - 2.0 * a) / 2.0;
    std::array<double, dim + 1> b_point;
    std::fill(b_point.begin(), b_point.begin() + dim - 1, a);
    b_point[dim - 1] = b;
    b_point[dim]     = b;

    b_weights.push_back(w);
    b_point_permutations.push_back(all_permutations(b_point));
  };

  // Apply a Barycentric permutation where three (or four) points
  // are different (since there are two inputs).
  auto process_point_3 = [&](const double a, const double b, const double w) {
    const double                c = 1.0 - (dim - 1.0) * a - b;
    std::array<double, dim + 1> b_point;
    std::fill(b_point.begin(), b_point.begin() + dim - 1, a);
    b_point[dim - 1] = b;
    b_point[dim]     = c;

    b_weights.push_back(w);
    b_point_permutations.push_back(all_permutations(b_point));
  };

  if (n_points_1D == 1)
    {
      b_point_permutations.push_back({centroid});
      b_weights.push_back(1.0);
    }
  else if (n_points_1D == 2)
    {
      // This is WV-4 in 2D and WV-3 in 3D
      if (dim == 2)
        {
          process_point_1(9.1576213509770743e-02, 1.0995174365532187e-01);
          process_point_1(4.4594849091596489e-01, 2.2338158967801147e-01);
        }
      else if (dim == 3)
        {
          process_point_1(3.281633025163817e-01, 1.362178425370874e-01);
          process_point_1(1.080472498984286e-01, 1.137821574629126e-01);
        }
    }
  else if (n_points_1D == 3)
    {
      // This is the WV-5 rule in both 2D and 3D
      if (dim == 2)
        {
          b_weights.push_back(0.225);
          b_point_permutations.push_back({centroid});

          process_point_1(1.0128650732345634e-01, 1.2593918054482714e-01);
          process_point_1(4.7014206410511511e-01, 1.3239415278850619e-01);
        }
      else if (dim == 3)
        {
          process_point_1(3.108859192633006e-01, 1.126879257180159e-01);
          process_point_1(9.273525031089125e-02, 7.349304311636196e-02);

          process_point_2(4.550370412564964e-02, 4.254602077708147e-02);
        }
    }
  else if (n_points_1D == 4)
    {
      // This is the WV-7 rule in both 2D and 3D
      if (dim == 2)
        {
          process_point_1(3.3730648554587850e-02, 1.6545050110792131e-02);
          process_point_1(4.7430969250471822e-01, 7.7086646185986069e-02);
          process_point_1(2.4157738259540357e-01, 1.2794417123015558e-01);
          process_point_3(4.7036644652595216e-02,
                          1.9868331479735168e-01,
                          5.5878732903199779e-02);
        }
      else if (dim == 3)
        {
          b_point_permutations.push_back({centroid});
          b_weights.push_back(9.548528946413085e-02);

          process_point_1(3.157011497782028e-01, 4.232958120996703e-02);
          process_point_2(5.048982259839635e-02, 3.189692783285758e-02);

          process_point_3(1.888338310260010e-01,
                          5.751716375870000e-01,
                          3.720713072833462e-02);
          process_point_3(2.126547254148314e-02,
                          8.108302410985486e-01,
                          8.110770829903342e-03);
        }
    }
  else if (n_points_1D == 5)
    {
      // This is the WV-9 rule in both 2D and 3D
      if (dim == 2)
        {
          b_point_permutations.push_back({centroid});
          b_weights.push_back(9.7135796282798836e-02);

          process_point_1(4.4729513394452691e-02, 2.5577675658698031e-02);
          process_point_1(4.8968251919873762e-01, 3.1334700227139071e-02);
          process_point_1(4.3708959149293664e-01, 7.7827541004774278e-02);
          process_point_1(1.8820353561903275e-01, 7.9647738927210249e-02);

          process_point_3(3.6838412054736258e-02,
                          2.2196298916076568e-01,
                          4.3283539377289376e-02);
        }
      else if (dim == 3)
        {
          b_point_permutations.push_back({centroid});
          b_weights.push_back(5.801054891248025e-02);

          process_point_1(6.198169755222693e-10, 6.431928175925639e-05);
          process_point_1(1.607745353952616e-01, 2.317333846242546e-02);
          process_point_1(3.222765218214210e-01, 2.956291233542929e-02);
          process_point_1(4.510891834541358e-02, 8.063979979616182e-03);

          process_point_2(1.122965460043761e-01, 3.813408010370246e-02);

          process_point_3(4.588714487524592e-01,
                          2.554579233041310e-03,
                          8.384422198298552e-03);
          process_point_3(3.377587068533860e-02,
                          7.183503264420745e-01,
                          1.023455935274533e-02);
          process_point_3(1.836413698099279e-01,
                          3.441591057817528e-02,
                          2.052491596798814e-02);
        }
    }
  else if (n_points_1D == 6)
    {
      // There is no WV-11 rule in 3D yet
      if (dim == 2)
        {
          b_point_permutations.push_back({centroid});
          b_weights.push_back(8.5761179732224219e-02);

          process_point_1(2.8485417614371900e-02, 1.0431870512894697e-02);
          process_point_1(4.9589190096589092e-01, 1.6606273054585369e-02);
          process_point_1(1.0263548271224643e-01, 3.8630759237019321e-02);
          process_point_1(4.3846592676435220e-01, 6.7316154079468296e-02);
          process_point_1(2.1021995670317828e-01, 7.0515684111716576e-02);

          process_point_3(7.3254276860644785e-03,
                          1.4932478865208237e-01,
                          1.0290289572953278e-02);
          process_point_3(4.6010500165429957e-02,
                          2.8958112563770588e-01,
                          4.0332476640500554e-02);
        }
      else if (dim == 3)
        {
          Assert(false, ExcNotImplemented());
        }
    }
  else
    {
      Assert(false, ExcNotImplemented());
    }

  Assert(b_point_permutations.size() == b_weights.size(), ExcInternalError());
  for (unsigned int permutation_n = 0; permutation_n < b_weights.size();
       ++permutation_n)
    {
      for (const std::array<double, dim + 1> &b_point :
           b_point_permutations[permutation_n])
        {
          const double volume = (dim == 2 ? 1.0 / 2.0 : 1.0 / 6.0);
          this->weights.emplace_back(volume * b_weights[permutation_n]);
          Point<dim> c_point;
          std::copy(b_point.begin(),
                    b_point.begin() + dim,
                    c_point.begin_raw());
          this->quadrature_points.emplace_back(c_point);
        }
    }
}



template <int dim>
QGaussWedge<dim>::QGaussWedge(const unsigned int n_points)
  : Quadrature<dim>()
{
  AssertDimension(dim, 3);

  QGaussSimplex<2> quad_tri(n_points);
  QGauss<1>        quad_line(n_points);

  for (unsigned int i = 0; i < quad_line.size(); ++i)
    for (unsigned int j = 0; j < quad_tri.size(); ++j)
      {
        this->quadrature_points.emplace_back(quad_tri.point(j)[0],
                                             quad_tri.point(j)[1],
                                             quad_line.point(i)[0]);
        this->weights.emplace_back(quad_tri.weight(j) * quad_line.weight(i));
      }

  AssertDimension(this->quadrature_points.size(), this->weights.size());
  Assert(this->quadrature_points.size() > 0,
         ExcMessage("No valid quadrature points!"));
}



template <int dim>
QGaussPyramid<dim>::QGaussPyramid(const unsigned int n_points_1D)
  : Quadrature<dim>()
{
  AssertDimension(dim, 3);

  if (n_points_1D == 1)
    {
      const double Q14 = 1.0 / 4.0;
      const double Q43 = 4.0 / 3.0;

      this->quadrature_points.emplace_back(0, 0, Q14);
      this->weights.emplace_back(Q43);
    }
  else if (n_points_1D == 2)
    {
      // clang-format off
        this->quadrature_points.emplace_back(-0.26318405556971, -0.26318405556971, 0.54415184401122);
        this->quadrature_points.emplace_back(-0.50661630334979, -0.50661630334979, 0.12251482265544);
        this->quadrature_points.emplace_back(-0.26318405556971, +0.26318405556971, 0.54415184401122);
        this->quadrature_points.emplace_back(-0.50661630334979, +0.50661630334979, 0.12251482265544);
        this->quadrature_points.emplace_back(+0.26318405556971, -0.26318405556971, 0.54415184401122);
        this->quadrature_points.emplace_back(+0.50661630334979, -0.50661630334979, 0.12251482265544);
        this->quadrature_points.emplace_back(+0.26318405556971, +0.26318405556971, 0.54415184401122);
        this->quadrature_points.emplace_back(+0.50661630334979, +0.50661630334979, 0.12251482265544);
      // clang-format on

      this->weights.emplace_back(0.10078588207983);
      this->weights.emplace_back(0.23254745125351);
      this->weights.emplace_back(0.10078588207983);
      this->weights.emplace_back(0.23254745125351);
      this->weights.emplace_back(0.10078588207983);
      this->weights.emplace_back(0.23254745125351);
      this->weights.emplace_back(0.10078588207983);
      this->weights.emplace_back(0.23254745125351);
    }

  AssertDimension(this->quadrature_points.size(), this->weights.size());
  Assert(this->quadrature_points.size() > 0,
         ExcMessage("No valid quadrature points!"));
}



// explicit specialization
// note that 1d formulae are specialized by implementation above
template class QGauss<2>;
template class QGaussLobatto<2>;
template class QMidpoint<2>;
template class QTrapezoid<2>;
template class QSimpson<2>;
template class QMilne<2>;
template class QWeddle<2>;

template class QGauss<3>;
template class QGaussLobatto<3>;
template class QMidpoint<3>;
template class QTrapezoid<3>;
template class QSimpson<3>;
template class QMilne<3>;
template class QWeddle<3>;

template class QSorted<1>;
template class QSorted<2>;
template class QSorted<3>;

template class QTelles<1>;
template class QTelles<2>;
template class QTelles<3>;

template class QGaussChebyshev<1>;
template class QGaussChebyshev<2>;
template class QGaussChebyshev<3>;

template class QGaussRadauChebyshev<1>;
template class QGaussRadauChebyshev<2>;
template class QGaussRadauChebyshev<3>;

template class QGaussLobattoChebyshev<1>;
template class QGaussLobattoChebyshev<2>;
template class QGaussLobattoChebyshev<3>;

template class QSimplex<1>;
template class QSimplex<2>;
template class QSimplex<3>;

template class QSplit<1>;
template class QSplit<2>;
template class QSplit<3>;

template class QGaussSimplex<0>;
template class QGaussSimplex<1>;
template class QGaussSimplex<2>;
template class QGaussSimplex<3>;
template class QGaussWedge<1>;
template class QGaussWedge<2>;
template class QGaussWedge<3>;
template class QGaussPyramid<1>;
template class QGaussPyramid<2>;
template class QGaussPyramid<3>;

template class QWitherdenVincentSimplex<1>;
template class QWitherdenVincentSimplex<2>;
template class QWitherdenVincentSimplex<3>;

DEAL_II_NAMESPACE_CLOSE
