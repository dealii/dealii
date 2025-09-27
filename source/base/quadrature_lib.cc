// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

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
        Utilities::pow(2, alpha + beta + 1) * gamma(alpha + q) *
        gamma(beta + q) / ((q - 1) * gamma(q) * gamma(alpha + beta + q + 1));
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


namespace internal
{
  namespace QGaussRadau
  {

    // Implements lookup table after affine transformation to [0,1].
    //
    // Analytical values for [-1,1] and n < 4 listed on
    // https://mathworld.wolfram.com/RadauQuadrature.html
    // Values for n > 3 calculated with the Julia Package
    // FastGaussQuadrature.jl
    // https://github.com/JuliaApproximation/FastGaussQuadrature.jl
    //
    std::vector<double>
    get_left_quadrature_points(const unsigned int n)
    {
      std::vector<double> q_points(n);
      switch (n)
        {
          case 1:
            q_points[0] = 0.;
            break;
          case 2:
            q_points[0] = 0.;
            q_points[1] = 2. / 3.;
            break;
          case 3:
            q_points[0] = 0.;
            q_points[1] = (6. - std::sqrt(6)) * 0.1;
            q_points[2] = (6. + std::sqrt(6)) * 0.1;
            break;

          case 4:
            q_points[0] = 0.000000000000000000;
            q_points[1] = 0.212340538239152943;
            q_points[2] = 0.590533135559265343;
            q_points[3] = 0.911412040487296071;
            break;
          case 5:
            q_points[0] = 0.000000000000000000;
            q_points[1] = 0.139759864343780571;
            q_points[2] = 0.416409567631083166;
            q_points[3] = 0.723156986361876197;
            q_points[4] = 0.942895803885482331;
            break;
          case 6:
            q_points[0] = 0.000000000000000000;
            q_points[1] = 0.098535085798826416;
            q_points[2] = 0.304535726646363913;
            q_points[3] = 0.562025189752613841;
            q_points[4] = 0.801986582126391845;
            q_points[5] = 0.960190142948531222;
            break;
          case 7:
            q_points[0] = 0.000000000000000000;
            q_points[1] = 0.073054328680258851;
            q_points[2] = 0.230766137969945495;
            q_points[3] = 0.441328481228449865;
            q_points[4] = 0.663015309718845702;
            q_points[5] = 0.851921400331515644;
            q_points[6] = 0.970683572840215114;
            break;
          case 8:
            q_points[0] = 0.000000000000000000;
            q_points[1] = 0.056262560536922135;
            q_points[2] = 0.180240691736892389;
            q_points[3] = 0.352624717113169672;
            q_points[4] = 0.547153626330555420;
            q_points[5] = 0.734210177215410598;
            q_points[6] = 0.885320946839095790;
            q_points[7] = 0.977520613561287499;
            break;
          default:
            DEAL_II_NOT_IMPLEMENTED();
            break;
        }
      return q_points;
    }

    std::vector<double>
    get_quadrature_points(const unsigned int                       n,
                          const ::dealii::QGaussRadau<1>::EndPoint end_point)
    {
      std::vector<double> left_points = get_left_quadrature_points(n);
      switch (end_point)
        {
          case ::dealii::QGaussRadau<1>::EndPoint::left:
            return left_points;
          case ::dealii::QGaussRadau<1>::EndPoint::right:
            {
              std::vector<double> points(n);
              for (unsigned int i = 0; i < n; ++i)
                {
                  points[n - i - 1] = 1. - left_points[i];
                }
              return points;
            }
          default:
            Assert(
              false,
              ExcMessage(
                "This constructor can only be called with either "
                "QGaussRadau::left or QGaussRadau::right as second argument."));
            return {};
        }
    }

    // Implements lookup table after affine transformation to [0,1].
    //
    // Analytical values for [-1,1] and n < 4 listed on
    // https://mathworld.wolfram.com/RadauQuadrature.html
    // Values for n > 3 calculated with the Julia Package
    // FastGaussQuadrature.jl
    // https://github.com/JuliaApproximation/FastGaussQuadrature.jl
    //
    std::vector<double>
    get_left_quadrature_weights(const unsigned int n)
    {
      std::vector<double> weights(n);
      switch (n)
        {
          case 1:
            weights[0] = 1.;
            break;
          case 2:
            weights[0] = 0.25;
            weights[1] = 0.75;
            break;
          case 3:
            weights[0] = 1. / 9.;
            weights[1] = (16. + std::sqrt(6)) / 36.;
            weights[2] = (16. - std::sqrt(6)) / 36.;
            break;
          case 4:
            weights[0] = 0.062500000000000000;
            weights[1] = 0.328844319980059696;
            weights[2] = 0.388193468843171852;
            weights[3] = 0.220462211176768369;
            break;
          case 5:
            weights[0] = 0.040000000000000001;
            weights[1] = 0.223103901083570894;
            weights[2] = 0.311826522975741427;
            weights[3] = 0.281356015149462124;
            weights[4] = 0.143713560791225797;
            break;
          case 6:
            weights[0] = 0.027777777777777776;
            weights[1] = 0.159820376610255471;
            weights[2] = 0.242693594234484888;
            weights[3] = 0.260463391594787597;
            weights[4] = 0.208450667155953895;
            weights[5] = 0.100794192626740456;
            break;
          case 7:
            weights[0] = 0.020408163265306121;
            weights[1] = 0.119613744612656100;
            weights[2] = 0.190474936822115581;
            weights[3] = 0.223554914507283209;
            weights[4] = 0.212351889502977870;
            weights[5] = 0.159102115733650767;
            weights[6] = 0.074494235556010341;
            break;
          case 8:
            weights[0] = 0.015625000000000000;
            weights[1] = 0.092679077401489660;
            weights[2] = 0.152065310323392683;
            weights[3] = 0.188258772694559262;
            weights[4] = 0.195786083726246729;
            weights[5] = 0.173507397817250691;
            weights[6] = 0.124823950664932445;
            weights[7] = 0.057254407372128648;
            break;

          default:
            DEAL_II_NOT_IMPLEMENTED();
            break;
        }
      return weights;
    }

    std::vector<double>
    get_quadrature_weights(const unsigned int                       n,
                           const ::dealii::QGaussRadau<1>::EndPoint end_point)
    {
      std::vector<double> left_weights = get_left_quadrature_weights(n);
      switch (end_point)
        {
          case ::dealii::QGaussRadau<1>::EndPoint::left:
            return left_weights;
          case ::dealii::QGaussRadau<1>::EndPoint::right:
            {
              std::vector<double> weights(n);
              for (unsigned int i = 0; i < n; ++i)
                {
                  weights[n - i - 1] = left_weights[i];
                }
              return weights;
            }
          default:
            Assert(false,
                   ExcMessage(
                     "This constructor can only be called with either "
                     "QGaussRadau::EndPoint::left or "
                     "QGaussRadau::EndPoint::right as second argument."));
            return {};
        }
    }
  } // namespace QGaussRadau
} // namespace internal

#ifndef DOXYGEN
template <>
QGaussRadau<1>::QGaussRadau(const unsigned int n, const EndPoint end_point)
  : Quadrature<1>(n)
  , end_point(end_point)
{
  Assert(n > 0, ExcMessage("Need at least one point for quadrature rules."));
  std::vector<double> p =
    internal::QGaussRadau::get_quadrature_points(n, end_point);
  std::vector<double> w =
    internal::QGaussRadau::get_quadrature_weights(n, end_point);

  for (unsigned int i = 0; i < this->size(); ++i)
    {
      this->quadrature_points[i] = dealii::Point<1>(p[i]);
      this->weights[i]           = w[i];
    }
}
#endif


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
        DEAL_II_NOT_IMPLEMENTED();
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
        DEAL_II_NOT_IMPLEMENTED();
        break;
    }

  return quadrature_weights;
}


template <>
QGaussLogR<1>::QGaussLogR(const unsigned int n,
                          const Point<1>    &origin,
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
  const QGaussLog<1> quad1(n, origin[0] != 0);
  const QGaussLog<1> quad2(n);
  const QGauss<1>    quad(n);

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
QGaussOneOverR<2>::quad_size(const Point<2> &singularity, const unsigned int n)
{
  const double eps     = 1e-8;
  bool         on_edge = false;
  for (unsigned int d = 0; d < 2; ++d)
    on_edge = on_edge || (std::abs(singularity[d]) < eps ||
                          std::abs(singularity[d] - 1.0) < eps);
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
                                  const Point<2>    &singularity,
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
  const QGauss<2> gauss(n);

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
  std::vector<double>   &ws  = this->weights;
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
QGaussRadau<dim>::QGaussRadau(const unsigned int n, EndPoint end_point)
  : Quadrature<dim>(
      QGaussRadau<1>(n, static_cast<QGaussRadau<1>::EndPoint>(end_point)))
  , end_point(end_point)
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
                      const Point<dim>    &singularity)
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
      gamma_bar = std::cbrt(eta_bar * eta_star + std::abs(eta_star)) +
                  std::cbrt(eta_bar * eta_star - std::abs(eta_star)) + eta_bar;
    }
  else
    {
      gamma_bar =
        (eta_bar * eta_star + std::abs(eta_star)) /
          std::abs(eta_bar * eta_star + std::abs(eta_star)) *
          std::cbrt(std::abs(eta_bar * eta_star + std::abs(eta_star))) +
        (eta_bar * eta_star - std::abs(eta_star)) /
          std::abs(eta_bar * eta_star - std::abs(eta_star)) *
          std::cbrt(std::abs(eta_bar * eta_star - std::abs(eta_star))) +
        eta_bar;
    }
  for (unsigned int q = 0; q < quadrature_points.size(); ++q)
    {
      double gamma = quadrature_points[q][0] * 2 - 1;
      double eta   = (Utilities::fixed_power<3>(gamma - gamma_bar) +
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
    get_quadrature_points(
      const unsigned int                                n,
      const ::dealii::QGaussRadauChebyshev<1>::EndPoint end_point)
    {
      std::vector<double> points(n);
      // n point quadrature: index from 0 to n-1
      for (unsigned short i = 0; i < n; ++i)
        // would be -cos(2i Pi/(2N+1))
        // put + Pi so we start from the smallest point
        // then map from [-1,1] to [0,1]
        switch (end_point)
          {
            case ::dealii::QGaussRadauChebyshev<1>::EndPoint::left:
              {
                points[i] =
                  1. / 2. *
                  (1. -
                   std::cos(numbers::PI *
                            (1 + 2 * double(i) / (2 * double(n - 1) + 1.))));
                break;
              }

            case ::dealii::QGaussRadauChebyshev<1>::EndPoint::right:
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
                  "QGaussRadauChebyshev::EndPoint::left or "
                  "QGaussRadauChebyshev::EndPoint:right as second argument."));
          }

      return points;
    }



    // Computes the weights of the quadrature formula.
    std::vector<double>
    get_quadrature_weights(
      const unsigned int                                n,
      const ::dealii::QGaussRadauChebyshev<1>::EndPoint end_point)
    {
      std::vector<double> weights(n);

      for (unsigned short i = 0; i < n; ++i)
        {
          // same weights as on [-1,1]
          weights[i] = 2. * numbers::PI / double(2 * (n - 1) + 1.);
          if (end_point == ::dealii::QGaussRadauChebyshev<1>::EndPoint::left &&
              i == 0)
            weights[i] /= 2.;
          else if (end_point ==
                     ::dealii::QGaussRadauChebyshev<1>::EndPoint::right &&
                   i == (n - 1))
            weights[i] /= 2.;
        }

      return weights;
    }
  } // namespace QGaussRadauChebyshev
} // namespace internal


template <>
QGaussRadauChebyshev<1>::QGaussRadauChebyshev(const unsigned int n,
                                              const EndPoint     end_point)
  : Quadrature<1>(n)
  , end_point(end_point)
{
  Assert(n > 0, ExcMessage("Need at least one point for quadrature rules."));
  std::vector<double> points =
    internal::QGaussRadauChebyshev::get_quadrature_points(n, end_point);
  std::vector<double> new_weights =
    internal::QGaussRadauChebyshev::get_quadrature_weights(n, end_point);

  for (unsigned int i = 0; i < this->size(); ++i)
    {
      this->quadrature_points[i] = Point<1>(points[i]);
      this->weights[i]           = new_weights[i];
    }
}


template <int dim>
QGaussRadauChebyshev<dim>::QGaussRadauChebyshev(const unsigned int n,
                                                EndPoint           end_point)
  : Quadrature<dim>(QGaussRadauChebyshev<1>(
      n,
      static_cast<QGaussRadauChebyshev<1>::EndPoint>(end_point)))
  , end_point(end_point)
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
template <int spacedim>
Quadrature<spacedim>
QSimplex<dim>::compute_affine_transformation(
  const std::array<Point<spacedim>, dim + 1> &vertices) const
{
  Assert(dim <= spacedim,
         ExcMessage("Invalid combination of dim and spacedim ."));
  DerivativeForm<1, spacedim, dim> Bt;
  for (unsigned int d = 0; d < dim; ++d)
    Bt[d] = vertices[d + 1] - vertices[0];

  const auto   B = Bt.transpose();
  const double J = std::abs(B.determinant());

  // if the determinant is zero, we return an empty quadrature
  if (J < 1e-12)
    return Quadrature<spacedim>();

  std::vector<Point<spacedim>> qp(this->size());
  std::vector<double>          w(this->size());

  for (unsigned int i = 0; i < this->size(); ++i)
    {
      qp[i] =
        Point<spacedim>(vertices[0] + apply_transformation(B, this->point(i)));
      w[i] = J * this->weight(i);
    }

  return Quadrature<spacedim>(qp, w);
}



template <int dim>
template <int spacedim>
Quadrature<spacedim>
QSimplex<dim>::mapped_quadrature(
  const std::vector<std::array<Point<spacedim>, dim + 1>> &simplices) const
{
  Assert(!(dim == 1 && spacedim == 1),
         ExcMessage("This function is not supposed to work in 1D-1d case."));
  Assert(dim <= spacedim,
         ExcMessage("Invalid combination of dim and spacedim ."));

  std::vector<Point<spacedim>> qp;
  std::vector<double>          ws;
  for (const auto &simplex : simplices)
    {
      const auto rule = this->compute_affine_transformation(simplex);
      std::transform(rule.get_points().begin(),
                     rule.get_points().end(),
                     std::back_inserter(qp),
                     [&](const Point<spacedim> &p) { return p; });
      std::transform(rule.get_weights().begin(),
                     rule.get_weights().end(),
                     std::back_inserter(ws),
                     [&](const double w) { return w; });
    }
  return Quadrature<spacedim>(qp, ws);
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
      const auto &q = base.point(i);
      const auto  w = base.weight(i);

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
      const auto &q = base.point(i);
      const auto  w = base.weight(i);

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
          // The Hillion 7 scheme, as communicated by quadpy
          //
          // See: Numerical Integration on a Triangle, International Journal for
          // Numerical Methods in Engineering, 1977
          const double Q12 = 1.0 / 2.0;
          this->quadrature_points.emplace_back(0.17855872826361643,
                                               0.1550510257216822);
          this->quadrature_points.emplace_back(0.07503111022260812,
                                               0.6449489742783178);
          this->quadrature_points.emplace_back(0.6663902460147014,
                                               0.1550510257216822);
          this->quadrature_points.emplace_back(0.28001991549907407,
                                               0.6449489742783178);

          this->weights.emplace_back(0.31804138174397717 * Q12);
          this->weights.emplace_back(0.18195861825602283 * Q12);
          this->weights.emplace_back(0.31804138174397717 * Q12);
          this->weights.emplace_back(0.18195861825602283 * Q12);
        }
      else if (n_points_1D == 3)
        {
          // The Hammer-Marlowe-Stroud 5 Scheme, as communicated by quadpy
          const double p0 = 2.0 / 7.0 - std::sqrt(15.0) / 21.0;
          const double p1 = 2.0 / 7.0 + std::sqrt(15.0) / 21.0;
          const double p2 = 3.0 / 7.0 - 2.0 * std::sqrt(15.0) / 21.0;
          const double p3 = 3.0 / 7.0 + 2.0 * std::sqrt(15.0) / 21.0;
          this->quadrature_points.emplace_back(1.0 / 3.0, 1.0 / 3.0);
          this->quadrature_points.emplace_back(p3, p0);
          this->quadrature_points.emplace_back(p0, p3);
          this->quadrature_points.emplace_back(p0, p0);
          this->quadrature_points.emplace_back(p2, p1);
          this->quadrature_points.emplace_back(p1, p2);
          this->quadrature_points.emplace_back(p1, p1);

          const double q12 = 0.5;
          const double w0  = 9.0 / 40.0;
          const double w1  = 31.0 / 240.0 - std::sqrt(15.0) / 1200.0;
          const double w2  = 31.0 / 240.0 + std::sqrt(15.0) / 1200.0;
          this->weights.emplace_back(q12 * w0);
          this->weights.emplace_back(q12 * w1);
          this->weights.emplace_back(q12 * w1);
          this->weights.emplace_back(q12 * w1);
          this->weights.emplace_back(q12 * w2);
          this->weights.emplace_back(q12 * w2);
          this->weights.emplace_back(q12 * w2);
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
      // The Xiao Gimbutas 03 scheme, as communicated by quadpy
      //
      // See: A numerical algorithm for the construction of efficient quadrature
      // rules in two and higher dimensions, Computers & Mathematics with
      // Applications, 2010
      else if (n_points_1D == 2)
        {
          const double Q16 = 1.0 / 6.0;
          this->weights.emplace_back(0.1223220027573451 * Q16);
          this->weights.emplace_back(0.1280664127107469 * Q16);
          this->weights.emplace_back(0.1325680271444452 * Q16);
          this->weights.emplace_back(0.1406244096604032 * Q16);
          this->weights.emplace_back(0.2244151669175574 * Q16);
          this->weights.emplace_back(0.2520039808095023 * Q16);

          this->quadrature_points.emplace_back(0.1620014916985245,
                                               0.1838503504920977,
                                               0.01271836631368145);
          this->quadrature_points.emplace_back(0.01090521221118924,
                                               0.2815238021235462,
                                               0.3621268299455338);
          this->quadrature_points.emplace_back(0.1901170024392839,
                                               0.01140332944455717,
                                               0.3586207204668839);
          this->quadrature_points.emplace_back(0.170816925164989,
                                               0.1528181430909273,
                                               0.6384932999617267);
          this->quadrature_points.emplace_back(0.1586851632274406,
                                               0.5856628056552158,
                                               0.1308471689520965);
          this->quadrature_points.emplace_back(0.5712260521491151,
                                               0.1469183900871696,
                                               0.1403728057942107);
        }
      // Past this point the best rules (positive weights, minimal number of
      // points) we have right now are the Witherden-Vincent ones
      else if (n_points_1D == 3)
        {
          Quadrature<dim>::operator=(
            QWitherdenVincentSimplex<dim>(n_points_1D));
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
  const unsigned int n_points_1D,
  const bool         use_odd_order)
  : QSimplex<dim>(Quadrature<dim>())
{
  Assert(1 <= dim && dim <= 3, ExcNotImplemented());
  // Just use Gauss in 1d: this is a high-order open rule so this is a
  // reasonable equivalent for generic programming.
  if (dim == 1)
    {
      Quadrature<dim>::operator=(QGauss<dim>(n_points_1D));
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
  // Equivalent to d3_aa and s31 in quadpy.
  auto process_point_1 = [&](const double a, const double w) {
    const double                b = 1.0 - dim * a;
    std::array<double, dim + 1> b_point;
    std::fill(b_point.begin(), b_point.begin() + dim, a);
    b_point[dim] = b;

    b_weights.push_back(w);
    b_point_permutations.push_back(all_permutations(b_point));
  };

  // Apply a Barycentric permutation where two points (in 3d) are different.
  // Equivalent to s22 in quadpy.
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
  // Equivalent to d3_ab and s211 in quadpy.
  auto process_point_3 = [&](const double a, const double b, const double w) {
    const double                c = 1.0 - (dim - 1.0) * a - b;
    std::array<double, dim + 1> b_point;
    std::fill(b_point.begin(), b_point.begin() + dim - 1, a);
    b_point[dim - 1] = b;
    b_point[dim]     = c;

    b_weights.push_back(w);
    b_point_permutations.push_back(all_permutations(b_point));
  };

  switch (n_points_1D)
    {
      case 1:
        switch (dim)
          {
            case 2:
              if (use_odd_order)
                {
                  // WV-1, 2d
                  b_point_permutations.push_back({centroid});
                  b_weights.push_back(1.0000000000000000e+00);
                }
              else
                {
                  // WV-2, 2d
                  process_point_1(1.6666666666666669e-01,
                                  3.3333333333333331e-01);
                }
              break;
            case 3:
              if (use_odd_order)
                {
                  // WV-1, 3d
                  b_point_permutations.push_back({centroid});
                  b_weights.push_back(1.0000000000000000e+00);
                }
              else
                {
                  // WV-2, 3d
                  process_point_1(1.3819660112501050e-01,
                                  2.5000000000000000e-01);
                }
              break;
            default:
              DEAL_II_NOT_IMPLEMENTED();
          }
        break;
      case 2:
        switch (dim)
          {
            case 2:
              // WV-4 in both cases (no WV-3 in 2d)
              process_point_1(9.1576213509770743e-02, 1.0995174365532187e-01);
              process_point_1(4.4594849091596489e-01, 2.2338158967801147e-01);
              break;
            case 3:
              if (use_odd_order)
                {
                  // WV-3, 3d
                  process_point_1(3.2816330251638171e-01,
                                  1.3621784253708741e-01);
                  process_point_1(1.0804724989842859e-01,
                                  1.1378215746291261e-01);
                }
              else
                {
                  // WV-5 (no WV-4 in 3d)
                  Quadrature<dim>::operator=(QWitherdenVincentSimplex<dim>(3));
                }
              break;
            default:
              DEAL_II_ASSERT_UNREACHABLE();
          }
        break;
      case 3:
        switch (dim)
          {
            case 2:
              if (use_odd_order)
                {
                  // WV-5, 2d
                  b_point_permutations.push_back({centroid});
                  b_weights.push_back(2.2500000000000001e-01);
                  process_point_1(1.0128650732345634e-01,
                                  1.2593918054482714e-01);
                  process_point_1(4.7014206410511511e-01,
                                  1.3239415278850619e-01);
                }
              else
                {
                  // WV-6, 2d
                  process_point_1(6.3089014491502227e-02,
                                  5.0844906370206819e-02);
                  process_point_1(2.4928674517091043e-01,
                                  1.1678627572637937e-01);
                  process_point_3(5.3145049844816938e-02,
                                  3.1035245103378439e-01,
                                  8.2851075618373571e-02);
                }
              break;
            case 3:
              if (use_odd_order)
                {
                  // WV-5, 3d
                  process_point_1(3.1088591926330061e-01,
                                  1.1268792571801590e-01);
                  process_point_1(9.2735250310891248e-02,
                                  7.3493043116361956e-02);
                  process_point_2(4.5503704125649642e-02,
                                  4.2546020777081472e-02);
                }
              else
                {
                  // WV-6, 3d
                  process_point_1(4.0673958534611372e-02,
                                  1.0077211055320640e-02);
                  process_point_1(3.2233789014227548e-01,
                                  5.5357181543654717e-02);
                  process_point_1(2.1460287125915201e-01,
                                  3.9922750258167487e-02);
                  process_point_3(6.3661001875017442e-02,
                                  6.0300566479164919e-01,
                                  4.8214285714285710e-02);
                }
              break;
            default:
              DEAL_II_ASSERT_UNREACHABLE();
          }
        break;
      case 4:
        switch (dim)
          {
            case 2:
              if (use_odd_order)
                {
                  // WV-7, 2d
                  process_point_1(3.3730648554587850e-02,
                                  1.6545050110792131e-02);
                  process_point_1(4.7430969250471822e-01,
                                  7.7086646185986069e-02);
                  process_point_1(2.4157738259540357e-01,
                                  1.2794417123015558e-01);
                  process_point_3(4.7036644652595216e-02,
                                  1.9868331479735168e-01,
                                  5.5878732903199779e-02);
                }
              else
                {
                  // WV-8, 2d
                  b_point_permutations.push_back({centroid});
                  b_weights.push_back(1.4431560767778717e-01);
                  process_point_1(5.0547228317030957e-02,
                                  3.2458497623198079e-02);
                  process_point_1(4.5929258829272313e-01,
                                  9.5091634267284619e-02);
                  process_point_1(1.7056930775176021e-01,
                                  1.0321737053471824e-01);
                  process_point_3(8.3947774099575878e-03,
                                  2.6311282963463811e-01,
                                  2.7230314174434993e-02);
                }
              break;
            case 3:
              if (use_odd_order)
                {
                  // WV-7, 3d
                  b_point_permutations.push_back({centroid});
                  b_weights.push_back(9.5485289464130846e-02);
                  process_point_1(3.1570114977820279e-01,
                                  4.2329581209967028e-02);
                  process_point_2(5.0489822598396350e-02,
                                  3.1896927832857580e-02);
                  process_point_3(1.8883383102600099e-01,
                                  5.7517163758699996e-01,
                                  3.7207130728334620e-02);
                  process_point_3(2.1265472541483140e-02,
                                  8.1083024109854862e-01,
                                  8.1107708299033420e-03);
                }
              else
                {
                  // WV-8, 3d
                  process_point_1(1.0795272496221089e-01,
                                  2.6426650908408830e-02);
                  process_point_1(1.8510948778258660e-01,
                                  5.2031747563738531e-02);
                  process_point_1(4.2316543684767283e-02,
                                  7.5252561535401989e-03);
                  process_point_1(3.1418170912403898e-01,
                                  4.1763782856934897e-02);
                  process_point_2(4.3559132858383021e-01,
                                  3.6280930261308818e-02);
                  process_point_3(2.1433930127130570e-02,
                                  7.1746406342630831e-01,
                                  7.1569028908444327e-03);
                  process_point_3(2.0413933387602909e-01,
                                  5.8379737830214440e-01,
                                  1.5453486150960340e-02);
                }
              break;
            default:
              DEAL_II_ASSERT_UNREACHABLE();
          }
        break;
      case 5:
        switch (dim)
          {
            case 2:
              if (use_odd_order)
                {
                  // WV-9, 2d
                  b_point_permutations.push_back({centroid});
                  b_weights.push_back(9.7135796282798836e-02);
                  process_point_1(4.4729513394452691e-02,
                                  2.5577675658698031e-02);
                  process_point_1(4.8968251919873762e-01,
                                  3.1334700227139071e-02);
                  process_point_1(4.3708959149293664e-01,
                                  7.7827541004774278e-02);
                  process_point_1(1.8820353561903275e-01,
                                  7.9647738927210249e-02);
                  process_point_3(3.6838412054736258e-02,
                                  2.2196298916076568e-01,
                                  4.3283539377289376e-02);
                }
              else
                {
                  // WV-10, 2d
                  b_point_permutations.push_back({centroid});
                  b_weights.push_back(8.1743329146285973e-02);
                  process_point_1(3.2055373216943517e-02,
                                  1.3352968813149567e-02);
                  process_point_1(1.4216110105656438e-01,
                                  4.5957963604744731e-02);
                  process_point_3(2.8367665339938453e-02,
                                  1.6370173373718250e-01,
                                  2.5297757707288385e-02);
                  process_point_3(2.9619889488729734e-02,
                                  3.6914678182781102e-01,
                                  3.4184648162959429e-02);
                  process_point_3(1.4813288578382056e-01,
                                  3.2181299528883545e-01,
                                  6.3904906396424044e-02);
                }
              break;
            case 3:
              if (use_odd_order)
                {
                  // WV-9, 3d
                  b_point_permutations.push_back({centroid});
                  b_weights.push_back(5.8010548912480253e-02);
                  process_point_1(6.1981697552226933e-10,
                                  6.4319281759256394e-05);
                  process_point_1(1.6077453539526160e-01,
                                  2.3173338462425461e-02);
                  process_point_1(3.2227652182142102e-01,
                                  2.9562912335429289e-02);
                  process_point_1(4.5108918345413578e-02,
                                  8.0639799796161822e-03);
                  process_point_2(1.1229654600437609e-01,
                                  3.8134080103702457e-02);
                  process_point_3(4.5887144875245922e-01,
                                  2.5545792330413102e-03,
                                  8.3844221982985519e-03);
                  process_point_3(3.3775870685338598e-02,
                                  7.1835032644207453e-01,
                                  1.0234559352745330e-02);
                  process_point_3(1.8364136980992790e-01,
                                  3.4415910578175279e-02,
                                  2.0524915967988139e-02);
                }
              else
                {
                  // WV-10, 3d
                  b_point_permutations.push_back({centroid});
                  b_weights.push_back(4.7399773556020743e-02);
                  process_point_1(3.1225006869518868e-01,
                                  2.6937059992268701e-02);
                  process_point_1(1.1430965385734609e-01,
                                  9.8691597167933822e-03);
                  process_point_3(4.1043073921896539e-01,
                                  1.6548602561961109e-01,
                                  1.1393881220195230e-02);
                  process_point_3(6.1380088247906528e-03,
                                  9.4298876734520487e-01,
                                  3.6194434433925362e-04);
                  process_point_3(1.2105018114558939e-01,
                                  4.7719037990428043e-01,
                                  2.5739731980456069e-02);
                  process_point_3(3.2779468216442620e-02,
                                  5.9425626948000698e-01,
                                  1.0135871679755789e-02);
                  process_point_3(3.2485281564823047e-02,
                                  8.0117728465834437e-01,
                                  6.5761472770359038e-03);
                  process_point_3(1.7497934218393901e-01,
                                  6.2807184547536599e-01,
                                  1.2907035798861989e-02);
                }
              break;
            default:
              DEAL_II_NOT_IMPLEMENTED();
          }
        break;
      case 6:
        // There is no WV-11 rule in 3d yet
        Assert(dim == 2, ExcNotImplemented());
        if (use_odd_order)
          {
            // WV-11, 2d
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
        else
          {
            // WV-12, 2d
            process_point_1(2.4646363436335583e-02, 7.9316425099736389e-03);
            process_point_1(4.8820375094554153e-01, 2.4266838081452032e-02);
            process_point_1(1.0925782765935427e-01, 2.8486052068877544e-02);
            process_point_1(4.4011164865859309e-01, 4.9918334928060942e-02);
            process_point_1(2.7146250701492608e-01, 6.2541213195902765e-02);
            process_point_3(2.1382490256170616e-02,
                            1.2727971723358933e-01,
                            1.5083677576511438e-02);
            process_point_3(2.3034156355267121e-02,
                            2.9165567973834094e-01,
                            2.1783585038607559e-02);
            process_point_3(1.1629601967792658e-01,
                            2.5545422863851736e-01,
                            4.3227363659414209e-02);
          }
        break;
      case 7:
        // There is no WV-13 rule in 3d yet
        Assert(dim == 2, ExcNotImplemented());
        if (use_odd_order)
          {
            // WV-13, 2d
            b_point_permutations.push_back({centroid});
            b_weights.push_back(6.7960036586831640e-02);
            process_point_1(2.1509681108843159e-02, 6.0523371035391717e-03);
            process_point_1(4.8907694645253935e-01, 2.3994401928894731e-02);
            process_point_1(4.2694141425980042e-01, 5.5601967530453329e-02);
            process_point_1(2.2137228629183292e-01, 5.8278485119199981e-02);
            process_point_3(5.1263891023823893e-03,
                            2.7251581777342970e-01,
                            9.5906810035432631e-03);
            process_point_3(2.4370186901093827e-02,
                            1.1092204280346341e-01,
                            1.4965401105165668e-02);
            process_point_3(8.7895483032197297e-02,
                            1.6359740106785048e-01,
                            2.4179039811593819e-02);
            process_point_3(6.8012243554206653e-02,
                            3.0844176089211778e-01,
                            3.4641276140848373e-02);
          }
        else
          {
            // WV-14, 2d
            process_point_1(1.9390961248701044e-02, 4.9234036024000819e-03);
            process_point_1(6.1799883090872587e-02, 1.4433699669776668e-02);
            process_point_1(4.8896391036217862e-01, 2.1883581369428889e-02);
            process_point_1(4.1764471934045394e-01, 3.2788353544125348e-02);
            process_point_1(1.7720553241254344e-01, 4.2162588736993016e-02);
            process_point_1(2.7347752830883865e-01, 5.1774104507291585e-02);
            process_point_3(1.2683309328720416e-03,
                            1.1897449769695684e-01,
                            5.0102288385006719e-03);
            process_point_3(1.4646950055654417e-02,
                            2.9837288213625779e-01,
                            1.4436308113533840e-02);
            process_point_3(5.7124757403647919e-02,
                            1.7226668782135557e-01,
                            2.4665753212563674e-02);
            process_point_3(9.2916249356971847e-02,
                            3.3686145979634496e-01,
                            3.8571510787060684e-02);
          }
        break;
      default:
        DEAL_II_NOT_IMPLEMENTED();
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
          for (int d = 0; d < dim; ++d)
            c_point[d] = b_point[d];
          this->quadrature_points.emplace_back(c_point);
        }
    }
}



namespace
{
  template <int dim>
  Quadrature<dim>
  setup_qiterated_1D(const Quadrature<dim> &, const unsigned int)
  {
    DEAL_II_ASSERT_UNREACHABLE();
    return Quadrature<dim>();
  }



  Quadrature<1>
  setup_qiterated_1D(const Quadrature<1> &base_quad,
                     const unsigned int   n_copies)
  {
    return QIterated<1>(base_quad, n_copies);
  }
} // namespace



template <int dim>
QIteratedSimplex<dim>::QIteratedSimplex(const Quadrature<dim> &base_quad,
                                        const unsigned int     n_copies)
{
  switch (dim)
    {
      case 1:
        static_cast<Quadrature<dim> &>(*this) =
          setup_qiterated_1D(base_quad, n_copies);
        break;
      case 2:
      case 3:
        {
          const auto n_refinements =
            static_cast<unsigned int>(std::round(std::log2(n_copies)));
          Assert((1u << n_refinements) == n_copies,
                 ExcMessage("The number of copies must be a power of 2."));
          Triangulation<dim> tria;
          const auto reference_cell = ReferenceCells::get_simplex<dim>();
          GridGenerator::reference_cell(tria, reference_cell);
          tria.refine_global(n_refinements);
          const Mapping<dim> &mapping =
            reference_cell.template get_default_linear_mapping<dim>();
          FE_Nothing<dim> fe(reference_cell);

          FEValues<dim>           fe_values(mapping,
                                  fe,
                                  base_quad,
                                  update_quadrature_points | update_JxW_values);
          std::vector<Point<dim>> points;
          std::vector<double>     weights;
          for (const auto &cell : tria.active_cell_iterators())
            {
              fe_values.reinit(cell);
              for (unsigned int qp = 0; qp < base_quad.size(); ++qp)
                {
                  points.push_back(fe_values.quadrature_point(qp));
                  weights.push_back(fe_values.JxW(qp));
                }
            }

          static_cast<Quadrature<dim> &>(*this) =
            Quadrature<dim>(points, weights);

          break;
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim>
QGaussWedge<dim>::QGaussWedge(const unsigned int n_points)
  : Quadrature<dim>()
{
  AssertDimension(dim, 3);

  const QGaussSimplex<2> quad_tri(n_points);
  const QGauss<1>        quad_line(n_points);

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
template class QGaussRadau<2>;
template class QGaussLobatto<2>;
template class QMidpoint<2>;
template class QTrapezoid<2>;
template class QSimpson<2>;
template class QMilne<2>;
template class QWeddle<2>;

template class QGauss<3>;
template class QGaussRadau<3>;
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

template class QIteratedSimplex<1>;
template class QIteratedSimplex<2>;
template class QIteratedSimplex<3>;

template class QSplit<1>;
template class QSplit<2>;
template class QSplit<3>;

template class QGaussSimplex<0>;
template class QGaussSimplex<1>;
template class QGaussSimplex<2>;
template class QGaussSimplex<3>;
template class QGaussWedge<0>;
template class QGaussWedge<1>;
template class QGaussWedge<2>;
template class QGaussWedge<3>;
template class QGaussPyramid<0>;
template class QGaussPyramid<1>;
template class QGaussPyramid<2>;
template class QGaussPyramid<3>;

template class QWitherdenVincentSimplex<1>;
template class QWitherdenVincentSimplex<2>;
template class QWitherdenVincentSimplex<3>;

#ifndef DOXYGEN
template Quadrature<1>
QSimplex<1>::compute_affine_transformation(
  const std::array<Point<1>, 1 + 1> &vertices) const;

template Quadrature<2>
QSimplex<1>::compute_affine_transformation(
  const std::array<Point<2>, 1 + 1> &vertices) const;

template Quadrature<2>
QSimplex<2>::compute_affine_transformation(
  const std::array<Point<2>, 2 + 1> &vertices) const;

template Quadrature<3>
QSimplex<1>::compute_affine_transformation(
  const std::array<Point<3>, 1 + 1> &vertices) const;

template Quadrature<3>
QSimplex<2>::compute_affine_transformation(
  const std::array<Point<3>, 2 + 1> &vertices) const;

template Quadrature<3>
QSimplex<3>::compute_affine_transformation(
  const std::array<Point<3>, 3 + 1> &vertices) const;

template Quadrature<2>
QSimplex<1>::mapped_quadrature(
  const std::vector<std::array<Point<2>, 1 + 1>> &simplices) const;

template Quadrature<3>
QSimplex<1>::mapped_quadrature(
  const std::vector<std::array<Point<3>, 1 + 1>> &simplices) const;

template Quadrature<2>
QSimplex<2>::mapped_quadrature(
  const std::vector<std::array<Point<2>, 2 + 1>> &simplices) const;

template Quadrature<3>
QSimplex<2>::mapped_quadrature(
  const std::vector<std::array<Point<3>, 2 + 1>> &simplices) const;

template Quadrature<3>
QSimplex<3>::mapped_quadrature(
  const std::vector<std::array<Point<3>, 3 + 1>> &simplices) const;
#endif

DEAL_II_NAMESPACE_CLOSE
