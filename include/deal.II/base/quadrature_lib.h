// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
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

#ifndef __deal2__quadrature_lib_h
#define __deal2__quadrature_lib_h


#include <deal.II/base/config.h>
#include <deal.II/base/quadrature.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup Quadrature */
/*@{*/

/**
 * The Gauss-Legendre family of quadrature rules for numerical
 * integration.
 *
 * The coefficients of these quadrature rules are computed by the function
 * described in <a
 * href="http://en.wikipedia.org/wiki/Numerical_Recipes">Numerical
 * Recipies</a>.
 *
 * @author Guido Kanschat, 2001
 */
template <int dim>
class QGauss : public Quadrature<dim>
{
public:
  /**
   * Generate a formula with
   * <tt>n</tt> quadrature points (in
   * each space direction), exact for
   * polynomials of degree
   * <tt>2n-1</tt>.
   */
  QGauss (const unsigned int n);
};


/**
 * The Gauss-Lobatto family of quadrature rules for numerical
 * integration.
 *
 * This modification of the Gauss quadrature uses the two interval end
 * points as well. Being exact for polynomials of degree <i>2n-3</i>,
 * this formula is suboptimal by two degrees.
 *
 * The quadrature points are interval end points plus the roots of
 * the derivative of the Legendre polynomial <i>P<sub>n-1</sub></i> of
 * degree <i>n-1</i>. The quadrature weights are
 * <i>2/(n(n-1)(P<sub>n-1</sub>(x<sub>i</sub>)<sup>2</sup>)</i>.
 *
 * @note This implementation has not been optimized concerning
 *   numerical stability and efficiency. It can be easily adapted
 *   to the general case of Gauss-Lobatto-Jacobi-Bouzitat quadrature
 *   with arbitrary parameters $\alpha$, $\beta$, of which
 *   the Gauss-Lobatto-Legendre quadrature ($\alpha = \beta = 0$)
 *   is a special case.
 *
 * @sa http://en.wikipedia.org/wiki/Handbook_of_Mathematical_Functions
 * @sa Karniadakis, G.E. and Sherwin, S.J.:
 *     Spectral/hp element methods for computational fluid dynamics.
 *     Oxford: Oxford University Press, 2005
 *
 * @author Guido Kanschat, 2005, 2006; F. Prill, 2006
 */
template <int dim>
class QGaussLobatto : public Quadrature<dim>
{
public:
  /**
   * Generate a formula with
   * <tt>n</tt> quadrature points
   * (in each space direction).
   */
  QGaussLobatto(const unsigned int n);

protected:
  /**
   * Compute Legendre-Gauss-Lobatto
   * quadrature points in the
   * interval $[-1, +1]$. They are
   * equal to the roots of the
   * corresponding Jacobi
   * polynomial (specified by @p
   * alpha, @p beta).  @p q is the
   * number of points.
   *
   * @return Vector containing nodes.
   */
  std::vector<long double>
  compute_quadrature_points (const unsigned int q,
                             const int alpha,
                             const int beta) const;

  /**
   * Compute Legendre-Gauss-Lobatto quadrature
   * weights.
   * The quadrature points and weights are
   * related to Jacobi polynomial specified
   * by @p alpha, @p beta.
   * @p x denotes the quadrature points.
   *
   * @return Vector containing weights.
   */
  std::vector<long double>
  compute_quadrature_weights (const std::vector<long double> &x,
                              const int alpha,
                              const int beta) const;

  /**
   * Evaluate a Jacobi polynomial
   * $ P^{\alpha, \beta}_n(x) $
   * specified by the parameters
   * @p alpha, @p beta, @p n.
   * Note: The Jacobi polynomials are
   * not orthonormal and defined on
   * the interval $[-1, +1]$.
   * @p x is the point of evaluation.
   */
  long double JacobiP(const long double x,
                      const int alpha,
                      const int beta,
                      const unsigned int n) const;

  /**
   * Evaluate the Gamma function
   * $ \Gamma(n) = (n-1)! $.
   * @param n  point of evaluation (integer).
   */
  long double gamma(const unsigned int n) const;
};



/**
 * The midpoint rule for numerical quadrature. This one-point formula is exact
 * for linear polynomials.
 */
template <int dim>
class QMidpoint : public Quadrature<dim>
{
public:
  QMidpoint ();
};


/**
 * The Simpson rule for numerical quadrature. This formula with 3 quadrature
 * points is exact for polynomials of degree 3.
 */
template <int dim>
class QSimpson : public Quadrature<dim>
{
public:
  QSimpson ();
};



/**
 * The trapezoidal rule for numerical quadrature. This formula with two
 * quadrature points is exact for linear polynomials.
 *
 * The class is poorly named since the proper name of the quadrature formula
 * is "trapezoidal rule", or sometimes also called the "trapezoid rule". The
 * misnomer results from the fact that its original authors's poor English
 * language skills led them to translate the name incorrectly from the German
 * "Trapezregel".
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class QTrapez : public Quadrature<dim>
{
public:
  QTrapez ();
};



/**
 * The Milne rule for numerical quadrature formula. The Milne rule is a
 * closed Newton-Cotes formula and is exact for polynomials of degree 5.
 *
 * @sa Stoer: Einführung in die Numerische Mathematik I, p. 102
 */
template <int dim>
class QMilne : public Quadrature<dim>
{
public:
  QMilne ();
};


/**
 * The Weddle rule for numerical quadrature. The Weddle rule is a
 * closed Newton-Cotes formula and is exact for polynomials of degree 7.
 *
 * @sa Stoer: Einführung in die Numerische Mathematik I, p. 102
 */
template <int dim>
class QWeddle : public Quadrature<dim>
{
public:
  QWeddle ();
};



/**
 * A class for Gauss quadrature with
 * logarithmic weighting function. This
 * formula is used to integrate $\ln|x|\;f(x)$ on the interval
 * $[0,1]$, where $f$ is a smooth function without
 * singularities. The collection of quadrature points and weights has
 * been obtained using <tt>Numerical Recipes</tt>.
 *
 * Notice that only the function $f(x)$ should be provided,
 * i.e., $\int_0^1 f(x) \ln|x| dx = \sum_{i=0}^N w_i f(q_i)$. Setting
 * the @p revert flag to true at construction time switches the weight
 * from $\ln|x|$ to $\ln|1-x|$.
 *
 * The weights and functions have been tabulated up to order 12.
 */
template <int dim>
class QGaussLog : public Quadrature<dim>
{
public:
  /**
   * Generate a formula with
   * <tt>n</tt> quadrature points
   */
  QGaussLog(const unsigned int n,
            const bool revert=false);

protected:
  /**
   * Sets the points of the
   * quadrature formula.
   */
  std::vector<double>
  set_quadrature_points(const unsigned int n) const;

  /**
   * Sets the weights of the
   * quadrature formula.
   */
  std::vector<double>
  set_quadrature_weights(const unsigned int n) const;

};




/**
 * A class for Gauss quadrature with arbitrary logarithmic weighting
 * function. This formula is used to to integrate
 * $\ln(|x-x_0|/\alpha)\;f(x)$ on the interval $[0,1]$,
 * where $f$ is a smooth function without singularities, and $x_0$ and
 * $\alpha$ are given at construction time, and are the location of the
 * singularity $x_0$ and an arbitrary scaling factor in the
 * singularity.
 *
 * You have to make sure that the point $x_0$ is not one of the Gauss
 * quadrature points of order $N$, otherwise an exception is thrown,
 * since the quadrature weights cannot be computed correctly.
 *
 * This quadrature formula is rather expensive, since it uses
 * internally two Gauss quadrature formulas of order n to integrate
 * the nonsingular part of the factor, and two GaussLog quadrature
 * formulas to integrate on the separate segments $[0,x_0]$ and
 * $[x_0,1]$. If the singularity is one of the extremes and the factor
 * alpha is 1, then this quadrature is the same as QGaussLog.
 *
 * The last argument from the constructor allows you to use this
 * quadrature rule in one of two possible ways:
 * \f[
 * \int_0^1 g(x) dx =
 * \int_0^1 f(x) \ln\left(\frac{|x-x_0|}{\alpha}\right) dx
 * = \sum_{i=0}^N w_i g(q_i) = \sum_{i=0}^N \bar{w}_i f(q_i)
 * \f]
 *
 * Which one of the two sets of weights is provided, can be selected
 * by the @p factor_out_singular_weight parameter. If it is false (the
 * default), then the $\bar{w}_i$ weigths are computed, and you should
 * provide only the smooth function $f(x)$, since the singularity is
 * included inside the quadrature. If the parameter is set to true,
 * then the singularity is factored out of the quadrature formula, and
 * you should provide a function $g(x)$, which should at least be
 * similar to $\ln(|x-x_0|/\alpha)$.
 *
 * Notice that this quadrature rule is worthless if you try to use it
 * for regular functions once you factored out the singularity.
 *
 * The weights and functions have been tabulated up to order 12.
 */
template <int dim>
class QGaussLogR : public Quadrature<dim>
{
public:
  /**
   * The constructor takes four arguments:
   * the order of the gauss formula on each
   * of the segments $[0,x_0]$ and
   * $[x_0,1]$, the actual location of the
   * singularity, the scale factor inside
   * the logarithmic function and a flag
   * that decides whether the singularity is
   * left inside the quadrature formula or
   * it is factored out, to be included in
   * the integrand.
   */
  QGaussLogR(const unsigned int n,
             const Point<dim> x0 = Point<dim>(),
             const double alpha = 1,
             const bool factor_out_singular_weight=false);

protected:
  /**
   * This is the length of interval
   * $(0,origin)$, or 1 if either of the two
   * extremes have been selected.
   */
  const double fraction;
};


/**
 * A class for Gauss quadrature with $1/R$ weighting function. This formula
 * can be used to to integrate $1/R \ f(x)$ on the reference
 * element $[0,1]^2$, where $f$ is a smooth function without
 * singularities, and $R$ is the distance from the point $x$ to the vertex
 * $\xi$, given at construction time by specifying its index. Notice that
 * this distance is evaluated in the reference element.
 *
 * This quadrature formula is obtained from two QGauss quadrature
 * formulas, upon transforming them into polar coordinate system
 * centered at the singularity, and then again into another reference
 * element. This allows for the singularity to be cancelled by part of
 * the Jacobian of the transformation, which contains $R$. In practice
 * the reference element is transformed into a triangle by collapsing
 * one of the sides adjacent to the singularity. The Jacobian of this
 * transformation contains $R$, which is removed before scaling the
 * original quadrature, and this process is repeated for the next half
 * element.
 *
 * Upon construction it is possible to specify whether we want the
 * singularity removed, or not. In other words, this quadrature can be
 * used to integrate $g(x) = 1/R\ f(x)$, or simply $f(x)$, with the $1/R$
 * factor already included in the quadrature weights.
 */
template <int dim>
class QGaussOneOverR : public Quadrature<dim>
{
public:
  /**
   * This constructor takes three arguments: the order of the Gauss
   * formula, the point of the reference element in which the
   * singularity is located, and whether we include the weighting
   * singular function inside the quadrature, or we leave it in the
   * user function to be integrated.
   *
   * Traditionally, quadrature formulas include their weighting
   * function, and the last argument is set to false by
   * default. There are cases, however, where this is undesirable
   * (for example when you only know that your singularity has the
   * same order of 1/R, but cannot be written exactly in this
   * way).
   *
   * In other words, you can use this function in either of
   * the following way, obtaining the same result:
   *
   * @code
   * QGaussOneOverR singular_quad(order, q_point, false);
   * // This will produce the integral of f(x)/R
   * for(unsigned int i=0; i<singular_quad.size(); ++i)
   *   integral += f(singular_quad.point(i))*singular_quad.weight(i);
   *
   * // And the same here
   * QGaussOneOverR singular_quad_noR(order, q_point, true);
   *
   * // This also will produce the integral of f(x)/R, but 1/R has to
   * // be specified.
   * for(unsigned int i=0; i<singular_quad.size(); ++i) {
   *   double R = (singular_quad_noR.point(i)-cell->vertex(vertex_id)).norm();
   *   integral += f(singular_quad_noR.point(i))*singular_quad_noR.weight(i)/R;
   * }
   * @endcode
   */
  QGaussOneOverR(const unsigned int n,
                 const Point<dim> singularity,
                 const bool factor_out_singular_weight=false);
  /**
   * The constructor takes three arguments: the order of the Gauss
   * formula, the index of the vertex where the singularity is
   * located, and whether we include the weighting singular function
   * inside the quadrature, or we leave it in the user function to
   * be integrated. Notice that this is a specialized version of the
   * previous constructor which works only for the vertices of the
   * quadrilateral.
   *
   * Traditionally, quadrature formulas include their weighting
   * function, and the last argument is set to false by
   * default. There are cases, however, where this is undesirable
   * (for example when you only know that your singularity has the
   * same order of 1/R, but cannot be written exactly in this
   * way).
   *
   * In other words, you can use this function in either of
   * the following way, obtaining the same result:
   *
   * @code
   * QGaussOneOverR singular_quad(order, vertex_id, false);
   * // This will produce the integral of f(x)/R
   * for(unsigned int i=0; i<singular_quad.size(); ++i)
   *   integral += f(singular_quad.point(i))*singular_quad.weight(i);
   *
   * // And the same here
   * QGaussOneOverR singular_quad_noR(order, vertex_id, true);
   *
   * // This also will produce the integral of f(x)/R, but 1/R has to
   * // be specified.
   * for(unsigned int i=0; i<singular_quad.size(); ++i) {
   *   double R = (singular_quad_noR.point(i)-cell->vertex(vertex_id)).norm();
   *   integral += f(singular_quad_noR.point(i))*singular_quad_noR.weight(i)/R;
   * }
   * @endcode
   */
  QGaussOneOverR(const unsigned int n,
                 const unsigned int vertex_index,
                 const bool factor_out_singular_weight=false);
private:
  /**
   * Given a quadrature point and a degree n, this function returns
   * the size of the singular quadrature rule, considering whether
   * the point is inside the cell, on an edge of the cell, or on a
   * corner of the cell.
   */
  static unsigned int quad_size(const Point<dim> singularity,
                                const unsigned int n);
};



/**
 * Sorted Quadrature. Given an arbitrary quadrature formula, this
 * class generates a quadrature formula where the quadrature points
 * are ordered according the weights, from those with smaller
 * corresponding weight, to those with higher corresponding weights.
 * This might be necessary, for example, when integrating high order
 * polynomials, since in these cases you might sum very big numbers
 * with very small numbers, and summation is not stable if the numbers
 * to sum are not close to each other.
 */
template <int dim>
class QSorted : public Quadrature<dim>
{
public:
  /**
   * The constructor takes an arbitrary quadrature formula.
   */
  QSorted (const Quadrature<dim>);

  /**
   * A rule to reorder pairs of points and weights.
   */
  bool operator()(const std::pair<double, Point<dim> > &a,
                  const std::pair<double, Point<dim> > &b);
};


/*@}*/

/* -------------- declaration of explicit specializations ------------- */

template <> QGauss<1>::QGauss (const unsigned int n);
template <> QGaussLobatto<1>::QGaussLobatto (const unsigned int n);
template <>
std::vector<long double> QGaussLobatto<1>::
compute_quadrature_points(const unsigned int, const int, const int) const;
template <>
std::vector<long double> QGaussLobatto<1>::
compute_quadrature_weights(const std::vector<long double> &, const int, const int) const;
template <>
long double QGaussLobatto<1>::
JacobiP(const long double, const int, const int, const unsigned int) const;
template <>
long double
QGaussLobatto<1>::gamma(const unsigned int n) const;

template <> std::vector<double> QGaussLog<1>::set_quadrature_points(const unsigned int) const;
template <> std::vector<double> QGaussLog<1>::set_quadrature_weights(const unsigned int) const;

template <> QMidpoint<1>::QMidpoint ();
template <> QTrapez<1>::QTrapez ();
template <> QSimpson<1>::QSimpson ();
template <> QMilne<1>::QMilne ();
template <> QWeddle<1>::QWeddle ();
template <> QGaussLog<1>::QGaussLog (const unsigned int n, const bool revert);
template <> QGaussLogR<1>::QGaussLogR (const unsigned int n, const Point<1> x0, const double alpha, const bool flag);
template <> QGaussOneOverR<2>::QGaussOneOverR (const unsigned int n, const unsigned int index, const bool flag);




DEAL_II_NAMESPACE_CLOSE

#endif
