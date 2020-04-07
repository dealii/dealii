// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2019 by the deal.II authors
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

#ifndef dealii_quadrature_lib_h
#define dealii_quadrature_lib_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup Quadrature */
/*@{*/

/**
 * The Gauss-Legendre family of quadrature rules for numerical integration.
 *
 * The coefficients of these quadrature rules are computed by the function
 * described in <a
 * href="http://en.wikipedia.org/wiki/Numerical_Recipes">Numerical
 * Recipes</a>.
 *
 * @author Guido Kanschat, 2001
 */
template <int dim>
class QGauss : public Quadrature<dim>
{
public:
  /**
   * Generate a formula with <tt>n</tt> quadrature points (in each space
   * direction), exact for polynomials of degree <tt>2n-1</tt>.
   */
  QGauss(const unsigned int n);
};


/**
 * The Gauss-Lobatto family of quadrature rules for numerical integration.
 *
 * This modification of the Gauss quadrature uses the two interval end points
 * as well. Being exact for polynomials of degree <i>2n-3</i>, this formula is
 * suboptimal by two degrees.
 *
 * The quadrature points are interval end points plus the roots of the
 * derivative of the Legendre polynomial <i>P<sub>n-1</sub></i> of degree
 * <i>n-1</i>. The quadrature weights are
 * <i>2/(n(n-1)(P<sub>n-1</sub>(x<sub>i</sub>)<sup>2</sup>)</i>.
 *
 * @note This implementation has not been optimized concerning numerical
 * stability and efficiency. It can be easily adapted to the general case of
 * Gauss-Lobatto-Jacobi-Bouzitat quadrature with arbitrary parameters
 * $\alpha$, $\beta$, of which the Gauss-Lobatto-Legendre quadrature ($\alpha
 * = \beta = 0$) is a special case.
 *
 * @sa http://en.wikipedia.org/wiki/Handbook_of_Mathematical_Functions @sa
 * Karniadakis, G.E. and Sherwin, S.J.: Spectral/hp element methods for
 * computational fluid dynamics. Oxford: Oxford University Press, 2005
 *
 * @author Guido Kanschat, 2005, 2006; F. Prill, 2006
 */
template <int dim>
class QGaussLobatto : public Quadrature<dim>
{
public:
  /**
   * Generate a formula with <tt>n</tt> quadrature points (in each space
   * direction).
   */
  QGaussLobatto(const unsigned int n);
};



/**
 * The midpoint rule for numerical quadrature. This one-point formula is exact
 * for linear polynomials.
 */
template <int dim>
class QMidpoint : public Quadrature<dim>
{
public:
  QMidpoint();
};


/**
 * The Simpson rule for numerical quadrature. This formula with 3 quadrature
 * points is exact for polynomials of degree 3.
 */
template <int dim>
class QSimpson : public Quadrature<dim>
{
public:
  QSimpson();
};



/**
 * The trapezoidal rule for numerical quadrature. This formula with two
 * quadrature points is exact for linear polynomials.
 *
 * The class is poorly named since the proper name of the quadrature formula
 * is "trapezoidal rule", or sometimes also called the "trapezoid rule". The
 * misnomer results from the fact that its original authors' poor English
 * language skills led them to translate the name incorrectly from the German
 * "Trapezregel".
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class QTrapez : public Quadrature<dim>
{
public:
  QTrapez();
};



/**
 * The Milne rule for numerical quadrature formula. The Milne rule is a closed
 * Newton-Cotes formula and is exact for polynomials of degree 5.
 *
 * @sa Stoer: Einführung in die Numerische Mathematik I, p. 102
 */
template <int dim>
class QMilne : public Quadrature<dim>
{
public:
  QMilne();
};


/**
 * The Weddle rule for numerical quadrature. The Weddle rule is a closed
 * Newton-Cotes formula and is exact for polynomials of degree 7.
 *
 * @sa Stoer: Einführung in die Numerische Mathematik I, p. 102
 */
template <int dim>
class QWeddle : public Quadrature<dim>
{
public:
  QWeddle();
};



/**
 * A class for Gauss quadrature with logarithmic weighting function. This
 * formula is used to integrate $\ln|x|\;f(x)$ on the interval $[0,1]$, where
 * $f$ is a smooth function without singularities. The collection of
 * quadrature points and weights has been obtained using <tt>Numerical
 * Recipes</tt>.
 *
 * Notice that only the function $f(x)$ should be provided, i.e., $\int_0^1
 * f(x) \ln|x| dx = \sum_{i=0}^N w_i f(q_i)$. Setting the @p revert flag to
 * true at construction time switches the weight from $\ln|x|$ to $\ln|1-x|$.
 *
 * The weights and functions have been tabulated up to order 12.
 */
template <int dim>
class QGaussLog : public Quadrature<dim>
{
public:
  /**
   * Generate a formula with <tt>n</tt> quadrature points
   */
  QGaussLog(const unsigned int n, const bool revert = false);

private:
  /**
   * Compute the points of the quadrature formula.
   */
  static std::vector<double>
  get_quadrature_points(const unsigned int n);

  /**
   * Compute the weights of the quadrature formula.
   */
  static std::vector<double>
  get_quadrature_weights(const unsigned int n);
};



/**
 * A class for Gauss quadrature with arbitrary logarithmic weighting function.
 * This formula is used to integrate $\ln(|x-x_0|/\alpha)\;f(x)$ on the
 * interval $[0,1]$, where $f$ is a smooth function without singularities, and
 * $x_0$ and $\alpha$ are given at construction time, and are the location of
 * the singularity $x_0$ and an arbitrary scaling factor in the singularity.
 *
 * You have to make sure that the point $x_0$ is not one of the Gauss
 * quadrature points of order $N$, otherwise an exception is thrown, since the
 * quadrature weights cannot be computed correctly.
 *
 * This quadrature formula is rather expensive, since it uses internally two
 * Gauss quadrature formulas of order n to integrate the nonsingular part of
 * the factor, and two GaussLog quadrature formulas to integrate on the
 * separate segments $[0,x_0]$ and $[x_0,1]$. If the singularity is one of the
 * extremes and the factor alpha is 1, then this quadrature is the same as
 * QGaussLog.
 *
 * The last argument from the constructor allows you to use this quadrature
 * rule in one of two possible ways: \f[ \int_0^1 g(x) dx = \int_0^1 f(x)
 * \ln\left(\frac{|x-x_0|}{\alpha}\right) dx = \sum_{i=0}^N w_i g(q_i) =
 * \sum_{i=0}^N \bar{w}_i f(q_i) \f]
 *
 * Which one of the two sets of weights is provided, can be selected by the @p
 * factor_out_singular_weight parameter. If it is false (the default), then
 * the $\bar{w}_i$ weights are computed, and you should provide only the
 * smooth function $f(x)$, since the singularity is included inside the
 * quadrature. If the parameter is set to true, then the singularity is
 * factored out of the quadrature formula, and you should provide a function
 * $g(x)$, which should at least be similar to $\ln(|x-x_0|/\alpha)$.
 *
 * Notice that this quadrature rule is worthless if you try to use it for
 * regular functions once you factored out the singularity.
 *
 * The weights and functions have been tabulated up to order 12.
 */
template <int dim>
class QGaussLogR : public Quadrature<dim>
{
public:
  /**
   * The constructor takes four arguments: the order of the Gauss formula on
   * each of the segments $[0,x_0]$ and $[x_0,1]$, the actual location of the
   * singularity, the scale factor inside the logarithmic function and a flag
   * that decides whether the singularity is left inside the quadrature
   * formula or it is factored out, to be included in the integrand.
   */
  QGaussLogR(const unsigned int n,
             const Point<dim>   x0                         = Point<dim>(),
             const double       alpha                      = 1,
             const bool         factor_out_singular_weight = false);

  /**
   * Move constructor. We cannot rely on the move constructor for `Quadrature`,
   * since it does not know about the additional member `fraction` of this
   * class.
   */
  QGaussLogR(QGaussLogR<dim> &&) noexcept = default;

protected:
  /**
   * This is the length of interval $(0,origin)$, or 1 if either of the two
   * extremes have been selected.
   */
  const double fraction;
};


/**
 * A class for Gauss quadrature with $1/R$ weighting function. This formula
 * can be used to integrate $1/R \ f(x)$ on the reference element $[0,1]^2$,
 * where $f$ is a smooth function without singularities, and $R$ is the
 * distance from the point $x$ to the vertex $\xi$, given at construction time
 * by specifying its index. Notice that this distance is evaluated in the
 * reference element.
 *
 * This quadrature formula is obtained from two QGauss quadrature formulas,
 * upon transforming them into polar coordinate system centered at the
 * singularity, and then again into another reference element. This allows for
 * the singularity to be cancelled by part of the Jacobian of the
 * transformation, which contains $R$. In practice the reference element is
 * transformed into a triangle by collapsing one of the sides adjacent to the
 * singularity. The Jacobian of this transformation contains $R$, which is
 * removed before scaling the original quadrature, and this process is
 * repeated for the next half element.
 *
 * Upon construction it is possible to specify whether we want the singularity
 * removed, or not. In other words, this quadrature can be used to integrate
 * $g(x) = 1/R\ f(x)$, or simply $f(x)$, with the $1/R$ factor already
 * included in the quadrature weights.
 */
template <int dim>
class QGaussOneOverR : public Quadrature<dim>
{
public:
  /**
   * This constructor takes three arguments: the order of the Gauss formula,
   * the point of the reference element in which the singularity is located,
   * and whether we include the weighting singular function inside the
   * quadrature, or we leave it in the user function to be integrated.
   *
   * Traditionally, quadrature formulas include their weighting function, and
   * the last argument is set to false by default. There are cases, however,
   * where this is undesirable (for example when you only know that your
   * singularity has the same order of 1/R, but cannot be written exactly in
   * this way).
   *
   * In other words, you can use this function in either of the following way,
   * obtaining the same result:
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
                 const Point<dim>   singularity,
                 const bool         factor_out_singular_weight = false);
  /**
   * The constructor takes three arguments: the order of the Gauss formula,
   * the index of the vertex where the singularity is located, and whether we
   * include the weighting singular function inside the quadrature, or we
   * leave it in the user function to be integrated. Notice that this is a
   * specialized version of the previous constructor which works only for the
   * vertices of the quadrilateral.
   *
   * Traditionally, quadrature formulas include their weighting function, and
   * the last argument is set to false by default. There are cases, however,
   * where this is undesirable (for example when you only know that your
   * singularity has the same order of 1/R, but cannot be written exactly in
   * this way).
   *
   * In other words, you can use this function in either of the following way,
   * obtaining the same result:
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
                 const bool         factor_out_singular_weight = false);

private:
  /**
   * Given a quadrature point and a degree n, this function returns the size
   * of the singular quadrature rule, considering whether the point is inside
   * the cell, on an edge of the cell, or on a corner of the cell.
   */
  static unsigned int
  quad_size(const Point<dim> singularity, const unsigned int n);
};



/**
 * Sorted Quadrature. Given an arbitrary quadrature formula, this class
 * generates a quadrature formula where the quadrature points are ordered
 * according the weights, from those with smaller corresponding weight, to
 * those with higher corresponding weights. This might be necessary, for
 * example, when integrating high order polynomials, since in these cases you
 * might sum very big numbers with very small numbers, and summation is not
 * stable if the numbers to sum are not close to each other.
 */
template <int dim>
class QSorted : public Quadrature<dim>
{
public:
  /**
   * The constructor takes an arbitrary quadrature formula @p quad and sorts
   * its points and weights according to ascending weights.
   */
  QSorted(const Quadrature<dim> &quad);

private:
  /**
   * A rule for std::sort to reorder pairs of points and weights.
   * @p a and @p b are indices into the weights array and the result will
   * be determined by comparing the weights.
   */
  bool
  compare_weights(const unsigned int a, const unsigned int b) const;
};

/**
 * Telles quadrature of arbitrary order.
 *
 * The coefficients of these quadrature rules are computed using a non linear
 * change of variables starting from a Gauss-Legendre quadrature formula. This
 * is done using a cubic polynomial, $n = a x^3 + b x^2 + c x + d$ in order to
 * integrate a singular integral, with singularity at a given point x_0.
 *
 * We start from a Gauss Quadrature Formula with arbitrary function. Then we
 * apply the cubic variable change. In the paper, J.C.F.Telles:A Self-Adaptive
 * Co-ordinate Transformation For Efficient Numerical Evaluation of General
 * Boundary Element Integrals. International Journal for Numerical Methods in
 * Engineering, vol 24, pages 959–973. year 1987, the author applies the
 * transformation on the reference cell $[-1, 1]$ getting
 * @f{align*}{
 * n(1) &= 1, \\ n(-1) &= -1, \\ \frac{dn}{dx} &= 0 \text{ at }
 * x = x_0, \\ \frac{d^2n}{dx^2} &= 0 \text{ at  } x = x_0
 * @f}
 * We get
 * @f{align*}{
 * a &= \frac{1}{q}, \\
 * b &= -3 \frac{\bar{\Gamma}}{q}, \\
 * c &= 3 \frac{\bar{\Gamma}}{q}, \\
 * d &= -b,
 * @f}
 * with
 * @f{align*}{
 * \eta^{*} &= \bar{\eta}^2 - 1, \\
 * \bar{\Gamma}  &= \sqrt[3]{\bar{\eta} \eta^{*} + |\eta^{*} | }
 *                  + \sqrt[3]{ \bar{\eta} \eta^{*} - |\eta^{*} | }
 *                  + \bar{\eta}, \\
 * q &= (\Gamma-\bar{\Gamma})^3 + \bar{\Gamma}
 *      \frac{\bar{\Gamma}^2+3}{1+3\bar{\Gamma}^2}
 * @f}
 * Since the library assumes $[0,1]$ as reference interval, we will map these
 * values on the proper reference interval in the implementation.
 *
 * This variable change can be used to integrate singular integrals. One
 * example is $f(x)/|x-x_0|$ on the reference interval $[0,1]$, where $x_0$ is
 * given at construction time, and is the location of the singularity $x_0$,
 * and $f(x)$ is a smooth non singular function.
 *
 * Singular quadrature formula are rather expensive, nevertheless Telles'
 * quadrature formula are much easier to compute with respect to other
 * singular integration techniques as Lachat-Watson.
 *
 * We have implemented the case for $dim = 1$. When we deal the case $dim >1$
 * we have computed the quadrature formula has a tensorial product of one
 * dimensional Telles' quadrature formulas considering the different
 * components of the singularity.
 *
 * The weights and functions for Gauss Legendre formula have been tabulated up
 * to order 12.
 *
 * @author Nicola Giuliani, Luca Heltai 2015
 */
template <int dim>
class QTelles : public Quadrature<dim>
{
public:
  /**
   * A constructor that takes a quadrature formula and a singular point as
   * argument. The quadrature formula will be mapped using Telles' rule. Make
   * sure that the order of the quadrature rule is appropriate for the
   * singularity in question.
   */
  QTelles(const Quadrature<1> &base_quad, const Point<dim> &singularity);
  /**
   * A variant of above constructor that takes as parameters the order @p n
   * and location of a singularity. A Gauss Legendre quadrature of order n
   * will be used
   */
  QTelles(const unsigned int n, const Point<dim> &singularity);
};

/**
 * Gauss-Chebyshev quadrature rules integrate the weighted product
 * $\int_{-1}^1 f(x) w(x) dx$ with weight given by: $w(x) = 1/\sqrt{1-x^2}$.
 * The nodes and weights are known analytically, and are exact for monomials
 * up to the order $2n-1$, where $n$ is the number of quadrature points. Here
 * we rescale the quadrature formula so that it is defined on the interval
 * $[0,1]$ instead of $[-1,1]$. So the quadrature formulas integrate exactly
 * the integral $\int_0^1 f(x) w(x) dx$ with the weight: $w(x) =
 * 1/\sqrt{x(1-x)}$. For details see: M. Abramowitz & I.A. Stegun: Handbook of
 * Mathematical Functions, par. 25.4.38
 *
 * @author Giuseppe Pitton, Luca Heltai 2015
 */
template <int dim>
class QGaussChebyshev : public Quadrature<dim>
{
public:
  /// Generate a formula with <tt>n</tt> quadrature points
  QGaussChebyshev(const unsigned int n);
};


/**
 * Gauss-Radau-Chebyshev quadrature rules integrate the weighted product
 * $\int_{-1}^1 f(x) w(x) dx$ with weight given by: $w(x) = 1/\sqrt{1-x^2}$
 * with the additional constraint that a quadrature point lies at one of the
 * two extrema of the interval. The nodes and weights are known analytically,
 * and are exact for monomials up to the order $2n-2$, where $n$ is the number
 * of quadrature points. Here we rescale the quadrature formula so that it is
 * defined on the interval $[0,1]$ instead of $[-1,1]$. So the quadrature
 * formulas integrate exactly the integral $\int_0^1 f(x) w(x) dx$ with the
 * weight: $w(x) = 1/\sqrt{x(1-x)}$. By default the quadrature is constructed
 * with the left endpoint as quadrature node, but the quadrature node can be
 * imposed at the right endpoint through the variable ep that can assume the
 * values left or right.
 *
 * @author Giuseppe Pitton, Luca Heltai 2015
 */
template <int dim>
class QGaussRadauChebyshev : public Quadrature<dim>
{
public:
  /* EndPoint is used to specify which of the two endpoints of the unit interval
   * is used also as quadrature point
   */
  enum EndPoint
  {
    /**
     * Left end point.
     */
    left,
    /**
     * Right end point.
     */
    right
  };
  /// Generate a formula with <tt>n</tt> quadrature points
  QGaussRadauChebyshev(const unsigned int n,
                       EndPoint           ep = QGaussRadauChebyshev::left);

  /**
   * Move constructor. We cannot rely on the move constructor for `Quadrature`,
   * since it does not know about the additional member `ep` of this class.
   */
  QGaussRadauChebyshev(QGaussRadauChebyshev<dim> &&) noexcept = default;

private:
  const EndPoint ep;
};

/**
 * Gauss-Lobatto-Chebyshev quadrature rules integrate the weighted product
 * $\int_{-1}^1 f(x) w(x) dx$ with weight given by: $w(x) = 1/\sqrt{1-x^2}$,
 * with the additional constraint that two of the quadrature points are
 * located at the endpoints of the quadrature interval. The nodes and weights
 * are known analytically, and are exact for monomials up to the order $2n-3$,
 * where $n$ is the number of quadrature points. Here we rescale the
 * quadrature formula so that it is defined on the interval $[0,1]$ instead of
 * $[-1,1]$. So the quadrature formulas integrate exactly the integral
 * $\int_0^1 f(x) w(x) dx$ with the weight: $w(x) = 1/\sqrt{x(1-x)}$. For
 * details see: M. Abramowitz & I.A. Stegun: Handbook of Mathematical
 * Functions, par. 25.4.40
 *
 * @author Giuseppe Pitton, Luca Heltai 2015
 */
template <int dim>
class QGaussLobattoChebyshev : public Quadrature<dim>
{
public:
  /// Generate a formula with <tt>n</tt> quadrature points
  QGaussLobattoChebyshev(const unsigned int n);
};

/**
 * Given an arbitrary quadrature formula, return one that chops the quadrature
 * points above the hyper-plane defined by $\sum_i x_i = 1$. In other words,
 * it extracts those quadrature points from the base formula that satisfy
 * $\sum_i (\mathbf x_q)_i \le 1+10^{-12}$."
 *
 * In general the resulting quadrature is not very useful, unless the
 * quadrature you started from has been constructed specifically to integrate
 * over triangles or tetrahedra. This class only ensures that the resulting
 * quadrature formula only has quadrature points in the reference simplex or on
 * its boundary.
 *
 * No transformation is applied to the weights, and the weights referring to
 * points that live outside the reference simplex are simply discarded.
 *
 * The main use of this quadrature formula is not to chop tensor product
 * quadratures. Ideally you should pass to this class a quadrature formula
 * constructed directly using points and weights in the reference simplex,
 * capable of integrating on triangles or tetrahedra.
 *
 * For finite elements based on quadrilaterals and hexahedra, a QSimplex
 * quadrature formula is not very useful on its own. This class is typically
 * used in conjunction with other classes, like QSplit, to patch the reference
 * element using several QSimplex quadrature formulas.
 *
 * Such quadrature formulas are useful to integrate functions with
 * singularities at certain points, or functions that present jumps along a
 * co-dimension one surface inside the reference element, like in the extended
 * finite element method (XFEM).
 *
 * @author Luca Heltai, 2017.
 */
template <int dim>
class QSimplex : public Quadrature<dim>
{
public:
  /**
   * Construct a quadrature that only contains the points that are in the lower
   * left reference simplex.
   *
   * @param[in] quad The input quadrature.
   */
  QSimplex(const Quadrature<dim> &quad);

  /**
   * Return an affine transformation of this quadrature, that can be used to
   * integrate on the simplex identified by `vertices`.
   *
   * Both the quadrature point locations and the weights are transformed, so
   * that you can effectively use the resulting quadrature to integrate on the
   * simplex.
   *
   * The transformation is defined as
   * \f[
   * x = v_0 + B \hat x
   * \f]
   * where the matrix $B$ is given by $B_{ij} = v[j][i]-v[0][i]$.
   *
   * The weights are scaled with the absolute value of the determinant of $B$,
   * that is $J \dealcoloneq |\text{det}(B)|$. If $J$ is zero, an empty
   * quadrature is returned. This may happen, in two dimensions, if the three
   * vertices are aligned, or in three dimensions if the four vertices are on
   * the same plane.
   *
   * @param[in] vertices The vertices of the simplex you wish to integrate on
   * @return A quadrature object that can be used to integrate on the simplex
   */
  Quadrature<dim>
  compute_affine_transformation(
    const std::array<Point<dim>, dim + 1> &vertices) const;
};

/**
 * A quadrature that implements a polar transformation from a square to a
 * triangle to integrate singularities in the origin of the reference simplex.
 * The quadrature is obtained through the following polar transformation:
 *
 * \f[
 *  \begin{pmatrix}
 *  x \\
 *  y
 *  \end{pmatrix}
 *  =
 * \begin{pmatrix}
 *  \frac{\hat x}{\sin(\theta)+\cos(\theta)} cos(\theta) \\
 *  \frac{\hat x}{\sin(\theta)+\cos(\theta)} sin(\theta)
 *  \end{pmatrix}
 *  \qquad \theta \dealcoloneq \frac\pi 2 \hat y
 * \f]
 *
 * @author Luca Heltai, 2017
 */
class QTrianglePolar : public QSimplex<2>
{
public:
  /**
   * Construct a QTrianglePolar quadrature, with different formulas in the
   * radial and angular directions.
   *
   * @param radial_quadrature Radial quadrature
   * @param angular_quadrature Angular quadrature
   */
  QTrianglePolar(const Quadrature<1> &radial_quadrature,
                 const Quadrature<1> &angular_quadrature);

  /**
   * Call the other constructor, with QGauss<1>(n) for both radial and
   * angular quadrature.
   *
   * @param n Order of QGauss quadrature
   */
  QTrianglePolar(const unsigned int n);
};

/**
 * A quadrature that implements the Duffy transformation from a square to a
 * triangle to integrate singularities in the origin of the reference
 * simplex.
 *
 * The Duffy transformation is defined as
 * \f[
 * \begin{pmatrix}
 * x\\
 * y
 * \end{pmatrix}
 * =
 * \begin{pmatrix}
 * \hat x^\beta (1-\hat y)\\
 * \hat x^\beta \hat y
 * \end{pmatrix}
 * \f]
 *
 * with determinant of the Jacobian equal to $J= \beta \hat x^{2\beta-1}$.
 * Such transformation maps the reference square $[0,1]\times[0,1]$ to the
 * reference simplex, by collapsing the left side of the square and squeezing
 * quadrature points towards the origin, and then shearing the resulting
 * triangle to the reference one. This transformation shows good convergence
 * properties when $\beta = 1$ with singularities of order $1/R$ in the origin,
 * but different $\beta$ values can be selected to increase convergence and/or
 * accuracy when higher order Gauss rules are used (see "Generalized Duffy
 * transformation for integrating vertex singularities", S. E. Mousavi, N.
 * Sukumar, Computational Mechanics 2009).
 *
 * When $\beta = 1$, this transformation is also known as the Lachat-Watson
 * transformation.
 *
 * @author Luca Heltai, Nicola Giuliani, 2017.
 */
class QDuffy : public QSimplex<2>
{
public:
  /**
   * Constructor that allows the specification of different quadrature rules
   * along the "radial" and "angular" directions.
   *
   * Since this quadrature is not based on a Polar change of coordinates, it
   * is not fully proper to talk about radial and angular directions. However,
   * the effect of the Duffy transformation is similar to a polar change
   * of coordinates, since the resulting quadrature points are aligned radially
   * with respect to the singularity.
   *
   * @param radial_quadrature Base quadrature to use in the radial direction
   * @param angular_quadrature Base quadrature to use in the angular direction
   * @param beta Exponent used in the transformation
   */
  QDuffy(const Quadrature<1> &radial_quadrature,
         const Quadrature<1> &angular_quadrature,
         const double         beta = 1.0);

  /**
   * Call the above constructor with QGauss<1>(n) quadrature formulas for
   * both the radial and angular quadratures.
   *
   * @param n Order of QGauss quadrature
   * @param beta Exponent used in the transformation
   */
  QDuffy(const unsigned int n, const double beta);
};

/**
 * A quadrature to use when the cell should be split into subregions to
 * integrate using one or more base quadratures.
 *
 * @author Luca Heltai, 2017.
 */
template <int dim>
class QSplit : public Quadrature<dim>
{
public:
  /**
   * Construct a quadrature formula by splitting the reference hyper cube into
   * the minimum number of simplices that have vertex zero coinciding with
   * @p split_point, and patch together affine transformations of the @p base
   * quadrature. The point @p split_point should be in the reference element,
   * and an exception is thrown if this is not the case.
   *
   * In two dimensions, the resulting quadrature formula will be composed of
   * two, three, or four triangular quadrature formulas if @p split_point
   * coincides with one of the vertices, if it lies on one of the edges, or if
   * it is internal to the reference element respectively.
   *
   * The same is true for the three dimensional case, with six, eight, ten, or
   * twelve tetrahedral quadrature formulas if @p split_point coincides with one
   * of the vertices, if it lies on one of the edges, on one of the faces, or
   * if it is internal to the reference element respectively.
   *
   * The resulting quadrature can be used, for example, to integrate functions
   * with integrable singularities at the split point, provided that you select
   * as base quadrature one that can integrate singular points on vertex zero
   * of the reference simplex.
   *
   * An example usage in dimension two is given by:
   * @code
   * const unsigned int order = 5;
   * QSplit<2> quad(QTrianglePolar(order), Point<2>(.3,.4));
   * @endcode
   *
   * The resulting quadrature will look like the following:
   * @image html split_quadrature.png ""
   *
   * @param base Base QSimplex quadrature to use
   * @param split_point Where to split the hyper cube
   */
  QSplit(const QSimplex<dim> &base, const Point<dim> &split_point);
};

/*@}*/

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN
template <>
QGauss<1>::QGauss(const unsigned int n);
template <>
QGaussLobatto<1>::QGaussLobatto(const unsigned int n);

template <>
std::vector<double>
QGaussLog<1>::get_quadrature_points(const unsigned int);
template <>
std::vector<double>
QGaussLog<1>::get_quadrature_weights(const unsigned int);

template <>
QMidpoint<1>::QMidpoint();
template <>
QTrapez<1>::QTrapez();
template <>
QSimpson<1>::QSimpson();
template <>
QMilne<1>::QMilne();
template <>
QWeddle<1>::QWeddle();
template <>
QGaussLog<1>::QGaussLog(const unsigned int n, const bool revert);
template <>
QGaussLogR<1>::QGaussLogR(const unsigned int n,
                          const Point<1>     x0,
                          const double       alpha,
                          const bool         flag);
template <>
QGaussOneOverR<2>::QGaussOneOverR(const unsigned int n,
                                  const unsigned int index,
                                  const bool         flag);
template <>
QTelles<1>::QTelles(const Quadrature<1> &base_quad,
                    const Point<1> &     singularity);
#endif // DOXYGEN



DEAL_II_NAMESPACE_CLOSE
#endif
