//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__quadrature_lib_h
#define __deal2__quadrature_lib_h


#include <base/config.h>
#include <base/quadrature.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup Quadrature */
/*@{*/

/**
 * Gauss-Legendre quadrature of arbitrary order.
 *
 * The coefficients of these quadrature rules are computed by the
 * function found in <tt>Numerical Recipies</tt>.
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
 * The Gauss-Lobatto quadrature rule.
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
 * Note: This implementation has not yet been optimized concerning
 *       numerical stability and efficiency. It can be easily adapted
 *       to the general case of Gauss-Lobatto-Jacobi-Bouzitat quadrature
 *       with arbitrary parameters <i>alpha</i>, <i>beta</i>, of which
 *       the Gauss-Lobatto-Legendre quadrature (<i>alpha = beta = 0</i>)
 *       is a special case.
 *
 * @sa http://en.wikipedia.org/wiki/Handbook_of_Mathematical_Functions 
 * @sa Karniadakis, G.E. and Sherwin, S.J.:
 *     Spectral/hp element methods for computational fluid dynamics. 
 *     Oxford: Oxford University Press, 2005 
 *
 * @author Guido Kanschat, 2005, 2006; F. Prill, 2006
 */
template<int dim>
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
				      * alpha, @p beta).  @p q is
				      * number of points.
				      *
				      * @return vector containing nodes.
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
				      * @return vector containing weights.
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
    unsigned int gamma(const unsigned int n) const;
};

  

/**
 * @deprecated Use QGauss for arbitrary order Gauss formulae instead!
 *
 *  2-Point-Gauss quadrature formula, exact for polynomials of degree 3.
 *
 *  Reference: Ward Cheney, David Kincaid: "Numerical Mathematics and Computing".
 *  For a comprehensive list of Gaussian quadrature formulae, see also:
 *  A. H. Strout, D. Secrest: "Gaussian Quadrature Formulas"
 */
template <int dim>
class QGauss2 : public Quadrature<dim>
{
  public:
    QGauss2 ();
};


/**
 * @deprecated Use QGauss for arbitrary order Gauss formulae instead!
 *
 *  3-Point-Gauss quadrature formula, exact for polynomials of degree 5.
 *
 *  Reference: Ward Cheney, David Kincaid: "Numerical Mathematics and Computing".
 *  For a comprehensive list of Gaussian quadrature formulae, see also:
 *  A. H. Strout, D. Secrest: "Gaussian Quadrature Formulas"
 */
template <int dim>
class QGauss3 : public Quadrature<dim>
{
  public:
    QGauss3 ();
};


/**
 * @deprecated Use QGauss for arbitrary order Gauss formulae instead!
 *
 * 4-Point-Gauss quadrature formula, exact for polynomials of degree 7.
 *
 *  Reference: Ward Cheney, David Kincaid: "Numerical Mathematics and Computing".
 *  For a comprehensive list of Gaussian quadrature formulae, see also:
 *  A. H. Strout, D. Secrest: "Gaussian Quadrature Formulas"
 */
template <int dim>
class QGauss4 : public Quadrature<dim>
{
  public:
    QGauss4 ();
};


/**
 * @deprecated Use QGauss for arbitrary order Gauss formulae instead!
 *
 *  5-Point-Gauss quadrature formula, exact for polynomials of degree 9.
 *
 *  Reference: Ward Cheney, David Kincaid: "Numerical Mathematics and Computing".
 *  For a comprehensive list of Gaussian quadrature formulae, see also:
 *  A. H. Strout, D. Secrest: "Gaussian Quadrature Formulas"
 */
template <int dim>
class QGauss5 : public Quadrature<dim>
{
  public:
    QGauss5 ();
};


/**
 * @deprecated Use QGauss for arbitrary order Gauss formulae instead!
 *
 *  6-Point-Gauss quadrature formula, exact for polynomials of degree 11.
 *  We have not found explicit
 *  representations of the zeros of the Legendre functions of sixth
 *  and higher degree. If anyone finds them, please replace the existing
 *  numbers by these expressions.
 *
 *  Reference: J. E. Akin: "Application and Implementation of Finite
 *  Element Methods"
 *  For a comprehensive list of Gaussian quadrature formulae, see also:
 *  A. H. Strout, D. Secrest: "Gaussian Quadrature Formulas"
 */
template <int dim>
class QGauss6 : public Quadrature<dim>
{
  public:
    QGauss6 ();
};


/**
 * @deprecated Use QGauss for arbitrary order Gauss formulae instead!
 *
 *  7-Point-Gauss quadrature formula, exact for polynomials of degree 13.
 *  We have not found explicit
 *  representations of the zeros of the Legendre functions of sixth
 *  and higher degree. If anyone finds them, please replace the existing
 *  numbers by these expressions.
 *
 *  Reference: J. E. Akin: "Application and Implementation of Finite
 *  Element Methods"
 *  For a comprehensive list of Gaussian quadrature formulae, see also:
 *  A. H. Strout, D. Secrest: "Gaussian Quadrature Formulas"
 */
template <int dim>
class QGauss7 : public Quadrature<dim>
{
  public:
    QGauss7 ();
};


/**
 * Midpoint quadrature rule, exact for linear polynomials.
 */
template <int dim>
class QMidpoint : public Quadrature<dim>
{
  public:
    QMidpoint ();
};


/**
 * Simpson quadrature rule, exact for polynomials of degree 3. 
 */
template <int dim>
class QSimpson : public Quadrature<dim>
{
  public:
    QSimpson ();
};


/**
 * Trapezoidal quadrature rule, exact for linear polynomials.
 */
template <int dim>
class QTrapez : public Quadrature<dim>
{
  public:
    QTrapez ();
};

/**
 * Milne-rule. Closed Newton-Cotes formula, exact for polynomials of degree 5.
 * See Stoer: Einführung in die Numerische Mathematik I, p. 102
 */
template <int dim>
class QMilne : public Quadrature<dim>
{
  public:
    QMilne ();
};


/**
 * Weddle-rule. Closed Newton-Cotes formula, exact for polynomials of degree 7.
 * See Stoer: Einführung in die Numerische Mathematik I, p. 102
 */
template <int dim>
class QWeddle : public Quadrature<dim>
{
  public:
    QWeddle ();
};



/**
 * Gauss Quadrature Formula with logarithmic weighting function. This
 * formula is used to to integrate <tt>ln|x|*f(x)</tt> on the interval
 * <tt>[0,1]</tt>, where f is a smooth function without
 * singularities. The collection of quadrature points and weights has
 * been obtained using <tt>Numerical Recipes</tt>.
 *
 * Notice that only the function <tt>f(x)</tt> should be provided,
 * i.e., $\int_0^1 f(x) ln|x| dx = \sum_{i=0}^N w_i f(q_i)$. Setting
 * the @p revert flag to true at construction time swithces the weight
 * from <tt>ln|x|</tt> to <tt>ln|1-x|</tt>.
 *
 * The weights and functions have been tabulated up to order 12.
 *
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
 * Gauss Quadrature Formula with arbitrary logarithmic weighting
 * function. This formula is used to to integrate
 * <tt>ln(|x-x0|/alpha)*f(x)</tt> on the interval <tt>[0,1]</tt>,
 * where f is a smooth function without singularities, and x0 and
 * alpha are given at construction time, and are the location of the
 * singularity <tt>x0</tt> and an arbitrary scaling factor in the
 * singularity.
 *
 * This quadrature formula is rather expensive, since it uses
 * internally two Gauss quadrature formulas of order n to integrate
 * the nonsingular part of the factor, and two GaussLog quadrature
 * formulas to integrate on the separate segments [0,x0] and
 * [x0,1]. If the singularity is one of the extremes and the factor
 * alpha is 1, then this quadrature is the same as QGaussLog.
 *
 * Notice again that only the function <tt>f(x)</tt> should be
 * provided, i.e., $\int_0^1 f(x) ln(|x-x0|/alpha) dx = \sum_{i=0}^N
 * w_i f(q_i)$. 
 *
 * The weights and functions have been tabulated up to order 12.
 *
 */
template<int dim>
class QGaussLogR : public Quadrature<dim> {
public:
    /** The constructor takes three arguments arguments: the order of
     * the gauss formula on each of the segments [0,x0] and [x0,1],
     * the actual location of the singularity and the scale factor for
     * the logarithmic function. */
    QGaussLogR(const unsigned int n, 
	       const Point<dim> x0 = Point<dim>(), 
	       const double alpha = 1);

protected:
    /** This is the length of interval (0,origin), or 1 if either of
     * the two extremes have been selected. */
    const double fraction;
};


/**
 * Gauss Quadrature Formula with 1/R weighting function. This formula
 * is used to to integrate <tt>1/(R)*f(x)</tt> on the reference
 * element <tt>[0,1]^2</tt>, where f is a smooth function without
 * singularities, and R is the distance from the point x to the vertex
 * xi, given at construction time by specifying its index.
 *
 * This quadrature formula is obtained from two QGauss quadrature
 * formula, upon transforming them into polar coordinate system
 * centered at the singularity, and then again into another reference
 * element. This allows for the singularity to be cancelled by part of
 * the Jacobian of the transformation, which contains R. In practice
 * the reference element is transformed into a triangle by collapsing
 * one of the side adjacent to the singularity. The Jacobian of this
 * transformation contains R, which is removed before scaling the
 * original quadrature, and this process is repeated for the next
 */
template<int dim>
class QGaussOneOverR : public Quadrature<dim> {
public:
    /** The constructor takes two arguments arguments: the order of
     * the gauss formula, and the index of the vertex where the
     * singularity is located. */
    QGaussOneOverR(const unsigned int n, 
		   const unsigned int vertex_index);
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
compute_quadrature_weights(const std::vector<long double>&, const int, const int) const;
template <>
long double QGaussLobatto<1>::
JacobiP(const long double, const int, const int, const unsigned int) const;
template <>
unsigned int 
QGaussLobatto<1>::gamma(const unsigned int n) const;

template <> std::vector<double> QGaussLog<1>::set_quadrature_points(const unsigned int) const;
template <> std::vector<double> QGaussLog<1>::set_quadrature_weights(const unsigned int) const;

template <> QGauss2<1>::QGauss2 ();
template <> QGauss3<1>::QGauss3 ();
template <> QGauss4<1>::QGauss4 ();
template <> QGauss5<1>::QGauss5 ();
template <> QGauss6<1>::QGauss6 ();
template <> QGauss7<1>::QGauss7 ();
template <> QMidpoint<1>::QMidpoint ();
template <> QTrapez<1>::QTrapez ();
template <> QSimpson<1>::QSimpson ();
template <> QMilne<1>::QMilne ();
template <> QWeddle<1>::QWeddle ();
template <> QGaussLog<1>::QGaussLog (const unsigned int n, const bool revert);
template <> QGaussLogR<1>::QGaussLogR (const unsigned int n, const Point<1> x0, const double alpha);
template <> QGaussOneOverR<2>::QGaussOneOverR (const unsigned int n, const unsigned int index);




DEAL_II_NAMESPACE_CLOSE

#endif
