//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
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
 * Gauss quadrature rule using the Gauss-Lobatto points.
 *
 * @note The quadrature weights are not implemented yet.
 *
 * @author Guido Kanschat, 2005
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

/*@}*/

/* -------------- declaration of explicit specializations ------------- */

template <> QGauss<1>::QGauss (const unsigned int n);
template <> QGaussLobatto<1>::QGaussLobatto (const unsigned int n);
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


#endif
