/*----------------------------   quadrature_lib.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __quadrature_lib_H
#define __quadrature_lib_H
/*----------------------------   quadrature_lib.h     ---------------------------*/


#include <base/quadrature.h>


/**
 *  2-Point-Gauss quadrature formula.
 *
 *  Reference: Ward Cheney, David Kincaid: Numerical Mathematics and Computing.
 */
template <int dim>
class QGauss2 : public Quadrature<dim>
{
  public:
    QGauss2 ();
};



/**
 *  3-Point-Gauss quadrature formula.
 *
 *  Reference: Ward Cheney, David Kincaid: Numerical Mathematics and Computing.
 */
template <int dim>
class QGauss3 : public Quadrature<dim>
{
  public:
    QGauss3 ();
};



/**
 * 4-Point-Gauss quadrature formula.
 *
 *  Reference: Ward Cheney, David Kincaid: Numerical Mathematics and Computing.
 */
template <int dim>
class QGauss4 : public Quadrature<dim>
{
  public:
    QGauss4 ();
};




/**
 *  5-Point-Gauss quadrature formula.
 *
 *  Reference: Ward Cheney, David Kincaid: Numerical Mathematics and Computing.
 */
template <int dim>
class QGauss5 : public Quadrature<dim>
{
  public:
    QGauss5 ();
};



/**
 *  6-Point-Gauss quadrature formula. I have not found explicite
 *  representations of the zeros of the Legendre functions of sixth
 *  and higher degree. If anyone finds them, please replace the existing
 *  numbers by these expressions.
 *
 *  Reference: J. E. Akin: Application and Implementation of Finite
 *  Element Methods
 */
template <int dim>
class QGauss6 : public Quadrature<dim>
{
  public:
    QGauss6 ();
};



/**
 *  7-Point-Gauss quadrature formula. I have not found explicite
 *  representations of the zeros of the Legendre functions of sixth
 *  and higher degree. If anyone finds them, please replace the existing
 *  numbers by these expressions.
 *
 *  Reference: J. E. Akin: Application and Implementation of Finite
 *  Element Methods
 */
template <int dim>
class QGauss7 : public Quadrature<dim>
{
  public:
    QGauss7 ();
};



/**
 *  8-Point-Gauss quadrature formula. I have not found explicite
 *  representations of the zeros of the Legendre functions of sixth
 *  and higher degree. If anyone finds them, please replace the existing
 *  numbers by these expressions.
 *
 *  Reference: J. E. Akin: Application and Implementation of Finite
 *  Element Methods
 */
template <int dim>
class QGauss8 : public Quadrature<dim>
{
  public:
    QGauss8 ();
};





/**
 * First order midpoint quadrature rule.
 * For compatibility, this rule may be accessed as #QGauss1#, too.
 */
template <int dim>
class QMidpoint : public Quadrature<dim>
{
  public:
    QMidpoint ();
};

#define QGauss1 QMidpoint

/**
 * Simpson quadrature rule.
 */
template <int dim>
class QSimpson : public Quadrature<dim>
{
  public:
    QSimpson ();
};



/**
 * Trapezoidal quadrature rule.
 */
template <int dim>
class QTrapez : public Quadrature<dim>
{
  public:
    QTrapez ();
};


/**
 * Iterated trapezoidal rule. The aim of this class is to provide a
 * low order formula, where the error constant can be tuned by
 * increasing the number of quadrature points. This is useful in
 * integrating non-differentiable functions on cells.
 *
 * For internal use, it may be worth to know that the points are
 * ordered in a fashion such that the last coordinate is the one which
 * runs fastest and then lexicographically from back to front.
 */
template <int dim>
class QIteratedTrapez :
  public Quadrature<dim>
{
public:
  QIteratedTrapez(const unsigned intervals);
};

/**
 * Iterated Simpson rule.
 * Like #QIteratedTrapez#, this class provides a lower order formula,
 * while the error constant can be tuned by choosing the number of sub-cells.
 */
template <int dim>
class QIteratedSimpson :
  public Quadrature<dim>
{
public:
  QIteratedSimpson(const unsigned intervals);
};

/*----------------------------   quadrature_lib.h     ---------------------------*/
/* end of #ifndef __quadrature_lib_H */
#endif
/*----------------------------   quadrature_lib.h     ---------------------------*/
