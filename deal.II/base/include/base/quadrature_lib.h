/*----------------------------   quadrature_lib.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __quadrature_lib_H
#define __quadrature_lib_H
/*----------------------------   quadrature_lib.h     ---------------------------*/


#include <fe/quadrature.h>


/**
 *  Second order Gauss quadrature formula.
 *
 *  Reference: Ward Cheney, David Kincaid: Numerical Mathematics and Computing.
 */
template <int dim>
class QGauss2 : public Quadrature<dim> {
  public:
    QGauss2 ();
};



/**
 *  Third order Gauss quadrature formula.
 *
 *  Reference: Ward Cheney, David Kincaid: Numerical Mathematics and Computing.
 */
template <int dim>
class QGauss3 : public Quadrature<dim> {
  public:
    QGauss3 ();
};



/**
 * Fourth order Gauss quadrature formula.
 *
 *  Reference: Ward Cheney, David Kincaid: Numerical Mathematics and Computing.
 */
template <int dim>
class QGauss4 : public Quadrature<dim> {
  public:
    QGauss4 ();
};




/**
 *  Fifth order Gauss quadrature formula.
 *
 *  Reference: Ward Cheney, David Kincaid: Numerical Mathematics and Computing.
 */
template <int dim>
class QGauss5 : public Quadrature<dim> {
  public:
    QGauss5 ();
};



/**
 *  Sixth order Gauss quadrature formula. I have not found explicite
 *  representations of the zeros of the Legendre functions of sixth
 *  and higher degree. If anyone finds them, please replace the existing
 *  numbers by these expressions.
 *
 *  Reference: J. E. Akin: Application and Implementation of Finite
 *  Element Methods
 */
template <int dim>
class QGauss6 : public Quadrature<dim> {
  public:
    QGauss6 ();
};



/**
 *  Seventh order Gauss quadrature formula. I have not found explicite
 *  representations of the zeros of the Legendre functions of sixth
 *  and higher degree. If anyone finds them, please replace the existing
 *  numbers by these expressions.
 *
 *  Reference: J. E. Akin: Application and Implementation of Finite
 *  Element Methods
 */
template <int dim>
class QGauss7 : public Quadrature<dim> {
  public:
    QGauss7 ();
};



/**
 *  Eighth order Gauss quadrature formula. I have not found explicite
 *  representations of the zeros of the Legendre functions of sixth
 *  and higher degree. If anyone finds them, please replace the existing
 *  numbers by these expressions.
 *
 *  Reference: J. E. Akin: Application and Implementation of Finite
 *  Element Methods
 */
template <int dim>
class QGauss8 : public Quadrature<dim> {
  public:
    QGauss8 ();
};





/**
 * First order midpoint quadrature rule.
 */
template <int dim>
class QMidpoint : public Quadrature<dim> {
  public:
    QMidpoint ();
};



/**
 * Simpson quadrature rule.
 */
template <int dim>
class QSimpson : public Quadrature<dim> {
  public:
    QSimpson ();
};



/**
 * Trapezoidal quadrature rule.
 */
template <int dim>
class QTrapez : public Quadrature<dim> {
  public:
    QTrapez ();
};




/*----------------------------   quadrature_lib.h     ---------------------------*/
/* end of #ifndef __quadrature_lib_H */
#endif
/*----------------------------   quadrature_lib.h     ---------------------------*/
