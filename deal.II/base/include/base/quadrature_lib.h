/*----------------------------   quadrature_lib.h     ---------------------------*/
/*      <Id:>                 */
#ifndef __quadrature_lib_H
#define __quadrature_lib_H
/*----------------------------   quadrature_lib.h     ---------------------------*/


#include <fe/quadrature.h>


template <int dim>
class QGauss2 : public Quadrature<dim> {
  public:
    QGauss2 ();
};



template <int dim>
class QGauss2x4 : public Quadrature<dim> {
  public:
    QGauss2x4 ();
};



template <int dim>
class QGauss4 : public Quadrature<dim> {
  public:
    QGauss4 ();
};



template <int dim>
class QGauss8 : public Quadrature<dim> {
  public:
    QGauss8 ();
};



template <int dim>
class QMidpoint : public Quadrature<dim> {
  public:
    QMidpoint ();
};



template <int dim>
class QSimpson : public Quadrature<dim> {
  public:
    QSimpson ();
};



template <int dim>
class QTrapez : public Quadrature<dim> {
  public:
    QTrapez ();
};




/*----------------------------   quadrature_lib.h     ---------------------------*/
/* end of #ifndef __quadrature_lib_H */
#endif
/*----------------------------   quadrature_lib.h     ---------------------------*/
