/*----------------------------   quadrature.h     ---------------------------*/
/*      <Id:>                 */
#ifndef __quadrature_H
#define __quadrature_H
/*----------------------------   quadrature.h     ---------------------------*/


#include <grid/point.h>
#include <vector.h>



/**
  Base class for quadrature formulae in arbitrary dimensions.
  */
template <int dim>
class Quadrature {
  public:
				     /**
				      * Number of quadrature points.
				      */
    const unsigned int n_quad_points;

				     /**
				      * Constructor.
				      */
    Quadrature (const unsigned int n_quad_points);
    
				     /**
				      * Return the #i#th quadrature point.
				      */
    const Point<dim> & quad_point (const unsigned int i) const;

				     /**
				      * Return the weight of the #i#th
				      * quadrature point.
				      */
    double weight (const unsigned int i) const;

  protected:
				     /**
				      * List of quadrature points. To be filled
				      * by the constructors of derived classes.
				      */
    vector<Point<dim> > quadrature_points;

				     /**
				      * List of weights of the quadrature points.
				      * To be filled by the constructors of
				      * derived classes.
				      */
    vector<double>      weights;
};





/*----------------------------   quadrature.h     ---------------------------*/
/* end of #ifndef __quadrature_H */
#endif
/*----------------------------   quadrature.h     ---------------------------*/
