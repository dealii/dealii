/*----------------------------   quadrature.h     ---------------------------*/
/*      $Id$                 */
#ifndef __quadrature_H
#define __quadrature_H
/*----------------------------   quadrature.h     ---------------------------*/


#include <grid/point.h>
#include <vector>



/**
 * Base class for quadrature formulae in arbitrary dimensions. This class
 * stores quadrature points and weights on the unit line [0,1], unit
 * square [0,1]x[0,1], etc. This information is used together with
 * objects of the \Ref{FiniteElement} class to compute the values stored
 * in the \Ref{FEValues} objects.
 *
 * There are a number of derived classes, denoting concrete integration
 * formulae. These are named by a prefixed #Q#, the name of the formula
 * (e.g. #Gauss#) and finally the order of integration. For example,
 * #QGauss2<dim># denotes a second order Gauss integration formula in
 * any dimension. Second order means that it integrates polynomials of
 * third order exact. In general, a formula of order #n# exactly
 * integrates polynomials of order #2n-1#.
 *
 * Most integration formulae in more than one space dimension are tensor
 * products of quadrature formulae in one space dimension, or more
 * generally the tensor product of a formula in #(dim-1)# dimensions and
 * one in one dimension. There is a special constructor to generate a
 * quadrature formula from two others.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class Quadrature {
  public:
				     /**
				      * Number of quadrature points.
				      */
    const unsigned int n_quadrature_points;

				     /**
				      * Constructor.
				      */
    Quadrature (const unsigned int n_quadrature_points);

				     /**
				      * Build this quadrature formula as the
				      * tensor product of a formula in a
				      * dimension one less than the present and
				      * a formula in one dimension.
				      */
    Quadrature (const Quadrature<dim-1> &,
		const Quadrature<1>     &);
    
				     /**
				      * Virtual destructor needed, since someone
				      * may want to use pointers to this class
				      * though the object pointed to is a derived
				      * class.
				      */
    virtual ~Quadrature ();
    
				     /**
				      * Return the #i#th quadrature point.
				      */
    const Point<dim> & quad_point (const unsigned int i) const;

				     /**
				      * Return a reference to the whole array of
				      * quadrature points.
				      */
    const vector<Point<dim> > & get_quad_points () const;
    
				     /**
				      * Return the weight of the #i#th
				      * quadrature point.
				      */
    double weight (const unsigned int i) const;

				     /**
				      * Return a reference to the whole array
				      * of weights.
				      */
    const vector<double> & get_weights () const;
    
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidIndex,
		    int, int,
		    << "The index " << arg1
		    << " is out of range, it should be less than " << arg2);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);
    
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




/**
 *  This class is a helper class to facilitate the usage of quadrature formulae
 *  on faces or subfaces of cells. It computes the locations of quadrature
 *  points on the unit cell from a quadrature object for a mannifold of
 *  one dimension less than that of the cell and the number of the face.
 *  For example, giving the Simpson rule in one dimension and using the
 *  #project_to_face# function with face number 1, the returned points will
 *  be $(1,0)$, $(1,0.5)$ and $(1,1)$. Note that faces have an orientation,
 *  so when projecting to face 3, you will get $(0,0)$, $(0,0.5)$ and $(0,1)$,
 *  which is in clockwise sense, while for face 1 the points were in
 *  counterclockwise sense.
 *
 *  For the projection to subfaces (i.e. to the children of a face of the
 *  unit cell), the same applies as above. Note the order in which the
 *  children of a face are numbered, which in two dimensions coincides
 *  with the orientation of the face.
 *  
 *  The different functions are grouped into a common class to avoid putting
 *  them into global namespace (and to make documentation easier, since
 *  presently the documentation tool can only handle classes, not global
 *  functions). However, since they have no local data, all functions are
 *  declared #static# and can be called without creating an object of this
 *  class.
 */
template <int dim>
class QProjector {
  public:
				     /**
				      * Compute the quadrature points on the
				      * cell if the given quadrature formula
				      * is used on face #face_no#. For further
				      * details, see the general doc for
				      * this class.
				      */
    static void project_to_face (const Quadrature<dim-1> &quadrature,
				 const unsigned int       face_no,
				 vector<Point<dim> >     &q_points);

    				     /**
				      * Compute the quadrature points on the
				      * cell if the given quadrature formula
				      * is used on face #face_no#, subface
				      * number #subface_no#. For further
				      * details, see the general doc for
				      * this class.
				      */
    static void project_to_subface (const Quadrature<dim-1> &quadrature,
				    const unsigned int       face_no,
				    const unsigned int       subface_no,
				    vector<Point<dim> >     &q_points);

				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidIndex,
		    int, int,
		    << "The index " << arg1
		    << " is out of range, it should be less than " << arg2);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);
};





/*----------------------------   quadrature.h     ---------------------------*/
/* end of #ifndef __quadrature_H */
#endif
/*----------------------------   quadrature.h     ---------------------------*/
