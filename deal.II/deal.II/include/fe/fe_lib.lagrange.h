/*----------------------------   fe_lib.h     ---------------------------*/
/*      $Id$                 */
#ifndef __fe_lib_H
#define __fe_lib_H
/*----------------------------   fe_lib.h     ---------------------------*/


#include <fe/fe.h>



/**
  Define a (bi-, tri-, etc)linear finite element in #dim# space dimensions,
  along with (bi-, tri-)linear (therefore isoparametric) transforms from the
  unit cell to the real cell.

  The linear, isoparametric mapping from a point $\vec \xi$ on the unit cell
  to a point $\vec x$ on the real cell is defined as
  $$ \vec x(\vec \xi)  = \sum_j {\vec p_j} N_j(\xi) $$
  where $\vec p_j$ is the vector to the $j$th corner point of the cell in
  real space and $N_j(\vec \xi)$ is the value of the basis function associated
  with the $j$th corner point, on the unit cell at point $\vec \xi$. The sum
  over $j$ runs over all corner points.
  */
template <int dim>
class FELinear : public FiniteElement<dim> {
  public:
				     /**
				      * Constructor
				      */
    FELinear ();

				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual Point<dim> shape_grad(const unsigned int i,
				  const Point<dim>& p) const;

				     /**
				      * Compute the Jacobian matrix and the
				      * quadrature points as well as the ansatz
				      * function locations on the real cell in
				      * real space from the given cell
				      * and the given quadrature points on the
				      * unit cell. The Jacobian matrix is to
				      * be computed at every quadrature point.
				      *
				      * Refer to the documentation of the
				      * \Ref{FEValues} class for a definition
				      * of the Jacobi matrix.
				      *
				      * For one dimensional elements, this
				      * function simply passes through to
				      * the one implemented in the base class.
				      * For two dimensional finite elements,
				      * these transformations are usually
				      * dependent on the actual finite element,
				      * which is expressed by the names
				      * sub- and isoparametric elements. This
				      * function is therefore not implemented
				      * by the FE<2> base class, but is made
				      * pure virtual.
				      *
				      * The function assumes that the fields
				      * already have the right number of
				      * elements.
				      */
    virtual void fill_fe_values (const Triangulation<dim>::cell_iterator &cell,
				 const vector<Point<dim> >               &unit_points,
				 vector<dFMatrix>    &jacobians,
				 const bool           compute_jacobians,
				 vector<Point<dim> > &ansatz_points,
				 const bool           compute_ansatz_points,
				 vector<Point<dim> > &q_points,
				 const bool           compute_q_points,
				 const Boundary<dim> &boundary) const;

				     /**
				      * Return the ansatz points this FE has
				      * on a face if a cell would have the
				      * given face as a side. Since we have no
				      * degrees of freedom on the faces for
				      * the linear ansatz, the ansatz points are
				      * simply the vertices of the face.
				      *
				      * The function assumes that the fields
				      * already have the right number of
				      * elements.
				      */
    virtual void face_ansatz_points (const Triangulation<dim>::face_iterator &face,
				     const Boundary<dim>  &boundary,
				     vector<Point<dim> >  &ansatz_points) const;
};




/**
  Define a (bi-, tri-, etc)quadratic finite element in #dim# space dimensions.
  In one space dimension, a linear (subparametric) mapping from the unit cell
  to the real cell is implemented.
  */
template <int dim>
class FEQuadratic : public FiniteElement<dim> {
  public:
				     /**
				      * Constructor
				      */
    FEQuadratic ();

				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual Point<dim> shape_grad(const unsigned int i,
				  const Point<dim>& p) const;

				     /**
				      * Compute the Jacobian matrix and the
				      * quadrature points as well as the ansatz
				      * function locations on the real cell in
				      * real space from the given cell
				      * and the given quadrature points on the
				      * unit cell. The Jacobian matrix is to
				      * be computed at every quadrature point.
				      *
				      * Refer to the documentation of the
				      * \Ref{FEValues} class for a definition
				      * of the Jacobi matrix.
				      *
				      * For one dimensional elements, this
				      * function simply passes through to
				      * the one implemented in the base class.
				      * For two dimensional finite elements,
				      * these transformations are usually
				      * dependent on the actual finite element,
				      * which is expressed by the names
				      * sub- and isoparametric elements. This
				      * function is therefore not implemented
				      * by the FE<2> base class, but is made
				      * pure virtual.
				      *
				      * The function assumes that the fields
				      * already have the right number of
				      * elements.
				      */
    virtual void fill_fe_values (const Triangulation<dim>::cell_iterator &cell,
				 const vector<Point<dim> >               &unit_points,
				 vector<dFMatrix>    &jacobians,
				 const bool           compute_jacobians,
				 vector<Point<dim> > &ansatz_points,
				 const bool           compute_ansatz_points,
				 vector<Point<dim> > &q_points,
				 const bool           compute_q_points,
				 const Boundary<dim> &boundary) const;
};




/**
  Define a (bi-, tri-, etc)cubic finite element in #dim# space dimensions.
  In one space dimension, a linear (subparametric) mapping from the unit cell
  to the real cell is implemented.
  */
template <int dim>
class FECubic : public FiniteElement<dim> {
  public:
				     /**
				      * Constructor
				      */
    FECubic ();

				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual Point<dim> shape_grad(const unsigned int i,
				  const Point<dim>& p) const;

				     /**
				      * Compute the Jacobian matrix and the
				      * quadrature points as well as the ansatz
				      * function locations on the real cell in
				      * real space from the given cell
				      * and the given quadrature points on the
				      * unit cell. The Jacobian matrix is to
				      * be computed at every quadrature point.
				      *
				      * Refer to the documentation of the
				      * \Ref{FEValues} class for a definition
				      * of the Jacobi matrix.
				      *
				      * For one dimensional elements, this
				      * function simply passes through to
				      * the one implemented in the base class.
				      * For two dimensional finite elements,
				      * these transformations are usually
				      * dependent on the actual finite element,
				      * which is expressed by the names
				      * sub- and isoparametric elements. This
				      * function is therefore not implemented
				      * by the FE<2> base class, but is made
				      * pure virtual.
				      *
				      * The function assumes that the fields
				      * already have the right number of
				      * elements.
				      */
    virtual void fill_fe_values (const Triangulation<dim>::cell_iterator &cell,
				 const vector<Point<dim> >               &unit_points,
				 vector<dFMatrix>    &jacobians,
				 const bool           compute_jacobians,
				 vector<Point<dim> > &ansatz_points,
				 const bool           compute_ansatz_points,
				 vector<Point<dim> > &q_points,
				 const bool           compute_q_points,
				 const Boundary<dim> &boundary) const;
};



/*----------------------------   fe_lib.h     ---------------------------*/
/* end of #ifndef __fe_lib_H */
#endif
/*----------------------------   fe_lib.h     ---------------------------*/
