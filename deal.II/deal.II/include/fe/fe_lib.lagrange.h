/*----------------------------   fe_lib.h     ---------------------------*/
/*      $Id$                 */
#ifndef __fe_lib_H
#define __fe_lib_H
/*----------------------------   fe_lib.h     ---------------------------*/


#include <fe/fe.h>



/**
 * Define a (bi-, tri-, etc)linear finite element in #dim# space dimensions,
 * along with (bi-, tri-)linear (therefore isoparametric) transforms from the
 * unit cell to the real cell.
 *
 * The linear, isoparametric mapping from a point $\vec \xi$ on the unit cell
 * to a point $\vec x$ on the real cell is defined as
 * $$ \vec x(\vec \xi)  = \sum_j {\vec p_j} N_j(\xi) $$
 * where $\vec p_j$ is the vector to the $j$th corner point of the cell in
 * real space and $N_j(\vec \xi)$ is the value of the basis function associated
 * with the $j$th corner point, on the unit cell at point $\vec \xi$. The sum
 * over $j$ runs over all corner points.
 *
 * @author Wolfgang Bangerth, 1998
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
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * For one dimensional elements, this
				      * function simply passes through to
				      * the one implemented in the base class.
				      * For higher dimensional finite elements
				      * we use linear mappings and therefore
				      * the boundary object is ignored since
				      * the boundary is approximated using
				      * piecewise straight boundary segments.
				      */
    virtual void fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
				 const vector<Point<dim> >            &unit_points,
				 vector<dFMatrix>    &jacobians,
				 const bool           compute_jacobians,
				 vector<Point<dim> > &ansatz_points,
				 const bool           compute_ansatz_points,
				 vector<Point<dim> > &q_points,
				 const bool           compute_q_points,
				 const Boundary<dim> &boundary) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_ansatz_points (const DoFHandler<dim>::cell_iterator &cell,
				    const Boundary<dim> &boundary,
				    vector<Point<dim> > &ansatz_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_ansatz_points (const DoFHandler<dim>::face_iterator &face,
					 const Boundary<dim> &boundary,
					 vector<Point<dim> > &ansatz_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * In two spatial dimensions, this function
				      * simply returns the length of the face.
				      */
    virtual void get_face_jacobians (const DoFHandler<dim>::face_iterator &face,
				     const Boundary<dim>         &boundary,
				     const vector<Point<dim-1> > &unit_points,
				     vector<double>      &face_jacobi_determinants) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * In two spatial dimensions, this function
				      * simply returns half the length of the
				      * whole face.
				      */
    virtual void get_subface_jacobians (const DoFHandler<dim>::face_iterator &face,
					const unsigned int           subface_no,
					const vector<Point<dim-1> > &unit_points,
					vector<double>      &face_jacobi_determinants) const;

				     /**
				      * Return the normal vectors to the
				      * face with number #face_no# of #cell#.
				      *
				      * For linear finite elements, this function
				      * is particularly simple since all normal
				      * vectors are equal and can easiliy be
				      * computed from the direction of the face
				      * without using the transformation (Jacobi)
				      * matrix, at least for two dimensions.
				      *
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int          face_no,
				     const Boundary<dim>         &boundary,
				     const vector<Point<dim-1> > &unit_points,
				     vector<Point<dim> >         &normal_vectors) const;    

				     /**
				      * Return the normal vectors to the
				      * subface with number #subface_no# of
				      * the face with number #face_no# of #cell#.
				      *
				      * For linear finite elements, this function
				      * is particularly simple since all normal
				      * vectors are equal and can easiliy be
				      * computed from the direction of the face
				      * without using the transformation (Jacobi)
				      * matrix, at least for two dimensions.
				      *
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int           face_no,
				     const unsigned int           subface_no,
				     const vector<Point<dim-1> > &unit_points,
				     vector<Point<dim> >         &normal_vectors) const;    

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;
};




/**
 * Define a (bi-, tri-, etc)quadratic finite element in #dim# space dimensions.
 * In one space dimension, a linear (subparametric) mapping from the unit cell
 * to the real cell is implemented.
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
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * For one dimensional elements, this
				      * function simply passes through to
				      * the one implemented in the base class.
				      */
    virtual void fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
				 const vector<Point<dim> >            &unit_points,
				 vector<dFMatrix>    &jacobians,
				 const bool           compute_jacobians,
				 vector<Point<dim> > &ansatz_points,
				 const bool           compute_ansatz_points,
				 vector<Point<dim> > &q_points,
				 const bool           compute_q_points,
				 const Boundary<dim> &boundary) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_ansatz_points (const DoFHandler<dim>::cell_iterator &cell,
				    const Boundary<dim> &boundary,
				    vector<Point<dim> > &ansatz_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_ansatz_points (const DoFHandler<dim>::face_iterator &face,
					 const Boundary<dim> &boundary,
					 vector<Point<dim> > &ansatz_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_jacobians (const DoFHandler<dim>::face_iterator &face,
				     const Boundary<dim>         &boundary,
				     const vector<Point<dim-1> > &unit_points,
				     vector<double>      &face_jacobi_determinants) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_subface_jacobians (const DoFHandler<dim>::face_iterator &face,
					const unsigned int           subface_no,
					const vector<Point<dim-1> > &unit_points,
					vector<double>      &face_jacobi_determinants) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int          face_no,
				     const Boundary<dim>         &boundary,
				     const vector<Point<dim-1> > &unit_points,
				     vector<Point<dim> >         &normal_vectors) const;    

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int           subface_no,
				     const unsigned int           face_no,
				     const vector<Point<dim-1> > &unit_points,
				     vector<Point<dim> >         &normal_vectors) const;    

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;
};




/**
 * Define a (bi-, tri-, etc)cubic finite element in #dim# space dimensions.
 * In one space dimension, a linear (subparametric) mapping from the unit cell
 * to the real cell is implemented.
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
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * For one dimensional elements, this
				      * function simply passes through to
				      * the one implemented in the base class.
				      */
    virtual void fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
				 const vector<Point<dim> >            &unit_points,
				 vector<dFMatrix>    &jacobians,
				 const bool           compute_jacobians,
				 vector<Point<dim> > &ansatz_points,
				 const bool           compute_ansatz_points,
				 vector<Point<dim> > &q_points,
				 const bool           compute_q_points,
				 const Boundary<dim> &boundary) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_ansatz_points (const DoFHandler<dim>::cell_iterator &cell,
				    const Boundary<dim> &boundary,
				    vector<Point<dim> > &ansatz_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_ansatz_points (const DoFHandler<dim>::face_iterator &face,
					 const Boundary<dim> &boundary,
					 vector<Point<dim> > &ansatz_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_jacobians (const DoFHandler<dim>::face_iterator &face,
				     const Boundary<dim>         &boundary,
				     const vector<Point<dim-1> > &unit_points,
				     vector<double>      &face_jacobi_determinants) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_subface_jacobians (const DoFHandler<dim>::face_iterator &face,
					const unsigned int           subface_no,
					const vector<Point<dim-1> > &unit_points,
					vector<double>      &face_jacobi_determinants) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int          face_no,
				     const Boundary<dim>         &boundary,
				     const vector<Point<dim-1> > &unit_points,
				     vector<Point<dim> >         &normal_vectors) const;    

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int           subface_no,
				     const unsigned int           face_no,
				     const vector<Point<dim-1> > &unit_points,
				     vector<Point<dim> >         &normal_vectors) const;    

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;
};



/*----------------------------   fe_lib.h     ---------------------------*/
/* end of #ifndef __fe_lib_H */
#endif
/*----------------------------   fe_lib.h     ---------------------------*/


