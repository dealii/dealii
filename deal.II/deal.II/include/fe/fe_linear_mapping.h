/*----------------------------   fe_linear_mapping.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __fe_linear_mapping_H
#define __fe_linear_mapping_H
/*----------------------------   fe_linear_mapping.h     ---------------------------*/


#include <fe/fe.h>




/**
 * Abstract base class for concrete finite elements which use a
 * (bi-,tri-)linear mapping from the unit cell to the real cell. Some
 * functions can be singled out from these elements and are collected
 * in this one.
 */
template <int dim>
class FELinearMapping : public FiniteElement<dim> {
  public:
				     /**
				      * Constructor. Simply passes through
				      * its arguments to the base class.
				      */
    FELinearMapping (const unsigned int dofs_per_vertex,
		     const unsigned int dofs_per_line,
		     const unsigned int dofs_per_quad=0) :
		    FiniteElement<dim> (dofs_per_vertex,
					dofs_per_line,
					dofs_per_quad,
					GeometryInfo<dim>::vertices_per_cell) {};

    				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      * Here, the (bi-)linear basis functions
				      * are meant, which are used for the
				      * computation of the transformation from
				      * unit cell to real space cell.
				      */
    virtual double shape_value_transform (const unsigned int i,
					  const Point<dim> &p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      * Here, the (bi-)linear basis functions
				      * are meant, which are used for the
				      * computation of the transformation from
				      * unit cell to real space cell.
				      */
    virtual Tensor<1,dim> shape_grad_transform (const unsigned int i,
						const Point<dim> &p) const;

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
				      *
				      * For one dimensional elements, this
				      * function simply passes through to
				      * the one implemented in the base class.
				      * For higher dimensional finite elements
				      * we use linear mappings and therefore
				      * the boundary object is ignored since
				      * the boundary is approximated using
				      * piecewise multilinear boundary segments.
				      */
    virtual void fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
				 const vector<Point<dim> >            &unit_points,
				 vector<Tensor<2,dim> >               &jacobians,
				 const bool           compute_jacobians,
				 vector<Point<dim> > &support_points,
				 const bool           compute_support_points,
				 vector<Point<dim> > &q_points,
				 const bool           compute_q_points,
				 const dFMatrix      &shape_values_transform,
				 const vector<vector<Tensor<1,dim> > > &shape_grad_transform,
				 const Boundary<dim> &boundary) const;
};




/*----------------------------   fe_linear_mapping.h     ---------------------------*/
/* end of #ifndef __fe_linear_mapping_H */
#endif
/*----------------------------   fe_linear_mapping.h     ---------------------------*/
