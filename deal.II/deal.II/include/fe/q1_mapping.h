/*----------------------------   q1_mapping.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __fe_q1_mapping_H
#define __fe_q1_mapping_H
/*----------------------------   q1_mapping.h     ---------------------------*/
// File has moved from fe_linear_mapping.h

#include <cmath>
#include <fe/fe.h>




/**
 * Implementation of Q1 transformation to the unit cell.
 * All finite element classes using a Q1 mapping of the grid cell to the
 * unit cell may be derived from this class. The grid transformation functions
 * are implemented here and do not have to be taken care of later.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FEQ1Mapping : public FiniteElement<dim>
{
  public:
				     /**
				      * Constructor. Simply passes through
				      * its arguments to the base class. For
				      * one space dimension, #dofs_per_quad#
				      * shall be zero; similarly, for one and
				      * two space dimensions, #dofs_per_hex#
				      * shall be zero.
				      */
    FEQ1Mapping (const unsigned int  dofs_per_vertex,
		 const unsigned int  dofs_per_line,
		 const unsigned int  dofs_per_quad,
		 const unsigned int  dofs_per_hex,
		 const unsigned int  n_components,
		 const vector<bool> &restriction_is_additive_flags);

    				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      * Here, the (bi-)linear basis functions
				      * are meant, which are used for the
				      * computation of the transformation from
				      * unit cell to real space cell.
				      */
    virtual double shape_value_transform (const unsigned int  i,
					  const Point<dim>   &p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      * Here, the (bi-)linear basis functions
				      * are meant, which are used for the
				      * computation of the transformation from
				      * unit cell to real space cell.
				      */
    virtual Tensor<1,dim> shape_grad_transform (const unsigned int  i,
						const Point<dim>   &p) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * In two spatial dimensions, this function
				      * simply returns the length of the face.
				      */
    virtual void get_face_jacobians (const DoFHandler<dim>::face_iterator &face,
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
				     const unsigned int           face_no,
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
				      * we use multilinear mappings.
				      */
    virtual void fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
				 const vector<Point<dim> >            &unit_points,
				 vector<Tensor<2,dim> >               &jacobians,
				 const bool              compute_jacobians,
				 vector<Tensor<3,dim> > &jacobians_grad,
				 const bool              compute_jacobians_grad,
				 vector<Point<dim> >    &support_points,
				 const bool              compute_support_points,
				 vector<Point<dim> >    &q_points,
				 const bool              compute_q_points,
				 const FullMatrix<double>         &shape_values_transform,
				 const vector<vector<Tensor<1,dim> > > &shape_grad_transform) const;

				     /**
				      * Compute the jacobian matrices
				      * of the mapping between unit
				      * and real cell at the given
				      * points on the unit cell.
				      */
    static void compute_jacobian_matrices (const DoFHandler<dim>::cell_iterator &cell,
					   const vector<Point<dim> >            &unit_points,
					   vector<Tensor<2,dim> >               &jacobians);

    				     /**
				      * Compute the gradients of the jacobian
				      * matrices of the mapping between unit
				      * and real cell at the given
				      * points on the unit cell.
				      */
    static void compute_jacobian_gradients (const DoFHandler<dim>::cell_iterator &cell,
					    const vector<Point<dim> >            &unit_points,
					    vector<Tensor<3,dim> >               &jacobians);


				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidData);
};




/*----------------------------   q1_mapping.h     ---------------------------*/
/* end of #ifndef __q1_mapping_H */
#endif
/*----------------------------   q1_mapping.h     ---------------------------*/
