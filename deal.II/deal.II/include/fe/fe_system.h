/*----------------------------   fe_lib.system.h     ---------------------------*/
/*      $Id$                 */
#ifndef __fe_system_H
#define __fe_system_H
/*----------------------------   fe_lib.system.h     ---------------------------*/


#include <fe/fe.h>



template <int dim>
class FESystem : public FiniteElement<dim> {
  public:

				     /**
				      * Constructor. Take a finite element type
				      * and the number of elements you want to
				      * group together using this class.
				      *
				      * In fact, the object #fe# is not used,
				      * apart from getting the number of dofs
				      * per vertex, line, etc for that finite
				      * element class. For this, it would have
				      * been possible to use the #get_fe_data#
				      * function that each element has to
				      * provide. The correct way to write
				      * this constructor would therefore have
				      * been to specify it without the first
				      * argument and let the user specify the
				      * desired finite element by an explicit
				      * template argument list, like this:
				      * #AnyClass::f<int>()#. However, #C++#
				      * does not allow this call sequence for
				      * constructors, so we have to use the
				      * way as shown here, to let the compiler
				      * deduce the template argument itself.
				      *
				      * Obviously, the tenplate finite element
				      * class needs to be of the same dimension
				      * as is this object.
				      */
    template <typename FE>
    FESystem (const FE &fe, const unsigned int n_elements);

				     /**
				      * Destructor.
				      */
    virtual ~FESystem ();

    				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      *
				      * For an element composed of #N#
				      * subelements, the first #N# shape
				      * functions refer to the zeroth shape
				      * function of the underlying object,
				      * the shape functions #N..2N-1# refer
				      * to the base shape function with
				      * number #1#, and so on. The #i# shape
				      * function therefore equals the
				      * #i/N# the shape function of the
				      * base object.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>  &p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      *
				      * For the ordering of shape functions
				      * refer to the #shape_value# function.
				      */
    virtual Tensor<1,dim> shape_grad(const unsigned int i,
				     const Point<dim>& p) const;

				     /**
				      * Return the tensor of second derivatives
				      * of the #i#th shape function at
				      * point #p# on the unit cell.
				      *
				      * For the ordering of shape functions
				      * refer to the #shape_value# function.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * For the ordering of shape functions
				      * refer to the #shape_value# function.
				      */
    virtual void get_unit_support_points (vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * For the ordering of shape functions
				      * refer to the #shape_value# function.
				      */
    virtual void get_support_points (const DoFHandler<dim>::cell_iterator &cell,
				     const Boundary<dim> &boundary,
				     vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  const Boundary<dim> &boundary,
					  vector<Point<dim> > &support_points) const;

    				     /**
				      * Fill the local mass matrix. The elements
				      * of this matrix are the integrals
				      * $\int_K \phi_i \phi_j dx$ over a given
				      * cell $K$. However, here only those
				      * elements of the matrix are set for which
				      * the shape functions $\phi_i$ and
				      * $\phi_j$ belong to the same subelement,
				      * i.e. the resulting matrix is a block
				      * matrix where each block is a diagonal
				      * matrix with diagonal values equal to
				      * the respective entry of the local mass
				      * matrix for the underlying finite element
				      * class. This definition of the mass
				      * matrix for systems of finite elements
				      * is consistent with the use of the matrix
				      * for the projection of initial values and
				      * the like, where the subelements are not
				      * coupled to each other. Also in most
				      * other cases you will not want the
				      * coupling terms to appear in the mass
				      * matrix.
				      *
				      * If the shape functions of this element
				      * were numbered such that the first
				      * numbers are for the shape functions of
				      * the first subelement, then those for
				      * the second subelement, and so on, then
				      * the mass matrix generated by this
				      * function would be a block diagonal
				      * matrix with each block being the mass
				      * matrix of the base finite element.
				      *
				      * Refer to the base class for more
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;

				     /**
				      * Return the value of the #i#th shape
				      * function of the transformation mapping
				      * from unit cell to real cell. Since
				      * the transform functions are not
				      * touched when clustering several finite
				      * element objects together using this
				      * class, this function simply passes down
				      * the call to the respective function of
				      * the underlying element.
				      */
    virtual double shape_value_transform (const unsigned int i,
					  const Point<dim> &p) const;

				     /**
				      * Same as above: return gradient of the
				      * #i#th shape function for the mapping
				      * from unit to real cell.
				      */
    virtual Tensor<1,dim> shape_grad_transform (const unsigned int i,
						const Point<dim> &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * Since this function is only about the
				      * mapping from unit to real cell, it
				      * is not affected by putting several
				      * equal elements together, so this
				      * function simply passes down to the
				      * underlying object.
				      */
    virtual void get_face_jacobians (const DoFHandler<dim>::face_iterator &face,
				     const Boundary<dim>         &boundary,
				     const vector<Point<dim-1> > &unit_points,
				     vector<double>      &face_jacobi_determinants) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * Since this function is only about the
				      * mapping from unit to real cell, it
				      * is not affected by putting several
				      * equal elements together, so this
				      * function simply passes down to the
				      * underlying object.
				      */
    virtual void get_subface_jacobians (const DoFHandler<dim>::face_iterator &face,
					const unsigned int           subface_no,
					const vector<Point<dim-1> > &unit_points,
					vector<double>      &face_jacobi_determinants) const;

				     /**
				      * Return the normal vectors to the
				      * face with number #face_no# of #cell#.
				      *
				      * Since this function is only about the
				      * mapping from unit to real cell, it
				      * is not affected by putting several
				      * equal elements together, so this
				      * function simply passes down to the
				      * underlying object.
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
				      * Since this function is only about the
				      * mapping from unit to real cell, it
				      * is not affected by putting several
				      * equal elements together, so this
				      * function simply passes down to the
				      * underlying object.
				      *
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int           face_no,
				     const unsigned int           subface_no,
				     const vector<Point<dim-1> > &unit_points,
				     vector<Point<dim> >         &normal_vectors) const;    
    
  private:
				     /**
				      * Pointer to an object of the underlying
				      * finite element class. This object is
				      * created by the constructor.
				      */
    const FiniteElement<dim> *const base_element;

				     /**
				      * Number of subelements of this object.
				      * Since these objects may have
				      * subobjects themselves, this may be
				      * smaller than the total number of finite
				      * elements composed into this structure.
				      */
    const unsigned int n_sub_elements;
    
				     /**
				      * Helper function used in the constructor:
				      * take a #FiniteElementData# object
				      * and return an object of the same type
				      * with the number of degrees of
				      * freedom per vertex, line, etc.
				      * multiplied by #n#. Don't touch the
				      * number of functions for the
				      * transformation from unit to real
				      * cell.
				      */
    static FiniteElementData<dim> multiply_dof_numbers (const FiniteElementData<dim> &fe_data,
							const unsigned int            N);
    
				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * #restriction# and #prolongation#
				      * matrices. Since the operation of this
				      * function can be done without explicit
				      * knowledge of the data type of the
				      * underlying finite element class, we
				      * don't want to have this function in
				      * the general template definition in
				      * the #.h# file.
				      */
    void initialize_matrices ();
};





/* ------------------------- template functions ------------------------- */

template <int dim>
template <typename FE>
FESystem<dim>::FESystem (const FE &fe, const unsigned int n_elements) :
		FiniteElement (multiply_dof_numbers(fe, n_elements)),
		base_element (new FE()),
		n_sub_elements (n_elements)
{
  base_element->subscribe ();
  initialize_matrices ();
};




/*----------------------------   fe_lib.system.h     ---------------------------*/
/* end of #ifndef __fe_system_H */
#endif
/*----------------------------   fe_lib.system.h     ---------------------------*/
