//----------------------------  fe.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe.h  ---------------------------
#ifndef __deal2__fe_h
#define __deal2__fe_h

#include <base/config.h>
#include <fe/fe_base.h>
#include <dofs/dof_handler.h>

template <int dim> class FEValuesData;
template <int dim> class FEValuesBase;
template <int dim> class FEValues;
template <int dim> class FEFaceValues;
template <int dim> class FESubfaceValues;
template <int dim> class FESystem;

/*!@addtogroup febase */
/*@{*/

/**
 * Common interface of all finite elements. Here, the functions to
 * fill the data fields of FEValues are declared. While
 * FiniteElementBase provides implementation of common
 * functionality, this class only serves as an abstract base class.
 *
 * The interface of this class is very restrictive. The reason is that
 * finite element values should be accessed only by use of
 * FEValues objects. These, together with @p FiniteElement are
 * responsible to provide an optimized implementation.
 *
 * This even holds for evaluating finite elements at their support
 * points (provided the element is based on Lagrangian interpolation):
 * first, it is necessary to construct a quadrature rule from the
 * support points. This is then fed into an object of class
 * FEValues. Even for evaluation on the unit cell, you will need
 * a triangulation containing that single cell.
 *
 * Basically, this class just declares the shape function and their
 * derivatives on the unit cell $[0,1]^d$, and the means to transform
 * them onto a given cell in physical space if provided by the
 * FEValues class with a Mapping object.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, Ralf Hartmann, 1998, 2000, 2001
 */
template <int dim>
class FiniteElement : public FiniteElementBase<dim>
{
  public:
                                     /**
                                      * Copy constructor. This one is declared
                                      * as a public constructor to avoid
                                      * certain compiler errors when a copy
                                      * constructor is required even if it is
                                      * not executed (for example when binding
                                      * a temporary object to a constant
                                      * reference). However, if you try to
                                      * actually call it, it will throw an
                                      * exception, since copying finite
                                      * element objects is not really
                                      * supported. If you want to copy such an
                                      * object, use the @p clone function.
                                      */
    FiniteElement (const FiniteElement &);
    
				     /**
				      * Constructor
				      */
    FiniteElement (const FiniteElementData<dim> &fe_data,
		   const std::vector<bool> &restriction_is_additive_flags,
		   const std::vector<std::vector<bool> > &nonzero_components);

				     /**
				      * Virtual destructor. Makes sure
				      * that pointers to this class
				      * are deleted properly.
				      */
    virtual ~FiniteElement ();
    
				     /**
				      * Number of base elements in a
				      * mixed discretization.
				      *
				      * Note that even for vector
				      * valued finite elements, the
				      * number of components needs not
				      * coincide with the number of
				      * base elements, since they may
				      * be reused. For example, if you
				      * create a FESystem with
				      * three identical finite element
				      * classes by using the
				      * constructor that takes one
				      * finite element and a
				      * multiplicity, then the number
				      * of base elements is still one,
				      * although the number of
				      * components of the finite
				      * element is equal to the
				      * multiplicity.
				      */
    virtual unsigned int n_base_elements () const = 0;
    
				     /**
				      * Access to base element
				      * objects. If the element is
				      * scalar, then
				      * @p{base_element(0)} is
				      * @p this.
				      */
    virtual
    const FiniteElement<dim> &
    base_element (const unsigned int index) const = 0;

                                     /**
                                      * This index denotes how often
                                      * the base element @p index is
                                      * used in a composed element. If
                                      * the element is scalar, then
                                      * the result is always equal to
                                      * one. See the documentation for
                                      * the @p n_base_elements
                                      * function for more details.
                                      */
    virtual
    unsigned int
    element_multiplicity (const unsigned int index) const = 0;
    
				     /**
				      * Check for non-zero values on a face.
				      *
				      * This function returns
				      * @p true, if the shape
				      * function @p shape_index has
				      * non-zero values on the face
				      * @p face_index.
				      */
    virtual bool has_support_on_face (const unsigned int shape_index,
				      const unsigned int face_index) const = 0;
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      *
				      * This function is made virtual,
				      * since finite element objects
				      * are usually accessed through
				      * pointers to their base class,
				      * rather than the class itself.
				      */
    virtual unsigned int memory_consumption () const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcBoundaryFaceUsed);
				     /**
				      * Exception
				      */
    DeclException0 (ExcJacobiDeterminantHasWrongSign);
				     /**
				      * Exception
				      */
    DeclException1 (ExcComputationNotUseful,
		    int,
		    << "The computation you required from this function is not "
		    << "feasible or not probable in the present dimension ("
		    << arg1 << ") because it would be prohibitively expensive.");

  protected:

				     /**
				      * Determine the values a finite
				      * element should compute on
				      * initialization of data for
				      * @p FEValues.
				      *
				      * Given a set of flags
				      * indicating what quantities are
				      * requested from a @p FEValues
				      * object, @p update_once and
				      * @p update_each compute which
				      * values must really be
				      * computed. Then, the
				      * @p{fill_*_values} functions
				      * are called with the result of
				      * these.
				      *
				      * Furthermore, values must be
				      * computed either on the unit
				      * cell or on the physical
				      * cell. For instance, the
				      * function values of @p FE_Q do
				      * only depend on the quadrature
				      * points on the unit
				      * cell. Therefore, this flags
				      * will be returned by
				      * @p update_once. The gradients
				      * require computation of the
				      * covariant transformation
				      * matrix. Therefore,
				      * @p update_covariant_transformation
				      * and @p update_gradients will
				      * be returned by
				      * @p update_each.
				      *
				      * For an example see the same
				      * function in the derived class
				      * @p FE_Q.
				      */
    virtual UpdateFlags update_once (const UpdateFlags flags) const = 0;
  
				     /**
				      * Complementary function for
				      * @p update_once.
				      *
				      * While @p update_once returns
				      * the values to be computed on
				      * the unit cell for yielding the
				      * required data, this function
				      * determines the values that
				      * must be recomputed on each
				      * cell.
				      *
				      * Refer to @p update_once for
				      * more details.
				      */
    virtual UpdateFlags update_each (const UpdateFlags flags) const = 0;
  
				     /**
				      * @p clone function instead of
				      * a copy constructor.
				      *
				      * This function is needed by the
				      * constructors of @p FESystem.
				      */
    virtual FiniteElement<dim> *clone() const = 0;
    
				     /**
				      * Prepare internal data
				      * structures and fill in values
				      * independent of the
				      * cell. Returns a pointer to an
				      * object of which the caller of
				      * this function then has to
				      * assume ownership (which
				      * includes destruction when it
				      * is no more needed).
				      */
    virtual typename Mapping<dim>::InternalDataBase*
    get_data (const UpdateFlags      flags,
	      const Mapping<dim>    &mapping,
	      const Quadrature<dim> &quadrature) const = 0;

				     /**
				      * Prepare internal data
				      * structure for transformation
				      * of faces and fill in values
				      * independent of the
				      * cell. Returns a pointer to an
				      * object of which the caller of
				      * this function then has to
				      * assume ownership (which
				      * includes destruction when it
				      * is no more needed).
				      */
    virtual typename Mapping<dim>::InternalDataBase*
    get_face_data (const UpdateFlags        flags,
		   const Mapping<dim>      &mapping,
		   const Quadrature<dim-1> &quadrature) const;

				     /**
				      * Prepare internal data
				      * structure for transformation
				      * of children of faces and fill
				      * in values independent of the
				      * cell. Returns a pointer to an
				      * object of which the caller of
				      * this function then has to
				      * assume ownership (which
				      * includes destruction when it
				      * is no more needed).
				      */
    virtual typename Mapping<dim>::InternalDataBase*
    get_subface_data (const UpdateFlags        flags,
		      const Mapping<dim>      &mapping,
		      const Quadrature<dim-1> &quadrature) const;

				     /**
				      * Fill the fields of
				      * FEValues. This function
				      * performs all the operations
				      * needed to compute the data of an
				      * @p FEValues object.
				      *
				      * The same function in
				      * @p mapping must have been
				      * called for the same cell first!
				      */				      
    virtual void
    fill_fe_values (const Mapping<dim>                   &mapping,
		    const typename DoFHandler<dim>::cell_iterator &cell,
		    const Quadrature<dim>                &quadrature,
		    typename Mapping<dim>::InternalDataBase       &mapping_internal,
		    typename Mapping<dim>::InternalDataBase       &fe_internal,
		    FEValuesData<dim>                    &data) const = 0;
    
				     /**
				      * Fill the fields of
				      * FEFaceValues. This function
				      * performs all the operations
				      * needed to compute the data of an
				      * @p FEFaceValues object.
				      *
				      * The same function in
				      * @p mapping must have been
				      * called for the same cell first!
				      */				      
    virtual void
    fill_fe_face_values (const Mapping<dim>                   &mapping,
			 const typename DoFHandler<dim>::cell_iterator &cell,
			 const unsigned int                    face_no,
			 const Quadrature<dim-1>              &quadrature,
			 typename Mapping<dim>::InternalDataBase       &mapping_internal,
			 typename Mapping<dim>::InternalDataBase       &fe_internal,
			 FEValuesData<dim>                    &data) const = 0;
    
				     /**
				      * Fill the fields of
				      * FESubfaceValues. This function
				      * performs all the operations
				      * needed to compute the data of an
				      * @p FESubfaceValues object.
				      *
				      * The same function in
				      * @p mapping must have been
				      * called for the same cell first!
				      */				      
    virtual void
    fill_fe_subface_values (const Mapping<dim>                   &mapping,
			    const typename DoFHandler<dim>::cell_iterator &cell,
			    const unsigned int                    face_no,
			    const unsigned int                    sub_no,
			    const Quadrature<dim-1>              &quadrature,
			    typename Mapping<dim>::InternalDataBase &mapping_internal,
			    typename Mapping<dim>::InternalDataBase &fe_internal,
			    FEValuesData<dim>                    &data) const = 0;

				     /**
				      * Declare some other classes as
				      * friends of this class.
				      */
    friend class FEValuesBase<dim>;
    friend class FEValues<dim>;
    friend class FEFaceValues<dim>;
    friend class FESubfaceValues<dim>;
    friend class FESystem<dim>;
};

/*@}*/

#endif
