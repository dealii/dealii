
//----------------------------  fe.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe.h  ---------------------------
#ifndef __deal2__fe_h
#define __deal2__fe_h

#include <fe/fe_base.h>
#include <dofs/dof_handler.h>

template <int dim> class FEValuesData;
template <int dim> class FEValues;
template <int dim> class FEFaceValues;
template <int dim> class FESubfaceValues;
template <int dim> class FESystem;
template <int dim> class MatrixCreator;


/**
 * Common interface of all finite elements. Here, the functions to
 * fill the data fields of @ref{FEValues} are declared. While
 * @ref{FiniteElementBase} provides implementation of common
 * functionality, this class only serves as an abstract base class.
 *
 * The interface of this class is very restrictive. The reason is that
 * finite element values should be accessed only by use of
 * @ref{FEValues} objects. These, together with @p{FiniteElement} are
 * responsible to provide an optimized implementation.
 *
 * This even holds for evaluating finite elements at their support
 * points (provided the element is based on Lagrangian interpolation):
 * first, it is necessary to construct a quadrature rule from the
 * support points. This is then fed into an object of class
 * @ref{FEValues}. Even for evaluation on the unit cell, you will need
 * a triangulation containing that single cell.
 * 
 * @author Wolfgang Bangerth, Guido Kanschat, Ralf Hartmann, 1998, 2000, 2001
 */
template <int dim>
class FiniteElement : public FiniteElementBase<dim>
{
  private:
				     /**
				      * Copy constructor prohibited.
				      */
    FiniteElement(const FESystem<dim>&);

  public:
				     /**
				      * Constructor
				      */
    FiniteElement (const FiniteElementData<dim> &fe_data,
		   const std::vector<bool> &restriction_is_additive_flags);

				     /**
				      * Virtual destructor. Makes sure
				      * that pointers to this class
				      * are deleted properly.
				      */
    virtual ~FiniteElement ();
    
				     /**
				      * Number of base elements in a mixed
				      * discretization. This function returns
				      * 1 for primitive elements.
				      */
    virtual unsigned int n_base_elements () const;
    
				     /**
				      * Access to base element
				      * objects.  By default,
				      * #base_element(0)# is #this#.
				      * This function is overloaded by
				      * system elements to allow
				      * access to the different
				      * components of mixed
				      * discretizations.
				      */
    virtual const FiniteElement<dim> & base_element (const unsigned int index) const;
    
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
				      * Compute flags for initial
				      * update only.
				      * @see{FEValuesBase}
				      */
    virtual UpdateFlags update_once (UpdateFlags flags) const = 0;
  
				     /**
				      * Compute flags for update on
				      * each cell.
				      * @see{FEValuesBase}
				      */
    virtual UpdateFlags update_each (UpdateFlags flags) const = 0;
  
				     /**
				      * @p{clone} function instead of
				      * a copy constructor.
				      *
				      * This function is needed by the
				      * constructors of @p{FESystem}.
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
				      * @ref{FEValues}. This function
				      * performs all the operations
				      * needed to compute the data of an
				      * @p{FEValues} object.
				      *
				      * The same function in
				      * @p{mapping} must have been
				      * called for the same cell first!
				      */				      
    virtual void
    fill_fe_values (const Mapping<dim>                   &mapping,
		    const DoFHandler<dim>::cell_iterator &cell,
		    const Quadrature<dim>                &quadrature,
		    Mapping<dim>::InternalDataBase       &mapping_internal,
		    Mapping<dim>::InternalDataBase       &fe_internal,
		    FEValuesData<dim>                    &data) const = 0;
    
				     /**
				      * Fill the fields of
				      * @ref{FEFaceValues}. This function
				      * performs all the operations
				      * needed to compute the data of an
				      * @p{FEFaceValues} object.
				      *
				      * The same function in
				      * @p{mapping} must have been
				      * called for the same cell first!
				      */				      
    virtual void
    fill_fe_face_values (const Mapping<dim>                   &mapping,
			 const DoFHandler<dim>::cell_iterator &cell,
			 const unsigned int                    face_no,
			 const Quadrature<dim-1>              &quadrature,
			 Mapping<dim>::InternalDataBase       &mapping_internal,
			 Mapping<dim>::InternalDataBase       &fe_internal,
			 FEValuesData<dim>                    &data) const = 0;
    
				     /**
				      * Fill the fields of
				      * @ref{FESubfaceValues}. This function
				      * performs all the operations
				      * needed to compute the data of an
				      * @p{FESubfaceValues} object.
				      *
				      * The same function in
				      * @p{mapping} must have been
				      * called for the same cell first!
				      */				      
    virtual void
    fill_fe_subface_values (const Mapping<dim>                   &mapping,
			    const DoFHandler<dim>::cell_iterator &cell,
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
    friend class FEValues<dim>;
    friend class FEFaceValues<dim>;
    friend class FESubfaceValues<dim>;
    friend class FESystem<dim>;
};


#endif
