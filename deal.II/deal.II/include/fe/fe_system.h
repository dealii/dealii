//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_system_h
#define __deal2__fe_system_h


/*----------------------------   fe_system.h     ---------------------------*/


#include <base/config.h>
#include <fe/fe.h>
#include <vector>
#include <utility>


/*!@addtogroup febase */
/*@{*/

/**
 * This class provides an interface to group several elements together
 * into one. To the outside world, the resulting object looks just
 * like a usual finite element object, which is composed of several
 * other finite elements that are possibly of different type.
 *
 * The overall numbering of degrees of freedom is as follows: for each
 * subobject (vertex, line, quad, or hex), the degrees of freedom are
 * numbered such that we run over all subelements first, before
 * turning for the next dof on this subobject or for the next
 * subobject. For example, for a element of three components in one
 * space dimension, the first two components being cubic lagrange
 * elements and the third being a quadratic lagrange element, the
 * ordering for the system <tt>s=(u,v,p)</tt> is:
 *
 * <ul>
 * <li> First vertex: <tt>u0, v0, p0 = s0, s1, s2</tt>
 * <li> Second vertex: <tt>u1, v1, p1 = s3, s4, s5</tt>
 * <li> First component on the line:
 *   <tt>u2, u3 = s4, s5</tt>
 * <li> Second component on the line:
 *   <tt>v2, v3 = s6, s7</tt>.
 * <li> Third component on the line:
 *   <tt>p2 = s8</tt>.
 * </ul>
 * Do not rely on this numbering in your application as these
 * internals might change in future. Rather use the functions
 * @p system_to_component_index and @p component_to_system_index,
 * instead.
 *
 * In the most cases, the composed element behaves as if it were a usual element
 * with more degrees of freedom. However the underlying structure is visible in
 * the restriction, prolongation and interface constraint matrices, which do not
 * couple the degrees of freedom of the subobject. E.g. the continuity requirement
 * is imposed for the shape functions of the subobjects separately; no requirement
 * exist between shape functions of different subobjects, i.e. in the above
 * example: on a hanging node, the respective value of the @p u velocity is only
 * coupled to @p u at the vertices and the line on the larger cell next to this
 * vertex, there is no interaction with @p v and @p w of this or the other cell.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2002, 2003, partial reimplementation Ralf Hartmann 2001.
 */
template <int dim>
class FESystem : public FiniteElement<dim>
{
  public:

				     /**
				      * Constructor. Take a finite element type
				      * and the number of elements you want to
				      * group together using this class.
				      *
				      * In fact, the object @p fe is not used,
				      * apart from getting the number of dofs
				      * per vertex, line, etc for that finite
				      * element class. The objects creates its
				      * own copy of the finite element object
				      * at construction time (but after
				      * the initialization of the base class
				      * @p FiniteElement, which is why we need
				      * a valid finite element object passed
				      * to the constructor).
				      *
				      * Obviously, the template finite element
				      * class needs to be of the same dimension
				      * as is this object.
				      */
    FESystem (const FiniteElement<dim> &fe, const unsigned int n_elements);

				     /** 
				      * Constructor for mixed
				      * discretizations with two
				      * base elements.
				      *
				      * See the other constructor.
				      */
    FESystem (const FiniteElement<dim> &fe1, const unsigned int n1,
	      const FiniteElement<dim> &fe2, const unsigned int n2);

				     /** 
				      * Constructor for mixed
				      * discretizations with three
				      * base elements.
				      *
				      * See the other constructor.
				      */
    FESystem (const FiniteElement<dim> &fe1, const unsigned int n1,
	      const FiniteElement<dim> &fe2, const unsigned int n2,
	      const FiniteElement<dim> &fe3, const unsigned int n3);

				     /**
				      * Destructor.
				      */
    virtual ~FESystem ();

				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. This element returns
				      * a string that is composed of
				      * the strings
				      * @p name1...@p nameN returned
				      * by the basis elements. From
				      * these, we create a sequence
				      * <tt>FESystem<dim>[name1^m1-name2^m2-...-nameN^mN]</tt>,
				      * where @p mi are the
				      * multiplicities of the basis
				      * elements. If a multiplicity is
				      * equal to one, then the
				      * superscript is omitted.
				      */
    virtual std::string get_name () const;

				     /**
				      * Return the value of the
				      * @p ith shape function at the
				      * point @p p.  @p p is a point
				      * on the reference element. Since
				      * this finite element is always
				      * vector-valued, we return the
				      * value of the only non-zero
				      * component of the vector value
				      * of this shape function. If the
				      * shape function has more than
				      * one non-zero component (which
				      * we refer to with the term
				      * non-primitive), then throw an
				      * exception of type
				      * @p ExcShapeFunctionNotPrimitive.
				      *
				      * An
				      * @p ExcUnitShapeValuesDoNotExist
				      * is thrown if the shape values
				      * of the @p FiniteElement
				      * (corresponding to the @p ith
				      * shape function) depend on the
				      * shape of the cell in real
				      * space.
				      */
    virtual double shape_value (const unsigned int i,
			        const Point<dim> &p) const;

				     /**
				      * Return the value of the
				      * @p componentth vector
				      * component of the @p ith shape
				      * function at the point
				      * @p p. See the
				      * FiniteElementBase base
				      * class for more information
				      * about the semantics of this
				      * function.
				      *
				      * Since this element is vector
				      * valued in general, it relays
				      * the computation of these
				      * values to the base elements.
				      */
    virtual double shape_value_component (const unsigned int i,
					  const Point<dim> &p,
					  const unsigned int component) const;
    
				     /**
				      * Return the gradient of the
				      * @p ith shape function at the
				      * point @p p. @p p is a point
				      * on the reference element, and
				      * likewise the gradient is the
				      * gradient on the unit cell with
				      * respect to unit cell
				      * coordinates. Since
				      * this finite element is always
				      * vector-valued, we return the
				      * value of the only non-zero
				      * component of the vector value
				      * of this shape function. If the
				      * shape function has more than
				      * one non-zero component (which
				      * we refer to with the term
				      * non-primitive), then throw an
				      * exception of type
				      * @p ExcShapeFunctionNotPrimitive.
				      *
				      * An
				      * @p ExcUnitShapeValuesDoNotExist
				      * is thrown if the shape values
				      * of the @p FiniteElement
				      * (corresponding to the @p ith
				      * shape function) depend on the
				      * shape of the cell in real
				      * space.
				      */
    virtual Tensor<1,dim> shape_grad (const unsigned int  i,
				      const Point<dim>   &p) const;

				     /**
				      * Return the gradient of the
				      * @p componentth vector
				      * component of the @p ith shape
				      * function at the point
				      * @p p. See the
				      * FiniteElementBase base
				      * class for more information
				      * about the semantics of this
				      * function.
				      *
				      * Since this element is vector
				      * valued in general, it relays
				      * the computation of these
				      * values to the base elements.
				      */
    virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
						const Point<dim> &p,
						const unsigned int component) const;

				     /**
				      * Return the tensor of second
				      * derivatives of the @p ith
				      * shape function at point @p p
				      * on the unit cell. The
				      * derivatives are derivatives on
				      * the unit cell with respect to
				      * unit cell coordinates. Since
				      * this finite element is always
				      * vector-valued, we return the
				      * value of the only non-zero
				      * component of the vector value
				      * of this shape function. If the
				      * shape function has more than
				      * one non-zero component (which
				      * we refer to with the term
				      * non-primitive), then throw an
				      * exception of type
				      * @p ExcShapeFunctionNotPrimitive.
				      *
				      * An
				      * @p ExcUnitShapeValuesDoNotExist
				      * is thrown if the shape values
				      * of the @p FiniteElement
				      * (corresponding to the @p ith
				      * shape function) depend on the
				      * shape of the cell in real
				      * space.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim> &p) const;

				     /**
				      * Return the second derivatives
				      * of the @p componentth vector
				      * component of the @p ith shape
				      * function at the point
				      * @p p. See the
				      * FiniteElementBase base
				      * class for more information
				      * about the semantics of this
				      * function.
				      *
				      * Since this element is vector
				      * valued in general, it relays
				      * the computation of these
				      * values to the base elements.
				      */
    virtual
    Tensor<2,dim>
    shape_grad_grad_component (const unsigned int i,
                               const Point<dim> &p,
                               const unsigned int component) const;
    
				     /**
				      * Return the matrix
				      * interpolating from the given
				      * finite element to the present
				      * one. The size of the matrix is
				      * then @p dofs_per_cell times
				      * <tt>source.dofs_per_cell</tt>.
				      *
				      * These matrices are available
				      * if source and destination
				      * element are both @p FESystem
				      * elements, have the same number
				      * of base elements with same
				      * element multiplicity, and if
				      * these base elements also
				      * implement their
				      * @p get_interpolation_matrix
				      * functions. Otherwise, an
				      * exception of type
				      * FiniteElementBase<dim>::ExcInterpolationNotImplemented
				      * is thrown.
				      */
    virtual void
    get_interpolation_matrix (const FiniteElementBase<dim> &source,
			      FullMatrix<double>           &matrix) const;

                                     /** 
				      * Number of different base
				      * elements of this object.
				      *
				      * Since these objects can have
				      * multiplicity and subobjects
				      * themselves, this may be
				      * smaller than the total number
				      * of finite elements composed
				      * into this structure.
				      */
    virtual unsigned int n_base_elements () const;

				     /**
				      * How often is a composing
				      * element used.
				      */
    virtual unsigned int element_multiplicity (const unsigned int index) const;

				     /**
				      * Access to a composing
				      * element. The index needs to be
				      * smaller than the number of
				      * base elements. Note that the
				      * number of base elements may in
				      * turn be smaller than the
				      * number of components of the
				      * system element, if the
				      * multiplicities are greater
				      * than one.
				      */
    virtual const FiniteElement<dim> &
    base_element (const unsigned int index) const;
    
				     /**
				      * Check for non-zero values on a
				      * face.
				      *
				      * This function returns
				      * @p true, if the shape
				      * function @p shape_index has
				      * non-zero values on the face
				      * @p face_index.
				      *
				      * Implementation of the
				      * interface in
				      * FiniteElement
				      */
    virtual bool has_support_on_face (const unsigned int shape_index,
				      const unsigned int face_index) const;

                                     /**
                                      * Implementation of the
                                      * respective function in the
                                      * base class.
                                      */
    virtual
    Point<dim>
    unit_support_point (const unsigned int index) const;
    
                                     /**
                                      * Implementation of the
                                      * respective function in the
                                      * base class.
                                      */
    virtual
    Point<dim-1>
    unit_face_support_point (const unsigned int index) const;
    
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
  
  protected:
				     /**
				      * Compute flags for initial
				      * update only.
				      */
    virtual UpdateFlags update_once (const UpdateFlags flags) const;
  
				     /**
				      * Compute flags for update on
				      * each cell.
				      */
    virtual UpdateFlags update_each (const UpdateFlags flags) const;

				     /**
				      * @p clone function instead of
				      * a copy constructor.
				      *
				      * This function is needed by the
				      * constructors of @p FESystem.
				      */
    virtual FiniteElement<dim> * clone() const;
  
				     /**
				      * Prepare internal data
				      * structures and fill in values
				      * independent of the cell.
				      */
    virtual typename Mapping<dim>::InternalDataBase*
    get_data (const UpdateFlags      update_flags,
	      const Mapping<dim>    &mapping,
	      const Quadrature<dim> &quadrature) const ;

				     /**
				      * Implementation of the same
				      * function in
				      * FiniteElement.
				      *
				      * Passes on control to
				      * @p compute_fill that does the
				      * work for all three
				      * <tt>fill_fe*_values</tt>
				      * functions.
				      */
    virtual void
    fill_fe_values (const Mapping<dim>                      &mapping,
		    const typename Triangulation<dim>::cell_iterator &cell,
		    const Quadrature<dim>                   &quadrature,
		    typename Mapping<dim>::InternalDataBase &mapping_data,
		    typename Mapping<dim>::InternalDataBase &fe_data,
		    FEValuesData<dim>                       &data) const;

				     /**
				      * Implementation of the same
				      * function in
				      * FiniteElement.
				      *
				      * Passes on control to
				      * @p compute_fill that does the
				      * work for all three
				      * <tt>fill_fe*_values</tt> functions.
				      */    
    virtual void
    fill_fe_face_values (const Mapping<dim>                   &mapping,
			 const typename Triangulation<dim>::cell_iterator &cell,
			 const unsigned int                    face_no,
			 const Quadrature<dim-1>              &quadrature,
			 typename Mapping<dim>::InternalDataBase      &mapping_data,
			 typename Mapping<dim>::InternalDataBase      &fe_data,
			 FEValuesData<dim>                    &data) const ;

    				     /**
				      * Implementation of the same
				      * function in
				      * FiniteElement.
				      *
				      * Passes on control to
				      * @p compute_fill that does the
				      * work for all three
				      * <tt>fill_fe*_values</tt> functions.
				      */
    virtual void
    fill_fe_subface_values (const Mapping<dim>                   &mapping,
			    const typename Triangulation<dim>::cell_iterator &cell,
			    const unsigned int                    face_no,
			    const unsigned int                    sub_no,
			    const Quadrature<dim-1>              &quadrature,
			    typename Mapping<dim>::InternalDataBase      &mapping_data,
			    typename Mapping<dim>::InternalDataBase      &fe_data,
			    FEValuesData<dim>                    &data) const ;
    

				     /**
				      * Do the work for the three
				      * <tt>fill_fe*_values</tt> functions.
				      * 
				      * Calls (among other things)
				      * <tt>fill_fe_([sub]face)_values</tt>
				      * of the base elements. Calls
				      * @p fill_fe_values if
				      * <tt>face_no==invalid_face_no</tt>
				      * and
				      * <tt>sub_no==invalid_face_no</tt>;
				      * calls @p fill_fe_face_values
				      * if
				      * <tt>face_no==invalid_face_no</tt>
				      * and
				      * <tt>sub_no!=invalid_face_no</tt>;
				      * and calls
				      * @p fill_fe_subface_values if
				      * <tt>face_no!=invalid_face_no</tt>
				      * and
				      * <tt>sub_no!=invalid_face_no</tt>.
				      */
    template <int dim_1>
    void compute_fill (const Mapping<dim>                   &mapping,
		       const typename Triangulation<dim>::cell_iterator &cell,
		       const unsigned int                    face_no,
		       const unsigned int                    sub_no,
		       const Quadrature<dim_1>              &quadrature,
		       typename Mapping<dim>::InternalDataBase &mapping_data,
		       typename Mapping<dim>::InternalDataBase &fe_data,
		       FEValuesData<dim>                    &data) const ;

  private:

				     /**
				      * Value to indicate that a given
				      * face or subface number is
				      * invalid.
				      */
    static const unsigned int invalid_face_number = deal_II_numbers::invalid_unsigned_int;
    
				     /**
				      * Pairs of multiplicity and
				      * element type.
				      */
    typedef std::pair<const FiniteElement<dim> *, unsigned int> ElementPair;
    
				     /**
				      * Pointer to underlying finite
				      * element classes.
				      *
				      * This object contains a pointer
				      * to each contributing element
				      * of a mixed discretization and
				      * its multiplicity. It is
				      * created by the constructor and
				      * constant afterwards.
				      */
    std::vector<ElementPair> base_elements;


				     /**
				      * Initialize the
				      * @p unit_support_points field
				      * of the FiniteElementBase
				      * class. Called from the
				      * constructor.
				      */
    void initialize_unit_support_points ();

				     /**
				      * Initialize the
				      * @p unit_face_support_points field
				      * of the FiniteElementBase
				      * class. Called from the
				      * constructor.
				      */
    void initialize_unit_face_support_points ();

				     /**
				      * Helper function used in the constructor:
				      * take a @p FiniteElementData object
				      * and return an object of the same type
				      * with the number of degrees of
				      * freedom per vertex, line, etc.
				      * multiplied by @p n. Don't touch the
				      * number of functions for the
				      * transformation from unit to real
				      * cell.
				      */
    static FiniteElementData<dim>
    multiply_dof_numbers (const FiniteElementData<dim> &fe_data,
			  const unsigned int            N);
    
				     /**
				      * Same as above for mixed elements
				      * with two different sub-elements.
				      */
    static FiniteElementData<dim>
    multiply_dof_numbers (const FiniteElementData<dim> &fe1,
			  const unsigned int            N1,
			  const FiniteElementData<dim> &fe2,
			  const unsigned int            N2);

				     /**
				      * Same as above for mixed elements
				      * with three different sub-elements.
				      */
    static FiniteElementData<dim>
    multiply_dof_numbers (const FiniteElementData<dim> &fe1,
			  const unsigned int            N1,
			  const FiniteElementData<dim> &fe2,
			  const unsigned int            N2,
			  const FiniteElementData<dim> &fe3,
			  const unsigned int            N3);


				     /**
				      * Helper function used in the
				      * constructor: takes a
				      * @p FiniteElement object and
				      * returns an boolean vector
				      * including the
				      * @p restriction_is_additive_flags
				      * of the mixed element
				      * consisting of @p N elements
				      * of the sub-element @p fe.
				      */
    static std::vector<bool>
    compute_restriction_is_additive_flags (const FiniteElement<dim> &fe,
					   const unsigned int        N);
    
				     /**
				      * Same as above for mixed elements
				      * with two different sub-elements.
				      */
    static std::vector<bool>
    compute_restriction_is_additive_flags (const FiniteElement<dim> &fe1,
					   const unsigned int        N1,
					   const FiniteElement<dim> &fe2,
					   const unsigned int        N2);

				     /**
				      * Same as above for mixed elements
				      * with three different sub-elements.
				      */
    static std::vector<bool>
    compute_restriction_is_additive_flags (const FiniteElement<dim> &fe1,
					   const unsigned int        N1,
					   const FiniteElement<dim> &fe2,
					   const unsigned int        N2,
					   const FiniteElement<dim> &fe3,
					   const unsigned int        N3);

				     /**
				      * Compute the named flags for a
				      * list of finite elements with
				      * multiplicities given in the
				      * second argument. This function
				      * is called from all the above
				      * functions.
				     */
    static std::vector<bool>
    compute_restriction_is_additive_flags (const std::vector<const FiniteElement<dim>*> &fes,
                                           const std::vector<unsigned int>              &multiplicities);
    
    
				     /**
				      * Compute the non-zero vector
				      * components of a composed
				      * finite element.
				      */
    static std::vector<std::vector<bool> >
    compute_nonzero_components (const FiniteElement<dim> &fe1,
				const unsigned int        N1);

				     /**
				      * Compute the non-zero vector
				      * components of a composed
				      * finite element.
				      */
    static std::vector<std::vector<bool> >
    compute_nonzero_components (const FiniteElement<dim> &fe1,
				const unsigned int        N1,
				const FiniteElement<dim> &fe2,
				const unsigned int        N2);

				     /**
				      * Compute the non-zero vector
				      * components of a composed
				      * finite element.
				      */
    static std::vector<std::vector<bool> >
    compute_nonzero_components (const FiniteElement<dim> &fe1,
				const unsigned int        N1,
				const FiniteElement<dim> &fe2,
				const unsigned int        N2,
				const FiniteElement<dim> &fe3,
				const unsigned int        N3);

				     /**
				      * Compute the nonzero components
				      * of a list of finite elements
				      * with multiplicities given in
				      * the second argument. This
				      * function is called from all
				      * the above functions.
				     */
    static std::vector<std::vector<bool> >
    compute_nonzero_components (const std::vector<const FiniteElement<dim>*> &fes,
				const std::vector<unsigned int>              &multiplicities);
    
				     /**
				      * This function is simply
				      * singled out of the
				      * constructors since there are
				      * several of them. It sets up
				      * the index table for the system
				      * as well as @p restriction and
				      * @p prolongation
				      * matrices.
				      */
    void initialize();

				     /**
				      * Used by @p initialize.
				      */
    void build_cell_tables();
    
				     /**
				      * Used by @p initialize.
				      */
    void build_face_tables();

				     /**
				      * Used by @p initialize.
				      */
    void build_interface_constraints ();

				     /**
				      * Usually: Fields of
				      * cell-independent data.
				      *
				      * However, here, this class does
				      * not itself store the data but
				      * only pointers to
				      * @p InternalData objects for
				      * each of the base elements.
				      */
    class InternalData : public FiniteElementBase<dim>::InternalDataBase
    {
      public:
					 /**
					  * Constructor. Is called by
					  * the @p get_data
					  * function. Sets the size of
					  * the @p base_fe_datas
					  * vector to
					  * @p n_base_elements.
					  */
	InternalData (const unsigned int n_base_elements);
	
					 /**
					  * Destructor. Deletes all
					  * @p InternalDatas whose
					  * pointers are stored by the
					  * @p base_fe_datas
					  * vector.
					  */
	~InternalData();

					 /**
					  * Flag indicating whether
					  * second derivatives shall
					  * be computed.
					  */
	bool compute_second_derivatives;
	
					 /**
					  * Gives write-access to the
					  * pointer to a
					  * @p InternalData of the
					  * @p base_noth base
					  * element.
					  */
	void set_fe_data(const unsigned int                        base_no,
			 typename FiniteElementBase<dim>::InternalDataBase *);

					 /**
					  * Gives read-access to the
					  * pointer to a
					  * @p InternalData of the
					  * @p base_noth base element.
					  */	
	typename FiniteElementBase<dim>::InternalDataBase &
	get_fe_data (const unsigned int base_no) const;


					 /**
					  * Gives write-access to the
					  * pointer to a
					  * @p FEValuesData for the
					  * @p base_noth base
					  * element.
					  */
	void set_fe_values_data (const unsigned int base_no,
				 FEValuesData<dim> *);

					 /**
					  * Gives read-access to the
					  * pointer to a
					  * @p FEValuesData for the
					  * @p base_noth base element.
					  */	
	FEValuesData<dim> & get_fe_values_data (const unsigned int base_no) const;

					 /**
					  * Deletes the
					  * @p FEValuesData the
					  * <tt>fe_datas[base_no]</tt>
					  * pointer is pointing
					  * to. Sets
					  * <tt>fe_datas[base_no]</tt> to
					  * zero.
					  *
					  * This function is used to
					  * delete @p FEValuesData
					  * that are needed only on
					  * the first cell but not any
					  * more afterwards.  This is
					  * the case for
					  * e.g. Lagrangian elements
					  * (see e.g. @p FE_Q
					  * classes).
					  */
	void delete_fe_values_data (const unsigned int base_no);

                                         /**
                                          * Set the @p first_cell
                                          * flag to @p false. Used by
                                          * the @p FEValues class to
                                          * indicate that we have
                                          * already done the work on
                                          * the first cell.
                                          *
                                          * In addition to calling the
                                          * respective function of the
                                          * base class, this function
                                          * also calls the functions
                                          * of the sub-data objects.
                                          */
        virtual void clear_first_cell ();
	
      private:
	
					 /**
					  * Pointers to
					  * @p InternalData objects
					  * for each of the base
					  * elements. They are
					  * accessed to by the
					  * @p set_ and
					  * @p get_fe_data functions.
					  *
					  * The size of this vector is
					  * set to @p n_base_elements
					  * by the InternalData
					  * constructor.  It is filled
					  * by the @p get_data
					  * function. Note that since
					  * the data for each instance
					  * of a base class is
					  * necessarily the same, we
					  * only need as many of these
					  * objects as there are base
					  * elements, irrespective of
					  * their multiplicity.
					  */
	typename std::vector<typename FiniteElementBase<dim>::InternalDataBase *> base_fe_datas;

					 /**
					  * Pointers to the
					  * @p FEValuesData objects
					  * that are given to the
					  * @p fill_fe_values
					  * function of the base
					  * elements. They are
					  * accessed to by the
					  * @p set_ and
					  * @p get_fe_values_data
					  * functions.
					  *
					  * The size of this vector is
					  * set to @p n_base_elements
					  * by the InternalData
					  * constructor.
					  */
	std::vector<FEValuesData<dim> *> base_fe_values_datas;
    };
};

/*@}*/

/* -------------- declaration of explicit specializations ------------- */


template <> void FESystem<1>::initialize_unit_face_support_points ();



/*----------------------------  fe_system.h  ---------------------------*/
#endif
/*----------------------------  fe_system.h  ---------------------------*/
