//----------------------------  fe_system.h  ---------------------------
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
//----------------------------  fe_system.h  ---------------------------
#ifndef __deal2__fe_system_h
#define __deal2__fe_system_h


/*----------------------------   fe_system.h     ---------------------------*/


#include <fe/fe.h>
#include <vector>
#include <utility>


/**
 * This class provides an interface to group several elements together
 * into one. To the outside world, the resulting object looks just
 * like a usual finite element object, which is composed of several
 * other finite elements that are possibly of different type.

 * The overall numbering of degrees of freedom is as follows: for each
 * subobject (vertex, line, quad, or hex), the degrees of freedom are
 * numbered such that we run over all subelements first, before
 * turning for the next dof on this subobject or for the next
 * subobject. For example, for a element of three components in one
 * space dimension, the first two components being cubic lagrange
 * elements and the third being a quadratic lagrange element, the
 * ordering for the system @p{s=(u,v,p)} is:
 *
 * @begin{itemize}
 * @item First vertex: @p{u0, v0, p0 = s0, s1, s2}
 * @item Second vertex: @p{u1, v1, p1 = s3, s4, s5}
 * @item First component on the line:
 *   @p{u2, u3 = s4, s5}
 * @item Second component on the line:
 *   @p{v2, v3 = s6, s7}.
 * @item Third component on the line:
 *   @p{p2 = s8}.
 * @end{itemize}
 * Do not rely on this numbering in your application as these
 * internals might change in future. Rather use the functions
 * @p{system_to_component_index} and @p{component_to_system_index},
 * instead.
 *
 * In the most cases, the composed element behaves as if it were a usual element
 * with more degrees of freedom. However the underlying structure is visible in
 * the restriction, prolongation and interface constraint matrices, which do not
 * couple the degrees of freedom of the subobject. E.g. the continuity requirement
 * is imposed for the shape functions of the subobjects separately; no requirement
 * exist between shape functions of different subobjects, i.e. in the above
 * example: on a hanging node, the respective value of the @p{u} velocity is only
 * coupled to @p{u} at the vertices and the line on the larger cell next to this
 * vertex, there is no interaction with @p{v} and @p{w} of this or the other cell.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, reimplementation Ralf Hartmann 2001.
 */
template <int dim>
class FESystem : public FiniteElement<dim>
{
				     /**
				      * Copy constructor prohibited.
				      */
    FESystem(const FESystem<dim>&);

  public:

				     /**
				      * Constructor. Take a finite element type
				      * and the number of elements you want to
				      * group together using this class.
				      *
				      * In fact, the object @p{fe} is not used,
				      * apart from getting the number of dofs
				      * per vertex, line, etc for that finite
				      * element class. The objects creates its
				      * own copy of the finite element object
				      * at construction time (but after
				      * the initialization of the base class
				      * @p{FiniteElement}, which is why we need
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
    virtual unsigned int n_base_elements() const;

				     /**
				      * How often is a composing element used.
				      *
				      */
    unsigned int element_multiplicity(unsigned int index) const;

				     /**
				      * Access to a composing element.
				      *
				      * If you assemble your system
				      * matrix, you usually will not
				      * want to have an FEValues object
				      * with a lot of equal entries. Ok,
				      * so initialize your FEValues with
				      * the @p{base_element} you get by
				      * this function. In a mixed
				      * discretization, you can choose
				      * the different base element types
				      * by index.
				      *
				      */
    virtual const FiniteElement<dim> & base_element(unsigned int index) const;

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
				      * Return the support points of
				      * the trial functions on the
				      * unit cell.
				      *
				      * The order of points in the
				      * array matches that returned by
				      * the @p{cell->get_dof_indices}
				      * function, but:
				      *
				      * If one of the base elements
				      * has no support points, then it
				      * makes no sense to define
				      * support points for the
				      * composed element, so return an
				      * empty array to demonstrate
				      * that fact.
				      */
    virtual void get_unit_support_points (typename std::vector<Point<dim> > &) const;    

				     /**
				      * Return the support points of
				      * the trial functions on the
				      * first face of the unit cell.
				      *
				      * The order of points in the
				      * array matches that returned by
				      * the @p{cell->get_dof_indices}
				      * function, but:
				      *
				      * If one of the base elements
				      * has no support points, then it
				      * makes no sense to define
				      * support points for the
				      * composed element, so return an
				      * empty array to demonstrate
				      * that fact.
				      */
    virtual void get_unit_face_support_points (typename std::vector<Point<dim-1> > &) const;    
  
  protected:
				     /**
				      * Compute flags for initial
				      * update only.
				      */
    virtual UpdateFlags update_once (UpdateFlags flags) const;
  
				     /**
				      * Compute flags for update on
				      * each cell.
				      */
    virtual UpdateFlags update_each (UpdateFlags flags) const;

				     /**
				      * @p{clone} function instead of
				      * a copy constructor.
				      *
				      * This function is needed by the
				      * constructors of @p{FESystem}.
				      */
    virtual FiniteElement<dim> *clone() const;
  
				     /**
				      * Prepare internal data
				      * structures and fill in values
				      * independent of the cell.
				      */
    virtual typename Mapping<dim>::InternalDataBase*
    get_data (const UpdateFlags,
	      const Mapping<dim>& mapping,
	      const Quadrature<dim>& quadrature) const ;

				     /**
				      * Implementation of the same
				      * function in
				      * @ref{FiniteElement}.
				      */
    virtual void
    fill_fe_values (const Mapping<dim>                   &mapping,
		    const typename DoFHandler<dim>::cell_iterator &cell,
		    const Quadrature<dim>                &quadrature,
		    Mapping<dim>::InternalDataBase      &mapping_data,
		    Mapping<dim>::InternalDataBase      &fe_data,
		    FEValuesData<dim>                    &data) const;

				     /**
				      * Implementation of the same
				      * function in
				      * @ref{FiniteElement}.
				      */    
    virtual void
    fill_fe_face_values (const Mapping<dim>                   &mapping,
			 const typename DoFHandler<dim>::cell_iterator &cell,
			 const unsigned int                    face_no,
			 const Quadrature<dim-1>              &quadrature,
			 typename Mapping<dim>::InternalDataBase      &mapping_data,
			 typename Mapping<dim>::InternalDataBase      &fe_data,
			 FEValuesData<dim>                    &data) const ;

    				     /**
				      * Implementation of the same
				      * function in
				      * @ref{FiniteElement}.
				      */
    virtual void
    fill_fe_subface_values (const Mapping<dim>                   &mapping,
			    const typename DoFHandler<dim>::cell_iterator &cell,
			    const unsigned int                    face_no,
			    const unsigned int                    sub_no,
			    const Quadrature<dim-1>              &quadrature,
			    typename Mapping<dim>::InternalDataBase      &mapping_data,
			    typename Mapping<dim>::InternalDataBase      &fe_data,
			    FEValuesData<dim>                    &data) const ;
    

				     /**
				      * Calls (among other things)
				      * @p{fill_fe_([sub]face)_values}
				      * of the base elements. Calls
				      * @p{fill_fe_values} if
				      * @p{face_no==-1} and
				      * @p{sub_no==-1}; calls
				      * @p{fill_fe_face_values} if
				      * @p{face_no==-1} and
				      * @p{sub_no!=-1}; and calls
				      * @p{fill_fe_subface_values} if 
				      * @p{face_no!=-1} and
				      * @p{sub_no!=-1}.
				      */
    template <int dim_1>
    void compute_fill (const Mapping<dim>                   &mapping,
		       const typename DoFHandler<dim>::cell_iterator &cell,
		       const unsigned int                    face_no,
		       const unsigned int                    sub_no,
		       const Quadrature<dim_1>              &quadrature,
		       typename Mapping<dim>::InternalDataBase      &mapping_data,
		       typename Mapping<dim>::InternalDataBase      &fe_data,
		       FEValuesData<dim>                    &data) const ;

  private:

				     /**
				      * Value to indicate that a given
				      * face or subface number is
				      * invalid.
				      */
    static const unsigned int invalid_face_number = static_cast<unsigned int>(-1);
    
				     /**
				      * Pairs of multiplicity and
				      * element type.
				      */
    typedef typename std::pair<const FiniteElement<dim> *, unsigned int> ElementPair;
    
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
    typename std::vector<ElementPair> base_elements;


				     /**
				      * Helper function used in the constructor:
				      * take a @p{FiniteElementData} object
				      * and return an object of the same type
				      * with the number of degrees of
				      * freedom per vertex, line, etc.
				      * multiplied by @p{n}. Don't touch the
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
				      * Helper function used in the constructor:
				      * takes a @p{FiniteElement} object
				      * and returns an boolean vector including
				      * the @p{restriction_is_additive_flags} of
				      * the mixed element consisting of @p{N}
				      * elements of the sub-element @p{fe}.
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
				      * This function is simply
				      * singled out of the
				      * constructors since there are
				      * several of them. It sets up
				      * the index table for the system
				      * as well as @p{restriction} and
				      * @p{prolongation}
				      * matrices.
				      */
    void initialize();

				     /**
				      * Used by @p{initialize}.
				      */
    void build_cell_table();
    
				     /**
				      * Used by @p{initialize}.
				      */
    void build_face_table();

				     /**
				      * Used by @p{initialize}.
				      */
    void build_interface_constraints ();

				     /**
				      * Usually: Fields of
				      * cell-independent data.
				      *
				      * But for @p{FESystem} this
				      * @p{InternalData} class does
				      * not itself store the data but
				      * only pointers to
				      * @p{InternalDatas} of the base
				      * elements.
				      */
    class InternalData : public FiniteElementBase<dim>::InternalDataBase
    {
      public:
					 /**
					  * Constructor. Is called by
					  * the @p{get_data}
					  * function. Sets the size of
					  * the @p{base_fe_datas}
					  * vector to
					  * @p{n_base_elements}.
					  */
	InternalData(const unsigned int n_base_elements);
	
					 /**
					  * Destructor. Deletes all
					  * @p{InternalDatas} whose
					  * pointers are stored by the
					  * @p{base_fe_datas}
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
					  * @p{InternalData} of the
					  * @p{base_no}th base
					  * element.
					  */
	void set_fe_data(unsigned int base_no,
			 typename FiniteElementBase<dim>::InternalDataBase *);

					 /**
					  * Gives read-access to the
					  * pointer to a
					  * @p{InternalData} of the
					  * @p{base_no}th base element.
					  */	
	typename FiniteElementBase<dim>::InternalDataBase &get_fe_data(unsigned int base_no) const;


					 /**
					  * Gives write-access to the
					  * pointer to a
					  * @p{FEValuesData} for the
					  * @p{base_no}th base
					  * element.
					  */
	void set_fe_values_data(unsigned int base_no,
				FEValuesData<dim> *);

					 /**
					  * Gives read-access to the
					  * pointer to a
					  * @p{FEValuesData} for the
					  * @p{base_no}th base element.
					  */	
	FEValuesData<dim> &get_fe_values_data(unsigned int base_no) const;

					 /**
					  * Deletes the
					  * @p{FEValuesData} the
					  * @p{fe_datas[base_no]}
					  * pointer is pointing
					  * to. Sets
					  * @p{fe_datas[base_no]} to
					  * zero.
					  *
					  * This function is used to
					  * delete @p{FEValuesData}
					  * that are needed only on
					  * the first cell but not any
					  * more afterwards.  This is
					  * the case for
					  * e.g. Lagrangian elements
					  * (see e.g. @p{FE_Q}
					  * classes).
					  */
	void delete_fe_values_data(unsigned int base_no);
	
      private:
	
					 /**
					  * Pointers to the
					  * @p{InternalDatas} of the
					  * base elements. They are
					  * accessed to by the
					  * @p{set_} and
					  * @p{get_fe_data}
					  * functions.
					  *
					  * The size of this vector is
					  * set to @p{n_base_elements}
					  * by the InternalData
					  * constructor.  It is
					  * filled by the @p{get_data}
					  * function.
					  */
	typename std::vector<FiniteElementBase<dim>::InternalDataBase *> base_fe_datas;

					 /**
					  * Pointers to the
					  * @p{FEValuesDatas}
					  * that are given to the
					  * @p{fill_fe_values}
					  * function of the base
					  * elements. They are
					  * accessed to by the
					  * @p{set_} and
					  * @p{get_fe_values_data}
					  * functions.
					  *
					  * The size of this vector is
					  * set to @p{n_base_elements}
					  * by the InternalData
					  * constructor.
					  */
	typename std::vector<FEValuesData<dim> *> base_fe_values_datas;
    };
};


/* ------------------------- inline functions ------------------------- */

template<int dim>
inline unsigned int
FESystem<dim>::n_base_elements() const
{
  return base_elements.size();
};



template<int dim>
inline unsigned int
FESystem<dim>::element_multiplicity (const unsigned int index) const
{
  Assert (index < base_elements.size(), 
	  ExcIndexRange(index, 0, base_elements.size()));
  return base_elements[index].second;
};



template <int dim>
inline const FiniteElement<dim> &
FESystem<dim>::base_element (const unsigned int index) const
{
  Assert (index < base_elements.size(), 
	  ExcIndexRange(index, 0, base_elements.size()));
  return *base_elements[index].first;
};





/*----------------------------  fe_system.h  ---------------------------*/
#endif
/*----------------------------  fe_system.h  ---------------------------*/
