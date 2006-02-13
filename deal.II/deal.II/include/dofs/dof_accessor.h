//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__dof_accessor_h
#define __deal2__dof_accessor_h


#include <base/config.h>
#include <grid/tria_accessor.h>
#include <dofs/dof_handler.h>
#include <dofs/hp_dof_handler.h>
#include <vector>

template <typename number> class FullMatrix;
template <typename number> class SparseMatrix;
template <typename number> class Vector;

template <int dim> class DoFHandler;
namespace hp
{
  template <int dim> class DoFHandler;
}
template <int dim> class FiniteElement;

template <int celldim, class DH> class DoFObjectAccessor;

template <int dim, typename Accessor> class TriaRawIterator;



// note: the file dof_accessor.templates.h is included at the end of
// this file.  this includes a lot of templates and thus makes
// compilation slower, but at the same time allows for more aggressive
// inlining and thus faster code.


/**
 * Define the base class for accessors to the degrees of
 * freedom. Accessors are used to, well, access the data that pertains
 * to edges, faces, and cells of a triangulation. The concept is
 * explained in more detail in connection to @ref Iterators.
 *
 * The template argument of this class refers to the type of DoF
 * handler we should work on. It can either be ::DoFHandler<dim> or
 * hp::DoFHandler<dim>. The space dimension of the object we are to
 * work on is automatically extracted from the DH template argument.
 *
 * @ingroup dofs
 * @ingroup Accessors 
 * @author Wolfgang Bangerth, 1998
 */
template <class DH>
class DoFAccessor
{
  public:
				     /**
				      * Constructor
				      */
    DoFAccessor ();

				     /**
				      * This should be the default constructor.
				      * We cast away the constness of the
				      * pointer which clearly is EVIL but
				      * we can't help without making all
				      * functions which could somehow use
				      * iterators (directly or indirectly) make
				      * non-const, even if they preserve
				      * constness.
				      */
    DoFAccessor (const DH *dof_handler);

				     /**
				      * Reset the DoF handler pointer.
				      */
    void set_dof_handler (DH *dh);

				     /**
				      * Return a handle on the
				      * DoFHandler object which we
				      * are using.
				      */
    const DH &
    get_dof_handler () const;

				     /**
				      * Copy operator.
				      */
    DoFAccessor<DH> &
    operator = (const DoFAccessor<DH> &da);

				     /**
				      * Exception for child classes
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcInvalidObject);
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcVectorNotEmpty);
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcVectorDoesNotMatch);
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcMatrixDoesNotMatch);
				     /**
				      * A function has been called for
				      * a cell which should be active,
				      * but is refined. @ref GlossActive
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcNotActive);
				     /**
				      * Exception
				      * 
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcCantCompareIterators);

  protected:    
				     /**
				      *  Compare for equality.            
				      */
    bool operator == (const DoFAccessor &) const;
	
				     /**
				      * Compare for inequality.
				      */
    bool operator != (const DoFAccessor &) const;

  protected:
				     /**
				      * Store the address of the DoFHandler object
				      * to be accessed.
				      */
    DH *dof_handler;  
};


/* -------------------------------------------------------------------------- */

/**
 * This is a switch class which only declares a @p typedef. It is meant to
 * determine which class a DoFObjectAccessor class is to be derived
 * from. By default, <tt>DoFObjectAccessor<celldim,dim></tt> derives from
 * the @p typedef in the general <tt>DoFObjectAccessor_Inheritance<celldim,dim></tt>
 * class, which is <tt>TriaObjectAccessor<celldim,dim></tt>,
 * but if <tt>celldim==dim</tt>, then the specialization <tt>DoFObjectAccessor_Inheritance<dim,dim></tt>
 * is used which declares its local type to be <tt>CellAccessor<dim></tt>. Therefore,
 * the inheritance is automatically chosen to be from CellAccessor if the
 * object under consideration has full dimension, i.e. constitutes a cell.
 *
 * @ingroup dofs
 * @ingroup Accessors 
 * @author Wolfgang Bangerth, 1999
 */
template <int celldim, int dim>
class DoFObjectAccessor_Inheritance 
{
  public:
				     /**
				      * Declaration of the @p typedef.
				      * See the full documentation for
				      * more information.
				      */
    typedef TriaObjectAccessor<celldim,dim> BaseClass;
};


/**
 * This is a switch class which only declares a @p typedef. It is meant to
 * determine which class a DoFObjectAccessor class is to be derived
 * from. By default, <tt>DoFObjectAccessor<celldim,dim></tt> derives from
 * the @p typedef in the general <tt>DoFObjectAccessor_Inheritance<celldim,dim></tt>
 * class, which is <tt>TriaObjectAccessor<celldim,dim></tt>,
 * but if <tt>celldim==dim</tt>, then the specialization <tt>DoFObjectAccessor_Inheritance<dim,dim></tt>
 * is used which declares its local type to be <tt>CellAccessor<dim></tt>. Therefore,
 * the inheritance is automatically chosen to be from CellAccessor if the
 * object under consideration has full dimension, i.e. constitutes a cell.
 *
 * @ingroup dofs
 * @ingroup Accessors 
 * @author Wolfgang Bangerth, 1999
 */
template <int dim>
class DoFObjectAccessor_Inheritance<dim,dim>
{
  public:
				     /**
				      * Declaration of the @p typedef.
				      * See the full documentation for
				      * more information.
				      */
    typedef CellAccessor<dim> BaseClass;
};


/* -------------------------------------------------------------------------- */


/**
 * Common template for access to the data on a line, quad, hex. Note
 * that this class is only here for documentation purposes, the actual
 * implementation of functions is in classes with specialized template
 * parameters. In this class here, we only collect all functions which
 * are in the specialized classes for simpler reference. Some
 * functions, however, might be missing in the specialized classes,
 * such as @p quad in the accessors for lines and quads, etc.
 *
 * This class follows mainly the route laid out by the accessor library
 * declared in the triangulation library (TriaAccessor). It enables
 * the user to access the degrees of freedom on the lines (there are similar
 * versions for the DoFs on quads, etc), where the dimension of the underlying
 * triangulation does not really matter (i.e. this accessor works with the
 * lines in 1D-, 2D-, etc dimensions).
 *
 *
 * <h3>Usage</h3>
 *
 * The DoFDimensionInfo classes inherited by the DoFHandler classes
 * declare typedefs to iterators using the accessors declared in this class
 * hierarchy tree. Usage is best to happen through these typedefs, since they
 * are more secure to changes in the class naming and template interface as well
 * as they provide easier typing (much less complicated names!).
 * 
 * 
 * <h3>Notes about the class hierarchy structure</h3>
 *
 * See the report on this subject, which is available from the general
 * documentation directory.
 *
 * @ingroup dofs
 * @ingroup Accessors 
 * @author Wolfgang Bangerth, 1998; Guido Kanschat, 1999
 */
template <int celldim, class DH>
class DoFObjectAccessor : public DoFAccessor<DH>,
			  public TriaObjectAccessor<celldim, DH::dimension>
{
  public:
				     /**
				      * Extract dimension from DH.
				      */
    static const unsigned int dim = DH::dimension;
    
				     /**
				      * Data type  passed by the iterator class.
				      */
    typedef DH AccessorData;

				     /**
				      * Declare base class as a local typedef
				      * for simpler access.
				      */
    typedef TriaObjectAccessor<celldim,dim> BaseClass;


				     /**
				      * Default constructor. Unused, thus
				      * not implemented.
				      */
    DoFObjectAccessor ();

    				     /**
				      * Constructor. The @p local_data
				      * argument is assumed to be a pointer
				      * to a DoFHandler object.
				      */
    DoFObjectAccessor (const Triangulation<dim> *tria,
		       const int                 level,
		       const int                 index,
		       const AccessorData       *local_data);
    
				     /**
				      * Index of the <i>i</i>th degree
				      * of freedom of this object.
				      */
    unsigned int dof_index (const unsigned int i) const;

    				     /**
				      * Set the index of the <i>i</i>th degree
				      * of freedom of this object to @p index.
				      */
    void set_dof_index (const unsigned int i,
			const int index) const;

				     /**
				      * Index of the <i>i</i> degree
				      * on the @p vertexth vertex.
				      */
    unsigned int vertex_dof_index (const unsigned int vertex,
				   const unsigned int i) const;

				     /**
				      * Set the index of the <i>i</i> degree
				      * on the @p vertexth vertex to @p index.
				      */
    void set_vertex_dof_index (const unsigned int vertex,
			       const unsigned int i,
			       const unsigned int index) const;

    				     /**
				      * Return the indices of the dofs of this
				      * quad in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * The vector has to have the
				      * right size before being passed
				      * to this function.
				      */
    void get_dof_indices (std::vector<unsigned int> &dof_indices) const;

    				     /**
				      * Return the values of the given vector
				      * restricted to the dofs of this
				      * cell in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * The vector has to have the
				      * right size before being passed
				      * to this function. This
				      * function is only callable for
				      * active cells.
				      *
				      * The input vector may be either a
				      * <tt>Vector<float></tt>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc vector if deal.II is compiled to
				      * support PETSc. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class InputVector, typename number>
    void get_dof_values (const InputVector &values,
			 Vector<number>    &local_values) const;

				     /**
				      * This function is the counterpart to
				      * get_dof_values(): it takes a vector
				      * of values for the degrees of freedom
				      * of the cell pointed to by this iterator
				      * and writes these values into the global
				      * data vector @p values. This function
				      * is only callable for active cells.
				      *
				      * Note that for continuous finite
				      * elements, calling this function affects
				      * the dof values on neighboring cells as
				      * well. It may also violate continuity
				      * requirements for hanging nodes, if
				      * neighboring cells are less refined than
				      * the present one. These requirements
				      * are not taken care of and must be
				      * enforced by the user afterwards.
				      *
				      * The vector has to have the
				      * right size before being passed
				      * to this function.
				      *
				      * The output vector may be either a
				      * Vector<float>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc vector if deal.II is compiled to
				      * support PETSc. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class OutputVector, typename number>
    void set_dof_values (const Vector<number> &local_values,
			 OutputVector         &values) const;
                         

                                     /**
                                      * Return the finite element that
                                      * is used on the cell pointed to
                                      * by this iterator. For non-hp
                                      * DoF handlers, this is of
                                      * course always the same
                                      * element, independent of the
                                      * cell we are presently on, but
                                      * for hp DoF handlers, this may
                                      * change from cell to cell.
                                      */
    const FiniteElement<DoFObjectAccessor<celldim,DH>::dim> &
    get_fe () const;
    
                                     /**
				      *  Returns the index inside the
				      *  hp::FECollection of the FiniteElement
				      *  used for this cell.
				      */
    unsigned int active_fe_index () const;

                                     /**
				      *  Sets the index of the FiniteElement used for
				      *  this cell.
				      */
    void set_active_fe_index (const unsigned int i);


    				     /**
				      *  Pointer to the <i>i</i>th line
				      *  bounding this Object.
				      */
    TriaIterator<DoFObjectAccessor<celldim,DH>::dim,DoFObjectAccessor<1, DH> >
    line (const unsigned int i) const;

    				     /**
				      *  Pointer to the <i>i</i>th quad
				      *  bounding this Object.
				      */
    TriaIterator<DoFObjectAccessor<celldim,DH>::dim,DoFObjectAccessor<2, DH> >
    quad (const unsigned int i) const;
    
				     /**
				      * <i>i</i>th child as a DoFObjectAccessor
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a hex accessor without
				      * access to the DoF data.
				      */
    TriaIterator<DoFObjectAccessor<celldim,DH>::dim,DoFObjectAccessor<celldim, DH> >
    child (const unsigned int) const;

				     /**
				      * Distribute a local (cell
				      * based) vector to a global one
				      * by mapping the local numbering
				      * of the degrees of freedom to
				      * the global one and entering
				      * the local values into the
				      * global vector.
				      *
				      * The elements are
				      * <em>added</em> up to the
				      * elements in the global vector,
				      * rather than just set, since
				      * this is usually what one
				      * wants.
				      */
    template <typename number, typename OutputVector>
    void
    distribute_local_to_global (const Vector<number> &local_source,
                                OutputVector         &global_destination) const;

				     /**
				      * This function does much the same as
				      * the
				      * <tt>distribute_local_to_global(Vector,Vector)</tt>
				      * function, but operates on matrices
				      * instead of vectors. The sparse matrix
				      * is supposed to have non-zero entry
				      * slots where required.
				      */
    template <typename number, typename OutputMatrix>
    void
    distribute_local_to_global (const FullMatrix<number> &local_source,
                                OutputMatrix             &global_destination) const;
    
				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const DoFObjectAccessor<celldim, DH> &a);

  protected:    
				     /**
				      *  Compare for equality.            
				      */
    bool operator == (const DoFObjectAccessor &) const;
	
				     /**
				      * Compare for inequality.
				      */
    bool operator != (const DoFObjectAccessor &) const;


                                     /**
                                      * Iterator classes need to be friends
                                      * because they need to access operator==
                                      * and operator!=.
                                      */
    template <int, typename> friend class TriaRawIterator;
};


/**
 * Closure class.
 * @ingroup dofs
 */
template <class DH>
class DoFObjectAccessor<0, DH> : public DoFAccessor<DH>,
				 public DoFObjectAccessor_Inheritance<0,DH::dimension>::BaseClass
{
  public:
    				     /**
				      * Extract dimension from DH.
				      */
    static const unsigned int dim = DH::dimension;

    typedef void AccessorData;

				     /**
				      * Constructor. Should never be called
				      * and thus throws an error.
				      */
    DoFObjectAccessor (Triangulation<dim> *,
		       const int,
		       const int,
		       const AccessorData *)
      {
	Assert (false, ExcInternalError());
      }
};


/**
 * Access to the degrees of freedom located on lines.
 * This class follows mainly the route laid out by the accessor library
 * declared in the triangulation library (TriaAccessor). It enables
 * the user to access the degrees of freedom on the lines (there are similar
 * versions for the DoFs on quads, etc), where the dimension of the underlying
 * triangulation does not really matter (i.e. this accessor works with the
 * lines in 1D-, 2D-, etc dimensions).
 *
 *
 * <h3>Usage</h3>
 *
 * The DoFDimensionInfo classes inherited by the DoFHandler classes
 * declare typedefs to iterators using the accessors declared in this class
 * hierarchy tree. Usage is best to happens through these typedefs, since they
 * are more secure to changes in the class naming and template interface as well
 * as they provide easier typing (much less complicated names!).
 * 
 * 
 * <h3>Notes about the class hierarchy structure</h3>
 *
 * Inheritance from <tt>DoFObjectAccessor_Inheritance<1,dim>::BaseClass</tt> yields
 * inheritance from <tt>CellAccessor<1></tt> if <tt>dim==1</tt> and from
 * <tt>TriaObjectAccessor<1,dim></tt> for all other @p dim values. Thus, an object
 * of this class shares all features of cells in one dimension, but behaves
 * like an ordinary line in all other cases.
 *
 * @ingroup dofs
 * @author Wolfgang Bangerth, 1998
 */
template <class DH>
class DoFObjectAccessor<1, DH> :
  public DoFAccessor<DH>,
  public DoFObjectAccessor_Inheritance<1,DH::dimension>::BaseClass
{
  public:
				     /**
				      * Extract dimension from DH.
				      */
    static const unsigned int dim = DH::dimension;
    
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef DH AccessorData;

				     /**
				      * Declare base class as a local typedef
				      * for simpler access.
				      */
    typedef typename DoFObjectAccessor_Inheritance<1,dim>::BaseClass BaseClass;

				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    DoFObjectAccessor ();
    
    				     /**
				      * Constructor. The @p local_data
				      * argument is assumed to be a pointer
				      * to a DoFHandler object.
				      */
    DoFObjectAccessor (const Triangulation<dim> *tria,
		       const int                 level,
		       const int                 index,
		       const AccessorData       *local_data);
    
				     /**
				      * Return the index of the @p ith degree
				      * of freedom of this line.
				      */
    unsigned int dof_index (const unsigned int i) const;

    				     /**
				      * Set the index of the @p ith degree
				      * of freedom of this line to @p index.
				      */
    void set_dof_index (const unsigned int i,
			const unsigned int index) const;

				     /**
				      * Return the index of the @p ith degree
				      * on the @p vertexth vertex.
				      */
    unsigned int vertex_dof_index (const unsigned int vertex,
				   const unsigned int i) const;

				     /**
				      * Set the index of the @p ith degree
				      * on the @p vertexth vertex to @p index.
				      */
    void set_vertex_dof_index (const unsigned int vertex,
			       const unsigned int i,
			       const unsigned int index) const;

    				     /**
				      * Return the indices of the dofs of this
				      * line in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, 
				      * dofs on line.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand.
				      */
    void get_dof_indices (std::vector<unsigned int> &dof_indices) const;

    				     /**
				      * Return the values of the given vector
				      * restricted to the dofs of this
				      * cell in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand. This
				      * function is only callable for active
				      * cells.
				      *
				      * The input vector may be either a
				      * Vector<float>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc vector if deal.II is compiled to
				      * support PETSc. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class InputVector, typename number>
    void get_dof_values (const InputVector &values,
			 Vector<number>    &local_values) const;

				     /**
				      * This function is the counterpart to
				      * get_dof_values(): it takes a vector
				      * of values for the degrees of freedom
				      * of the cell pointed to by this iterator
				      * and writes these values into the global
				      * data vector @p values. This function
				      * is only callable for active cells.
				      *
				      * Note that for continuous finite
				      * elements, calling this function affects
				      * the dof values on neighboring cells as
				      * well. It may also violate continuity
				      * requirements for hanging nodes, if
				      * neighboring cells are less refined than
				      * the present one. These requirements
				      * are not taken care of and must be
				      * enforced by the user afterwards.
				      *
				      * It is assumed that both vectors already
				      * have the right size beforehand.
				      *
				      * The output vector may be either a
				      * Vector<float>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc vector if deal.II is compiled to
				      * support PETSc. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class OutputVector, typename number>
    void set_dof_values (const Vector<number> &local_values,
			 OutputVector         &values) const;


                                     /**
                                      * Return the finite element that
                                      * is used on the cell pointed to
                                      * by this iterator. For non-hp
                                      * DoF handlers, this is of
                                      * course always the same
                                      * element, independent of the
                                      * cell we are presently on, but
                                      * for hp DoF handlers, this may
                                      * change from cell to cell.
                                      */
    const FiniteElement<DoFObjectAccessor<1,DH>::dim> &
    get_fe () const;
    
                                     /**
				      *  Returns the index inside the
				      *  hp::FECollection of the FiniteElement
				      *  used for this cell.
				      */
    unsigned int active_fe_index () const;

                                     /**
				      *  Sets the index of the FiniteElement used for
				      *  this cell.
				      */
    void set_active_fe_index (const unsigned int i);

				     /**
				      * Return the @p ith child as a DoF line
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a line accessor without
				      * access to the DoF data.
				      */
    TriaIterator<DoFObjectAccessor<1,DH>::dim,DoFObjectAccessor<1,DH> >
    child (const unsigned int) const;

				     /**
				      * Distribute a local (cell based) vector
				      * to a global one by mapping the local
				      * numbering of the degrees of freedom to
				      * the global one and entering the local
				      * values into the global vector.
				      *
				      * The elements are
				      * <em>added</em> up to the
				      * elements in the global vector,
				      * rather than just set, since
				      * this is usually what one
				      * wants.
				      */
    template <typename number, typename OutputVector>
    void
    distribute_local_to_global (const Vector<number> &local_source,
                                OutputVector         &global_destination) const;

				     /**
				      * This function does much the same as
				      * the
				      * <tt>distribute_local_to_global(Vector,Vector)</tt>
				      * function, but operates on matrices
				      * instead of vectors. The sparse matrix
				      * is supposed to have non-zero entry
				      * slots where required.
				      */
    template <typename number, typename OutputMatrix>
    void
    distribute_local_to_global (const FullMatrix<number> &local_source,
                                OutputMatrix             &global_destination) const;
    
				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const DoFObjectAccessor<1,DH> &a);

  protected:    
				     /**
				      *  Compare for equality.            
				      */
    bool operator == (const DoFObjectAccessor<1,DH> &) const;
	
				     /**
				      * Compare for inequality.
				      */
    bool operator != (const DoFObjectAccessor<1,DH> &) const;


                                     /**
                                      * Iterator classes need to be friends
                                      * because they need to access operator==
                                      * and operator!=.
                                      */
    template <int, typename> friend class TriaRawIterator;
};


/**
 * Grant access to the degrees of freedom located on quads.
 *
 * @ingroup dofs
 */
template <class DH>
class DoFObjectAccessor<2, DH> :
  public DoFAccessor<DH>,
  public DoFObjectAccessor_Inheritance<2,DH::dimension>::BaseClass
{
  public:
				     /**
				      * Extract dimension from DH.
				      */
    static const unsigned int dim = DH::dimension;

				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef DH AccessorData;

				     /**
				      * Declare base class as a local typedef
				      * for simpler access.
				      */
    typedef typename DoFObjectAccessor_Inheritance<2,dim>::BaseClass BaseClass;
    
				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    DoFObjectAccessor ();

    				     /**
				      * Constructor. The @p local_data
				      * argument is assumed to be a pointer
				      * to a DoFHandler object.
				      */
    DoFObjectAccessor (const Triangulation<dim> *tria,
		       const int                 level,
		       const int                 index,
		       const AccessorData       *local_data);
    
				     /**
				      * Return the index of the @p ith degree
				      * of freedom of this quad.
				      */
    unsigned int dof_index (const unsigned int i) const;

    				     /**
				      * Set the index of the @p ith degree
				      * of freedom of this quad to @p index.
				      */
    void set_dof_index (const unsigned int i,
			const unsigned int index) const;

				     /**
				      * Return the index of the @p ith degree
				      * on the @p vertexth vertex.
				      */
    unsigned int vertex_dof_index (const unsigned int vertex,
				   const unsigned int i) const;

				     /**
				      * Set the index of the @p ith degree
				      * on the @p vertexth vertex to @p index.
				      */
    void set_vertex_dof_index (const unsigned int vertex,
			       const unsigned int i,
			       const unsigned int index) const;

    				     /**
				      * Return the indices of the dofs of this
				      * quad in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand.
				      */
    void get_dof_indices (std::vector<unsigned int> &dof_indices) const;

    				     /**
				      * Return the values of the given vector
				      * restricted to the dofs of this
				      * cell in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand. This
				      * function is only callable for active
				      * cells.
				      *
				      * The input vector may be either a
				      * Vector<float>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc vector if deal.II is compiled to
				      * support PETSc. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class InputVector, typename number>
    void get_dof_values (const InputVector &values,
			 Vector<number>    &local_values) const;

				     /**
				      * This function is the counterpart to
				      * get_dof_values(): it takes a vector
				      * of values for the degrees of freedom
				      * of the cell pointed to by this iterator
				      * and writes these values into the global
				      * data vector @p values. This function
				      * is only callable for active cells.
				      *
				      * Note that for continuous finite
				      * elements, calling this function affects
				      * the dof values on neighboring cells as
				      * well. It may also violate continuity
				      * requirements for hanging nodes, if
				      * neighboring cells are less refined than
				      * the present one. These requirements
				      * are not taken care of and must be
				      * enforced by the user afterwards.
				      *
				      * It is assumed that both vectors already
				      * have the right size beforehand.
				      *
				      * The output vector may be either a
				      * Vector<float>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc vector if deal.II is compiled to
				      * support PETSc. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class OutputVector, typename number>
    void set_dof_values (const Vector<number> &local_values,
			 OutputVector         &values) const;


                                     /**
                                      * Return the finite element that
                                      * is used on the cell pointed to
                                      * by this iterator. For non-hp
                                      * DoF handlers, this is of
                                      * course always the same
                                      * element, independent of the
                                      * cell we are presently on, but
                                      * for hp DoF handlers, this may
                                      * change from cell to cell.
                                      */
    const FiniteElement<DoFObjectAccessor<2,DH>::dim> &
    get_fe () const;
    
                                     /**
				      *  Returns the index inside the
				      *  hp::FECollection of the FiniteElement
				      *  used for this cell.
				      */
    unsigned int active_fe_index () const;

                                     /**
				      *  Sets the index of the FiniteElement used for
				      *  this cell.
				      */
    void set_active_fe_index (const unsigned int i);

    				     /**
				      *  Return a pointer to the @p ith line
				      *  bounding this @p Quad.
				      */
    TriaIterator<DoFObjectAccessor<2,DH>::dim,DoFObjectAccessor<1, DH> >
    line (const unsigned int i) const;
    
				     /**
				      * Return the @p ith child as a DoF quad
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a quad accessor without
				      * access to the DoF data.
				      */
    TriaIterator<DoFObjectAccessor<2,DH>::dim,DoFObjectAccessor<2, DH> >
    child (const unsigned int) const;

				     /**
				      * Distribute a local (cell based) vector
				      * to a global one by mapping the local
				      * numbering of the degrees of freedom to
				      * the global one and entering the local
				      * values into the global vector.
				      *
				      * The elements are
				      * <em>added</em> up to the
				      * elements in the global vector,
				      * rather than just set, since
				      * this is usually what one
				      * wants.
				      */
    template <typename number, typename OutputVector>
    void
    distribute_local_to_global (const Vector<number> &local_source,
                                OutputVector         &global_destination) const;

				     /**
				      * This function does much the same as
				      * the
				      * <tt>distribute_local_to_global(Vector,Vector)</tt>
				      * function, but operates on matrices
				      * instead of vectors. The sparse matrix
				      * is supposed to have non-zero entry
				      * slots where required.
				      */
    template <typename number, typename OutputMatrix>
    void
    distribute_local_to_global (const FullMatrix<number> &local_source,
                                OutputMatrix             &global_destination) const;
    
				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const DoFObjectAccessor<2, DH> &a);

  protected:    
				     /**
				      *  Compare for equality.            
				      */
    bool operator == (const DoFObjectAccessor<2,DH> &) const;
	
				     /**
				      * Compare for inequality.
				      */
    bool operator != (const DoFObjectAccessor<2,DH> &) const;


                                     /**
                                      * Iterator classes need to be friends
                                      * because they need to access operator==
                                      * and operator!=.
                                      */
    template <int, typename> friend class TriaRawIterator;
};


/**
 * Grant access to the degrees of freedom located on hexes.
 *
 * @ingroup dofs
 */
template <class DH>
class DoFObjectAccessor<3, DH> :
  public DoFAccessor<DH>,
  public DoFObjectAccessor_Inheritance<3,DH::dimension>::BaseClass
{
  public:
				     /**
				      * Extract dimension from DH.
				      */
    static const unsigned int dim = DH::dimension;
    
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef DH AccessorData;

				     /**
				      * Declare base class as a local typedef
				      * for simpler access.
				      */
    typedef typename DoFObjectAccessor_Inheritance<3,dim>::BaseClass BaseClass;
    
				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    DoFObjectAccessor ();

    				     /**
				      * Constructor. The @p local_data
				      * argument is assumed to be a pointer
				      * to a DoFHandler object.
				      */
    DoFObjectAccessor (const Triangulation<dim> *tria,
		       const int                 level,
		       const int                 index,
		       const AccessorData       *local_data);
    
				     /**
				      * Return the index of the @p ith degree
				      * of freedom of this hex.
				      */
    unsigned int dof_index (const unsigned int i) const;

    				     /**
				      * Set the index of the @p ith degree
				      * of freedom of this hex to @p index.
				      */
    void set_dof_index (const unsigned int i,
			const unsigned int index) const;

				     /**
				      * Return the index of the @p ith degree
				      * on the @p vertexth vertex.
				      */
    unsigned int vertex_dof_index (const unsigned int vertex,
				   const unsigned int i) const;

				     /**
				      * Set the index of the @p ith degree
				      * on the @p vertexth vertex to @p index.
				      */
    void set_vertex_dof_index (const unsigned int vertex,
			       const unsigned int i,
			       const unsigned int index) const;

    				     /**
				      * Return the indices of the dofs of this
				      * hex in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand.
				      */
    void get_dof_indices (std::vector<unsigned int> &dof_indices) const;

    				     /**
				      * Return the values of the given vector
				      * restricted to the dofs of this
				      * cell in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand. This
				      * function is only callable for active
				      * cells.
				      *
				      * The input vector may be either a
				      * Vector<float>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc vector if deal.II is compiled to
				      * support PETSc. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class InputVector, typename number>
    void get_dof_values (const InputVector &values,
			 Vector<number>    &local_values) const;

				     /**
				      * This function is the counterpart to
				      * get_dof_values(): it takes a vector
				      * of values for the degrees of freedom
				      * of the cell pointed to by this iterator
				      * and writes these values into the global
				      * data vector @p values. This function
				      * is only callable for active cells.
				      *
				      * Note that for continuous finite
				      * elements, calling this function affects
				      * the dof values on neighboring cells as
				      * well. It may also violate continuity
				      * requirements for hanging nodes, if
				      * neighboring cells are less refined than
				      * the present one. These requirements
				      * are not taken care of and must be
				      * enforced by the user afterwards.
				      *
				      * It is assumed that both vectors already
				      * have the right size beforehand.
				      *
				      * The output vector may be either a
				      * Vector<float>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc vector if deal.II is compiled to
				      * support PETSc. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class OutputVector, typename number>
    void set_dof_values (const Vector<number> &local_values,
			 OutputVector         &values) const;


                                     /**
                                      * Return the finite element that
                                      * is used on the cell pointed to
                                      * by this iterator. For non-hp
                                      * DoF handlers, this is of
                                      * course always the same
                                      * element, independent of the
                                      * cell we are presently on, but
                                      * for hp DoF handlers, this may
                                      * change from cell to cell.
                                      */
    const FiniteElement<DoFObjectAccessor<3, DH>::dim> &
    get_fe () const;

                                     /**
				      *  Returns the index inside the
				      *  hp::FECollection of the FiniteElement
				      *  used for this cell.
				      */
    unsigned int active_fe_index () const;

                                     /**
				      *  Sets the index of the FiniteElement used for
				      *  this cell.
				      */
    void set_active_fe_index (const unsigned int i);

    				     /**
				      *  Return a pointer to the @p ith line
				      *  bounding this @p Hex.
				      */
    TriaIterator<DoFObjectAccessor<3, DH>::dim,DoFObjectAccessor<1, DH> >
    line (const unsigned int i) const;

    				     /**
				      *  Return a pointer to the @p ith quad
				      *  bounding this @p Hex.
				      */
    TriaIterator<DoFObjectAccessor<3, DH>::dim,DoFObjectAccessor<2, DH> >
    quad (const unsigned int i) const;
    
				     /**
				      * Return the @p ith child as a DoF hex
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a hex accessor without
				      * access to the DoF data.
				      */
    TriaIterator<DoFObjectAccessor<3, DH>::dim,DoFObjectAccessor<3, DH> >
    child (const unsigned int) const;

				     /**
				      * Distribute a local (cell
				      * based) vector to a global one
				      * by mapping the local numbering
				      * of the degrees of freedom to
				      * the global one and entering
				      * the local values into the
				      * global vector.
				      *
				      * The elements are
				      * <em>added</em> up to the
				      * elements in the global vector,
				      * rather than just set, since
				      * this is usually what one
				      * wants.
				      */
    template <typename number, typename OutputVector>
    void
    distribute_local_to_global (const Vector<number> &local_source,
                                OutputVector         &global_destination) const;

				     /**
				      * This function does much the same as
				      * the
				      * <tt>distribute_local_to_global(Vector,Vector)</tt>
				      * function, but operates on matrices
				      * instead of vectors. The sparse matrix
				      * is supposed to have non-zero entry
				      * slots where required.
				      */
    template <typename number, typename OutputMatrix>
    void
    distribute_local_to_global (const FullMatrix<number> &local_source,
                                OutputMatrix             &global_destination) const;
    
				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const DoFObjectAccessor<3, DH> &a);

  protected:    
				     /**
				      *  Compare for equality.            
				      */
    bool operator == (const DoFObjectAccessor<3,DH> &) const;
	
				     /**
				      * Compare for inequality.
				      */
    bool operator != (const DoFObjectAccessor<3,DH> &) const;


                                     /**
                                      * Iterator classes need to be friends
                                      * because they need to access operator==
                                      * and operator!=.
                                      */
    template <int, typename> friend class TriaRawIterator;
};



/**
 * Grant access to the degrees of freedom on a cell. In fact, since all
 * access to the degrees of freedom has been enabled by the classes
 * <tt>DoFObjectAccessor<1, 1></tt> and <tt>DoFObjectAccessor<2, 2></tt> for the space dimension
 * one and two, respectively, this class only collects the pieces
 * together by deriving from the appropriate <tt>DoF*Accessor</tt> and the
 * right <tt>CellAccessor<dim></tt> and finally adding two functions which give
 * access to the neighbors and children as DoFCellAccessor objects
 * rather than CellAccessor objects (the latter function was inherited
 * from the <tt>CellAccessor<dim></tt> class).
 *
 * Note that since for the class we derive from, i.e. <tt>DoFObjectAccessor<dim,dim></tt>,
 * the two template parameters are equal, the base class is actually derived from
 * CellAccessor, which makes the functions of this class available to the
 * DoFCellAccessor class as well.
 *
 * @ingroup dofs
 * @author Wolfgang Bangerth, 1998
 */
template <class DH>
class DoFCellAccessor :  public DoFObjectAccessor<DH::dimension,DH>
{
  public:
				     /**
				      * Extract dimension from DH.
				      */
    static const unsigned int dim = DH::dimension;
    
				     /**
				      * Declare the data type that
				      * this accessor class expects to
				      * get passed from the iterator
				      * classes.
				      */
    typedef typename DoFObjectAccessor<dim,DH>::AccessorData AccessorData;
    
    				     /**
				      * Constructor
				      */
    DoFCellAccessor (const Triangulation<dim> *tria,
		     const int                 level,
		     const int                 index,
		     const AccessorData       *local_data);

				     /**
				      * Return the @p ith neighbor as
				      * a DoF cell iterator. This
				      * function is needed since the
				      * neighbor function of the base
				      * class returns a cell accessor
				      * without access to the DoF
				      * data.
				      */
    TriaIterator<dim,DoFCellAccessor<DH> >
    neighbor (const unsigned int) const;

    				     /**
				      * Return the @p ith child as a
				      * DoF cell iterator. This
				      * function is needed since the
				      * child function of the base
				      * class returns a cell accessor
				      * without access to the DoF
				      * data.
				      */
    TriaIterator<dim,DoFCellAccessor<DH> >
    child (const unsigned int) const;

    				     /**
				      * Return an iterator to the @p ith face
				      * of this cell.
				      *
				      * This function is not implemented in
				      * 1D, and maps to DoFObjectAccessor<2,
				      * dim>::line in 2D.
				      */
    TriaIterator<dim, DoFObjectAccessor<dim-1, DH> >
    face (const unsigned int i) const;

                                     /**
				      * Return the result of the
				      * @p neighbor_child_on_subface
				      * function of the base class,
				      * but convert it so that one can
				      * also access the DoF data (the
				      * function in the base class
				      * only returns an iterator with
				      * access to the triangulation
				      * data).
				      */
    TriaIterator<dim,DoFCellAccessor<DH> >
    neighbor_child_on_subface (const unsigned int face_no,
                               const unsigned int subface_no) const;

                                     /**
				      * Return the interpolation of
				      * the given finite element
				      * function to the present
				      * cell. In the simplest case,
				      * the cell is a terminal one,
				      * i.e. has no children; then,
				      * the returned value is the
				      * vector of nodal values on that
				      * cell. You could then as well
				      * get the desired values through
				      * the @p get_dof_values
				      * function In the other case,
				      * when the cell has children, we
				      * use the restriction matrices
				      * provided by the finite element
				      * class to compute the
				      * interpolation from the
				      * children to the present cell.
				      *
				      * It is assumed that both
				      * vectors already have the right
				      * size beforehand. This function
				      * assumes the existence of an
				      * interpolation from child cells
				      * to the mother cell, denoted by
				      * the restriction matrices of
				      * the finite element class.
				      * Futhermore, this interpolation
				      * should be possible for each
				      * child alone, i.e.  it should
				      * be possible to compute the
				      * restriction by writing values
			              * obtained from each child
				      * directly into the output
				      * vector, without, for example,
				      * computing an average over all
				      * children.  These properties,
				      * however, do not exist for all
				      * elements; an example is the
				      * DG(0) element, which
				      * represents piecewise constant
				      * elements: for these the
				      * restriction to mother cell
				      * could be the average of the
				      * values of the children, maybe
				      * weighted by the measure of
				      * each child. It is not yet
				      * decided what the this function
				      * does in these cases.
				      *
				      * Unlike the get_dof_values()
				      * function, this function is
				      * associated to cells rather
				      * than to lines, quads, and
				      * hexes, since interpolation is
				      * presently only provided for
				      * cells by the finite element
				      * objects.
				      */
    template <class InputVector, typename number>
    void get_interpolated_dof_values (const InputVector &values,
				      Vector<number>    &interpolated_values) const;

				     /**
				      * This, again, is the
				      * counterpart to
				      * get_interpolated_dof_values():
				      * you specify the dof values on
				      * a cell and these are
				      * interpolated to the children
				      * of the present cell and set on
				      * the terminal cells.
				      *
				      * In principle, it works as
				      * follows: if the cell pointed
				      * to by this object is terminal,
				      * then the dof values are set in
				      * the global data vector by
				      * calling the set_dof_values()
				      * function; otherwise, the
				      * values are prolonged to each
				      * of the children and this
				      * function is called for each of
				      * them.
				      *
				      * Using the
				      * get_interpolated_dof_values()
				      * and this function, you can
				      * compute the interpolation of a
				      * finite element function to a
				      * coarser grid by first getting
				      * the interpolated solution on a
				      * cell of the coarse grid and
				      * afterwards redistributing it
				      * using this function.
				      *
				      * Note that for continuous
				      * finite elements, calling this
				      * function affects the dof
				      * values on neighboring cells as
				      * well. It may also violate
				      * continuity requirements for
				      * hanging nodes, if neighboring
				      * cells are less refined than
				      * the present one, or if their
				      * children are less refined than
				      * the children of this
				      * cell. These requirements are
				      * not taken care of and must be
				      * enforced by the user
				      * afterwards.
				      *
				      * It is assumed that both
				      * vectors already have the right
				      * size beforehand. This function
				      * relies on the existence of a
				      * natural interpolation property
				      * of finite element spaces of a
				      * cell to its children, denoted
				      * by the prolongation matrices
				      * of finite element classes. For
				      * some elements, the spaces on
				      * coarse and fine grids are not
				      * nested, in which case the
				      * interpolation to a child is
				      * not the identity; refer to the
				      * documentation of the
				      * respective finite element
				      * class for a description of
				      * what the prolongation matrices
				      * represent in this case.
				      *
				      * Unlike the set_dof_values()
				      * function, this function is
				      * associated to cells rather
				      * than to lines, quads, and
				      * hexes, since interpolation is
				      * presently only provided for
				      * cells by the finite element
				      * objects.
				      *
				      * The output vector may be either a
				      * Vector<float>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc vector if deal.II is compiled to
				      * support PETSc. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class OutputVector, typename number>
    void set_dof_values_by_interpolation (const Vector<number> &local_values,
					  OutputVector         &values) const;
    
};


/* -------------- declaration of explicit specializations ------------- */

template <>
TriaIterator<1, DoFObjectAccessor<0,DoFHandler<1> > >
DoFCellAccessor<DoFHandler<1> >::face (const unsigned int) const;
template <>
TriaIterator<2, DoFObjectAccessor<1,DoFHandler<2> > >
DoFCellAccessor<DoFHandler<2> >::face (const unsigned int i) const;
template <>
TriaIterator<3, DoFObjectAccessor<2,DoFHandler<3> > >
DoFCellAccessor<DoFHandler<3> >::face (const unsigned int i) const;

template <>
const FiniteElement<1> &
DoFObjectAccessor<1,hp::DoFHandler<1> >::get_fe () const;
template <>
const FiniteElement<2> &
DoFObjectAccessor<1,hp::DoFHandler<2> >::get_fe () const;
template <>
const FiniteElement<3> &
DoFObjectAccessor<1,hp::DoFHandler<3> >::get_fe () const;
template <>
const FiniteElement<2> &
DoFObjectAccessor<2,hp::DoFHandler<2> >::get_fe () const;
template <>
const FiniteElement<3> &
DoFObjectAccessor<2,hp::DoFHandler<3> >::get_fe () const;
template <>
const FiniteElement<3> &
DoFObjectAccessor<3,hp::DoFHandler<3> >::get_fe () const;

template <>
void DoFObjectAccessor<1,hp::DoFHandler<1> >::set_dof_index (const unsigned int,
							     const unsigned int) const;
template <>
void DoFObjectAccessor<1,hp::DoFHandler<2> >::set_dof_index (const unsigned int,
							     const unsigned int) const;
template <>
void DoFObjectAccessor<1,hp::DoFHandler<3> >::set_dof_index (const unsigned int,
							     const unsigned int) const;
template <>
void DoFObjectAccessor<2,hp::DoFHandler<2> >::set_dof_index (const unsigned int,
							     const unsigned int) const;
template <>
void DoFObjectAccessor<2,hp::DoFHandler<3> >::set_dof_index (const unsigned int,
							     const unsigned int) const;
template <>
void DoFObjectAccessor<3,hp::DoFHandler<3> >::set_dof_index (const unsigned int,
							     const unsigned int) const;


template <>
unsigned int
DoFObjectAccessor<1,hp::DoFHandler<1> >::active_fe_index () const;
template <>
unsigned int
DoFObjectAccessor<1,hp::DoFHandler<2> >::active_fe_index () const;
template <>
unsigned int
DoFObjectAccessor<1,hp::DoFHandler<3> >::active_fe_index () const;
template <>
unsigned int
DoFObjectAccessor<2,hp::DoFHandler<2> >::active_fe_index () const;
template <>
unsigned int
DoFObjectAccessor<2,hp::DoFHandler<3> >::active_fe_index () const;
template <>
unsigned int
DoFObjectAccessor<3,hp::DoFHandler<3> >::active_fe_index () const;

template <>
void
DoFObjectAccessor<1,hp::DoFHandler<1> >::set_active_fe_index (const unsigned int i);
template <>
void
DoFObjectAccessor<1,hp::DoFHandler<2> >::set_active_fe_index (const unsigned int i);
template <>
void
DoFObjectAccessor<1,hp::DoFHandler<3> >::set_active_fe_index (const unsigned int i);
template <>
void
DoFObjectAccessor<2,hp::DoFHandler<2> >::set_active_fe_index (const unsigned int i);
template <>
void
DoFObjectAccessor<2,hp::DoFHandler<3> >::set_active_fe_index (const unsigned int i);
template <>
void
DoFObjectAccessor<3,hp::DoFHandler<3> >::set_active_fe_index (const unsigned int i);


// include more templates
#include "dof_accessor.templates.h"


/*----------------------------   dof_iterator.h     ---------------------------*/

#endif
/*----------------------------   dof_iterator.h     ---------------------------*/
