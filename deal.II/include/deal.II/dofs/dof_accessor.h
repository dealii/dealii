//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__dof_accessor_h
#define __deal2__dof_accessor_h


#include <deal.II/base/config.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

template <typename number> class FullMatrix;
template <typename number> class SparseMatrix;
template <typename number> class Vector;
class ConstraintMatrix;

template <typename Accessor> class TriaRawIterator;

namespace internal
{
  namespace DoFCellAccessor
  {
    struct Implementation;
  }

  namespace DoFHandler
  {
    struct Implementation;
    namespace Policy
    {
      struct Implementation;
    }
  }

  namespace hp
  {
    namespace DoFHandler
    {
      struct Implementation;
    }
  }
}

// note: the file dof_accessor.templates.h is included at the end of
// this file.  this includes a lot of templates and thus makes
// compilation slower, but at the same time allows for more aggressive
// inlining and thus faster code.


namespace internal
{
  namespace DoFAccessor
  {
/**
 * This is a switch class which only declares a @p typedef. It is meant to
 * determine which class a DoFAccessor class is to be derived from. By
 * default, <tt>DoFAccessor@<structdim,dim,spacedim@></tt> derives from the
 * typedef in the general
 * <tt>Inheritance@<structdim,dim,spacedim@></tt> class, which is
 * <tt>TriaAccessor@<structdim,dim,spacedim@></tt>, but if
 * <tt>structdim==dim</tt>, then the specialization
 * <tt>Inheritance@<dim,dim,spacedim@></tt> is used which declares
 * its local type to be <tt>CellAccessor@<dim,spacedim@></tt>. Therefore, the
 * inheritance is automatically chosen to be from CellAccessor if the object
 * under consideration has full dimension, i.e. constitutes a cell.
 *
 * @ingroup dofs
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1999
 */
    template <int structdim, int dim, int spacedim>
    struct Inheritance
    {
					 /**
					  * Declaration of the @p typedef.
					  * See the full documentation for
					  * more information.
					  */
	typedef dealii::TriaAccessor<structdim,dim,spacedim> BaseClass;
    };


/**
 * This is the specialization of the general template used for the case where
 * an object has full dimension, i.e. is a cell. See the general template for
 * more details.
 */
    template <int dim, int spacedim>
    struct Inheritance<dim,dim,spacedim>
    {
					 /**
					  * Declaration of the @p typedef.
					  * See the full documentation for
					  * more information.
					  */
	typedef dealii::CellAccessor<dim,spacedim> BaseClass;
    };
  }
}


/* -------------------------------------------------------------------------- */



/**
 * A class that gives access to the degrees of freedom stored in a DoFHandler
 * or hp::DoFHandler object. Accessors are used to, well, access the data that
 * pertains to edges, faces, and cells of a triangulation. The concept is
 * explained in more detail in connection to @ref Iterators.
 *
 * This class follows mainly the route laid out by the accessor library
 * declared in the triangulation library (TriaAccessor). It enables the user
 * to access the degrees of freedom on lines, quads, or hexes. The first
 * template argument of this class determines the dimensionality of the object
 * under consideration: 1 for lines, 2 for quads, and 3 for hexes. From the
 * second template argument we can deduce the dimensionality of the
 * triangulation to which this object belongs as well as the dimensionality of
 * the space in which it is embedded. The second argument also denotes the
 * type of DoF handler we should work on. It can either be
 * ::DoFHandler or hp::DoFHandler.
 *
 * Depending on whether the structural dimension of the object
 * accessed equals the dimension on which the DoF handler object
 * operates, this class is derived from CellAccessor or
 * TriaAccessor. This means that, for example accessors to quads
 * in 2d have access to all the mesh aspects of cells, whereas
 * accessors to quads in 3d can only access things that make sense for
 * faces.
 *
 *
 * <h3>Usage</h3>
 *
 * Usage is best to happen through the typedefs to the various kinds
 * of iterators provided by the DoFHandler and hp::DoFHandler classes,
 * since they are more secure to changes in the class naming and
 * template interface as well as providing easier typing (much less
 * complicated names!).
 *
 *
 * <h3>Inheritance</h3>
 *
 * If the structural dimension given by the first template argument
 * equals the dimension of the DoFHandler (given as the second
 * template argument), then we are obviously dealing with cells,
 * rather than lower-dimensional objects. In that case, inheritance is
 * from CellAccessor, to provide access to all the cell specific
 * information afforded by that class. Otherwise, i.e. for
 * lower-dimensional objects, inheritance is from TriaAccessor.
 *
 * There is a DoFCellAccessor class that provides the
 * equivalent to the CellAccessor class.
 *
 * @ingroup dofs
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998, 2006, 2008
 */
template <int structdim, class DH>
class DoFAccessor : public internal::DoFAccessor::Inheritance<structdim, DH::dimension, DH::space_dimension>::BaseClass
{
  public:

				     /**
				      * A static variable that allows users of
				      * this class to discover the value of
				      * the second template argument.
				      */
    static const unsigned int dimension=DH::dimension;

				     /**
				      * A static variable that allows users of
				      * this class to discover the value of
				      * the third template argument.
				      */
    static const unsigned int space_dimension=DH::space_dimension;

				     /**
				      * Declare a typedef to the base
				      * class to make accessing some
				      * of the exception classes
				      * simpler.
				      */
    typedef
    typename internal::DoFAccessor::Inheritance<structdim, dimension, space_dimension>::BaseClass
    BaseClass;

				     /**
				      * Data type passed by the iterator class.
				      */
    typedef DH AccessorData;

				     /**
				      * @name Constructors
				      */
				     /**
				      * @{
				      */

				     /**
				      * Default constructor. Provides
				      * an accessor that can't be
				      * used.
				      */
    DoFAccessor ();

				     /**
				      * Constructor
				      */
    DoFAccessor (const Triangulation<DH::dimension,DH::space_dimension> *tria,
		 const int                 level,
		 const int                 index,
		 const DH                 *local_data);

				     /**
				      * Conversion constructor. This
				      * constructor exists to make certain
				      * constructs simpler to write in
				      * dimension independent code. For
				      * example, it allows assigning a face
				      * iterator to a line iterator, an
				      * operation that is useful in 2d but
				      * doesn't make any sense in 3d. The
				      * constructor here exists for the
				      * purpose of making the code conform to
				      * C++ but it will unconditionally abort;
				      * in other words, assigning a face
				      * iterator to a line iterator is better
				      * put into an if-statement that checks
				      * that the dimension is two, and assign
				      * to a quad iterator in 3d (an operator
				      * that, without this constructor would
				      * be illegal if we happen to compile for
				      * 2d).
				      */
    template <int structdim2, int dim2, int spacedim2>
    DoFAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &);

				     /**
				      * Another conversion operator
				      * between objects that don't
				      * make sense, just like the
				      * previous one.
				      */
    template <int dim2, class DH2>
    DoFAccessor (const DoFAccessor<dim2, DH2> &);

				     /**
				      * @}
				      */

				     /**
				      * Return a handle on the
				      * DoFHandler object which we
				      * are using.
				      */
    const DH &
    get_dof_handler () const;

				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const DoFAccessor<structdim, DH> &a);

				     /**
				      * Copy operator used by the
				      * iterator class. Keeps the
				      * previously set dof handler,
				      * but sets the object
				      * coordinates of the TriaAccessor.
				      */
    void copy_from (const TriaAccessorBase<structdim, DH::dimension, DH::space_dimension> &da);

				     /**
				      * Return an iterator pointing to
				      * the the parent.
				      */
    TriaIterator<DoFAccessor<structdim,DH> >
    parent () const;

				     /**
				      *  @name Accessing sub-objects
				      */
				     /**
				      * @{
				      */

				     /**
				      * Return an iterator pointing to
				      * the the @p c-th child.
				      */
    TriaIterator<DoFAccessor<structdim,DH> >
    child (const unsigned int c) const;

				     /**
				      * Pointer to the @p ith line
				      * bounding this object.
				      */
    typename internal::DoFHandler::Iterators<DH>::line_iterator
    line (const unsigned int i) const;

    				     /**
				      * Pointer to the @p ith quad
				      * bounding this object.
				      */
    typename internal::DoFHandler::Iterators<DH>::quad_iterator
    quad (const unsigned int i) const;

				     /**
				      * @}
				      */

				     /**
				      *  @name Accessing the DoF indices of this object
				      */
				     /**
				      * @{
				      */

    				     /**
				      * Return the indices of the dofs of this
				      * object in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * The vector has to have the
				      * right size before being passed
				      * to this function.
				      *
				      * This function is most often
				      * used on active objects (edges,
				      * faces, cells). It can be used
				      * on non-active objects as well
				      * (i.e. objects that have
				      * children), but only if the
				      * finite element under
				      * consideration has degrees of
				      * freedom exclusively on
				      * vertices. Otherwise, the
				      * function doesn't make much
				      * sense, since for example
				      * inactive edges do not have
				      * degrees of freedom associated
				      * with them at all.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      *
				      * For cells, there is only a
				      * single possible finite element
				      * index (namely the one for that
				      * cell, returned by
				      * <code>cell-@>active_fe_index</code>. Consequently,
				      * the derived DoFCellAccessor
				      * class has an overloaded
				      * version of this function that
				      * calls the present function
				      * with
				      * <code>cell-@>active_fe_index</code>
				      * as last argument.
				      */
    void get_dof_indices (std::vector<unsigned int> &dof_indices,
			  const unsigned int fe_index = DH::default_fe_index) const;

				     /**
				      * Global DoF index of the <i>i</i>
				      * degree associated with the @p vertexth
				      * vertex of the present cell.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      */
    unsigned int vertex_dof_index (const unsigned int vertex,
				   const unsigned int i,
				   const unsigned int fe_index = DH::default_fe_index) const;

				     /**
				      * Index of the <i>i</i>th degree
				      * of freedom of this object.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      */
    unsigned int dof_index (const unsigned int i,
			    const unsigned int fe_index = DH::default_fe_index) const;

				     /**
				      * @}
				      */

				     /**
				      *  @name Accessing the finite element associated with this object
				      */
				     /**
				      * @{
				      */

                                     /**
                                      * Return the number of finite
                                      * elements that are active on a
                                      * given object.
                                      *
                                      * For non-hp DoFHandler objects,
                                      * the answer is of course always
                                      * one. However, for
                                      * hp::DoFHandler objects, this
                                      * isn't the case: If this is a
                                      * cell, the answer is of course
                                      * one. If it is a face, the
                                      * answer may be one or two,
                                      * depending on whether the two
                                      * adjacent cells use the same
                                      * finite element or not. If it
                                      * is an edge in 3d, the possible
                                      * return value may be one or any
                                      * other value larger than that.
                                      */
    unsigned int
    n_active_fe_indices () const;

				     /**
				      * Return the @p n-th active fe
				      * index on this object. For
				      * cells and all non-hp objects,
				      * there is only a single active
				      * fe index, so the argument must
				      * be equal to zero. For
				      * lower-dimensional hp objects,
				      * there are
				      * n_active_fe_indices() active
				      * finite elements, and this
				      * function can be queried for
				      * their indices.
				      */
    unsigned int
    nth_active_fe_index (const unsigned int n) const;

				     /**
				      * Return true if the finite
				      * element with given index is
				      * active on the present
				      * object. For non-hp DoF
				      * accessors, this is of course
				      * the case only if @p fe_index
				      * equals zero. For cells, it is
				      * the case if @p fe_index equals
				      * active_fe_index() of this
				      * cell. For faces and other
				      * lower-dimensional objects,
				      * there may be more than one @p
				      * fe_index that are active on
				      * any given object (see
				      * n_active_fe_indices()).
				      */
    bool
    fe_index_is_active (const unsigned int fe_index) const;

				     /**
				      * Return a reference to the finite
				      * element used on this object with the
				      * given @p fe_index. @p fe_index must be
				      * used on this object,
				      * i.e. <code>fe_index_is_active(fe_index)</code>
				      * must return true.
				      */
    const FiniteElement<DH::dimension,DH::space_dimension> &
    get_fe (const unsigned int fe_index) const;

				     /**
				      * @}
				      */

                                     /**
				      * Exceptions for child classes
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
				      * Store the address of the DoFHandler object
				      * to be accessed.
				      */
    DH *dof_handler;

				     /**
				      *  Compare for equality.
				      */
    bool operator == (const DoFAccessor &) const;

				     /**
				      * Compare for inequality.
				      */
    bool operator != (const DoFAccessor &) const;

				     /**
				      * Reset the DoF handler pointer.
				      */
    void set_dof_handler (DH *dh);

    				     /**
				      * Set the index of the
				      * <i>i</i>th degree of freedom
				      * of this object to @p index.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      */
    void set_dof_index (const unsigned int i,
			const unsigned int index,
			const unsigned int fe_index = DH::default_fe_index) const;

				     /**
				      * Set the global index of the <i>i</i>
				      * degree on the @p vertex-th vertex of
				      * the present cell to @p index.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      */
    void set_vertex_dof_index (const unsigned int vertex,
			       const unsigned int i,
			       const unsigned int index,
			       const unsigned int fe_index = DH::default_fe_index) const;

                                     /**
                                      * Iterator classes need to be friends
                                      * because they need to access operator==
                                      * and operator!=.
                                      */
    template <typename> friend class TriaRawIterator;


  private:
    				     /**
				      *  Copy operator. This is normally used
				      *  in a context like <tt>iterator a,b;
				      *  *a=*b;</tt>. Presumably, the intent
				      *  here is to copy the object pointed to
				      *  by @p b to the object pointed to by
				      *  @p a. However, the result of
				      *  dereferencing an iterator is not an
				      *  object but an accessor; consequently,
				      *  this operation is not useful for
				      *  iterators on triangulations. We
				      *  declare this function here private,
				      *  thus it may not be used from outside.
				      *  Furthermore it is not implemented and
				      *  will give a linker error if used
				      *  anyway.
				      */
    DoFAccessor<structdim,DH> &
    operator = (const DoFAccessor<structdim,DH> &da);

				     /**
				      * Make the DoFHandler class a friend so
				      * that it can call the set_xxx()
				      * functions.
				      */
    template <int dim, int spacedim> friend class DoFHandler;
    template <int dim, int spacedim> friend class hp::DoFHandler;

    friend struct internal::DoFHandler::Policy::Implementation;
    friend struct internal::DoFHandler::Implementation;
    friend struct internal::hp::DoFHandler::Implementation;
    friend struct internal::DoFCellAccessor::Implementation;
};



/**
 * Specialization of the general DoFAccessor class template for the
 * case of zero-dimensional objects (a vertex) that are the face of a
 * one-dimensional cell in spacedim space dimensions. Since vertices
 * function differently than general faces, this class does a few
 * things differently than the general template, but the interface
 * should look the same.
 *
 * @author Wolfgang Bangerth, 2010
 */
template <template <int, int> class DH, int spacedim>
class DoFAccessor<0,DH<1,spacedim> > : public TriaAccessor<0,1,spacedim>
{
  public:

				     /**
				      * A static variable that allows users of
				      * this class to discover the value of
				      * the second template argument.
				      */
    static const unsigned int dimension=1;

				     /**
				      * A static variable that allows users of
				      * this class to discover the value of
				      * the third template argument.
				      */
    static const unsigned int space_dimension=spacedim;

				     /**
				      * Declare a typedef to the base
				      * class to make accessing some
				      * of the exception classes
				      * simpler.
				      */
    typedef TriaAccessor<0,1,spacedim> BaseClass;

				     /**
				      * Data type passed by the iterator class.
				      */
    typedef DH<1,spacedim> AccessorData;

				     /**
				      * @name Constructors
				      */
				     /**
				      * @{
				      */

				     /**
				      * Default constructor. Provides
				      * an accessor that can't be
				      * used.
				      */
    DoFAccessor ();

				     /**
				      * Constructor to be used if the
				      * object here refers to a vertex
				      * of a one-dimensional
				      * triangulation, i.e. a face of
				      * the triangulation.
				      *
				      * Since there is no mapping from
				      * vertices to cells, an accessor
				      * object for a point has no way
				      * to figure out whether it is at
				      * the boundary of the domain or
				      * not. Consequently, the second
				      * argument must be passed by the
				      * object that generates this
				      * accessor -- e.g. a 1d cell
				      * that can figure out whether
				      * its left or right vertex are
				      * at the boundary.
				      *
				      * The third argument is the
				      * global index of the vertex we
				      * point to.
				      *
				      * The fourth argument is a
				      * pointer to the DoFHandler
				      * object.
				      *
				      * This iterator can only be
				      * called for one-dimensional
				      * triangulations.
				      */
    DoFAccessor (const Triangulation<1,spacedim> * tria,
		 const typename TriaAccessor<0,1,spacedim>::VertexKind vertex_kind,
		 const unsigned int    vertex_index,
		 const DH<1,spacedim> * dof_handler);

				     /**
				      * Constructor. This constructor
				      * exists in order to maintain
				      * interface compatibility with
				      * the other accessor
				      * classes. However, it doesn't
				      * do anything useful here and so
				      * may not actually be called.
				      */
    DoFAccessor (const Triangulation<1,spacedim> *,
		 const int = 0,
		 const int = 0,
		 const DH<1,spacedim> * = 0);

				     /**
				      * Conversion constructor. This
				      * constructor exists to make certain
				      * constructs simpler to write in
				      * dimension independent code. For
				      * example, it allows assigning a face
				      * iterator to a line iterator, an
				      * operation that is useful in 2d but
				      * doesn't make any sense in 3d. The
				      * constructor here exists for the
				      * purpose of making the code conform to
				      * C++ but it will unconditionally abort;
				      * in other words, assigning a face
				      * iterator to a line iterator is better
				      * put into an if-statement that checks
				      * that the dimension is two, and assign
				      * to a quad iterator in 3d (an operator
				      * that, without this constructor would
				      * be illegal if we happen to compile for
				      * 2d).
				      */
    template <int structdim2, int dim2, int spacedim2>
    DoFAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &);

				     /**
				      * Another conversion operator
				      * between objects that don't
				      * make sense, just like the
				      * previous one.
				      */
    template <int dim2, class DH2>
    DoFAccessor (const DoFAccessor<dim2, DH2> &);

				     /**
				      * @}
				      */

				     /**
				      * Return a handle on the
				      * DoFHandler object which we
				      * are using.
				      */
    const DH<1,spacedim> &
    get_dof_handler () const;

				     /**
				      * Copy operator.
				      */
    DoFAccessor<0,DH<1,spacedim> > &
    operator = (const DoFAccessor<0,DH<1,spacedim> > &da);

				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const DoFAccessor<0, DH<1,spacedim> > &a);

				     /**
				      * Copy operator used by the
				      * iterator class. Keeps the
				      * previously set dof handler,
				      * but sets the object
				      * coordinates of the TriaAccessor.
				      */
    void copy_from (const TriaAccessorBase<0, 1, spacedim> &da);

				     /**
				      * Return an iterator pointing to
				      * the the parent.
				      */
    TriaIterator<DoFAccessor<0,DH<1,spacedim> > >
    parent () const;

				     /**
				      *  @name Accessing sub-objects
				      */
				     /**
				      * @{
				      */

				     /**
				      * Return an iterator pointing to
				      * the the @p c-th child.
				      */
    TriaIterator<DoFAccessor<0,DH<1,spacedim> > >
    child (const unsigned int c) const;

				     /**
				      * Pointer to the @p ith line
				      * bounding this object.
				      */
    typename internal::DoFHandler::Iterators<DH<1,spacedim> >::line_iterator
    line (const unsigned int i) const;

    				     /**
				      * Pointer to the @p ith quad
				      * bounding this object.
				      */
    typename internal::DoFHandler::Iterators<DH<1,spacedim> >::quad_iterator
    quad (const unsigned int i) const;

				     /**
				      * @}
				      */

				     /**
				      *  @name Accessing the DoF indices of this object
				      */
				     /**
				      * @{
				      */

    				     /**
				      * Return the indices of the dofs of this
				      * object in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * The vector has to have the
				      * right size before being passed
				      * to this function.
				      *
				      * This function is most often
				      * used on active objects (edges,
				      * faces, cells). It can be used
				      * on non-active objects as well
				      * (i.e. objects that have
				      * children), but only if the
				      * finite element under
				      * consideration has degrees of
				      * freedom exclusively on
				      * vertices. Otherwise, the
				      * function doesn't make much
				      * sense, since for example
				      * inactive edges do not have
				      * degrees of freedom associated
				      * with them at all.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      *
				      * For cells, there is only a
				      * single possible finite element
				      * index (namely the one for that
				      * cell, returned by
				      * <code>cell-@>active_fe_index</code>. Consequently,
				      * the derived DoFCellAccessor
				      * class has an overloaded
				      * version of this function that
				      * calls the present function
				      * with
				      * <code>cell-@>active_fe_index</code>
				      * as last argument.
				      */
    void get_dof_indices (std::vector<unsigned int> &dof_indices,
			  const unsigned int fe_index = AccessorData::default_fe_index) const;

				     /**
				      * Global DoF index of the <i>i</i>
				      * degree associated with the @p vertexth
				      * vertex of the present cell.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      */
    unsigned int vertex_dof_index (const unsigned int vertex,
				   const unsigned int i,
				   const unsigned int fe_index = AccessorData::default_fe_index) const;

				     /**
				      * Index of the <i>i</i>th degree
				      * of freedom of this object.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      */
    unsigned int dof_index (const unsigned int i,
			    const unsigned int fe_index = AccessorData::default_fe_index) const;

				     /**
				      * @}
				      */

				     /**
				      *  @name Accessing the finite element associated with this object
				      */
				     /**
				      * @{
				      */

                                     /**
                                      * Return the number of finite
                                      * elements that are active on a
                                      * given object.
                                      *
                                      * For non-hp DoFHandler objects,
                                      * the answer is of course always
                                      * one. However, for
                                      * hp::DoFHandler objects, this
                                      * isn't the case: If this is a
                                      * cell, the answer is of course
                                      * one. If it is a face, the
                                      * answer may be one or two,
                                      * depending on whether the two
                                      * adjacent cells use the same
                                      * finite element or not. If it
                                      * is an edge in 3d, the possible
                                      * return value may be one or any
                                      * other value larger than that.
                                      */
    unsigned int
    n_active_fe_indices () const;

				     /**
				      * Return the @p n-th active fe
				      * index on this object. For
				      * cells and all non-hp objects,
				      * there is only a single active
				      * fe index, so the argument must
				      * be equal to zero. For
				      * lower-dimensional hp objects,
				      * there are
				      * n_active_fe_indices() active
				      * finite elements, and this
				      * function can be queried for
				      * their indices.
				      */
    unsigned int
    nth_active_fe_index (const unsigned int n) const;

				     /**
				      * Return true if the finite
				      * element with given index is
				      * active on the present
				      * object. For non-hp DoF
				      * accessors, this is of course
				      * the case only if @p fe_index
				      * equals zero. For cells, it is
				      * the case if @p fe_index equals
				      * active_fe_index() of this
				      * cell. For faces and other
				      * lower-dimensional objects,
				      * there may be more than one @p
				      * fe_index that are active on
				      * any given object (see
				      * n_active_fe_indices()).
				      */
    bool
    fe_index_is_active (const unsigned int fe_index) const;

				     /**
				      * Return a reference to the finite
				      * element used on this object with the
				      * given @p fe_index. @p fe_index must be
				      * used on this object,
				      * i.e. <code>fe_index_is_active(fe_index)</code>
				      * must return true.
				      */
    const FiniteElement<DH<1,spacedim>::dimension,DH<1,spacedim>::space_dimension> &
    get_fe (const unsigned int fe_index) const;

				     /**
				      * @}
				      */

                                     /**
				      * Exceptions for child classes
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
				      * Store the address of the DoFHandler object
				      * to be accessed.
				      */
    DH<1,spacedim> *dof_handler;

				     /**
				      *  Compare for equality.
				      */
    bool operator == (const DoFAccessor &) const;

				     /**
				      * Compare for inequality.
				      */
    bool operator != (const DoFAccessor &) const;

				     /**
				      * Reset the DoF handler pointer.
				      */
    void set_dof_handler (DH<1,spacedim> *dh);

    				     /**
				      * Set the index of the
				      * <i>i</i>th degree of freedom
				      * of this object to @p index.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      */
    void set_dof_index (const unsigned int i,
			const unsigned int index,
			const unsigned int fe_index = AccessorData::default_fe_index) const;

				     /**
				      * Set the global index of the <i>i</i>
				      * degree on the @p vertex-th vertex of
				      * the present cell to @p index.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      */
    void set_vertex_dof_index (const unsigned int vertex,
			       const unsigned int i,
			       const unsigned int index,
			       const unsigned int fe_index = AccessorData::default_fe_index) const;

                                     /**
                                      * Iterator classes need to be friends
                                      * because they need to access operator==
                                      * and operator!=.
                                      */
    template <typename> friend class TriaRawIterator;


				     /**
				      * Make the DoFHandler class a friend so
				      * that it can call the set_xxx()
				      * functions.
				      */
    template <int, int> friend class DoFHandler;
    template <int, int> friend class hp::DoFHandler;

    friend struct internal::DoFHandler::Policy::Implementation;
    friend struct internal::DoFHandler::Implementation;
    friend struct internal::hp::DoFHandler::Implementation;
    friend struct internal::DoFCellAccessor::Implementation;
};


/* -------------------------------------------------------------------------- */


/**
 * Grant access to the degrees of freedom on a cell.
 *
 * Note that since for the class we derive from, i.e. <tt>DoFAccessor<dim></tt>,
 * the two template parameters are equal, the base class is actually derived from
 * CellAccessor, which makes the functions of this class available to the
 * DoFCellAccessor class as well.
 *
 * @ingroup dofs
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998
 */
template <class DH>
class DoFCellAccessor :  public DoFAccessor<DH::dimension,DH>
{
  public:
				     /**
				      * Extract dimension from DH.
				      */
    static const unsigned int dim = DH::dimension;

				     /**
				      * Extract space dimension from DH.
				      */
    static const unsigned int spacedim = DH::space_dimension;

				     /**
				      * Declare the data type that
				      * this accessor class expects to
				      * get passed from the iterator
				      * classes.
				      */
    typedef typename DoFAccessor<DH::dimension,DH>::AccessorData AccessorData;

				     /**
				      * Declare a typedef to the base
				      * class to make accessing some
				      * of the exception classes
				      * simpler.
				      */
    typedef DoFAccessor<DH::dimension,DH> BaseClass;

    				     /**
				      * Define the type of the
				      * container this is part of.
				      */
    typedef DH Container;

				     /**
				      * @name Constructors
				      */
				     /**
				      * @{
				      */

    				     /**
				      * Constructor
				      */
    DoFCellAccessor (const Triangulation<DH::dimension,DH::space_dimension> *tria,
		     const int                 level,
		     const int                 index,
		     const AccessorData       *local_data);

				     /**
				      * Conversion constructor. This
				      * constructor exists to make certain
				      * constructs simpler to write in
				      * dimension independent code. For
				      * example, it allows assigning a face
				      * iterator to a line iterator, an
				      * operation that is useful in 2d but
				      * doesn't make any sense in 3d. The
				      * constructor here exists for the
				      * purpose of making the code conform to
				      * C++ but it will unconditionally abort;
				      * in other words, assigning a face
				      * iterator to a line iterator is better
				      * put into an if-statement that checks
				      * that the dimension is two, and assign
				      * to a quad iterator in 3d (an operator
				      * that, without this constructor would
				      * be illegal if we happen to compile for
				      * 2d).
				      */
    template <int structdim2, int dim2, int spacedim2>
    DoFCellAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &);

				     /**
				      * Another conversion operator
				      * between objects that don't
				      * make sense, just like the
				      * previous one.
				      */
    template <int dim2, class DH2>
    DoFCellAccessor (const DoFAccessor<dim2, DH2> &);

				     /**
				      * @}
				      */

    				 /**
				      * Return the parent as a DoF
				      * cell iterator. This
				      * function is needed since the
				      * parent function of the base
				      * class returns a cell accessor
				      * without access to the DoF
				      * data.
				      */
    typename internal::DoFHandler::Iterators<DH>::cell_iterator
    parent () const;

				     /**
				      *  @name Accessing sub-objects and neighbors
				      */
				     /**
				      * @{
				      */

				     /**
				      * Return the @p ith neighbor as
				      * a DoF cell iterator. This
				      * function is needed since the
				      * neighbor function of the base
				      * class returns a cell accessor
				      * without access to the DoF
				      * data.
				      */
    typename internal::DoFHandler::Iterators<DH>::cell_iterator
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
    typename internal::DoFHandler::Iterators<DH>::cell_iterator
    child (const unsigned int) const;

    				     /**
				      * Return an iterator to the @p ith face
				      * of this cell.
				      *
				      * This function is not implemented in
				      * 1D, and maps to DoFAccessor::line
				      * in 2D.
				      */
    typename internal::DoFHandler::Iterators<DH>::face_iterator
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
    typename internal::DoFHandler::Iterators<DH>::cell_iterator
    neighbor_child_on_subface (const unsigned int face_no,
                               const unsigned int subface_no) const;

				     /**
				      * @}
				      */

				     /**
				      *  @name Extracting values from global vectors
				      */
				     /**
				      * @{
				      */

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
				      * The input vector may be either
				      * a <tt>Vector<float></tt>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc or Trilinos vector if
				      * deal.II is compiled to support
				      * these libraries. It is in the
				      * responsibility of the caller
				      * to assure that the types of
				      * the numbers stored in input
				      * and output vectors are
				      * compatible and with similar
				      * accuracy.
				      */
    template <class InputVector, typename number>
    void get_dof_values (const InputVector &values,
			 Vector<number>    &local_values) const;

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
				      * The input vector may be either
				      * a <tt>Vector<float></tt>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc or Trilinos vector if
				      * deal.II is compiled to support
				      * these libraries. It is in the
				      * responsibility of the caller
				      * to assure that the types of
				      * the numbers stored in input
				      * and output vectors are
				      * compatible and with similar
				      * accuracy.
				      */
    template <class InputVector, typename ForwardIterator>
    void get_dof_values (const InputVector &values,
			 ForwardIterator    local_values_begin,
			 ForwardIterator    local_values_end) const;

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
				      * BlockVector<double>, or a PETSc or
				      * Trilinos vector if deal.II is
				      * compiled to support these
				      * libraries. It is in the
				      * responsibility of the caller to
				      * assure that the types of the numbers
				      * stored in input and output vectors
				      * are compatible and with similar
				      * accuracy. The ConstraintMatrix
				      * passed as an argument to this
				      * function makes sure that constraints
				      * are correctly distributed when the
				      * dof values are calculated.
				      */
    template <class InputVector, typename ForwardIterator>
    void get_dof_values (const ConstraintMatrix &constraints,
			 const InputVector      &values,
			 ForwardIterator         local_values_begin,
			 ForwardIterator         local_values_end) const;

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
				      * support these libraries. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class OutputVector, typename number>
    void set_dof_values (const Vector<number> &local_values,
			 OutputVector         &values) const;

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
				      * size beforehand.
				      *
				      * Unlike the get_dof_values()
				      * function, this function works
				      * on cells rather than to lines,
				      * quads, and hexes, since
				      * interpolation is presently
				      * only provided for cells by the
				      * finite element classes.
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
				      * support these libraries. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class OutputVector, typename number>
    void set_dof_values_by_interpolation (const Vector<number> &local_values,
					  OutputVector         &values) const;

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
				      * Distribute a local (cell based)
				      * vector in iterator format to a
				      * global one by mapping the local
				      * numbering of the degrees of freedom
				      * to the global one and entering the
				      * local values into the global vector.
				      *
				      * The elements are <em>added</em> up
				      * to the elements in the global
				      * vector, rather than just set, since
				      * this is usually what one wants.
				      */
    template <typename ForwardIterator, typename OutputVector>
    void
    distribute_local_to_global (ForwardIterator   local_source_begin,
				ForwardIterator   local_source_end,
				OutputVector     &global_destination) const;

				     /**
				      * Distribute a local (cell based)
				      * vector in iterator format to a
				      * global one by mapping the local
				      * numbering of the degrees of freedom
				      * to the global one and entering the
				      * local values into the global vector.
				      *
				      * The elements are <em>added</em> up
				      * to the elements in the global
				      * vector, rather than just set, since
				      * this is usually what one
				      * wants. Moreover, the
				      * ConstraintMatrix passed to this
				      * function makes sure that also
				      * constraints are eliminated in this
				      * process.
				      */
    template <typename ForwardIterator, typename OutputVector>
    void
    distribute_local_to_global (const ConstraintMatrix &constraints,
                                ForwardIterator         local_source_begin,
				ForwardIterator         local_source_end,
				OutputVector           &global_destination) const;

				     /**
				      * This function does much the
				      * same as the
				      * <tt>distribute_local_to_global(Vector,Vector)</tt>
				      * function, but operates on
				      * matrices instead of
				      * vectors. If the matrix type is
				      * a sparse matrix then it is
				      * supposed to have non-zero
				      * entry slots where required.
				      */
    template <typename number, typename OutputMatrix>
    void
    distribute_local_to_global (const FullMatrix<number> &local_source,
                                OutputMatrix             &global_destination) const;

				     /**
				      * This function does what the two
				      * <tt>distribute_local_to_global</tt>
				      * functions with vector and matrix
				      * argument do, but all at once.
				      */
    template <typename number, typename OutputMatrix, typename OutputVector>
    void
    distribute_local_to_global (const FullMatrix<number> &local_matrix,
				const Vector<number>     &local_vector,
                                OutputMatrix             &global_matrix,
				OutputVector             &global_vector) const;

				     /**
				      * @}
				      */

				     /**
				      *  @name Accessing the DoF indices of this object
				      */
				     /**
				      * @{
				      */

    				     /**
				      * Return the indices of the dofs of this
				      * quad in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand.
				      *
				      * This function reimplements the
				      * same function in the base
				      * class. The functions in the
				      * base classes are available for
				      * all geometric objects,
				      * i.e. even in 3d they can be
				      * used to access the dof indices
				      * of edges, for example. On the
				      * other hand, the most common
				      * case is clearly the use on
				      * cells, which is why we cache
				      * the array for each cell, but
				      * not edge. To retrieve the
				      * cached values, rather than
				      * collect the necessary
				      * information every time, this
				      * function overwrites the one in
				      * the base class.
				      *
				      * This function is most often
				      * used on active objects (edges,
				      * faces, cells). It can be used
				      * on non-active objects as well
				      * (i.e. objects that have
				      * children), but only if the
				      * finite element under
				      * consideration has degrees of
				      * freedom exclusively on
				      * vertices. Otherwise, the
				      * function doesn't make much
				      * sense, since for example
				      * inactive edges do not have
				      * degrees of freedom associated
				      * with them at all.
				      */
    void get_dof_indices (std::vector<unsigned int> &dof_indices) const;

				     /**
				      * @}
				      */

				     /**
				      *  @name Accessing the finite element associated with this object
				      */
				     /**
				      * @{
				      */

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
    const FiniteElement<DH::dimension,DH::space_dimension> &
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
				      * @}
				      */

				     /**
				      * Set the DoF indices of this
				      * cell to the given values. This
				      * function bypasses the DoF
				      * cache, if one exists for the
				      * given DoF handler class.
				      */
    void set_dof_indices (const std::vector<unsigned int> &dof_indices);

				     /**
				      * Update the cache in which we
				      * store the dof indices of this
				      * cell.
				      */
    void update_cell_dof_indices_cache () const;

  private:
    				     /**
				      *  Copy operator. This is normally used
				      *  in a context like <tt>iterator a,b;
				      *  *a=*b;</tt>. Presumably, the intent
				      *  here is to copy the object pointed to
				      *  by @p b to the object pointed to by
				      *  @p a. However, the result of
				      *  dereferencing an iterator is not an
				      *  object but an accessor; consequently,
				      *  this operation is not useful for
				      *  iterators on triangulations. We
				      *  declare this function here private,
				      *  thus it may not be used from outside.
				      *  Furthermore it is not implemented and
				      *  will give a linker error if used
				      *  anyway.
				      */
    DoFCellAccessor<DH> &
    operator = (const DoFCellAccessor<DH> &da);

				     /**
				      * Make the DoFHandler class a
				      * friend so that it can call the
				      * update_cell_dof_indices_cache()
				      * function
				      */
    template <int dim, int spacedim> friend class DoFHandler;
    friend struct internal::DoFCellAccessor::Implementation;
};


DEAL_II_NAMESPACE_CLOSE

// include more templates
#include "dof_accessor.templates.h"


/*----------------------------   dof_iterator.h     ---------------------------*/

#endif
/*----------------------------   dof_iterator.h     ---------------------------*/
