//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__tria_accessor_h
#define __deal2__tria_accessor_h


#include <base/config.h>
#include <base/exceptions.h>
#include <grid/tria_iterator_base.h>


namespace std
{
  template<class T1, class T2>
  struct pair;
}

DEAL_II_NAMESPACE_OPEN

template <int dim> class Point;

template <int dim> class Triangulation;
template <int dim, typename Accessor> class TriaRawIterator;
template <int dim, typename Accessor> class TriaIterator;
template <int dim, typename Accessor> class TriaActiveIterator;


namespace internal
{
  namespace Triangulation
  {
    template <int dim> class TriaObject;
    template <typename G> class TriaObjects;
  }
}
template <int celldim, int dim> class TriaObjectAccessor;
template <int dim>              class TriaObjectAccessor<0, dim>;
template <int dim>              class TriaObjectAccessor<1, dim>;
template <int dim>              class TriaObjectAccessor<2, dim>;
template <int dim>              class TriaObjectAccessor<3, dim>;


// note: the file tria_accessor.templates.h is included at the end of
// this file.  this includes a lot of templates. originally, this was
// only done in debug mode, but led to cyclic reduction problems and
// so is now on by default.


/**
 * A namespace that contains exception classes used by the accessor classes.
 */
namespace TriaAccessorExceptions
{
				   /**
				    * @ingroup Exceptions
				    */
  DeclException0 (ExcCellNotUsed);
				   /**
				    * @ingroup Exceptions
				    */
  DeclException1 (ExcInvalidNeighbor,
		  int,
		  << "Neighbor indices must be between 0 and 2^dim-1, but "
		  << "yours was " << arg1);
				   /**
				    * @ingroup Exceptions
				    */
  DeclException0 (ExcCellHasNoChildren);
				   /**
				    * @ingroup Exceptions
				    */
  DeclException0 (ExcUnusedCellAsChild);
				   /**
				    * @ingroup Exceptions
				    */
  DeclException1 (ExcCantSetChildren,
		  int,
		  << "You can only set the child index if the cell has no "
		  << "children, or clear it. The given "
		  << "index was " << arg1 << " (-1 means: clear children)");
				   /**
				    * @ingroup Exceptions
				    */
  DeclException0 (ExcUnusedCellAsNeighbor);
				   /**
				    * @ingroup Exceptions
				    */
  DeclException0 (ExcUncaughtCase);
				   /**
				    * @ingroup Exceptions
				    */
  DeclException0 (ExcDereferenceInvalidObject);
				   /**
				    * @ingroup Exceptions
				    */
  DeclException0 (ExcCantCompareIterators);
				   /**
				    * @ingroup Exceptions
				    */
  DeclException0 (ExcNeighborIsCoarser);
				   /**
				    * @ingroup Exceptions
				    */
  DeclException0 (ExcNeighborIsNotCoarser);
				   /**
				    * @ingroup Exceptions
				    */
  DeclException0 (ExcFacesHaveNoLevel);
}


/**
 * A base class for the accessor classes used by TriaRawIterator and
 * derived classes.
 *
 * This class offers only the basic functionality required by the
 * iterators (stores the necessary data members, offers comparison
 * operators and the like), but has no functionality to actually
 * dereference data. This is done in the derived classes.
 *
 * This template is used for faces and edges, which have no level in
 * the triangulation hierarchy associated with them. There exists a
 * partial specialization of the current class template where
 * <tt>structdim</tt> (the dimensionality of the object represented,
 * for example 1 for a line) equals <tt>dim</tt> (the dimensionality
 * of the space the object lives in, for example 3 if we solve
 * three-dimensional problems).
 *
 * Some internals of this class are discussed in @ref IteratorAccessorInternals .
 *
 * @ingroup grid
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998
 */
template <int structdim, int dim>
class TriaAccessor
{
  protected:
				     /**
				      * Declare the data type that
				      * this accessor class expects to
				      * get passed from the iterator
				      * classes. Since the pure
				      * triangulation iterators need
				      * no additional data, this data
				      * type is @p void.
				      */
    typedef void AccessorData;

				     /**
				      *  Constructor. Protected, thus
				      *  only callable from friend
				      *  classes.
				      */
    TriaAccessor (const Triangulation<dim> *parent =  0,
		  const int                 level  = -1,
		  const int                 index  = -1,
		  const AccessorData       *       =  0);

				     /**
				      *  Copy operator. Since this is
				      *  only called from iterators,
				      *  do not return anything, since
				      *  the iterator will return
				      *  itself.
				      *
				      *  This method is protected,
				      *  since it is only to be called
				      *  from the iterator class.
				      */
    void copy_from (const TriaAccessor &);

				     /**
				      *  Copy operator. This is normally
				      *  used in a context like
				      *  <tt>iterator a,b;  *a=*b;</tt>. Since
				      *  the meaning is to copy the object
				      *  pointed to by @p b to the object
				      *  pointed to by @p a and since
				      *  accessors are not real but
				      *  virtual objects, this operation
				      *  is not useful for iterators on
				      *  triangulations. We declare this
				      *  function here private, thus it may
				      *  not be used from outside.
				      *  Furthermore it is not implemented
				      *  and will give a linker error if
				      *  used anyway.
				      */
    void operator = (const TriaAccessor *);
    
				     /**
				      *  Same as above.
				      */
    TriaAccessor & operator = (const TriaAccessor &);

  protected:
    
				     /**
				      *  Compare for equality.            
				      */
    bool operator == (const TriaAccessor &) const;
	
				     /**
				      * Compare for inequality.
				      */
    bool operator != (const TriaAccessor &) const;
    
  public:
				     /**
				      * Data type to be used for passing
				      * parameters from iterators to the
				      * accessor classes in a unified
				      * way, no matter what the type of
				      * number of these parameters is.
				      */
    typedef void * LocalData;
    
				     /**@ name Iterator address and state
				      */
				     /*@{*/
				     /**
				      *  Return the level the element
				      *  pointed to belongs to.
				      *  This is only valid for cells.
				      */
    static int level ();
    
				     /**
				      *  Return the index of the
				      *  element presently pointed to
				      *  on the present level.
				      */
    int index () const;
    
				     /**
				      *  Return the state of the
				      *  iterator.  For the different
				      *  states an accessor can be in,
				      *  refer to the
				      *  TriaRawIterator
				      *  documentation.
				      */
    IteratorState::IteratorStates state () const;

				     /**
				      * Return a pointer to the
				      * triangulation which the object
				      * pointed to by this class
				      * belongs to.
				      */
    const Triangulation<dim> & get_triangulation () const;
    
				     /*@}*/
  protected:
				     /**
				      *  Used to store the index of
				      *  the element presently pointed
				      *  to on the level presentl
				      *  used.
				      */
    int present_index;
    
				     /**
				      *  Pointer to the triangulation
				      *  which we act on.
				      */
    const Triangulation<dim> *tria;

  private:
    
				     /**
				      * Dimension of the TriaObject
				      * this accessor gives access to.
				      */
    static const unsigned int objectdim = dim;

    template <int anydim, typename Accessor> friend class TriaRawIterator;
    template <int anydim, typename Accessor> friend class TriaIterator;
    template <int anydim, typename Accessor> friend class TriaActiveIterator;
};



/**
 * A base class for the accessor classes used by TriaRawIterator and
 * derived classes.
 *
 * This class offers only the basic functionality required by the
 * iterators (stores the necessary data members, offers comparison
 * operators and the like), but has no functionality to actually
 * dereference data. This is done in the derived classes.
 *
 * This template is a partial specialization of the general template
 * used for cells only, i.e. for the case where <tt>structdim</tt>
 * (the dimensionality of the object represented, for example 1 for a
 * line) equals <tt>dim</tt> (the dimensionality of the space the
 * object lives in, for example 3 if we solve three-dimensional
 * problems). The difference is that cells have a level in the
 * triangulation hierarchy associated with them, whereas faces and
 * edges do not.
 *
 * Some internals of this class are discussed in @ref IteratorAccessorInternals .
 *
 * @ingroup grid
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class TriaAccessor<dim,dim>
{
  protected:
				     /**
				      * Declare the data type that
				      * this accessor class expects to
				      * get passed from the iterator
				      * classes. Since the pure
				      * triangulation iterators need
				      * no additional data, this data
				      * type is @p void.
				      */
    typedef void AccessorData;

				     /**
				      *  Constructor. Protected, thus
				      *  only callable from friend
				      *  classes.
				      */
    TriaAccessor (const Triangulation<dim> *parent =  0,
		  const int                 level  = -1,
		  const int                 index  = -1,
		  const AccessorData       *       =  0);

				     /**
				      *  Copy operator. Since this is
				      *  only called from iterators,
				      *  do not return anything, since
				      *  the iterator will return
				      *  itself.
				      *
				      *  This method is protected,
				      *  since it is only to be called
				      *  from the iterator class.
				      */
    void copy_from (const TriaAccessor &);

				     /**
				      *  Copy operator. This is normally
				      *  used in a context like
				      *  <tt>iterator a,b;  *a=*b;</tt>. Since
				      *  the meaning is to copy the object
				      *  pointed to by @p b to the object
				      *  pointed to by @p a and since
				      *  accessors are not real but
				      *  virtual objects, this operation
				      *  is not useful for iterators on
				      *  triangulations. We declare this
				      *  function here private, thus it may
				      *  not be used from outside.
				      *  Furthermore it is not implemented
				      *  and will give a linker error if
				      *  used anyway.
				      */
    void operator = (const TriaAccessor *);
    
				     /**
				      *  Same as above.
				      */
    TriaAccessor & operator = (const TriaAccessor &);

  protected:
    
				     /**
				      *  Compare for equality.            
				      */
    bool operator == (const TriaAccessor &) const;
	
				     /**
				      * Compare for inequality.
				      */
    bool operator != (const TriaAccessor &) const;
    
  public:
				     /**
				      * Data type to be used for passing
				      * parameters from iterators to the
				      * accessor classes in a unified
				      * way, no matter what the type of
				      * number of these parameters is.
				      */
    typedef void * LocalData;
    
				     /**@ name Iterator address and state
				      */
				     /*@{*/
				     /**
				      *  Return the level the element
				      *  pointed to belongs to.
				      *  This is only valid for cells.
				      */
    int level () const;
    
				     /**
				      *  Return the index of the
				      *  element presently pointed to
				      *  on the present level.
				      */
    int index () const;
    
				     /**
				      *  Return the state of the
				      *  iterator.  For the different
				      *  states an accessor can be in,
				      *  refer to the
				      *  TriaRawIterator
				      *  documentation.
				      */
    IteratorState::IteratorStates state () const;

				     /**
				      * Return a pointer to the
				      * triangulation which the object
				      * pointed to by this class
				      * belongs to.
				      */
    const Triangulation<dim> & get_triangulation () const;
    
				     /*@}*/
  protected:
				     /**
				      *  Used to store the level
				      *  presently pointed to. For
				      *  faces this should always be
				      *  0.
				      */
    int present_level;
    
				     /**
				      *  Used to store the index of
				      *  the element presently pointed
				      *  to on the level presentl
				      *  used.
				      */
    int present_index;
    
				     /**
				      *  Pointer to the triangulation
				      *  which we act on.
				      */
    const Triangulation<dim> *tria;

  private:
    
				     /**
				      * Dimension of the TriaObject
				      * this accessor gives access to.
				      */
    static const unsigned int objectdim = dim;

    template <int anydim, typename Accessor> friend class TriaRawIterator;
    template <int anydim, typename Accessor> friend class TriaIterator;
    template <int anydim, typename Accessor> friend class TriaActiveIterator;
};




/**
 * Common template for line, quad, hex accessors.  According to
 * @p celldim, objects of this class represent lines, quadrilaterals,
 * or hexahedra in @p dim space dimensions. Concrete implementations
 * are done for specialized @p celldim template parameter. For easier
 * access, we nevertheless document all functions of the specialized
 * classes here as well. However, they are not implemented.
 *
 * @ingroup grid
 * @ingroup Accessors
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
template <int celldim, int dim>
class TriaObjectAccessor :  public TriaAccessor<celldim,dim>
{
  public:
				     /**
				      * Propagate typedef from
				      * base class to this class.
				      */
    typedef typename TriaAccessor<celldim,dim>::AccessorData AccessorData;

				     /**
				      * Constructor.
				      * By default, an illegal
				      * accessor is constructed.
				      */
    TriaObjectAccessor (const Triangulation<dim> *parent     =  0,
			const int                 level      = -1,
			const int                 index      = -1,
			const AccessorData       *local_data =  0);

				     /**
				      * Copy the data of a line. Only
				      * implemented for
				      * <tt>celldim==1</tt>.
				      */

    void set (const internal::Triangulation::TriaObject<1>&) const;
    
				     /**
				      * Copy the data of the given
				      * quad. Only implemented for
				      * <tt>celldim==2</tt>.
				      */
    void set (const internal::Triangulation::TriaObject<2>&) const;
    
				     /**
				      * Copy the data of the given
				      * hex. Only implemented for
				      * <tt>celldim==3</tt>.
				      */
    void set (const internal::Triangulation::TriaObject<3>&) const;

				     /**
				      *  Index of vertex. The convention
				      *  regarding the numbering of vertices
				      *  is laid down in the documentation of
				      *  the GeometryInfo class.
				      *
				      *  Note that the returned value is only
				      *  the index of the geometrical
				      *  vertex. It has nothing to do with
				      *  possible degrees of freedom
				      *  associated with it. For this, see the
				      *  @p DoFAccessor::vertex_dof_index
				      *  functions.
				      */ 
    int vertex_index (const unsigned int i) const;

    				     /**
				      *  Return a reference to the
				      *  @p ith vertex.
				      */
    Point<dim> & vertex (const unsigned int i) const;

				     /**
				      *  Pointer to the @p ith line
				      *  bounding this object.
				      *
				      * Implemented only for <tt>celldim>1</tt>.
				      */
    TriaIterator<dim,TriaObjectAccessor<1, dim> >
    line (const unsigned int i) const;

				     /**
				      * Line index of the @p ith
				      * line. The level is naturally
				      * the same as that of the
				      * object.
				      *
				      * Implemented only for <tt>celldim>1</tt>.
				      */
    unsigned int line_index (const unsigned int i) const;
    
    				     /**
				      *  Pointer to the @p ith quad
				      *  bounding this object.
				      *
				      * Implemented only for <tt>celldim>2</tt>.
				      */
    TriaIterator<dim,TriaObjectAccessor<2, dim> >
    quad (const unsigned int i) const;

				     /**
				      * Quad index of the @p ith
				      * quad. The level is naturally
				      * the same as that of the object.
				      *
				      * Implemented only for <tt>celldim>2</tt>.
				      */
    unsigned int quad_index (const unsigned int i) const;

				     /**
				      *  Test for the element being
				      *  used or not.  The return
				      *  value is @p true for all
				      *  iterators that are either
				      *  normal iterators or active
				      *  iterators, only raw iterators
				      *  can return @p false. Since
				      *  raw iterators are only used
				      *  in the interiors of the
				      *  library, you will not usually
				      *  need this function.
				      */
    bool used () const;

				     /**
				      * @name User data
				      */
				     /*@{*/
				     /**
				      *  Set the @p used flag. You
				      *  should know quite exactly
				      *  what you are doing of you
				      *  touch this function. It is
				      *  exclusively for internal use
				      *  in the library.
				      */
    void set_used_flag () const;

				     /**
				      *  Clear the @p used flag. You
				      *  should know quite exactly
				      *  what you are doing of you
				      *  touch this function. It is
				      *  exclusively for internal use
				      *  in the library.
				      */
    void clear_used_flag () const;

    				     /**
				      *  Read the user flag.
				      */
    bool user_flag_set () const;

				     /**
				      *  Set the user flag.
				      */
    void set_user_flag () const;

				     /**
				      *  Clear the user flag.
				      */
    void clear_user_flag () const;

				     /**
				      *  Set the user flag for this
				      * and all descendants.
				      */
    void recursively_set_user_flag () const;

    				     /**
				      *  Clear the user flag for this
				      * and all descendants. 
				      */
    void recursively_clear_user_flag () const;

				     /**
				      * Reset the user data to zero,
				      * independent if pointer or index.
				      */
    void clear_user_data () const;
    
				     /**
				      * Set the user pointer
				      * to @p p.
				      */
    void set_user_pointer (void *p) const;

				     /**
				      * Reset the user pointer
				      * to a @p NULL pointer.
				      */
    void clear_user_pointer () const;
    
				     /**
				      * Access the value of the user
				      * pointer. It is in the
				      * responsibility of the user to
				      * make sure that the pointer
				      * points to something
				      * useful. You should use the new
				      * style cast operator to
				      * maintain a minimum of
				      * typesafety, e.g.
				      * <tt>A *a=static_cast<A*>(cell->user_pointer());</tt>.
				      */
    void * user_pointer () const;

				     /**
				      * Set the user pointer of this
				      * object and all its children to
				      * the given value. This is
				      * useful for example if all
				      * cells of a certain subdomain,
				      * or all faces of a certain part
				      * of the boundary should have
				      * user pointers pointing to
				      * objects describing this part
				      * of the domain or boundary.
				      *
				      * Note that the user pointer is
				      * not inherited under mesh
				      * refinement, so after mesh
				      * refinement there might be
				      * cells or faces that don't have
				      * user pointers pointing to the
				      * describing object. In this
				      * case, simply loop over all the
				      * elements of the coarsest level
				      * that has this information, and
				      * use this function to
				      * recursively set the user
				      * pointer of all finer levels of
				      * the triangulation.
				      */
    void recursively_set_user_pointer (void *p) const;
    
				     /**
				      * Clear the user pointer of this
				      * object and all of its
				      * descendants. The same holds as
				      * said for the
				      * recursively_set_user_pointer()
				      * function.
				      */
    void recursively_clear_user_pointer () const;
    
				     /**
				      * Set the user index
				      * to @p p.
				      */
    void set_user_index (unsigned int p) const;
    
				     /**
				      * Reset the user index to 0.
				      */
    void clear_user_index () const;
    
				     /**
				      * Access the value of the user
				      * index.
				      */
    unsigned int user_index () const;

				     /**
				      * Set the user index of this
				      * object and all its children.
				      *
				      * Note that the user index is
				      * not inherited under mesh
				      * refinement, so after mesh
				      * refinement there might be
				      * cells or faces that don't have
				      * the expected user indices. In
				      * this case, simply loop over
				      * all the elements of the
				      * coarsest level that has this
				      * information, and use this
				      * function to recursively set
				      * the user index of all finer
				      * levels of the triangulation.
				      */
    void recursively_set_user_index (unsigned int p) const;

				     /**
				      * Clear the user index of this
				      * object and all of its
				      * descendants. The same holds as
				      * said for the
				      * recursively_set_user_index()
				      * function.
				      */
    void recursively_clear_user_index () const;
				     /*@}*/
    
				     /**
				      *  Pointer to the @p ith
				      *  child.
				      */
    TriaIterator<dim,TriaObjectAccessor<celldim,dim> >
    child (const unsigned int i) const;

				     /**
				      *  Index of the @p ith child.
				      *  The level of the child is one
				      *  higher than that of the
				      *  present cell, if the children
				      *  of a cell are accessed. The
				      *  children of faces have no level.
				      *  If the child does not exist, -1
				      *  is returned.
				      */
    int child_index (const unsigned int i) const;

				     /**
				      *  Set the child field. Since we
				      *  only store the index of the
				      *  first child (the others
				      *  follow directly) only one
				      *  child index is to be
				      *  given. The level of the child
				      *  is one level up of the level
				      *  of the cell to which this
				      *  iterator points.
				      */
    void set_children (const int index) const;
	
				     /**
				      *  Clear the child field,
				      *  i.e. set it to a value which
				      *  indicates that this cell has
				      *  no children.
				      */
    void clear_children () const;
    
				     /**
				      *  Test whether the object has
				      *  children.
				      */
    bool has_children () const;

				     /**
				      * Return the number of immediate
				      * children of this object.
				      *
				      * An assertion assures that this
				      * function can only be called
				      * for TriaObjects with
				      * <tt>has_children()==true</tt>.
				      */
    unsigned int n_children() const;

				     /**
				      * Number of times that this
				      * object is refined. Note that
				      * not all its children are
				      * refined that often (which is
				      * why we prepend @p max_), the
				      * returned number is rather the
				      * maximum number of refinement
				      * in any branch of children of
				      * this object.
				      *
				      * For example, if this object is
				      * refined, and one of its children is
				      * refined exactly one more time, then
				      * <tt>max_refinement_depth</tt> should
				      * return 2.
				      *
				      * If this object is not refined (i.e. it
				      * is active), then the return value is
				      * zero.
				      */
    unsigned int max_refinement_depth () const;    
    
				     /**
				      * Boundary indicator of this
				      * object.  If the return value
				      * is 255, then this line is in
				      * the interior of the domain.
				      *
				      * Since boundary data is only
				      * useful for structures with a
				      * dimension less than the
				      * dimension of a cell, this
				      * function issues an error if
				      * <tt>dim<4</tt>.
				      */
    unsigned char boundary_indicator () const;

				     /**
				      * Set the boundary indicator.
				      * The same applies as for the
				      * <tt>boundary_indicator()</tt>
				      * function.
				      *
				      * Caution: Never set the
				      * boundary indicator to 255,
				      * unless you exactly know what
				      * you are doing!  This value is
				      * reserved for another purpose
				      * and algorithms may not work if
				      * boundary cells have this
				      * boundary indicator or if
				      * interior cells have boundary
				      * indicators other than 255.
				      */
    void set_boundary_indicator (const unsigned char) const;

				     /**
				      * Return whether this object is
				      * at the boundary. This is
				      * checked via the boundary
				      * indicator field, which is
				      * always 255 if the object is in
				      * the interior of the
				      * domain. Obviously, the use of
				      * this function is only possible
				      * for <tt>dim>celldim</tt>; however,
				      * for <tt>dim==celldim</tt>, an
				      * object is a cell and the
				      * CellAccessor class
				      * offers another possibility to
				      * determine whether a cell is at
				      * the boundary or not.
				      */
    bool at_boundary () const;

				     /**
				      * Diameter of the object.
				      *
				      * The diameter of an object is
				      * computed to be the largest
				      * diagonal. This is not
				      * necessarily the true diameter,
				      * but completely sufficient for
				      * computations.
				      */
    double diameter () const;

    				     /**
				      * Center of the object. The
				      * center of an object is defined
				      * to be the average of the
				      * vertices, which is also where
				      * the (dim-)linear mapping
				      * places the midpoint of the
				      * unit cell in real space.
				      * However, this may not be the
				      * barycenter of the object.
				      */
    Point<dim> center () const;

				     /**
				      * Barycenter of the object.
				      */
    Point<dim> barycenter () const;

				     /**
				      * Volume of the object.  Here,
				      * the volume is defined to be
				      * confined by the (dim-)linear
				      * mapping of the unit cell.  No
				      * information about the boundary
				      * is used. If a more
				      * sophisticated computation is
				      * needed, try the volume of an
				      * appropriate finite element
				      * class.
				      */
    double measure () const;

				     /**
				      * Number of active descendants.
				      * This function only counts the
				      * number of active descendants,
				      * i.e. the number of descendants
				      * which are not further
				      * refined. Thus, if all of the
				      * eight children of a hex are
				      * further refined exactly once,
				      * the returned number will be
				      * 64, not 80.
				      *
				      * If the present cell is not refined,
				      * one is returned.
				      */
    unsigned int number_of_children () const;

                                     /**
                                      * Return whether the face with
                                      * index @p face has its normal
                                      * pointing in the standard
                                      * direction (@p true) or
                                      * whether it is the opposite
                                      * (@p false). Which is the
                                      * standard direction is
                                      * documented with the
                                      * GeometryInfo class. In
                                      * 1d and 2d, this is always
                                      * @p true, but in 3d it may be
                                      * different, see the respective
                                      * discussion in the
                                      * documentation of the
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool face_orientation (const unsigned int face) const;

                                     /**
                                      * Return whether the face with index @p
                                      * face is rotated by 180 degrees (@p true)
                                      * or or not (@p false). In 1d and 2d, this
                                      * is always @p false, but in 3d it may be
                                      * different, see the respective discussion
                                      * in the documentation of the
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool face_flip (const unsigned int face) const;

                                     /**
                                      * Return whether the face with index @p
                                      * face is rotated by 90 degrees (@p true)
                                      * or or not (@p false). In 1d and 2d, this
                                      * is always @p false, but in 3d it may be
                                      * different, see the respective discussion
                                      * in the documentation of the
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool face_rotation (const unsigned int face) const;
    
                                     /**
                                      * Return whether the line with index @p
                                      * line is oriented in standard
                                      * direction. @p true indicates, that the
                                      * line is oriented from vertex 0 to vertex
                                      * 1, whereas it is the other way around
                                      * otherwise. In 1d and 2d, this is always
                                      * @p true, but in 3d it may be different,
                                      * see the respective discussion in the
                                      * documentation of the 
                                      * GeometryInfo classe.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool line_orientation (const unsigned int line) const;
    
                                     /**
                                      * Set whether the quad with
                                      * index @p face has its normal
                                      * pointing in the standard
                                      * direction (@p true) or
                                      * whether it is the opposite
                                      * (@p false). Which is the
                                      * standard direction is
                                      * documented with the
                                      * GeometryInfo class.
                                      *
                                      * This function is only for
                                      * internal use in the
                                      * library. Setting this flag to
                                      * any other value than the one
                                      * that the triangulation has
                                      * already set is bound to bring
                                      * you desaster.
				      *
				      * It is only possible to set the
				      * face_orientation of cells in 3d
				      * (i.e. <code>celldim==3 &&
				      * dim==3</code>).
                                      */
    void set_face_orientation (const unsigned int face,
                               const bool         orientation) const;

				     /**
                                      * Set the flag indicating, what
                                      * <code>face_flip()</code> will
                                      * return.
				      *
				      * It is only possible to set the
				      * face_flip of cells in 3d
				      * (i.e. <code>celldim==3 &&
				      * dim==3</code>).
                                      */
    void set_face_flip (const unsigned int face,
			const bool         flip) const;

				     /**
                                      * Set the flag indicating, what
                                      * <code>face_rotation()</code> will
                                      * return.
				      *
				      * It is only possible to set the
				      * face_rotation of cells in 3d
				      * (i.e. <code>celldim==3 &&
				      * dim==3</code>).
                                      */
    void set_face_rotation (const unsigned int face,
			    const bool         rotation) const;

                                     /**
                                      * Set the flag indicating, what
                                      * <code>line_orientation()</code> will
                                      * return.
				      *
				      * It is only possible to set the
				      * line_orientation of faces in 3d
				      * (i.e. <code>celldim==2 &&
				      * dim==3</code>).
                                      */
    void set_line_orientation (const unsigned int line,
                               const bool         orientation) const;

  private:
    				     /**
				      *  Copy operator. This is normally
				      *  used in a context like
				      *  <tt>iterator a,b;  *a=*b;</tt>. Since
				      *  the meaning is to copy the object
				      *  pointed to by @p b to the object
				      *  pointed to by @p a and since
				      *  accessors are not real but
				      *  virtual objects, this operation
				      *  is not useful for iterators on
				      *  triangulations. We declare this
				      *  function here private, thus it may
				      *  not be used from outside.
				      *  Furthermore it is not implemented
				      *  and will give a linker error if
				      *  used anyway.
				      */
    void operator = (const TriaObjectAccessor<celldim, dim> &);

  protected:
				     /**@name Advancement of iterators*/
				     /*@{*/
				     /**
				      *  This operator advances the
				      *  iterator to the next element.
				      *
				      *  For cells only:
				      *  The next element is next on
				      *  this level if there are
				      *  more. If the present element
				      *  is the last on this level,
				      *  the first on the next level
				      *  is accessed.
				      */
    void operator ++ ();

    				     /**
				      *  This operator moves the
				      *  iterator to the previous
				      *  element.
				      *
				      *  For cells only:
				      *  The previous element is
				      *  previous on this level if
				      *  <tt>index>0</tt>. If the present
				      *  element is the first on this
				      *  level, the last on the
				      *  previous level is accessed.
				      */
    void operator -- ();

  private:
    
				     /**
				      * Dimension of the TriaObject
				      * this accessor gives access to.
				      */
    static const unsigned int objectdim = celldim;
				     /*@}*/

				     /**
				      * Declare some friends.
				      */
    template <int anydim, typename AnyAccessor> friend class TriaRawIterator;
};



/**
 * Closure class to stop induction of classes. Should never be called
 * and thus produces an error when created.
 *
 * @ingroup grid
 */
template<int dim>
class TriaObjectAccessor<0, dim> : public TriaAccessor<0,dim>
{
  public:
				     /**
				      * Propagate typedef from
				      * base class to this class.
				      */
    typedef typename TriaAccessor<0,dim>::AccessorData AccessorData;

				     /**
				      * Constructor. Should never be
				      * called and thus produces an
				      * error.
				      */
    TriaObjectAccessor (const Triangulation<dim> *parent     =  0,
			const int                 level      = -1,
			const int                 index      = -1,
			const AccessorData       *local_data =  0)
                    :
		    TriaAccessor<0,dim> (parent, level, index, local_data)
      {
	Assert (false, ExcInternalError());
      }
};



/**
 *   Accessor to dereference the data of lines. This accessor is used
 *   to point to lines in @p dim space dimensions. There is a derived
 *   class for lines in one space dimension, in which case a line is
 *   also a cell and thus has much more functionality than in other
 *   dimensions.
 *
 * @ingroup grid
 *   @author Wolfgang Bangerth, 1998
 */
template <int dim>
class TriaObjectAccessor<1, dim> :  public TriaAccessor<1,dim>
{
  public:
				     /**
				      * Propagate typedef from
				      * base class to this class.
				      */
    typedef typename TriaAccessor<1,dim>::AccessorData AccessorData;

				     /**
				      *  Constructor.
				      */
    TriaObjectAccessor (const Triangulation<dim> *parent     =  0,
			const int                 level      = -1,
			const int                 index      = -1,
			const AccessorData       *local_data =  0);

				     /**
				      *  Copy the data of the given
				      *  line.
				      */
    void set (const internal::Triangulation::TriaObject<1> &l) const;

				     /**
				      *  Return the index of vertex
				      *  <tt>i=0,1</tt> of a line.
				      *
				      *  Note that the returned value is only
				      *  the index of the geometrical
				      *  vertex. It has nothing to do with
				      *  possible degrees of freedom
				      *  associated with it. For this, see the
				      *  @p DoFAccessor::vertex_dof_index
				      *  functions.
				      */ 
    int vertex_index (const unsigned int i) const;

				     /**
				      *  Return a reference to the
				      *  @p ith vertex.
				      */
    Point<dim> & vertex (const unsigned int i) const;

				     /**
				      *  Test for the element being
				      *  used or not.  The return
				      *  value is @p true for all
				      *  iterators that are either
				      *  normal iterators or active
				      *  iterators, only raw iterators
				      *  can return @p false. Since
				      *  raw iterators are only used
				      *  in the interiors of the
				      *  library, you will not usually
				      *  need this function.
				      */
    bool used () const;

				     /**
				      *  Set the @p used flag. Only
				      *  for internal use in the
				      *  library.
				      */
    void set_used_flag () const;

				     /**
				      *  Clear the @p used flag. Only
				      *  for internal use in the
				      *  library.
				      */
    void clear_used_flag () const;

    				     /**
				      *  Return whether the user flag
				      *  is set or not.
				      */
    bool user_flag_set () const;

				     /**
				      *  Flag the user flag for this cell.
				      */
    void set_user_flag () const;

				     /**
				      *  Clear the user flag.
				      */
    void clear_user_flag () const;

				     /**
				      *  Set the user flag of this
				      *  object and of all its
				      *  children and their children,
				      *  etc.
				      */
    void recursively_set_user_flag () const;

    				     /**
				      *  Clear the user flag of this
				      *  object and of all its
				      *  children and their children,
				      *  etc.
				      */
    void recursively_clear_user_flag () const;

				     /**
				      * Reset the user data to zero,
				      * independent if pointer or index.
				      */
    void clear_user_data () const;
    
				     /**
				      * Set the user pointer of this
				      * line to @p p.
				      */
    void set_user_pointer (void *p) const;

				     /**
				      * Set the user pointer of this
				      * object and all its children to
				      * the given value. This is
				      * useful for example if all
				      * cells of a certain subdomain,
				      * or all faces of a certain part
				      * of the boundary should have
				      * user pointers pointing to
				      * objects describing this part
				      * of the domain or boundary.
				      *
				      * Note that the user pointer is
				      * not inherited under mesh
				      * refinement, so after mesh
				      * refinement there might be
				      * cells or faces that don't have
				      * user pointers pointing to the
				      * describing object. In this
				      * case, simply loop over all the
				      * elements of the coarsest level
				      * that has this information, and
				      * use this function to
				      * recursively set the user
				      * pointer of all finer levels of
				      * the triangulation.
				      */
    void recursively_set_user_pointer (void *p) const;

				     /**
				      * Clear the user pointer of this
				      * object and all of its
				      * descendants. The same holds as
				      * said for the
				      * recursively_set_user_pointer()
				      * function.
				      */
    void recursively_clear_user_pointer () const;

				     /**
				      * Reset the user pointer of this
				      * line to a @p NULL pointer.
				      */
    void clear_user_pointer () const;

				     /**
				      * Access the value of the user
				      * pointer of this line. It is in
				      * the responsibility of the user
				      * to make sure that the pointer
				      * points to something
				      * useful. You should use the new
				      * style cast operator to
				      * maintain a minimum of
				      * typesafety, e.g.
				      * <tt>A *a=static_cast<A*>(cell->user_pointer());</tt>.
				      */
    void * user_pointer () const;
    
				     /**
				      * Set the user index
				      * to @p p.
				      */
    void set_user_index (unsigned int p) const;
    
				     /**
				      * Reset the user index to 0.
				      */
    void clear_user_index () const;
    
				     /**
				      * Access the value of the user
				      * index.
				      */
    unsigned int user_index () const;

				     /**
				      * Set the user index of this
				      * object and all its children.
				      *
				      * Note that the user index is
				      * not inherited under mesh
				      * refinement, so after mesh
				      * refinement there might be
				      * cells or faces that don't have
				      * the expected user indices. In
				      * this case, simply loop over
				      * all the elements of the
				      * coarsest level that has this
				      * information, and use this
				      * function to recursively set
				      * the user index of all finer
				      * levels of the triangulation.
				      */
    void recursively_set_user_index (unsigned int p) const;

				     /**
				      * Clear the user index of this
				      * object and all of its
				      * descendants. The same holds as
				      * said for the
				      * recursively_set_user_index()
				      * function.
				      */
    void recursively_clear_user_index () const;

				     /**
				      *  Return a pointer to the
				      *  @p ith child.
				      */
    TriaIterator<dim,TriaObjectAccessor<1, dim> >
    child (const unsigned int i) const;

				     /**
				      *  Return the index of the
				      *  @p ith child.  The level of
				      *  the child is one higher than
				      *  that of the present cell if
				      *  @p dim=1 and lines are cells.
				      *  Otherwise the line and its children
				      *  have no level.
				      *  If the child does not exist, -1
				      *  is returned.
				      */
    int child_index (const unsigned int i) const;

				     /**
				      *  Set the child field. Since we
				      *  only store the index of the
				      *  first child (the others
				      *  follow directly) only one
				      *  child index is to be
				      *  given. The level of the child
				      *  is one level up of the level
				      *  of the cell to which this
				      *  iterator points if
				      *  @p dim=1 and lines are cells.
				      *  Otherwise the line and its children
				      *  have no level.
				      */
    void set_children (const int index) const;
	
				     /**
				      *  Clear the child field,
				      *  i.e. set it to a value which
				      *  indicates that this cell has
				      *  no children.
				      */
    void clear_children () const;

				     /**
				      *  Test whether the line has
				      *  children.
				      */
    bool has_children () const;

				     /**
				      * Return the number of immediate
				      * children of this object.
				      *
				      * An assertion assures that this
				      * function can only be called
				      * for TriaObjects with
				      * <tt>has_children()==true</tt>.
				      */
    unsigned int n_children() const;

				     /**
				      * Return the number of times
				      * that this cell is
				      * refined. Note that not all its
				      * children are refined that
				      * often (which is why we prepend
				      * @p max_), the returned number
				      * is rather the maximum number
				      * of refinement in any branch of
				      * children of this object.
				      *
				      * For example, if this object is
				      * refined, and one of its children is
				      * refined exactly one more time, then
				      * <tt>max_refinement_depth</tt> should
				      * return 2.
				      *
				      * If this object is not refined (i.e. it
				      * is active), then the return value is
				      * zero.
				      */
    unsigned int max_refinement_depth () const;
    
				     /**
				      * Return the boundary indicator
				      * of this line. Since boundary
				      * data is only useful for
				      * structures with a dimension
				      * less than the dimension of a
				      * cell, this function issues an
				      * error if <tt>dim<2</tt>.
				      *
				      * If the return value is 255,
				      * then this line is in the
				      * interior of the domain.
				      */
    unsigned char boundary_indicator () const;

				     /**
				      * Set the boundary indicator of
				      * this line.  The same applies
				      * as for the
				      * <tt>boundary_indicator()</tt>
				      * function.
				      *
				      * You should be careful with
				      * this function and especially
				      * never try to set the boundary
				      * indicator to 255, unless you
				      * exactly know what you are
				      * doing, since this value is
				      * reserved for another purpose
				      * and algorithms may not work if
				      * boundary cells have this
				      * boundary indicator or if
				      * interior cells have boundary
				      * indicators other than 255.
				      */
    void set_boundary_indicator (const unsigned char) const;

				     /**
				      * Return whether this line is at
				      * the boundary. This is checked
				      * via the the boundary indicator
				      * field, which is always 255 if
				      * the line is in the interior of
				      * the domain. Obviously, this is
				      * only possible for <tt>dim>1</tt>;
				      * however, for <tt>dim==1</tt>, a
				      * line is a cell and the
				      * CellAccessor class
				      * offers another possibility to
				      * determine whether a cell is at
				      * the boundary or not.
				      */
    bool at_boundary () const;

    				     /**
				      * Return the length of the
				      * line. If the line describes
				      * part of the boundary (e.g. if
				      * it is face to a cell in 2D)
				      * and is not a straight one, ask
				      * the finite element class for
				      * the correct length!
				      */
    double diameter () const;

				     /**
				      * Return the center of the
				      * line. This is the average of
				      * the two vertices, which is the
				      * obvious definition for
				      * straight lines. However, if
				      * you use higher order mappings
				      * from the unit cell to the real
				      * cell (in more than one space
				      * dimension), the bounding lines
				      * may not necessarily be
				      * straight.  In that case ask
				      * the finite element class for
				      * the correct place of the
				      * midpoint of the line in real
				      * space.
				      */
    Point<dim> center () const;

				     /**
				      * Return the barycenter of the
				      * line, which is the
				      * midpoint. The same applies as
				      * for the @p center function
				      * with regard to lines at the
				      * boundary.
				      */
    Point<dim> barycenter () const;
    
				     /**
				      * Return the length of the line.
				      * The same applies as for the
				      * @p center function with
				      * regard to lines at the
				      * boundary.
				      */
    double measure () const;

				     /**
				      * Compute and return the number
				      * of children of this
				      * line. Actually, this function
				      * only counts the number of
				      * active children, i.e. the
				      * number if lines which are not
				      * further refined. Thus, if both
				      * of the two children of a line
				      * are further refined exactly
				      * once, the returned number will
				      * be four, not six.
				      *
				      * If the present cell is not
				      * refined, one is returned.
				      */
    unsigned int number_of_children () const;
    
                                     /**
                                      * Return whether the face with
                                      * index @p face has its normal
                                      * pointing in the standard
                                      * direction (@p true) or
                                      * whether it is the opposite
                                      * (@p false). Which is the
                                      * standard direction is
                                      * documented with the
                                      * GeometryInfo class. In
                                      * 1d and 2d, this is always
                                      * @p true, but in 3d it may be
                                      * different, see the respective
                                      * discussion in the
                                      * documentation of the
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool face_orientation (const unsigned int face) const;

                                     /**
                                      * Return whether the face with index @p
                                      * face is rotated by 180 degrees (@p true)
                                      * or or not (@p false). In 1d and 2d, this
                                      * is always @p false, but in 3d it may be
                                      * different, see the respective discussion
                                      * in the documentation of the
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool face_flip (const unsigned int face) const;

                                     /**
                                      * Return whether the face with index @p
                                      * face is rotated by 90 degrees (@p true)
                                      * or or not (@p false). In 1d and 2d, this
                                      * is always @p false, but in 3d it may be
                                      * different, see the respective discussion
                                      * in the documentation of the
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool face_rotation (const unsigned int face) const;
    
                                     /**
                                      * Return whether the line with index @p
                                      * line is oriented in standard
                                      * direction. @p true indicates, that the
                                      * line is oriented from vertex 0 to vertex
                                      * 1, whereas it is the other way around
                                      * otherwise. In 1d and 2d, this is always
                                      * @p true, but in 3d it may be different,
                                      * see the respective discussion in the
                                      * documentation of the 
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool line_orientation (const unsigned int line) const;
    
                                     /**
                                      * Set the flag indicating, what
                                      * <code>line_orientation()</code> will
                                      * return.
				      *
				      * It is only possible to set the
				      * line_orientation of faces in 3d
				      * (i.e. <code>celldim==2 &&
				      * dim==3</code>).
                                      */
    void set_line_orientation (const unsigned int line,
                               const bool         orientation) const;
    
  private:
    				     /**
				      *  Copy operator. This is normally
				      *  used in a context like
				      *  <tt>iterator a,b;  *a=*b;</tt>. Since
				      *  the meaning is to copy the object
				      *  pointed to by @p b to the object
				      *  pointed to by @p a and since
				      *  accessors are not real but
				      *  virtual objects, this operation
				      *  is not useful for iterators on
				      *  triangulations. We declare this
				      *  function here private, thus it may
				      *  not be used from outside.
				      *  Furthermore it is not implemented
				      *  and will give a linker error if
				      *  used anyway.
				      */
    void operator = (const TriaObjectAccessor<1, dim> &);

  protected:
    
				     /**@name Advancement of iterators*/
				     /*@{*/
				     /**
				      *  This operator advances the
				      *  iterator to the next element.
				      *
				      *  For @p dim=1 only:
				      *  The next element is next on
				      *  this level if there are
				      *  more. If the present element
				      *  is the last on this level,
				      *  the first on the next level
				      *  is accessed.
				      */
    void operator ++ ();

    				     /**
				      *  This operator moves the
				      *  iterator to the previous
				      *  element.
				      *
				      *  For @p dim=1 only:
				      *  The previous element is
				      *  previous on this level if
				      *  <tt>index>0</tt>. If the present
				      *  element is the first on this
				      *  level, the last on the
				      *  previous level is accessed.
				      */
    void operator -- ();

				     /**
				      * Access to the other objects of
				      * a Triangulation with same
				      * dimension. This function is
				      * specialized for dim==1 because
				      * lines are stored in
				      * TriaLevel<1> for dim==1 but in
				      * TriaFace<dim> for dim!=1. All
				      * accessor functions using this
				      * function can be implemented
				      * for all dim without
				      * specialization.
				      */
    internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<1> > & objects() const;

  private:
    
				     /**
				      * Dimension of the TriaObject
				      * this accessor gives access to.
				      */
    static const unsigned int objectdim = 1;

				     /*@}*/

				     /**
				      * Declare some friends.
				      */
    template <int anydim, typename AnyAccessor> friend class TriaRawIterator;
};



/**
 *   Accessor to dereference the data of quads. This accessor is used
 *   to point to quads in @p dim space dimensions (only <tt>dim>=2</tt>
 *   seems reasonable to me). There is a derived class for quads in
 *   two space dimension, in which case a quad is also a cell and thus
 *   has much more functionality than in other dimensions.
 *
 * @ingroup grid
 *   @author Wolfgang Bangerth, 1998
 */
template <int dim>
class TriaObjectAccessor<2, dim> :  public TriaAccessor<2,dim>
{
  public:
				     /**
				      * Propagate typedef from base
				      * class to this class.
				      */
    typedef typename TriaAccessor<2,dim>::AccessorData AccessorData;

				     /**
				      *  Constructor.
				      */
    TriaObjectAccessor (const Triangulation<dim> *parent     =  0,
			const int                 level      = -1,
			const int                 index      = -1,
			const AccessorData       *local_data =  0);

				     /**
				      *  Copy the data of the given quad.
				      */
    void set (const internal::Triangulation::TriaObject<2> &q) const;

				     /**
				      * Return index of a vertex of a
				      * quad in the internal data
				      * structure of a Triangulation.
				      *
				      * For local numbering of the
				      * vertices, see the page on
				      * Geometry.
				      *
				      * @attention The value returned
				      * is only the index of the
				      * geometrical vertex. It has
				      * nothing to do with possible
				      * degrees of freedom associated
				      * with it. For this, see the
				      * DoFAccessor::vertex_dof_index
				      * functions.
				      */ 
    int vertex_index (const unsigned int i) const;

    				     /**
				      *  Return a reference to the
				      *  @p ith vertex.
				      */
    Point<dim> & vertex (const unsigned int i) const;

				     /**
				      *  Return a pointer to the
				      *  @p ith line bounding this
				      *  Quad.
				      */
    TriaIterator<dim,TriaObjectAccessor<1, dim> >
    line (const unsigned int i) const;

				     /**
				      * Return the line index of the
				      * @p ith side (a line). The
				      * line has no level.
				      */
    unsigned int line_index (const unsigned int i) const;
    
				     /**
				      *  Test for the element being
				      *  used or not.  The return
				      *  value is @p true for all
				      *  iterators that are either
				      *  normal iterators or active
				      *  iterators, only raw iterators
				      *  can return @p false. Since
				      *  raw iterators are only used
				      *  in the interiors of the
				      *  library, you will not usually
				      *  need this function.
				      */
    bool used () const;

				     /**
				      *  Set the @p used flag. You
				      *  should know quite exactly
				      *  what you are doing of you
				      *  touch this function. It is
				      *  exclusively for internal use
				      *  in the library.
				      */
    void set_used_flag () const;

				     /**
				      *  Clear the @p used flag. You
				      *  should know quite exactly
				      *  what you are doing of you
				      *  touch this function. It is
				      *  exclusively for internal use
				      *  in the library.
				      */
    void clear_used_flag () const;

    				     /**
				      *  Return whether the user flag
				      *  is set or not.
				      */
    bool user_flag_set () const;

				     /**
				      *  Flag the user flag for this
				      *  cell.
				      */
    void set_user_flag () const;

				     /**
				      *  Clear the user flag.
				      */
    void clear_user_flag () const;

				     /**
				      *  Set the user flag of this
				      *  object and of all its
				      *  children and their children,
				      *  etc.
				      */
    void recursively_set_user_flag () const;

    				     /**
				      *  Clear the user flag of this
				      *  object and of all its
				      *  children and their children,
				      *  etc.
				      */
    void recursively_clear_user_flag () const;

				     /**
				      * Reset the user data to zero,
				      * independent if pointer or index.
				      */
    void clear_user_data () const;
    
				     /**
				      * Set the user pointer of this
				      * quad to @p p.
				      */
    void set_user_pointer (void *p) const;

				     /**
				      * Reset the user pointer of this
				      * quad to a @p NULL pointer.
				      */
    void clear_user_pointer () const;

				     /**
				      * Access the value of the user
				      * pointer of this quad. It is in
				      * the responsibility of the user
				      * to make sure that the pointer
				      * points to something
				      * useful. You should use the new
				      * style cast operator to
				      * maintain a minimum of
				      * typesafety, e.g.
				      * <tt>A *a=static_cast<A*>(cell->user_pointer());</tt>.
				      */
    void * user_pointer () const;

				     /**
				      * Set the user pointer of this
				      * object and all its children to
				      * the given value. This is
				      * useful for example if all
				      * cells of a certain subdomain,
				      * or all faces of a certain part
				      * of the boundary should have
				      * user pointers pointing to
				      * objects describing this part
				      * of the domain or boundary.
				      *
				      * Note that the user pointer is
				      * not inherited under mesh
				      * refinement, so after mesh
				      * refinement there might be
				      * cells or faces that don't have
				      * user pointers pointing to the
				      * describing object. In this
				      * case, simply loop over all the
				      * elements of the coarsest level
				      * that has this information, and
				      * use this function to
				      * recursively set the user
				      * pointer of all finer levels of
				      * the triangulation.
				      */
    void recursively_set_user_pointer (void *p) const;

				     /**
				      * Clear the user pointer of this
				      * object and all of its
				      * descendants. The same holds as
				      * said for the
				      * recursively_set_user_pointer()
				      * function.
				      */
    void recursively_clear_user_pointer () const;

				     /**
				      * Set the user index
				      * to @p p.
				      */
    void set_user_index (unsigned int p) const;
    
				     /**
				      * Reset the user index to 0.
				      */
    void clear_user_index () const;
    
				     /**
				      * Access the value of the user
				      * index.
				      */
    unsigned int user_index () const;

				     /**
				      * Set the user index of this
				      * object and all its children.
				      *
				      * Note that the user index is
				      * not inherited under mesh
				      * refinement, so after mesh
				      * refinement there might be
				      * cells or faces that don't have
				      * the expected user indices. In
				      * this case, simply loop over
				      * all the elements of the
				      * coarsest level that has this
				      * information, and use this
				      * function to recursively set
				      * the user index of all finer
				      * levels of the triangulation.
				      */
    void recursively_set_user_index (unsigned int p) const;

				     /**
				      * Clear the user index of this
				      * object and all of its
				      * descendants. The same holds as
				      * said for the
				      * recursively_set_user_index()
				      * function.
				      */
    void recursively_clear_user_index () const;
				     /**
				      *  Return a pointer to the @p ith
				      *  child.
				      */
    TriaIterator<dim,TriaObjectAccessor<2, dim> > child (const unsigned int i) const;

				     /**
				      *  Return the index of the
				      *  @p ith child.  The level of
				      *  the child is one higher than
				      *  that of the present cell if
				      *  @p dim=2 and quads are cells.
				      *  Otherwise the quad and its children
				      *  have no level.
				      *  If the child does not exist, -1
				      *  is returned.
				      */
    int child_index (const unsigned int i) const;

				     /**
				      *  Set the child field. Since we
				      *  only store the index of the
				      *  first child (the others
				      *  follow directly) only one
				      *  child index is to be
				      *  given. The level of the child
				      *  is one level up of the level
				      *  of the cell to which this
				      *  iterator points.
				      */
    void set_children (const int index) const;
	
				     /**
				      *  Clear the child field,
				      *  i.e. set it to a value which
				      *  indicates that this cell has
				      *  no children.
				      */
    void clear_children () const;
    
				     /**
				      *  Test whether the quad has
				      *  children.
				      */
    bool has_children () const;

				     /**
				      * Return the number of immediate
				      * children of this object.
				      *
				      * An assertion assures that this
				      * function can only be called
				      * for TriaObjects with
				      * <tt>has_children()==true</tt>.
				      */
    unsigned int n_children() const;

				     /**
				      * Return the number of times
				      * that this cell is
				      * refined. Note that not all its
				      * children are refined that
				      * often (which is why we prepend
				      * @p max_), the returned number
				      * is rather the maximum number
				      * of refinement in any branch of
				      * children of this object.
				      *
				      * For example, if this object is
				      * refined, and one of its children is
				      * refined exactly one more time, then
				      * <tt>max_refinement_depth</tt> should
				      * return 2.
				      *
				      * If this object is not refined (i.e. it
				      * is active), then the return value is
				      * zero.
				      */
    unsigned int max_refinement_depth () const;
    
				     /**
				      * Return the boundary indicator
				      * of this quad. Since boundary
				      * data is only useful for
				      * structures with a dimension
				      * less than the dimension of a
				      * cell, this function issues an
				      * error if <tt>dim<3</tt>.
				      *
				      * If the return value is 255,
				      * then this quad is in the
				      * interior of the domain.
				      */
    unsigned char boundary_indicator () const;

				     /**
				      * Set the boundary indicator of
				      * this quad.  The same applies
				      * as for the
				      * <tt>boundary_indicator()</tt>
				      * function.
				      *
				      * You should be careful with
				      * this function and especially
				      * never try to set the boundary
				      * indicator to 255, unless you
				      * exactly know what you are
				      * doing, since this value is
				      * reserved for another purpose
				      * and algorithms may not work if
				      * boundary cells have this
				      * boundary indicator or if
				      * interior cells have boundary
				      * indicators other than 255.
				      */
    void set_boundary_indicator (const unsigned char) const;

				     /**
				      * Return whether this quad is at
				      * the boundary. This is checked
				      * via the the boundary indicator
				      * field, which is always 255 if
				      * the quad is in the interior of
				      * the domain. Obviously, this
				      * function is only useful for
				      * <tt>dim>2</tt>; however, for
				      * <tt>dim==2</tt>, a quad is a cell
				      * and the CellAccessor
				      * class offers another
				      * possibility to determine
				      * whether a cell is at the
				      * boundary or not.
				      */
    bool at_boundary () const;
    
				     /**
				      * Return the diameter of the
				      * quad. If the quad describes
				      * part of the boundary (e.g. if
				      * it is face to a cell in 3D)
				      * and is not a plane one, ask
				      * the finite element class for
				      * the correct length!
				      *
				      * The diameter of a quad is
				      * computed to be the larger of
				      * the two diagonals. You should
				      * absolutely be clear about the
				      * fact that this definitely is
				      * not the diameter of all
				      * quadrilaterals; however it may
				      * serve as an approximation and
				      * is exact in many cases,
				      * especially if the
				      * quadrilateral is not too much
				      * distorted.
				      */
    double diameter () const;

    				     /**
				      * Return the center of the
				      * quad. The center of a quad is
				      * defined to be the average of
				      * the four vertices, which is
				      * also the point where the
				      * bilinear mapping places the
				      * midpoint of the unit quad in
				      * real space.  However, this may
				      * not be the point of the
				      * barycenter of the quad.
				      *
				      * Also note that if you use
				      * higher order mappings from the
				      * unit cell to the real cell (in
				      * more than two space
				      * dimension), the bounding quads
				      * may not necessarily be
				      * straight.  In that case ask
				      * the finite element class for
				      * the correct place of the
				      * midpoint of the quad in real
				      * space.
				      */
    Point<dim> center () const;

				     /**
				      * Return the barycenter of the
				      * qaud. The same applies as for
				      * the @p center function with
				      * regard to quads at the
				      * boundary.
				      */
    Point<dim> barycenter () const;

				     /**
				      * Return the area of the
				      * quad. With area, we mean the
				      * area included by the straight
				      * lines connecting the four
				      * vertices, i.e. no information
				      * about the boundary is used; if
				      * you want other than bilinearly
				      * mapped unit quadrilaterals,
				      * ask the appropriate finite
				      * element class for the area of
				      * this quad.
				      */
    double measure () const;

				     /**
				      * Compute and return the number
				      * of children of this
				      * quad. Actually, this function
				      * only counts the number of
				      * active children, i.e. the
				      * number if quads which are not
				      * further refined. Thus, if all
				      * of the four children of a quad
				      * are further refined exactly
				      * once, the returned number will
				      * be 16, not 20.
				      *
				      * If the present cell is not
				      * refined, one is returned.
				      */
    unsigned int number_of_children () const;

                                     /**
                                      * Return whether the face with
                                      * index @p face has its normal
                                      * pointing in the standard
                                      * direction (@p true) or
                                      * whether it is the opposite
                                      * (@p false). Which is the
                                      * standard direction is
                                      * documented with the
                                      * GeometryInfo class. In
                                      * 1d and 2d, this is always
                                      * @p true, but in 3d it may be
                                      * different, see the respective
                                      * discussion in the
                                      * documentation of the
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool face_orientation (const unsigned int face) const;

                                     /**
                                      * Return whether the face with index @p
                                      * face is rotated by 180 degrees (@p true)
                                      * or or not (@p false). In 1d and 2d, this
                                      * is always @p false, but in 3d it may be
                                      * different, see the respective discussion
                                      * in the documentation of the
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool face_flip (const unsigned int face) const;

                                     /**
                                      * Return whether the face with index @p
                                      * face is rotated by 90 degrees (@p true)
                                      * or or not (@p false). In 1d and 2d, this
                                      * is always @p false, but in 3d it may be
                                      * different, see the respective discussion
                                      * in the documentation of the
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool face_rotation (const unsigned int face) const;
    
                                     /**
                                      * Return whether the line with index @p
                                      * line is oriented in standard
                                      * direction. @p true indicates, that the
                                      * line is oriented from vertex 0 to vertex
                                      * 1, whereas it is the other way around
                                      * otherwise. In 1d and 2d, this is always
                                      * @p true, but in 3d it may be different,
                                      * see the respective discussion in the
                                      * documentation of the 
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool line_orientation (const unsigned int line) const;
    
                                     /**
                                      * Set the flag indicating, what
                                      * <code>line_orientation()</code> will
                                      * return.
				      *
				      * It is only possible to set the
				      * line_orientation of faces in 3d
				      * (i.e. <code>celldim==2 &&
				      * dim==3</code>).
                                      */
    void set_line_orientation (const unsigned int line,
                               const bool         orientation) const;

  private:
    				     /**
				      * Copy operator. This is
				      * normally used in a context
				      * like <tt>iterator a,b;
				      * *a=*b;</tt>. Since the meaning
				      * is to copy the object pointed
				      * to by @p b to the object
				      * pointed to by @p a and since
				      * accessors are not real but
				      * virtual objects, this
				      * operation is not useful for
				      * iterators on
				      * triangulations. We declare
				      * this function here private,
				      * thus it may not be used from
				      * outside.  Furthermore it is
				      * not implemented and will give
				      * a linker error if used anyway.
				      */
    void operator = (const TriaObjectAccessor<2, dim> &);

  protected:
				     /**@name Advancement of iterators*/
				     /*@{*/
				     /**
				      *  This operator advances the
				      *  iterator to the next element.
				      *
				      *  For @p dim=2 only:
				      *  The next element is next on
				      *  this level if there are
				      *  more. If the present element
				      *  is the last on this level,
				      *  the first on the next level
				      *  is accessed.
				      */
    void operator ++ ();

    				     /**
				      *  This operator moves the
				      *  iterator to the previous
				      *  element.
				      *
				      *  For @p dim=2 only:
				      *  The previous element is
				      *  previous on this level if
				      *  <tt>index>0</tt>. If the present
				      *  element is the first on this
				      *  level, the last on the
				      *  previous level is accessed.
				      */
    void operator -- ();

				     /**
				      * Access to the other objects of
				      * a Triangulation of same
				      * dimension. This function is
				      * specialized for dim==2 and
				      * dim==3 because quads are
				      * stored in TriaLevel<2> for
				      * dim==2 but in TriaFace<3> for
				      * dim==3. All accessor functions
				      * using this function can be
				      * implemented for all dimensions
				      * without specialization.
				      */
    internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<2> > & objects() const;

  private:
    
				     /**
				      * Dimension of the TriaObject
				      * this accessor gives access to.
				      */
    static const unsigned int objectdim = 2;
				     /*@}*/

				     /**
				      * Declare some friends.
				      */
    template <int anydim, typename AnyAccessor> friend class TriaRawIterator;
};



/**
 *   Accessor to dereference the data of hexahedra. This accessor is
 *   used to point to hexs in @p dim space dimensions (only
 *   <tt>dim>=3</tt> seems reasonable to me). There is a derived class for
 *   hexs in three space dimension, in which case a hex is also a cell
 *   and thus has much more functionality than in other dimensions.
 *
 * @ingroup grid
 *   @author Wolfgang Bangerth, 1998
 */
template <int dim>
class TriaObjectAccessor<3, dim> :  public TriaAccessor<3,dim>
{
  public:
				     /**
				      * Propagate typedef from base
				      * class to this class.
				      */
    typedef typename TriaAccessor<3,dim>::AccessorData AccessorData;

				     /**
				      *  Constructor.
				      */
    TriaObjectAccessor (const Triangulation<dim> *parent     =  0,
			const int                 level      = -1,
			const int                 index      = -1,
			const AccessorData       *local_data =  0);

				     /**
				      *  Copy the data of the given
				      *  hex.
				      */
    void set (const internal::Triangulation::TriaObject<3> &h) const;
    
				     /**
				      *  Return index of a vertex of a hex in the internal structures of Triangulation.
				      *
				      * For local numbering of the
				      * vertices, see the page on
				      * Geometry.
				      *
				      *
				      *  Note that the returned value is only
				      *  the index of the geometrical
				      *  vertex. It has nothing to do with
				      *  possible degrees of freedom
				      *  associated with it. For this, see the
				      *  @p DoFAccessor::vertex_dof_index
				      *  functions.
				      */ 
    int vertex_index (const unsigned int i) const;

    				     /**
				      *  Return a reference to the
				      *  @p ith vertex.
				      */
    Point<dim> & vertex (const unsigned int i) const;

				     /**
				      *  Return a pointer to the
				      *  @p ith line bounding this
				      *  @p Hex.
				      */
    TriaIterator<dim,TriaObjectAccessor<1, dim> >
    line (const unsigned int i) const;

				     /**
				      * Return the line index of the
				      * @p ith line. The line has no
				      * level.
				      */
    unsigned int line_index (const unsigned int i) const;
    
    				     /**
				      *  Return a pointer to the
				      *  @p ith quad bounding this
				      *  @p Hex.
				      */
    TriaIterator<dim,TriaObjectAccessor<2, dim> >
    quad (const unsigned int i) const;

				     /**
				      * Return the quad index of the
				      * @p ith side (a quad). The quad
				      * has no level.
				      */
    unsigned int quad_index (const unsigned int i) const;

				     /**
				      *  Test for the element being
				      *  used or not.  The return
				      *  value is @p true for all
				      *  iterators that are either
				      *  normal iterators or active
				      *  iterators, only raw iterators
				      *  can return @p false. Since
				      *  raw iterators are only used
				      *  in the interiors of the
				      *  library, you will not usually
				      *  need this function.
				      */
    bool used () const;

				     /**
				      *  Set the @p used flag. You
				      *  should know quite exactly
				      *  what you are doing of you
				      *  touch this function. It is
				      *  exclusively for internal use
				      *  in the library.
				      */
    void set_used_flag () const;

				     /**
				      *  Clear the @p used flag. You
				      *  should know quite exactly
				      *  what you are doing of you
				      *  touch this function. It is
				      *  exclusively for internal use
				      *  in the library.
				      */
    void clear_used_flag () const;

    				     /**
				      *  Return whether the user flag
				      *  is set or not.
				      */
    bool user_flag_set () const;

				     /**
				      *  Flag the user flag for this
				      *  cell.
				      */
    void set_user_flag () const;

				     /**
				      *  Clear the user flag.
				      */
    void clear_user_flag () const;

				     /**
				      *  Set the user flag of this
				      *  object and of all its
				      *  children and their children,
				      *  etc.
				      */
    void recursively_set_user_flag () const;

    				     /**
				      *  Clear the user flag of this
				      *  object and of all its
				      *  children and their children,
				      *  etc.
				      */
    void recursively_clear_user_flag () const;

				     /**
				      * Reset the user data to zero,
				      * independent if pointer or index.
				      */
    void clear_user_data () const;
    
				     /**
				      * Set the user pointer of this
				      * hex to @p p.
				      */
    void set_user_pointer (void *p) const;

				     /**
				      * Reset the user pointer of this
				      * hex to a @p NULL pointer.
				      */
    void clear_user_pointer () const;

				     /**
				      * Access the value of the user
				      * pointer of this hex. It is in
				      * the responsibility of the user
				      * to make sure that the pointer
				      * points to something
				      * useful. You should use the new
				      * style cast operator to
				      * maintain a minimum of
				      * typesafety, e.g.
				      * <tt>A *a=static_cast<A*>(cell->user_pointer());</tt>.
				      */
    void * user_pointer () const;

				     /**
				      * Set the user pointer of this
				      * object and all its children to
				      * the given value. This is
				      * useful for example if all
				      * cells of a certain subdomain,
				      * or all faces of a certain part
				      * of the boundary should have
				      * user pointers pointing to
				      * objects describing this part
				      * of the domain or boundary.
				      *
				      * Note that the user pointer is
				      * not inherited under mesh
				      * refinement, so after mesh
				      * refinement there might be
				      * cells or faces that don't have
				      * user pointers pointing to the
				      * describing object. In this
				      * case, simply loop over all the
				      * elements of the coarsest level
				      * that has this information, and
				      * use this function to
				      * recursively set the user
				      * pointer of all finer levels of
				      * the triangulation.
				      */
    void recursively_set_user_pointer (void *p) const;

				     /**
				      * Clear the user pointer of this
				      * object and all of its
				      * descendants. The same holds as
				      * said for the
				      * recursively_set_user_pointer()
				      * function.
				      */
    void recursively_clear_user_pointer () const;

				     /**
				      * Set the user index
				      * to @p p.
				      */
    void set_user_index (unsigned int p) const;
    
				     /**
				      * Reset the user index to 0.
				      */
    void clear_user_index () const;
    
				     /**
				      * Access the value of the user
				      * index.
				      */
    unsigned int user_index () const;

				     /**
				      * Set the user index of this
				      * object and all its children.
				      *
				      * Note that the user index is
				      * not inherited under mesh
				      * refinement, so after mesh
				      * refinement there might be
				      * cells or faces that don't have
				      * the expected user indices. In
				      * this case, simply loop over
				      * all the elements of the
				      * coarsest level that has this
				      * information, and use this
				      * function to recursively set
				      * the user index of all finer
				      * levels of the triangulation.
				      */
    void recursively_set_user_index (unsigned int p) const;

				     /**
				      * Clear the user index of this
				      * object and all of its
				      * descendants. The same holds as
				      * said for the
				      * recursively_set_user_index()
				      * function.
				      */
    void recursively_clear_user_index () const;
    
				     /**
				      *  Return a pointer to the
				      *  @p ith child.
				      */
    TriaIterator<dim,TriaObjectAccessor<3, dim> >
    child (const unsigned int i) const;

				     /**
				      *  Return the index of the
				      *  @p ith child.  The level of
				      *  the child is one higher than
				      *  that of the present cell.  If
				      *  the child does not exist, -1
				      *  is returned.
				      */
    int child_index (const unsigned int i) const;

				     /**
				      *  Set the child field. Since we
				      *  only store the index of the
				      *  first child (the others
				      *  follow directly) only one
				      *  child index is to be
				      *  given. The level of the child
				      *  is one level up of the level
				      *  of the cell to which this
				      *  iterator points.
				      */
    void set_children (const int index) const;
	
				     /**
				      *  Clear the child field,
				      *  i.e. set it to a value which
				      *  indicates that this cell has
				      *  no children.
				      */
    void clear_children () const;
    
				     /**
				      *  Test whether the hex has
				      *  children.
				      */
    bool has_children () const;

				     /**
				      * Return the number of immediate
				      * children of this object.
				      *
				      * An assertion assures that this
				      * function can only be called
				      * for TriaObjects with
				      * <tt>has_children()==true</tt>.
				      */
    unsigned int n_children() const;

				     /**
				      * Return the number of times
				      * that this cell is
				      * refined. Note that not all its
				      * children are refined that
				      * often (which is why we prepend
				      * @p max_), the returned number
				      * is rather the maximum number
				      * of refinement in any branch of
				      * children of this object.
				      *
				      * For example, if this object is
				      * refined, and one of its children is
				      * refined exactly one more time, then
				      * <tt>max_refinement_depth</tt> should
				      * return 2.
				      *
				      * If this object is not refined (i.e. it
				      * is active), then the return value is
				      * zero.
				      */
    unsigned int max_refinement_depth () const;    
    
				     /**
				      * Return the boundary indicator
				      * of this hex. Since boundary
				      * data is only useful for
				      * structures with a dimension
				      * less than the dimension of a
				      * cell, this function issues an
				      * error if <tt>dim<4</tt>.
				      *
				      * If the return value is 255,
				      * then this line is in the
				      * interior of the domain.
				      */
    unsigned char boundary_indicator () const;

				     /**
				      * Set the boundary indicator of
				      * this hex.  The same applies as
				      * for the
				      * <tt>boundary_indicator()</tt>
				      * function.
				      *
				      * You should be careful with
				      * this function and especially
				      * never try to set the boundary
				      * indicator to 255, unless you
				      * exactly know what you are
				      * doing, since this value is
				      * reserved for another purpose
				      * and algorithms may not work if
				      * boundary cells have this
				      * boundary indicator or if
				      * interior cells have boundary
				      * indicators other than 255.
				      */
    void set_boundary_indicator (const unsigned char) const;

				     /**
				      * Return whether this hex is at
				      * the boundary. This is checked
				      * via the boundary indicator
				      * field, which is always 255 if
				      * the hex is in the interior of
				      * the domain. Obviously, the use
				      * of this function is only
				      * possible for <tt>dim>3</tt>;
				      * however, for <tt>dim==3</tt>, a hex
				      * is a cell and the
				      * CellAccessor class
				      * offers another possibility to
				      * determine whether a cell is at
				      * the boundary or not.
				      */
    bool at_boundary () const;

				     /**
				      * Return the diameter of the hex.
				      *
				      * The diameter of a hex is
				      * computed to be the largest
				      * diagonal. You should
				      * absolutely be clear about the
				      * fact that this definitely is
				      * not the diameter of all
				      * hexahedra; however it may
				      * serve as an approximation and
				      * is exact in many cases,
				      * especially if the hexahedron
				      * is not too much distorted.
				      */
    double diameter () const;

    				     /**
				      * Return the center of the
				      * hex. The center of a hex is
				      * defined to be the average of
				      * the vertices, which is also
				      * the point where the trilinear
				      * mapping places the midpoint of
				      * the unit hex in real space.
				      * However, this may not be the
				      * point of the barycenter of the
				      * hex.
				      */
    Point<dim> center () const;

				     /**
				      * Return the barycenter of the
				      * hex.
				      */
    Point<dim> barycenter () const;

				     /**
				      * Return the volume of the
				      * hex. With volume, we mean the
				      * area included by the
				      * hexahedron if its faces are
				      * supposed to be derived by a
				      * trilinear mapping from the
				      * unit cell, only using the
				      * location of the vertices.
				      * Therefore, no information
				      * about the boundary is used; if
				      * you want other than
				      * trilinearly mapped unit
				      * hexahedra, ask the appropriate
				      * finite element class for the
				      * volume.
				      */
    double measure () const;

				     /**
				      * Compute and return the number
				      * of children of this
				      * hex. Actually, this function
				      * only counts the number of
				      * active children, i.e. the
				      * number if hexs which are not
				      * further refined. Thus, if all
				      * of the eight children of a hex
				      * are further refined exactly
				      * once, the returned number will
				      * be 64, not 80.
				      *
				      * If the present cell is not
				      * refined, one is returned.
				      */
    unsigned int number_of_children () const;

                                     /**
                                      * Return whether the face with
                                      * index @p face has its normal
                                      * pointing in the standard
                                      * direction (@p true) or
                                      * whether it is the opposite
                                      * (@p false). Which is the
                                      * standard direction is
                                      * documented with the
                                      * GeometryInfo class.  In
                                      * 1d and 2d, this is always
                                      * @p true, but in 3d it may be
                                      * different, see the respective
                                      * discussion in the
                                      * documentation of the
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool face_orientation (const unsigned int face) const;
                                     /**
                                      * Return whether the face with index @p
                                      * face is rotated by 180 degrees (@p true)
                                      * or or not (@p false). In 1d and 2d, this
                                      * is always @p false, but in 3d it may be
                                      * different, see the respective discussion
                                      * in the documentation of the
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool face_flip (const unsigned int face) const;

                                     /**
                                      * Return whether the face with index @p
                                      * face is rotated by 90 degrees (@p true)
                                      * or or not (@p false). In 1d and 2d, this
                                      * is always @p false, but in 3d it may be
                                      * different, see the respective discussion
                                      * in the documentation of the
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool face_rotation (const unsigned int face) const;
    
                                     /**
                                      * Return whether the line with index @p
                                      * line is oriented in standard
                                      * direction. @p true indicates, that the
                                      * line is oriented from vertex 0 to vertex
                                      * 1, whereas it is the other way around
                                      * otherwise. In 1d and 2d, this is always
                                      * @p true, but in 3d it may be different,
                                      * see the respective discussion in the
                                      * documentation of the 
                                      * GeometryInfo class.
                                      *
                                      * This function is really only
                                      * for internal use in the
                                      * library unless you absolutely
                                      * know what this is all about.
                                      */
    bool line_orientation (const unsigned int line) const;
    
                                     /**
                                      * Set whether the quad with
                                      * index @p face has its normal
                                      * pointing in the standard
                                      * direction (@p true) or
                                      * whether it is the opposite
                                      * (@p false). Which is the
                                      * standard direction is
                                      * documented with the
                                      * GeometryInfo class.
                                      *
                                      * This function is only for
                                      * internal use in the
                                      * library. Setting this flag to
                                      * any other value than the one
                                      * that the triangulation has
                                      * already set is bound to bring
                                      * you desaster.
                                      */
    void set_face_orientation (const unsigned int face,
                               const bool         orientation) const;

				     /**
                                      * Set the flag indicating, what
                                      * <code>face_flip()</code> will
                                      * return.
				      *
				      * It is only possible to set the
				      * face_orientation of cells in 3d
				      * (i.e. <code>celldim==3 &&
				      * dim==3</code>).
                                      */
    void set_face_flip (const unsigned int face,
			const bool         flip) const;

				     /**
                                      * Set the flag indicating, what
                                      * <code>face_rotation()</code> will
                                      * return.
				      *
				      * It is only possible to set the
				      * face_orientation of cells in 3d
				      * (i.e. <code>celldim==3 &&
				      * dim==3</code>).
                                      */
    void set_face_rotation (const unsigned int face,
			    const bool         rotation) const;
    
  private:
    				     /**
				      *  Copy operator. This is normally
				      *  used in a context like
				      *  <tt>iterator a,b;  *a=*b;</tt>. Since
				      *  the meaning is to copy the
				      *  object pointed to by @p b to
				      *  the object pointed to by
				      *  @p a and since accessors are
				      *  not real but virtual objects,
				      *  this operation is not useful
				      *  for iterators on
				      *  triangulations. We declare
				      *  this function here private,
				      *  thus it may not be used from
				      *  outside.  Furthermore it is
				      *  not implemented and will give
				      *  a linker error if used
				      *  anyway.
				      */
    void operator = (const TriaObjectAccessor<3, dim> &);

  protected:
				     /**@name Advancement of iterators*/
				     /*@{*/
				     /**
				      *  This operator advances the
				      *  iterator to the next element.
				      *
				      *  The next element is next on
				      *  this level if there are
				      *  more. If the present element
				      *  is the last on this level,
				      *  the first on the next level
				      *  is accessed.
				      */
    void operator ++ ();

    				     /**
				      *  This operator moves the
				      *  iterator to the previous
				      *  element.
				      *
				      *  The previous element is
				      *  previous on this level if
				      *  <tt>index>0</tt>. If the present
				      *  element is the first on this
				      *  level, the last on the
				      *  previous level is accessed.
				      */
    void operator -- ();

  private:
    
				     /**
				      * Dimension of the TriaObject
				      * this accessor gives access to.
				      */
    static const unsigned int objectdim = 3;
				     /*@}*/

				     /**
				      * Declare some friends.
				      */
    template <int anydim, typename AnyAccessor> friend class TriaRawIterator;
};



/**
 * This class allows access to a cell: a line in one dimension, a quad
 * in two dimension, etc.
 *
 * The following refers to any space dimension:
 * 
 * This class allows access to a <tt>cell</tt>, which is a line in 1D
 * and a quad in 2D. Cells have more functionality than lines or quads
 * by themselves, for example they can be flagged for refinement, they
 * have neighbors, they have the possibility to check whether they are
 * at the boundary etc. This class offers access to all this data.
 *
 * @ingroup grid
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998, 1999, 2000
 */
template <int dim>
class CellAccessor :  public TriaObjectAccessor<dim,dim>
{
  public:
				     /**
				      * Propagate the AccessorData type
				      * into the present class.
				      */
    typedef typename TriaObjectAccessor<dim,dim>::AccessorData AccessorData;
    
				     /**
				      *  Constructor.
				      */
    CellAccessor (const Triangulation<dim> *parent     =  0,
		  const int                 level      = -1,
		  const int                 index      = -1,
		  const AccessorData       *local_data =  0);

				     /**
				      *  Return a pointer to the
				      *  @p ith neighbor.  If the
				      *  neighbor does not exist, an
				      *  invalid iterator is returned.
				      *  
				      *  <b>Note</b> (cf. TriaLevel<0>):
				      *  The neighbor of a cell has at most the
				      *  same level as this cell, i.e. it may
				      *  or may not be refined.
				      */
    TriaIterator<dim,CellAccessor<dim> >
    neighbor (const unsigned int i) const;

				     /**
				      *  Return the index of the
				      *  @p ith neighbor.  If the
				      *  neighbor does not exist, its
				      *  index is -1.
				      */
    int neighbor_index (const unsigned int i) const;

    				     /**
				      *  Return the level of the
				      *  @p ith neighbor.  If the
				      *  neighbor does not exist, its
				      *  level is -1.
				      */
    int neighbor_level (const unsigned int i) const;

				     /**
				      *  Set the neighbor @p i of
				      *  this cell to the cell pointed
				      *  to by @p pointer.  This line
				      *  must be used.
				      */
    void set_neighbor (const unsigned int i,
		       const TriaIterator<dim,CellAccessor<dim> > &pointer) const;

				     /**
				      * Return the how-many'th
				      * neighbor this cell is of
				      * <tt>cell->neighbor(neighbor)</tt>,
				      * i.e. return the number @p n
				      * such that
				      * <tt>cell->neighbor(neighbor)->neighbor(n)==cell</tt>. This
				      * function is the right one if
				      * you want to know how to get
				      * back from a neighbor to the
				      * present cell.
				      *
				      * Note that this operation is
				      * only useful if the neighbor is
				      * not on a coarser level than
				      * the present cell
				      * (i.e. <tt>cell->neighbor(neighbor)->level()</tt>
				      * needs to be equal to
				      * <tt>cell->level()</tt>. Use the
				      * @p neighbor_of_coarser_neighbor
				      * function in that case.
				      */
    unsigned int neighbor_of_neighbor (const unsigned int neighbor) const;
    
				     /**
				      * This function is a
				      * generalization of the
				      * @p neighbor_of_neighbor
				      * function for the case of a
				      * coarser neighbor. It returns a
				      * pair of numbers, face_no and
				      * subface_no, with the following
				      * property:
				      * <tt>cell->neighbor(neighbor)->face(face_no)->child(subface_no)==cell->face(neighbor)</tt>.
				      *
				      * This function is impossible
				      * for <tt>dim==1</tt>.
				      */
    std::pair<unsigned int, unsigned int>
    neighbor_of_coarser_neighbor (const unsigned int neighbor) const;
    
				     /**
				      *  Return whether the @p ith
				      *  vertex or face (depending on
				      *  the dimension) is part of the
				      *  boundary. This is true, if
				      *  the @p ith neighbor does not
				      *  exist.
				      */
    bool at_boundary (const unsigned int i) const;

				     /**
				      *  Return whether the cell is at
				      *  the boundary. Being at the
				      *  boundary is defined by one
				      *  face being on the
				      *  boundary. Note that this does
				      *  not catch cases where only one
				      *  vertex of a quad or of a hex
				      *  is at the boundary, or where
				      *  only one line of a hex is at
				      *  the boundary while the
				      *  interiors of all faces are in
				      *  the interior of the
				      *  domain. For the latter case,
				      *  the @p has_boundary_lines
				      *  function is the right one to
				      *  ask.
				      */
    bool at_boundary () const;

				     /**
				      * This is a slight variation to
				      * the @p at_boundary function:
				      * for 1 and 2 space dimensions,
				      * it is equivalent, for three
				      * space dimensions it returns
				      * whether at least one of the 12
				      * lines of the hexahedron is at
				      * a boundary. This, of course,
				      * includes the case where a
				      * whole face is at the boundary,
				      * but also some other cases.
				      */
    bool has_boundary_lines () const;

				     /**
				      *  Return whether the refinement
				      *  flag is set or not.
				      */
    bool refine_flag_set () const;

				     /**
				      *  Flag the cell pointed to fo
				      *  refinement. This function is
				      *  only allowed for active
				      *  cells.
				      */
    void set_refine_flag () const;

				     /**
				      *  Clear the refinement flag.
				      */
    void clear_refine_flag () const;

				     /**
				      *  Return whether the coarsen flag
				      *  is set or not.
				      */
    bool coarsen_flag_set () const;

				     /**
				      *  Flag the cell pointed to for
				      *  coarsening. This function is
				      *  only allowed for active
				      *  cells.
				      */
    void set_coarsen_flag () const;

				     /**
				      *  Clear the coarsen flag.
				      */
    void clear_coarsen_flag () const;

				     /**
				      *  Return a pointer to the
				      *  @p ith child. Overloaded
				      *  version which returns a more
				      *  reasonable iterator class.
				      */
    TriaIterator<dim,CellAccessor<dim> >
    child (const unsigned int i) const;

				     /**
				      * Return an iterator to the
				      * @p ith face of this cell.
				      *
				      * This function is not
				      * implemented in 1D, and maps to
				      * QuadAccessor::line in 2D.
				      */
    TriaIterator<dim,TriaObjectAccessor<dim-1, dim> >
    face (const unsigned int i) const;

				     /**
				      * Return the (global) index of the
				      * @p ith face of this cell.
				      *
				      * This function is not
				      * implemented in 1D, and maps to
				      * line_index in 2D and quad_index
				      * in 3D.
				      */
    unsigned int
    face_index (const unsigned int i) const;

                                     /**
                                      * Return an iterator to that
                                      * cell that neighbors the
                                      * present cell on the given face
                                      * and subface number.
                                      *
                                      * To succeed, the present cell
                                      * must not be further refined,
                                      * and the neighbor on the given
                                      * face must be further refined
                                      * exactly once; the returned
                                      * cell is then a child of that
                                      * neighbor.
                                      *
                                      * The function may not be called
                                      * in 1d, since there we have no
                                      * subfaces.  The implementation
                                      * of this function is rather
                                      * straightforward in 2d, by
                                      * first determining which face
                                      * of the neighbor cell the
                                      * present cell is bordering on
                                      * (this is what the
                                      * @p neighbor_of_neighbor
                                      * function does), and then
                                      * asking
                                      * @p GeometryInfo::child_cell_on_subface
                                      * for the index of the
                                      * child.
                                      *
                                      * However, the function is more
                                      * complicated in 3d, since there faces may
                                      * have more than one orientation, and we
                                      * have to use @p face_orientation, @p
                                      * face_flip and @p face_rotation for both
                                      * this and the neighbor cell to figure out
                                      * which cell we want to have.
                                      *
                                      * This can lead to surprising
                                      * results: if we are sitting on
                                      * a cell and are asking for a
                                      * cell behind subface
                                      * <tt>sf</tt>, then this means
                                      * that we are considering the
                                      * subface for the face in the
                                      * natural direction for the
                                      * present cell. However, if the
                                      * face as seen from this cell
                                      * has
                                      * <tt>face_orientation()==false</tt>,
                                      * then the child of the face
                                      * that separates the present
                                      * cell from the neighboring
                                      * cell's child is not
                                      * necessarily the @p sf-th child
                                      * of the face of this cell. This
                                      * is so because the @p
                                      * subface_no parameter to this
                                      * function corresponds to the
                                      * subface with respect to the
                                      * intrinsic ordering of the
                                      * present cell, whereas children
                                      * of face iterators are computed
                                      * with respect to the intrinsic
                                      * ordering of faces; these two
                                      * orderings are only identical
                                      * if the face orientation is @p
                                      * true, and reversed otherwise.
				      *
				      * Similarly, effects of
				      * <tt>face_flip()==true</tt> and
				      * <tt>face_rotation()==true()</tt>, both
				      * of which indicate a non-standard face
				      * have to be considered.
                                      *
                                      * Fortunately, this is only very rarely
                                      * of concern, since one is usually only
                                      * concerned either in the exact number
                                      * of a subface, or in the cell that is
                                      * behind it, not in both at the same
                                      * time.
                                      */
    TriaIterator<dim,CellAccessor<dim> >
    neighbor_child_on_subface (const unsigned int face_no,
                               const unsigned int subface_no) const;
    
				     /**
				      * Return the material id of this
				      * cell.
				      *
				      * For a typical use of this
				      * function, see the @ref step_28
				      * "step-28" tutorial program.
				      */
    unsigned char material_id () const;

				     /**
				      * Set the material id of this
				      * cell.
				      *
				      * For a typical use of this
				      * function, see the @ref step_28
				      * "step-28" tutorial program.
				      */
    void set_material_id (const unsigned char new_material_id) const;

				     /**
				      * Set the material id of this
				      * cell and all its children (and
				      * grand-children, and so on) to
				      * the given value.
				      */
    void recursively_set_material_id (const unsigned char new_material_id) const;

				     /**
				      * Return the subdomain id of
				      * this cell.
				      */
    unsigned int subdomain_id () const;

				     /**
				      * Set the subdomain id of this
				      * cell.
				      */
    void set_subdomain_id (const unsigned int new_subdomain_id) const;

				     /**
				      * Test whether the cell has children
				      * (this is the criterion for activity
				      * of a cell).
				      */
    bool active () const;

				     /**
				      * Test whether the point @p p
				      * is inside this cell. Points on
				      * the boundary are counted as
				      * being inside the cell.
				      *
				      * Note that this function
				      * assumes that the mapping
				      * between unit cell and real
				      * cell is (bi-, tri-)linear,
				      * i.e. that faces in 2d and
				      * edges in 3d are straight
				      * lines. If you have higher
				      * order transformations, results
				      * may be different as to whether
				      * a point is in- or outside the
				      * cell in real space.
				      */
    bool point_inside (const Point<dim> &p) const;
    
    
				     /**
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcRefineCellNotActive);
				     /**
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcCellFlaggedForRefinement);
				     /**
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcCellFlaggedForCoarsening);
    
  private:
    
    				     /**
				      *  Copy operator. This is
				      *  normally used in a context
				      *  like <tt>iterator a,b;
				      *  *a=*b;</tt>. Since the meaning is
				      *  to copy the object pointed to
				      *  by @p b to the object
				      *  pointed to by @p a and since
				      *  accessors are not real but
				      *  virtual objects, this
				      *  operation is not useful for
				      *  iterators on
				      *  triangulations. We declare
				      *  this function here private,
				      *  thus it may not be used from
				      *  outside.  Furthermore it is
				      *  not implemented and will give
				      *  a linker error if used
				      *  anyway.
				      */
    void operator = (const CellAccessor<dim> &);
};



/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <> Point<2> TriaObjectAccessor<2, 2>::barycenter () const;
template <> Point<3> TriaObjectAccessor<2, 3>::barycenter () const;
template <> Point<3> TriaObjectAccessor<3, 3>::barycenter () const;
template <> internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<1> >&
TriaObjectAccessor<1, 1>::objects() const;
template <> internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<2> >&
TriaObjectAccessor<2, 2>::objects() const;
template <> bool CellAccessor<1>::at_boundary () const;
template <> bool CellAccessor<1>::point_inside (const Point<1> &p) const;
template <> bool CellAccessor<2>::at_boundary () const;
template <> bool CellAccessor<2>::point_inside (const Point<2> &p) const;
template <> bool CellAccessor<3>::at_boundary () const;
template <> bool CellAccessor<3>::point_inside (const Point<3> &) const;

template <> bool CellAccessor<1>::has_boundary_lines () const;

template <>
TriaIterator<1,CellAccessor<1> >
CellAccessor<1>::neighbor_child_on_subface (const unsigned int,
                                            const unsigned int) const;
template <>
TriaIterator<2,CellAccessor<2> >
CellAccessor<2>::neighbor_child_on_subface (const unsigned int,
                                            const unsigned int) const;
template <>
TriaIterator<3,CellAccessor<3> >
CellAccessor<3>::neighbor_child_on_subface (const unsigned int,
                                            const unsigned int) const;

template <> double TriaObjectAccessor<2, 2>::measure () const;
template <> double TriaObjectAccessor<2, 3>::measure () const;
template <> double TriaObjectAccessor<3, 3>::measure () const;


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

// include more templates in debug and optimized mode
#  include "tria_accessor.templates.h"

#endif
