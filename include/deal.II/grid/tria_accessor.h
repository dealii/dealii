// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__tria_accessor_h
#define __deal2__tria_accessor_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/tria_iterator_base.h>
#include <deal.II/grid/tria_iterator_selector.h>
#include <deal.II/grid/cell_id.h>

#include <utility>


DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class Triangulation;
template <typename Accessor> class TriaRawIterator;
template <typename Accessor> class TriaIterator;
template <typename Accessor> class TriaActiveIterator;

template <int dim, int spacedim> class Boundary;
template <int dim, int spacedim> class Manifold;


namespace internal
{
  namespace Triangulation
  {
    template <int dim> class TriaObject;
    template <typename G> class TriaObjects;
    struct Implementation;
  }

  namespace TriaAccessor
  {
    struct Implementation;

    /**
     * Implementation of a type with
     * which to store the level of an
     * accessor object. We only need
     * it for the case that
     * <tt>structdim ==
     * dim</tt>. Otherwise, an empty
     * object is sufficient.
     */
    template <int structdim, int dim> struct PresentLevelType
    {
      struct type
      {
        /**
         * Default constructor.
         */
        type ()
        {}

        /**
         * Dummy
         * constructor. Only
         * level zero is allowed.
         */
        type (const int level)
        {
          Assert (level == 0, ExcInternalError());
          (void)level; // removes -Wunused-parameter warning in optimized mode
        }

        /**
         * Dummy conversion
         * operator. Returns
         * level zero.
         */
        operator int () const
        {
          return 0;
        }

        void operator ++ () const
        {
          Assert (false, ExcInternalError());
        }

        void operator -- () const
        {
          Assert (false, ExcInternalError());
        }
      };
    };


    /**
     * Implementation of a type with
     * which to store the level of an
     * accessor object. We only need
     * it for the case that
     * <tt>structdim ==
     * dim</tt>. Otherwise, an empty
     * object is sufficient.
     */
    template <int dim> struct PresentLevelType<dim,dim>
    {
      typedef int type;
    };

  }
}
template <int structdim, int dim, int spacedim> class TriaAccessor;
template <int dim, int spacedim>                class TriaAccessor<0, dim, spacedim>;
template <int spacedim>                         class TriaAccessor<0, 1, spacedim>;

// note: the file tria_accessor.templates.h is included at the end of
// this file.  this includes a lot of templates. originally, this was
// only done in debug mode, but led to cyclic reduction problems and
// so is now on by default.


/**
 * A namespace that contains exception classes used by the accessor classes.
 */
namespace TriaAccessorExceptions
{
//TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException0 (ExcCellNotUsed);
  /**
   * The cell is not an @ref
   * GlossActive "active" cell, but
   * it already has children. Some
   * operations, like setting
   * refinement flags or accessing
   * degrees of freedom are only
   * possible on active cells.
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcCellNotActive);
  /**
   * Trying to access the children of
   * a cell which is in fact active.
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcCellHasNoChildren);
  /**
   * Trying to access the parent of
   * a cell which is in the coarsest
   * level of the triangulation.
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcCellHasNoParent);
//TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException0 (ExcUnusedCellAsChild);
//TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException1 (ExcCantSetChildren,
                  int,
                  << "You can only set the child index if the cell has no "
                  << "children, or clear it. The given "
                  << "index was " << arg1 << " (-1 means: clear children)");
//TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException0 (ExcUnusedCellAsNeighbor);
//TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException0 (ExcUncaughtCase);
//TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException0 (ExcDereferenceInvalidObject);
//TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException0 (ExcCantCompareIterators);
//TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException0 (ExcNeighborIsCoarser);
//TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException0 (ExcNeighborIsNotCoarser);
  /**
   * You are trying to access the
   * level of a face, but faces have
   * no inherent level. The level of
   * a face can only be determined by
   * the level of an adjacent face,
   * which in turn implies that a
   * face can have several levels.
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcFacesHaveNoLevel);
//TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException1 (ExcSetOnlyEvenChildren,
                  int,
                  << "You can only set the child index of an even numbered child."
                  << "The number of the child given was " << arg1 << ".");
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
 * In the implementation, the behavior of this class differs between
 * the cases where <tt>structdim==dim</tt> (cells of a mesh) and
 * <tt>structdim&lt;dim</tt> (faces and edges). For the latter,
 * #present_level is always equal to zero and the constructors may not
 * receive a positive value there. For cells, any level is
 * possible, but only those within the range of the levels of the
 * Triangulation are reasonable. Furthermore, the function objects()
 * returns either the container with all cells on the same level or
 * the container with all objects of this dimension (<tt>structdim&lt;dim</tt>).
 *
 * Some internals of this class are discussed in @ref IteratorAccessorInternals .
 *
 * @ingroup grid
 * @ingroup Accessors
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 2010
 */
template <int structdim, int dim, int spacedim=dim>
class TriaAccessorBase
{
public:
  /**
   *  Dimension of the space the
   *  object represented by this
   *  accessor lives in. For
   *  example, if this accessor
   *  represents a quad that is
   *  part of a two-dimensional
   *  surface in four-dimensional
   *  space, then this value is
   *  four.
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * Dimensionality of the object
   * that the thing represented by
   * this accessopr is part of. For
   * example, if this accessor
   * represents a line that is part
   * of a hexahedron, then this
   * value will be three.
   */
  static const unsigned int dimension = dim;

  /**
   * Dimensionality of the current
   * object represented by this
   * accessor. For example, if it
   * is line (irrespective of
   * whether it is part of a quad
   * or hex, and what dimension we
   * are in), then this value
   * equals 1.
   */
  static const unsigned int structure_dimension = structdim;

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
  TriaAccessorBase (const Triangulation<dim,spacedim> *parent =  0,
                    const int                 level  = -1,
                    const int                 index  = -1,
                    const AccessorData             * =  0);

  /**
   *  Copy constructor. Creates an
   *  object with exactly the same data.
   */
  TriaAccessorBase (const TriaAccessorBase &);

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
  void copy_from (const TriaAccessorBase &);

  /**
   *  Copy operator. Creates an
   *  object with exactly the same data.
   */
  TriaAccessorBase &operator = (const TriaAccessorBase &);

  /**
   * Ordering of accessors. If #structure_dimension is less than
   * #dimension, we simply compare the index of such an object. If
   * #structure_dimension equals #dimension, we compare the level()
   * first, and the index() only if levels are equal.
   */
  bool operator < (const TriaAccessorBase &other) const;

protected:
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
  void operator = (const TriaAccessorBase *);

  /**
   *  Compare for equality.
   */
  bool operator == (const TriaAccessorBase &) const;

  /**
   * Compare for inequality.
   */
  bool operator != (const TriaAccessorBase &) const;

  /**
   * @name Advancement of iterators
   */
  /**
   * @{
   */
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
   * @}
   */

  /**
   * Access to the other objects of
   * a Triangulation with same
   * dimension.
   */
  dealii::internal::Triangulation::TriaObjects<dealii::internal::Triangulation::TriaObject<structdim> > &
  objects () const;

public:
  /**
   * Data type to be used for passing
   * parameters from iterators to the
   * accessor classes in a unified
   * way, no matter what the type of
   * number of these parameters is.
   */
  typedef void *LocalData;

  /**
   * @name Iterator address and state
   */
  /**
   * @{
   */

  /**
   * For cells, this function returns the level within the mesh hierarchy at
   * which this cell is located. For all other objects, the function returns
   * zero.
   *
   * @note Within a Triangulation object, cells are uniquely identified by a
   * pair <code>(level, index)</code> where the former is the cell's
   * refinement level and the latter is the index of the cell within
   * this refinement level (the former being what this function
   * returns). Consequently, there may be multiple cells on different
   * refinement levels but with the same index within their level.
   * Contrary to this, if the current object corresponds to a face or
   * edge, then the object is uniquely identified solely by its index
   * as faces and edges do not have a refinement level. For these objects,
   * the current function always returns zero as the level.
   */
  int level () const;

  /**
   * Return the index of the
   * element presently pointed to
   * on the present level.
   *
   * Within a Triangulation object, cells are uniquely identified by a
   * pair <code>(level, index)</code> where the former is the cell's
   * refinement level and the latter is the index of the cell within
   * this refinement level (the latter being what this function
   * returns). Consequently, there may be multiple cells on different
   * refinement levels but with the same index within their level.
   * Contrary to this, if the current object corresponds to a face or
   * edge, then the object is uniquely identified solely by its index
   * as faces and edges do not have a refinement level.
   *
   * @note The indices objects returned by this function are not a
   * contiguous set of numbers on each level: going from cell to cell,
   * some of the indices in a level may be unused.
   *
   * @note If the triangulation is actually of type
   * parallel::distributed::Triangulation then the indices are
   * relatively only to that part of the distributed triangulation
   * that is stored on the current processor. In other words, cells
   * living in the partitions of the triangulation stored on different
   * processors may have the same index even if they refer to the same
   * cell, and the may have different indices even if they do refer to
   * the same cell (e.g., if a cell is owned by one processor but is a
   * ghost cell on another).
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
  const Triangulation<dim,spacedim> &get_triangulation () const;

  /**
   * @}
   */
protected:
  /**
   * The level if this is a cell
   * (<tt>structdim==dim</tt>). Else,
   * contains zero.
   */
  typename dealii::internal::TriaAccessor::PresentLevelType<structdim,dim>::type present_level;

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
  const Triangulation<dim,spacedim> *tria;

private:

  template <typename Accessor> friend class TriaRawIterator;
  template <typename Accessor> friend class TriaIterator;
  template <typename Accessor> friend class TriaActiveIterator;
};



/**
 * A class that represents accessor objects to iterators that don't
 * make sense such as quad iterators in on 1d meshes.  This class can
 * not be used to create objects (it will in fact throw an exception
 * if this should ever be attempted but it sometimes allows code to be
 * written in a simpler way in a dimension independent way. For
 * example, it allows to write code that works on quad iterators that
 * is dimension independent because quad iterators (with the current
 * class) exist and are syntactically correct. You can not expect,
 * however, to ever generate one of these iterators, meaning you need
 * to expect to wrap the code block in which you use quad iterators
 * into something like <code>if (dim@>1)</code> -- which makes eminent
 * sense anyway.
 *
 * This class provides the minimal interface necessary for Accessor
 * classes to interact with Iterator classes. However, this is only
 * for syntactic correctness, none of the functions do anything but
 * generate errors.
 *
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 2008
 */
template <int structdim, int dim, int spacedim=dim>
class InvalidAccessor :  public TriaAccessorBase<structdim,dim,spacedim>
{
public:
  /**
   * Propagate typedef from
   * base class to this class.
   */
  typedef typename TriaAccessorBase<structdim,dim,spacedim>::AccessorData AccessorData;

  /**
   * Constructor.  This class is
   * used for iterators that make
   * sense in a given dimension,
   * for example quads for 1d
   * meshes. Consequently, while
   * the creation of such objects
   * is syntactically valid, they
   * make no semantic sense, and we
   * generate an exception when
   * such an object is actually
   * generated.
   */
  InvalidAccessor (const Triangulation<dim,spacedim> *parent     =  0,
                   const int                 level      = -1,
                   const int                 index      = -1,
                   const AccessorData       *local_data =  0);

  /**
   * Copy constructor.  This class
   * is used for iterators that
   * make sense in a given
   * dimension, for example quads
   * for 1d meshes. Consequently,
   * while the creation of such
   * objects is syntactically
   * valid, they make no semantic
   * sense, and we generate an
   * exception when such an object
   * is actually generated.
   */
  InvalidAccessor (const InvalidAccessor &);

  /**
   * Conversion from other
   * accessors to the current
   * invalid one. This of course
   * also leads to a run-time
   * error.
   */
  template <typename OtherAccessor>
  InvalidAccessor (const OtherAccessor &);

  /**
   * Dummy copy operation.
   */
  void copy_from (const InvalidAccessor &);

  /**
   * Dummy comparison operators.
   */
  bool operator == (const InvalidAccessor &) const;
  bool operator != (const InvalidAccessor &) const;

  /**
   * Dummy operators to make things
   * compile. Does nothing.
   */
  void operator ++ () const;
  void operator -- () const;

  /**
   * Dummy function representing
   * whether the accessor points to
   * a used or an unused object.
   */
  bool used () const;

  /**
   * Dummy function representing
   * whether the accessor points to
   * an object that has children.
   */
  bool has_children () const;
};



/**
 * A class that provides access to objects in a triangulation such as
 * its vertices, sub-objects, children, geometric information, etc.
 * This class represents objects of dimension <code>structdim</code>
 * (i.e. 1 for lines, 2 for quads, 3 for hexes) in a triangulation of
 * dimensionality <code>dim</code> (i.e. 1 for a triangulation of
 * lines, 2 for a triangulation of quads, and 3 for a triangulation of
 * hexes) that is embedded in a space of dimensionality
 * <code>spacedim</code> (for <code>spacedim==dim</code> the
 * triangulation represents a domain in $R^{dim}$, for
 * <code>spacedim@>dim</code> the triangulation is of a manifold
 * embedded in a higher dimensional space).
 *
 * @ingroup Accessors
 * @author Wolfgang Bangerth and others, 1998, 2000, 2008
 */
template <int structdim, int dim, int spacedim>
class TriaAccessor : public TriaAccessorBase<structdim, dim, spacedim>
{
public:
  /**
   * Propagate typedef from
   * base class to this class.
   */
  typedef
  typename TriaAccessorBase<structdim,dim,spacedim>::AccessorData
  AccessorData;

  /**
   *  Constructor.
   */
  TriaAccessor (const Triangulation<dim,spacedim> *parent     =  0,
                const int                 level      = -1,
                const int                 index      = -1,
                const AccessorData       *local_data =  0);

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
  TriaAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &);

  /**
   * Another conversion operator
   * between objects that don't
   * make sense, just like the
   * previous one.
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor (const TriaAccessor<structdim2,dim2,spacedim2> &);

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
   *  @name Accessing sub-objects
   */
  /**
   * @{
   */

  /**
   *  Return the global index of i-th
   *  vertex of the current object. The
   *  convention regarding the numbering of
   *  vertices is laid down in the
   *  documentation of the GeometryInfo
   *  class.
   *
   *  Note that the returned value is only
   *  the index of the geometrical
   *  vertex. It has nothing to do with
   *  possible degrees of freedom
   *  associated with it. For this, see the
   *  @p DoFAccessor::vertex_dof_index
   *  functions.
   *
   *  @note Despite the name, the index returned here is only
   *  global in the sense that it is specific to a particular
   *  Triangulation object or, in the case the triangulation is
   *  actually of type parallel::distributed::Triangulation,
   *  specific to that part of the distributed triangulation stored
   *  on the current processor.
   */
  unsigned int vertex_index (const unsigned int i) const;

  /**
   *  Return a reference to the
   *  @p ith vertex. The reference is not const, i.e., it is possible
   *  to call this function on the left hand side of an assignment,
   *  thereby moving the vertex of a cell within the triangulation. Of
   *  course, doing so requires that you ensure that the new location
   *  of the vertex remains useful -- for example, avoiding inverted
   *  or otherwise distorted (see also @ref GlossDistorted "this glossary entry").
   *
   *  @note When a cell is refined, its children inherit the position
   *  of the vertex positions of those vertices they share with the mother
   *  cell (plus the locations of the new vertices on edges, faces, and
   *  cell interiors that are created for the new child cells). If the
   *  vertex of a cell is moved, this implies that its children will also
   *  use these new locations. On the other hand, imagine a 2d situation
   *  where you have one cell that is refined (with four children) and then
   *  you move the central vertex connecting all four children. If you
   *  coarsen these four children again to the mother cell, then the
   *  location of the moved vertex is lost and if, in a later step, you
   *  refine the mother cell again, the then again new vertex will be
   *  placed again at the same position as the first time around -- i.e.,
   *  not at the location you had previously moved it to.
   *
   *  @note The behavior described above is relevant if
   *  you have a parallel::distributed::Triangulation object. There,
   *  refining a mesh always involves a re-partitioning. In other words,
   *  vertices of locally owned cells (see
   *  @ref GlossLocallyOwnedCell "this glossary entry")
   *  that you may have moved to a different location on one processor
   *  may be moved to a different processor upon mesh refinement (even
   *  if these particular cells were not refined) which will re-create
   *  their position based on the position of the coarse cells they
   *  previously had, not based on the position these vertices had on
   *  the processor that previously owned them. In other words, in
   *  parallel computations, you will probably have to move nodes
   *  explicitly after every mesh refinement because vertex positions
   *  may or may not be preserved across the re-partitioning that
   *  accompanies mesh refinement.
   */
  Point<spacedim> &vertex (const unsigned int i) const;

  /**
   * Pointer to the @p ith line
   * bounding this object.
   */
  typename dealii::internal::Triangulation::Iterators<dim,spacedim>::line_iterator
  line (const unsigned int i) const;

  /**
   * Line index of the @p ith
   * line bounding this object.
   *
   * Implemented only for
   * <tt>structdim>1</tt>,
   * otherwise an exception
   * generated.
   */
  unsigned int line_index (const unsigned int i) const;

  /**
   * Pointer to the @p ith quad
   * bounding this object.
   */
  typename dealii::internal::Triangulation::Iterators<dim,spacedim>::quad_iterator
  quad (const unsigned int i) const;

  /**
   * Quad index of the @p ith
   * quad bounding this object.
   *
   * Implemented only for
   * <tt>structdim>2</tt>,
   * otherwise an exception
   * generated.
   */
  unsigned int quad_index (const unsigned int i) const;
  /**
   * @}
   */

  /**
   *  @name Orientation of sub-objects
   */
  /**
   * @{
   */

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
   * @}
   */

  /**
   *  @name Accessing children
   */
  /**
   * @{
   */

  /**
   *  Test whether the object has
   *  children.
   */
  bool has_children () const;

  /**
   * Return the number of immediate
   * children of this object. The
   * number of children of an
   * unrefined cell is zero.
   */
  unsigned int n_children() const;

  /**
   * Compute and return the number
   * of active descendants of this
   * objects. For example, if all
   * of the eight children of a hex
   * are further refined
   * isotropically exactly once,
   * the returned number will be
   * 64, not 80.
   *
   * If the present cell is not
   * refined, one is returned.
   *
   * If one considers a triangulation as a
   * forest where the root of each tree are
   * the coarse mesh cells and nodes have
   * descendents (the children of a cell),
   * then this function returns the number
   * of terminal nodes in the sub-tree
   * originating from the current object;
   * consequently, if the current object is
   * not further refined, the answer is
   * one.
   */
  unsigned int number_of_children () const;

  /**
   * Return the number of times
   * that this object is
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
   *  Return an iterator to the @p ith
   *  child.
   */
  TriaIterator<TriaAccessor<structdim,dim,spacedim> >
  child (const unsigned int i) const;

  /**
   *  Return an iterator to that object
   *  that is identical to the ith child
   *  for isotropic refinement. If the
   *  current object is refined
   *  isotropically, then the returned
   *  object is the ith child. If the
   *  current object is refined
   *  anisotropically, the returned child
   *  may in fact be a grandchild of the
   *  object, or may not exist at all (in
   *  which case an exception is
   *  generated).
   */
  TriaIterator<TriaAccessor<structdim,dim,spacedim> >
  isotropic_child (const unsigned int i) const;

  /**
   * Return the RefinementCase
   * of this cell.
   */
  RefinementCase<structdim> refinement_case () const;

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
   *  Index of the @p ith isotropic child.
   *  See the isotropic_child() function
   *  for a definition of this concept.  If
   *  the child does not exist, -1 is
   *  returned.
   */
  int isotropic_child_index (const unsigned int i) const;
  /**
   * @}
   */

  /**
   *  @name Dealing with boundary indicators
   */
  /**
   * @{
   */

  /**
   * Boundary indicator of this
   * object.
   *
   * If the return value is the special
   * value numbers::internal_face_boundary_id,
   * then this object is in the
   * interior of the domain.
   *
   * @see @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  types::boundary_id boundary_indicator () const;

  /**
   * Set the boundary indicator.
   * The same applies as for the
   * <tt>boundary_indicator()</tt>
   * function.
   *
   * Note that it only sets the
   * boundary object of the current
   * object itself, not the
   * indicators of the ones that
   * bound it. For example, in 3d,
   * if this function is called on
   * a face, then the boundary
   * indicator of the 4 edges that
   * bound the face remain
   * unchanged. If you want to set
   * the boundary indicators of
   * face and edges at the same
   * time, use the
   * set_all_boundary_indicators()
   * function. You can see the result of not using the correct function in the
   * results section of step-49.
   *
   * @warning You should never set the
   * boundary indicator of an interior face
   * (a face not at the boundary of the
   * domain), or set set the boundary
   * indicator of an exterior face to
   * numbers::internal_face_boundary_id
   * (this value is reserved for another
   * purpose). Algorithms may not work or
   * produce very confusing results if
   * boundary cells have a boundary
   * indicator of numbers::internal_face_boundary_id
   * or if interior cells have boundary
   * indicators other than numbers::internal_face_boundary_id.
   * Unfortunately, the current object
   * has no means of finding out whether it
   * really is at the boundary of the
   * domain and so cannot determine whether
   * the value you are trying to set makes
   * sense under the current circumstances.
   *
   * @ingroup boundary
   *
   * @see @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  void set_boundary_indicator (const types::boundary_id) const;

  /**
   * Do as set_boundary_indicator()
   * but also set the boundary
   * indicators of the objects that
   * bound the current object. For
   * example, in 3d, if
   * set_boundary_indicator() is
   * called on a face, then the
   * boundary indicator of the 4
   * edges that bound the face
   * remain unchanged. In contrast, if you call the current function,
   * the boundary indicators
   * of face and edges are all set to the given value.
   *
   * This function is useful if you set boundary indicators of faces
   * in 3d (in 2d, the function does the same as set_boundary_indicator())
   * and you do so because you want a curved boundary object to
   * represent the part of the boundary that corresponds to the
   * current face. In that case, the Triangulation class needs to figure
   * out where to put new vertices upon mesh refinement, and higher order
   * Mapping objects also need to figure out where new interpolation points
   * for a curved boundary approximation should be. In either case, the
   * two classes first determine where interpolation points on the edges
   * of a boundary face should be, asking the boundary object, before
   * asking the boundary object for the interpolation points corresponding
   * to the interior of the boundary face. For this to work properly, it is
   * not sufficient to have set the boundary indicator for the face alone,
   * but you also need to set the boundary indicators of the edges that
   * bound the face. This function does all of this at once. You can see
   * the result of not using the correct function in the
   * results section of step-49.
   *
   * @ingroup boundary
   *
   * @see @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  void set_all_boundary_indicators (const types::boundary_id) const;

  /**
   * Return whether this object is at the
   * boundary. Obviously, the use of this
   * function is only possible for
   * <tt>dim@>structdim</tt>; however, for
   * <tt>dim==structdim</tt>, an object is a
   * cell and the CellAccessor class offers
   * another possibility to determine
   * whether a cell is at the boundary or
   * not.
   */
  bool at_boundary () const;

  /**
   * Return a constant reference to the manifold object used for this
   * object. This function exists for backward compatibility and calls
   * get_manifold() internally.
   */
  const Boundary<dim,spacedim> &get_boundary () const;

  /**
   * Return a constant reference to the manifold object used for this object.
   *
   * As explained in
   * @ref boundary "Boundary and manifold description for triangulations",
   * the process involved in finding the appropriate manifold description
   * involves querying both the manifold or boundary indicators. See there
   * for more information.
   */
  const Manifold<dim,spacedim> &get_manifold () const;

  /**
   * @}
   */

  /**
   *  @name Dealing with manifold indicators
   */
  /**
   * @{
   */

  /**
   * Return the manifold indicator of this
   * object.
   *
   * If the return value is the special value
   * numbers::flat_manifold_id, then this object is associated with a
   * standard Cartesian Manifold Description.
   *
   * @see @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  types::manifold_id manifold_id () const;

  /**
   * Set the manifold indicator.  The same applies as for the
   * <tt>manifold_id()</tt> function.
   *
   * Note that it only sets the manifold object of the current object
   * itself, not the indicators of the ones that bound it, nor of its
   * children. For example, in 3d, if this function is called on a
   * face, then the manifold indicator of the 4 edges that bound the
   * face remain unchanged. If you want to set the manifold indicators
   * of face, edges and all children at the same time, use the
   * set_all_manifold_ids() function.
   *
   *
   * @ingroup manifold
   *
   * @see @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  void set_manifold_id (const types::manifold_id) const;

  /**
   * Do as set_manifold_id()
   * but also set the manifold
   * indicators of the objects that
   * bound the current object. For
   * example, in 3d, if
   * set_manifold_id() is
   * called on a face, then the
   * manifold indicator of the 4
   * edges that bound the face
   * remain unchanged. On the other
   * hand, the manifold indicators
   * of face and edges are all set
   * at the same time using the
   * current function.
   *
   * @ingroup manifold
   *
   * @see @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  void set_all_manifold_ids (const types::manifold_id) const;

  /**
   * @}
   */


  /**
   * @name User data
   */
  /**
   * @{
   */
  /**
   *  Read the user flag.
  *  See @ref GlossUserFlags for more information.
   */
  bool user_flag_set () const;

  /**
   *  Set the user flag.
  *  See @ref GlossUserFlags for more information.
   */
  void set_user_flag () const;

  /**
   *  Clear the user flag.
  *  See @ref GlossUserFlags for more information.
   */
  void clear_user_flag () const;

  /**
   *  Set the user flag for this
   * and all descendants.
  *  See @ref GlossUserFlags for more information.
   */
  void recursively_set_user_flag () const;

  /**
   *  Clear the user flag for this
   * and all descendants.
  *  See @ref GlossUserFlags for more information.
   */
  void recursively_clear_user_flag () const;

  /**
   * Reset the user data to zero,
   * independent if pointer or index.
  * See @ref GlossUserData for more information.
   */
  void clear_user_data () const;

  /**
   * Set the user pointer
   * to @p p.
   *
   * @note User pointers and user
   * indices are mutually
   * exclusive. Therefore, you can
   * only use one of them, unless
   * you call
   * Triangulation::clear_user_data()
   * in between.
  *
  * See @ref GlossUserData for more information.
   */
  void set_user_pointer (void *p) const;

  /**
   * Reset the user pointer
   * to a @p NULL pointer.
  * See @ref GlossUserData for more information.
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
   *
   * @note User pointers and user
   * indices are mutually
   * exclusive. Therefore, you can
   * only use one of them, unless
   * you call
   * Triangulation::clear_user_data()
   * in between.
   * <tt>A *a=static_cast<A*>(cell->user_pointer());</tt>.
  *
  * See @ref GlossUserData for more information.
   */
  void *user_pointer () const;

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
   *
   * @note User pointers and user
   * indices are mutually
   * exclusive. Therefore, you can
   * only use one of them, unless
   * you call
   * Triangulation::clear_user_data()
   * in between.
  *
  * See @ref GlossUserData for more information.
   */
  void recursively_set_user_pointer (void *p) const;

  /**
   * Clear the user pointer of this
   * object and all of its
   * descendants. The same holds as
   * said for the
   * recursively_set_user_pointer()
   * function.
  * See @ref GlossUserData for more information.
   */
  void recursively_clear_user_pointer () const;

  /**
   * Set the user index
   * to @p p.
   *
   * @note User pointers and user
   * indices are mutually
   * exclusive. Therefore, you can
   * only use one of them, unless
   * you call
   * Triangulation::clear_user_data()
   * in between.
  * See @ref GlossUserData for more information.
   */
  void set_user_index (const unsigned int p) const;

  /**
   * Reset the user index to 0.
  * See @ref GlossUserData for more information.
   */
  void clear_user_index () const;

  /**
   * Access the value of the user
   * index.
   *
   * @note User pointers and user
   * indices are mutually
   * exclusive. Therefore, you can
   * only use one of them, unless
   * you call
   * Triangulation::clear_user_data()
   * in between.
  *
  * See @ref GlossUserData for more information.
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
   *
   * @note User pointers and user
   * indices are mutually
   * exclusive. Therefore, you can
   * only use one of them, unless
   * you call
   * Triangulation::clear_user_data()
   * in between.
  *
  * See @ref GlossUserData for more information.
   */
  void recursively_set_user_index (const unsigned int p) const;

  /**
   * Clear the user index of this
   * object and all of its
   * descendants. The same holds as
   * said for the
   * recursively_set_user_index()
   * function.
  *
  * See @ref GlossUserData for more information.
   */
  void recursively_clear_user_index () const;
  /**
   * @}
   */

  /**
   * @name Geometric information about an object
   */
  /**
   * @{
   */

  /**
   * Diameter of the object.
   *
   * The diameter of an object is computed
   * to be the largest diagonal. This is
   * not necessarily the true diameter for
   * objects that may use higher order
   * mappings, but completely sufficient
   * for most computations.
   */
  double diameter () const;

  /**
   * Length of an object in the direction
   * of the given axis, specified in the
   * local coordinate system. See the
   * documentation of GeometryInfo for the
   * meaning and enumeration of the local
   * axes.
   *
   * Note that the "length" of an object
   * can be interpreted in a variety of
   * ways. Here, we choose it as the
   * maximal length of any of the edges of
   * the object that are parallel to the
   * chosen axis on the reference cell.
   */
  double extent_in_direction (const unsigned int axis) const;

  /**
   * Returns the minimal distance between
   * any two vertices.
   */
  double minimum_vertex_distance () const;

  /**
   * Returns a point belonging to the Manifold<dim,spacedim> where
   * this object lives, given its parametric coordinates on the
   * reference #structdim cell. This function queries the underlying
   * manifold object, and can be used to obtain the exact geometrical
   * location of arbitrary points on this object.
   *
   * Notice that the argument @p coordinates are the coordinates on
   * the <emph>reference cell</emph>, given in reference
   * coordinates. In other words, the argument provides a weighting
   * between the different vertices. For example, for lines, calling
   * this function with argument Point<1>(.5), is equivalent to asking
   * the line for its center.
   */
  Point<spacedim> intermediate_point(const Point<structdim> &coordinates) const;

  /**
   * Center of the object. The center of an object is defined to be
   * the average of the locations of the vertices. If required, the
   * user may ask this function to return the average of the point
   * according to the underlyinging Manifold object, by setting to
   * true the optional parameter @p respect_manifold.
   *
   * When the geometry of a TriaAccessor is not flat, or when part of
   * the bounding objects of this TriaAccessor are not flat, the
   * result given by the TriaAccessor::center() function may not be
   * accurate enough, even when parameter @p respect_manifold is set
   * to true. If you find this to be case, than you can further refine
   * the computation of the center by setting to true the second
   * additional parameter @p use_laplace_transformation, which will
   * force this function to compute the location of the center by
   * solving a linear elasticity problem with Dirichlet boundary
   * conditions set to the location of the bounding vertices and the
   * centers of the bounding lines and quads.
   */
  Point<spacedim> center (const bool respect_manifold=false,
                          const bool use_laplace_transformation=false) const;

  /**
   * Barycenter of the object.
   */
  Point<spacedim> barycenter () const;

  /**
   * Compute the dim-dimensional measure of the object. For a
   * dim-dimensional cell in dim-dimensional space, this equals its
   * volume. On the other hand, for a 2d cell in 3d space, or if the
   * current object pointed to is a 2d face of a 3d cell in 3d space,
   * then the function computes the area the object occupies. For a
   * one-dimensional object, return its length.
   *
   * The function only computes the measure of cells, faces or edges
   * assumed to be represented by (bi-/tri-)linear mappings. In other
   * words, it only takes into account the locations of the vertices
   * that bound the current object but not how the interior of the
   * object may actually be mapped. In most simple cases, this is
   * exactly what you want. However, for objects that are not
   * "straight", e.g. 2d cells embedded in 3d space as part of a
   * triangulation of a curved domain, two-dimensional faces of 3d
   * cells that are not just parallelograms, or for faces that are at
   * the boundary of a domain that is not just bounded by straight
   * line segments or planes, this function only computes the
   * dim-dimensional measure of a (bi-/tri-)linear interpolation of
   * the "real" object as defined by the manifold or boundary object
   * describing the real geometry of the object in question. If you
   * want to consider the "real" geometry, you will need to compute
   * the measure by integrating a function equal to one over the
   * object, which after applying quadrature equals the summing the
   * JxW values returned by the FEValues or FEFaceValues object you
   * will want to use for the integral.
   */
  double measure () const;

  /**
   * Return true if the current object is a
   * translation of the given argument.
   *
   * @note For the purpose of a
   * triangulation, cells, faces, etc are
   * only characterized by their
   * vertices. The current function
   * therefore only compares the locations
   * of vertices. For many practical
   * applications, however, it is not only
   * the vertices that determine whether
   * one cell is a translation of another,
   * but also how the cell is mapped from
   * the reference cell to its location in
   * real space. For example, if we are
   * using higher order mappings, then not
   * only do the vertices have to be
   * translations of each other, but also
   * the points along edges. In these
   * questions, therefore, it would be
   * appropriate to ask the mapping, not
   * the current function, whether two
   * objects are translations of each
   * other.
   */
  bool
  is_translation_of (const TriaIterator<TriaAccessor<structdim,dim,spacedim> > &o) const;

  /**
   * @}
   */


private:
  /**
   *  Copy the data of the given
   *  object into the internal data
   *  structures of a
   *  triangulation.
   */
  void set (const dealii::internal::Triangulation::TriaObject<structdim> &o) const;

  /**
   * Set the flag indicating, what
   * <code>line_orientation()</code> will
   * return.
   *
   * It is only possible to set the
   * line_orientation of faces in 3d
   * (i.e. <code>structdim==2 &&
   * dim==3</code>).
   */
  void set_line_orientation (const unsigned int line,
                             const bool         orientation) const;

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
   * (i.e. <code>structdim==3 &&
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
   * (i.e. <code>structdim==3 &&
   * dim==3</code>).
   */
  void set_face_rotation (const unsigned int face,
                          const bool         rotation) const;

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
   * Set the @p RefinementCase<dim> this
   * TriaObject is refined with.
   * Not defined for
   * <tt>structdim=1</tt> as lines
   * are always refined resulting
   * in 2 children lines (isotropic
   * refinement).
   *
   * You should know quite exactly
   * what you are doing if you
   * touch this function. It is
   * exclusively for internal use
   * in the library.
   */
  void set_refinement_case (const RefinementCase<structdim> &ref_case) const;

  /**
   * Clear the RefinementCase<dim> of
   * this TriaObject, i.e. reset it
   * to RefinementCase<dim>::no_refinement.
   *
   * You should know quite exactly
   * what you are doing if you
   * touch this function. It is
   * exclusively for internal use
   * in the library.
   */
  void clear_refinement_case () const;

  /**
   *  Set the index of the ith
   *  child. Since the children
   *  come at least in pairs, we
   *  need to store the index of
   *  only every second child,
   *  i.e. of the even numbered
   *  children. Make sure, that the
   *  index of child i=0 is set
   *  first. Calling this function
   *  for odd numbered children is
   *  not allowed.
   */
  void set_children (const unsigned int i, const int index) const;

  /**
   *  Clear the child field,
   *  i.e. set it to a value which
   *  indicates that this cell has
   *  no children.
   */
  void clear_children () const;

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
  void operator = (const TriaAccessor &);

  template <int, int> friend class Triangulation;

  friend struct dealii::internal::Triangulation::Implementation;
  friend struct dealii::internal::TriaAccessor::Implementation;
};






/**
 * Closure class to stop induction of classes. Should never be called
 * and thus produces an error when created.
 *
 * @ingroup grid
 */
template<int dim, int spacedim>
class TriaAccessor<0, dim, spacedim>
{
private:
  /**
   * Constructor. Made private to
   * make sure that this class
   * can't be used.
   */
  TriaAccessor ();
};



/**
 * A class that represents an access to a face in 1d -- i.e. to a
 * point. This is not a full fledged access from which you can build
 * an iterator: for example, you can't iterate from one such point to
 * the next. Point also don't have children, and they don't have
 * neighbors.
 *
 * @author Wolfgang Bangerth, 2010
 */
template <int spacedim>
class TriaAccessor<0, 1, spacedim>
{
public:
  /**
   *  Dimension of the space the
   *  object represented by this
   *  accessor lives in. For
   *  example, if this accessor
   *  represents a quad that is
   *  part of a two-dimensional
   *  surface in four-dimensional
   *  space, then this value is
   *  four.
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * Dimensionality of the object
   * that the thing represented by
   * this accessopr is part of. For
   * example, if this accessor
   * represents a line that is part
   * of a hexahedron, then this
   * value will be three.
   */
  static const unsigned int dimension = 1;

  /**
   * Dimensionality of the current
   * object represented by this
   * accessor. For example, if it
   * is line (irrespective of
   * whether it is part of a quad
   * or hex, and what dimension we
   * are in), then this value
   * equals 1.
   */
  static const unsigned int structure_dimension = 0;

  /**
   * Pointer to internal data.
   */
  typedef void AccessorData;

  /**
   * Whether the vertex represented
   * here is at the left end of the
   * domain, the right end, or in
   * the interior.
   */
  enum VertexKind
  {
    left_vertex,
    interior_vertex,
    right_vertex
  };

  /**
   * Constructor.
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
   */
  TriaAccessor (const Triangulation<1,spacedim> *tria,
                const VertexKind      vertex_kind,
                const unsigned int    vertex_index);

  /**
   * Constructor. This constructor
   * exists in order to maintain
   * interface compatibility with
   * the other accessor
   * classes. However, it doesn't
   * do anything useful here and so
   * may not actually be called.
   */
  TriaAccessor (const Triangulation<1,spacedim> *tria = 0,
                const int = 0,
                const int = 0,
                const AccessorData * = 0);

  /**
   * Constructor. Should never be
   * called and thus produces an
   * error.
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor (const TriaAccessor<structdim2,dim2,spacedim2> &);

  /**
   * Constructor. Should never be
   * called and thus produces an
   * error.
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &);

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
   *  Return the state of the
   *  iterator. Since an iterator
   *  to points can not be
   *  incremented or decremented,
   *  its state remains constant,
   *  and in particular equal to
   *  IteratorState::valid.
   */
  static IteratorState::IteratorStates state ();

  /**
   * Level of this object. Vertices
   * have no level, so this
   * function always returns zero.
   */
  static int level ();

  /**
   * Index of this object. Returns
   * the global index of the vertex
   * this object points to.
   */
  int index () const;

  /**
   * @name Advancement of iterators
   */
  /**
   * @{
   */
  /**
   *  This operator advances the
   *  iterator to the next
   *  element. For points, this
   *  operation is not defined, so
   *  you can't iterate over point
   *  iterators.
   */
  void operator ++ () const;

  /**
   *  This operator moves the
   *  iterator to the previous
   *  element. For points, this
   *  operation is not defined, so
   *  you can't iterate over point
   *  iterators.
   */
  void operator -- () const;
  /**
   *  Compare for equality.
   */
  bool operator == (const TriaAccessor &) const;

  /**
   * Compare for inequality.
   */
  bool operator != (const TriaAccessor &) const;

  /**
   * @}
   */

  /**
   *  @name Accessing sub-objects
   */
  /**
   * @{
   */

  /**
   *  Return the global index of
   *  i-th vertex of the current
   *  object. If i is zero, this
   *  returns the index of the
   *  current point to which this
   *  object refers. Otherwise, it
   *  throws an exception.
   *
   *  Note that the returned value is only
   *  the index of the geometrical
   *  vertex. It has nothing to do with
   *  possible degrees of freedom
   *  associated with it. For this, see the
   *  @p DoFAccessor::vertex_dof_index
   *  functions.
   *
   *  @note Despite the name, the index returned here is only
   *  global in the sense that it is specific to a particular
   *  Triangulation object or, in the case the triangulation is
   *  actually of type parallel::distributed::Triangulation,
   *  specific to that part of the distributed triangulation stored
   *  on the current processor.
   */
  unsigned int vertex_index (const unsigned int i = 0) const;

  /**
   *  Return a reference to the
   *  @p ith vertex. If i is zero, this
   *  returns a reference to the
   *  current point to which this
   *  object refers. Otherwise, it
   *  throws an exception.
   */
  Point<spacedim> &vertex (const unsigned int i = 0) const;

  /**
   * Return the center of this object,
   * which of course co-incides with the
   * location of the vertex this object
   * refers to.
   */
  Point<spacedim> center () const;

  /**
   * Pointer to the @p ith line
   * bounding this object. Will
   * point to an invalid object.
   */
  typename dealii::internal::Triangulation::Iterators<1,spacedim>::line_iterator
  static line (const unsigned int);

  /**
   * Line index of the @p ith
   * line bounding this object.
   *
   * Implemented only for
   * <tt>structdim>1</tt>,
   * otherwise an exception
   * generated.
   */
  static unsigned int line_index (const unsigned int i);

  /**
   * Pointer to the @p ith quad
   * bounding this object.
   */
  static
  typename dealii::internal::Triangulation::Iterators<1,spacedim>::quad_iterator
  quad (const unsigned int i);

  /**
   * Quad index of the @p ith
   * quad bounding this object.
   *
   * Implemented only for
   * <tt>structdim>2</tt>,
   * otherwise an exception
   * generated.
   */
  static unsigned int quad_index (const unsigned int i);

  /**
   * @}
   */


  /**
   * Return whether this point is
   * at the boundary of the
   * one-dimensional triangulation
   * we deal with here.
   */
  bool at_boundary () const;

  /**
   * Return the boundary indicator of this object. The convention for one
   * dimensional triangulations is that left end vertices (of each line
   * segment from which the triangulation may be constructed) have boundary
   * indicator zero, and right end vertices have boundary indicator one,
   * unless explicitly set differently.
   *
   * If the return value is the special value
   * numbers::internal_face_boundary_id, then this object is in the interior
   * of the domain.
   *
   * @see @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  types::boundary_id boundary_indicator () const;

  /**
   * Return the manifold indicator of this
   * object.
   *
   * @see @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  types::manifold_id manifold_id () const;


  /**
   *  @name Orientation of sub-objects
   */
  /**
   * @{
   */

  /**
   * @brief Always return false
   */
  static bool face_orientation (const unsigned int face);

  /**
   * @brief Always return false
   */
  static bool face_flip (const unsigned int face);

  /**
   * @brief Always return false
   */
  static bool face_rotation (const unsigned int face);

  /**
   * @brief Always return false
   */
  static bool line_orientation (const unsigned int line);

  /**
   * @}
   */

  /**
   *  @name Accessing children
   */
  /**
   * @{
   */

  /**
   *  Test whether the object has
   *  children. Always false.
   */
  static bool has_children ();

  /**
   * Return the number of immediate
   * children of this object.This
   * is always zero in dimension 0.
   */
  static unsigned int n_children();

  /**
   * Compute and return the number
   * of active descendants of this
   * objects. Always zero.
   */
  static unsigned int number_of_children ();

  /**
   * Return the number of times
   * that this object is
   * refined. Always 0.
   */
  static unsigned int max_refinement_depth ();

  /**
   * @brief Return an invalid object
   */
  static
  TriaIterator<TriaAccessor<0,1,spacedim> >
  child (const unsigned int);

  /**
   * @brief Return an invalid object
   */
  static
  TriaIterator<TriaAccessor<0,1,spacedim> >
  isotropic_child (const unsigned int);

  /**
   * Always return no refinement.
   */
  static
  RefinementCase<0> refinement_case ();

  /**
   * @brief Returns -1
   */
  static
  int child_index (const unsigned int i);

  /**
   * @brief Returns -1
   */
  static
  int isotropic_child_index (const unsigned int i);
  /**
   * @}
   */

  /**
   *  @name Dealing with boundary indicators
   */
  /**
   * @{
   */

  /**
   * Set the boundary indicator.
   * The same applies as for the
   * <tt>boundary_indicator()</tt>
   * function.
   *
   * @warning You should never set the
   * boundary indicator of an interior face
   * (a face not at the boundary of the
   * domain), or set set the boundary
   * indicator of an exterior face to
   *  numbers::internal_face_boundary_id
   * (this value is reserved for another
   * purpose). Algorithms may not work or
   * produce very confusing results if
   * boundary cells have a boundary
   * indicator of numbers::internal_face_boundary_id
   * or if interior cells have boundary
   * indicators other than numbers::internal_face_boundary_id.
   * Unfortunately, the current object
   * has no means of finding out whether it
   * really is at the boundary of the
   * domain and so cannot determine whether
   * the value you are trying to set makes
   * sense under the current circumstances.
   *
   * @ingroup boundary
   *
   * @see @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  void
  set_boundary_indicator (const types::boundary_id);

  /**
   * Set the manifold indicator of this vertex. This does nothing so far since
   * manifolds are only used to refine and map objects, but vertices are not
   * refined and the mapping is trivial. This function is here only to allow
   * dimension independent programming.
   */
  void
  set_manifold_id (const types::manifold_id);

  /**
   * Set the boundary indicator of this object and all of its
   * lower-dimensional sub-objects.  Since this object only represents a
   * single vertex, there are no lower-dimensional obejct and this function is
   * equivalent to calling set_boundary_indicator with the same argument.
   *
   * @ingroup boundary
   *
   * @see @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  void
  set_all_boundary_indicators (const types::boundary_id);

  /**
   * Set the manifold indicator of this object and all of its
   * lower-dimensional sub-objects.  Since this object only represents a
   * single vertex, there are no lower-dimensional obejct and this function is
   * equivalent to calling set_manifold_indicator with the same argument.
   *
   * @ingroup manifold
   *
   * @see @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */

  void
  set_all_manifold_ids (const types::manifold_id);
  /**
   * @}
   */

  /**
   * Return whether the vertex
   * pointed to here is used.
   */
  bool used () const;

protected:
  /**
   * Pointer to the triangulation
   * we operate on.
   */
  const Triangulation<1,spacedim> *tria;

  /**
   * Whether this is a left end,
   * right end, or interior
   * vertex. This information is
   * provided by the cell at the
   * time of creation.
   */
  VertexKind      vertex_kind;

  /**
   * The global vertex index of the
   * vertex this object corresponds
   * to.
   */
  unsigned int    global_vertex_index;
};




/**
 * This class allows access to a cell: a line in one dimension, a quad
 * in two dimension, etc.
 *
 * The following refers to any dimension:
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
template <int dim, int spacedim=dim>
class CellAccessor :  public TriaAccessor<dim,dim,spacedim>
{
public:
  /**
   * Propagate the AccessorData type
   * into the present class.
   */
  typedef typename TriaAccessor<dim,dim,spacedim>::AccessorData AccessorData;

  /**
   * Define the type of the
   * container this is part of.
   */
  typedef Triangulation<dim, spacedim> Container;

  /**
   *  @name Constructors
   */
  /**
   * @{
   */

  /**
   *  Constructor.
   */
  CellAccessor (const Triangulation<dim,spacedim> *parent     =  0,
                const int                 level      = -1,
                const int                 index      = -1,
                const AccessorData       *local_data =  0);

  /**
   * Copy constructor.
   */
  CellAccessor (const TriaAccessor<dim,dim,spacedim> &cell_accessor);

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
  CellAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &);

  /**
   * Another conversion operator
   * between objects that don't
   * make sense, just like the
   * previous one.
   */
  template <int structdim2, int dim2, int spacedim2>
  CellAccessor (const TriaAccessor<structdim2,dim2,spacedim2> &);

  /**
   * @}
   */

  /**
   *  @name Accessing sub-objects and neighbors
   */
  /**
   * @{
   */

  /**
   *  Return a pointer to the
   *  @p ith child. Overloaded
   *  version which returns a more
   *  reasonable iterator class.
   */
  TriaIterator<CellAccessor<dim, spacedim> >
  child (const unsigned int i) const;

  /**
   * Return an iterator to the
   * @p ith face of this cell.
   */
  TriaIterator<TriaAccessor<dim-1,dim,spacedim> >
  face (const unsigned int i) const;

  /**
   * Return the (global) index of the
   * @p ith face of this cell.
   *
   *  @note Despite the name, the index returned here is only
   *  global in the sense that it is specific to a particular
   *  Triangulation object or, in the case the triangulation is
   *  actually of type parallel::distributed::Triangulation,
   *  specific to that part of the distributed triangulation stored
   *  on the current processor.
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
   * However, the situation is more
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
   * subface_no on a cell corresponds to the
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
   * Fortunately, this is only very rarely of
   * concern, since usually one simply wishes
   * to loop over all finer neighbors at a
   * given face of an active cell. Only in
   * the process of refinement of a
   * Triangulation we want to set neighbor
   * information for both our child cells and
   * the neighbor's children. Since we can
   * respect orientation of faces from our
   * current cell in that case, we do NOT
   * respect face_orientation, face_flip and
   * face_rotation of the present cell within
   * this function, i.e. the returned
   * neighbor's child is behind subface @p
   * subface concerning the intrinsic
   * ordering of the given face.
   */
  TriaIterator<CellAccessor<dim, spacedim> >
  neighbor_child_on_subface (const unsigned int face_no,
                             const unsigned int subface_no) const;

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
  TriaIterator<CellAccessor<dim, spacedim>  >
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
   * Return the how-many'th
   * neighbor this cell is of
   * <tt>cell->neighbor(neighbor)</tt>,
   * i.e. return the @p face_no
   * such that
   * <tt>cell->neighbor(neighbor)->neighbor(face_no)==cell</tt>. This
   * function is the right one if
   * you want to know how to get
   * back from a neighbor to the
   * present cell.
   *
   * Note that this operation is
   * only useful if the neighbor is
   * not coarser than the present
   * cell. If the neighbor is
   * coarser this function throws
   * an exception. Use the @p
   * neighbor_of_coarser_neighbor
   * function in that case.
   */
  unsigned int neighbor_of_neighbor (const unsigned int neighbor) const;

  /**
   * Return, whether the neighbor
   * is coarser then the present
   * cell. This is important in
   * case of ansiotropic
   * refinement where this
   * information does not depend on
   * the levels of the cells.
   *
   * Note, that in an anisotropic
   * setting, a cell can only be
   * coarser than another one at a
   * given face, not on a general
   * basis. The face of the finer
   * cell is contained in the
   * corresponding face of the
   * coarser cell, the finer face
   * is either a child or a
   * grandchild of the coarser
   * face.
   */
  bool neighbor_is_coarser (const unsigned int neighbor) const;

  /**
   * This function is a generalization of the
   * @p neighbor_of_neighbor function for the
   * case of a coarser neighbor. It returns a
   * pair of numbers, face_no and subface_no,
   * with the following property, if the
   * neighbor is not refined:
   * <tt>cell->neighbor(neighbor)->neighbor_child_on_subface(face_no,subface_no)==cell</tt>.
   * In 3D, a coarser neighbor can still be
   * refined. In that case subface_no denotes the child index of the neighbors face that relates to our face:
   * <tt>cell->neighbor(neighbor)->face(face_no)->child(subface_no)==cell->face(neighbor)</tt>.
   * This case in 3d and how it can happen
   * is discussed in the introduction of the
   * step-30 tutorial program.
   *
   * This function is impossible
   * for <tt>dim==1</tt>.
   */
  std::pair<unsigned int, unsigned int>
  neighbor_of_coarser_neighbor (const unsigned int neighbor) const;

  /**
   * This function is a generalization of the
   * @p neighbor_of_neighbor and the @p
   * neighbor_of_coarser_neighbor
   * functions. It checks whether the
   * neighbor is coarser or not and calls the
   * respective function. In both cases, only
   * the face_no is returned.
   */
  unsigned int neighbor_face_no (const unsigned int neighbor) const;

  /**
   * @}
   */

  /**
   *  @name Dealing with boundary indicators
   */
  /**
   * @{
   */

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
   * for 1 and 2 dimensions,
   * it is equivalent, for three
   * dimensions it returns
   * whether at least one of the 12
   * lines of the hexahedron is at
   * a boundary. This, of course,
   * includes the case where a
   * whole face is at the boundary,
   * but also some other cases.
   */
  bool has_boundary_lines () const;
  /**
   * @}
   */

  /**
   *  @name Dealing with refinement indicators
   */
  /**
   * @{
   */

  /**
   * Return the @p RefinementCase<dim> this cell was flagged to be refined
   * with.  The return value of this function can be compared to a bool to
   * check if this cell is flagged for any kind of refinement. For example, if
   * you have previously called cell->set_refine_flag() for a cell, then you
   * will enter the 'if' block in the following snippet:
   *
   * @code
   * if (cell->refine_flag_set())
   * {
   *   // yes, this cell is marked for refinement.
   * }
   * @endcode
   */
  RefinementCase<dim> refine_flag_set () const;

  /**
   *  Flag the cell pointed to for
   *  refinement. This function is
   *  only allowed for active
   *  cells. Keeping the default value for @p ref_case will mark this cell for
   *  isotropic refinement.
   */
  void set_refine_flag (const RefinementCase<dim> ref_case = RefinementCase<dim>::isotropic_refinement) const;

  /**
   *  Clear the refinement flag.
   */
  void clear_refine_flag () const;

  /**
   *  Modify the refinement flag of the cell
   *  to ensure (at least) the given
   *  refinement case @p face_refinement_case at
   *  face <tt>face_no</tt>, taking into
   *  account orientation, flip and rotation
   *  of the face. Return, whether the
   *  refinement flag had to be
   *  modified. This function is only allowed
   *  for active cells.
   */
  bool flag_for_face_refinement (const unsigned int face_no,
                                 const RefinementCase<dim-1> &face_refinement_case=RefinementCase<dim-1>::isotropic_refinement) const;

  /**
   *  Modify the refinement flag of the cell
   *  to ensure that line <tt>face_no</tt>
   *  will be refined. Return, whether the
   *  refinement flag had to be
   *  modified. This function is only allowed
   *  for active cells.
   */
  bool flag_for_line_refinement (const unsigned int line_no) const;

  /**
   * Return the SubfaceCase of face
   * <tt>face_no</tt>. Note that this is not
   * identical to asking
   * <tt>cell->face(face_no)->refinement_case()</tt>
   * since the latter returns a RefinementCase<dim-1>
   * and thus only considers one
   * (anisotropic) refinement, whereas this
   * function considers the complete
   * refinement situation including possible
   * refinement of the face's children. This
   * function may only be called for active
   * cells in 2d and 3d.
   */
  dealii::internal::SubfaceCase<dim> subface_case(const unsigned int face_no) const;

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
   * @}
   */

  /**
   *  @name Dealing with material indicators
   */
  /**
   * @{
   */

  /**
   * Return the material id of this
   * cell.
   *
   * For a typical use of this
   * function, see the @ref step_28
   * "step-28" tutorial program.
   *
   * See the @ref GlossMaterialId
   * "glossary" for more
   * information.
   */
  types::material_id material_id () const;

  /**
   * Set the material id of this
   * cell.
   *
   * For a typical use of this
   * function, see the @ref step_28
   * "step-28" tutorial program.
   *
   * See the @ref GlossMaterialId
   * "glossary" for more
   * information.
   */
  void set_material_id (const types::material_id new_material_id) const;

  /**
   * Set the material id of this
   * cell and all its children (and
   * grand-children, and so on) to
   * the given value.
   *
   * See the @ref GlossMaterialId
   * "glossary" for more
   * information.
   */
  void recursively_set_material_id (const types::material_id new_material_id) const;
  /**
   * @}
   */

  /**
   *  @name Dealing with subdomain indicators
   */
  /**
   * @{
   */

  /**
   * Return the subdomain id of
   * this cell.
   *
   * See the @ref GlossSubdomainId
   * "glossary" for more
   * information.
   *
   * @note The subdomain of a cell is a property only defined
   * for active cells, i.e., cells that are not further
   * refined. Consequently, you can only call this function if
   * the cell it refers to has no children. For multigrid
   * methods in parallel, it is also important to know which
   * processor owns non-active cells, and for this you can call
   * level_subdomain_id().
   */
  types::subdomain_id subdomain_id () const;

  /**
   * Set the subdomain id of this
   * cell.
   *
   * See the @ref GlossSubdomainId
   * "glossary" for more
   * information. This function
   * should not be called if you
   * use a
   * parallel::distributed::Triangulation
   * object.
   *
   * @note The subdomain of a cell is a property only defined
   * for active cells, i.e., cells that are not further
   * refined. Consequently, you can only call this function if
   * the cell it refers to has no children. For multigrid
   * methods in parallel, it is also important to know which
   * processor owns non-active cells, and for this you can call
   * level_subdomain_id().
   */
  void set_subdomain_id (const types::subdomain_id new_subdomain_id) const;

  /**
   * Get the level subdomain id of this cell. This is used for parallel multigrid.
   */
  types::subdomain_id level_subdomain_id () const;

  /**
   * Set the level subdomain id of this cell. This is used for parallel multigrid.
   */
  void set_level_subdomain_id (const types::subdomain_id new_level_subdomain_id) const;


  /**
   * Set the subdomain id of this cell (if it is active) or all its terminal
   * children (and grand-children, and so on, as long as they have no children
   * of their own) to the given value. Since the subdomain id is a concept
   * that is only defined for cells that are active (i.e., have no children
   * of their own), this function only sets the subdomain ids for all
   * children and grand children of this cell that are actually active,
   * skipping intermediate child cells.
   *
   * See the @ref GlossSubdomainId
   * "glossary" for more
   * information. This function
   * should not be called if you
   * use a
   * parallel::distributed::Triangulation
   * object since there the subdomain id is implicitly defined by which
   * processor you're on.
   */
  void recursively_set_subdomain_id (const types::subdomain_id new_subdomain_id) const;
  /**
   * @}
   */

  /**
   *  @name Dealing with codim 1 cell orientation
   */
  /**
   * @{
   */

  /**
   * Return the orientation of
   * this cell.
   *
   * For the meaning of this flag, see
   * @ref GlossDirectionFlag .
   */
  bool direction_flag () const;

  /**
   *  Index of the parent of this cell within the level of the triangulation
   *  to which the parent cell belongs. The level of the parent is of course
   *  one lower than that of the present cell. If the parent does not exist
   *  (i.e., if the object is at the coarsest level of the mesh hierarchy), an
   *  exception is generated.
   */
  int parent_index () const;

  /**
   *  Return an iterator to the
   *  parent. If the
   *  parent does not exist (i.e., if the object is at the coarsest level of
   *  the mesh hierarchy), an exception is generated.
   */
  TriaIterator<CellAccessor<dim,spacedim> >
  parent () const;

  /**
   * @}
   */

  /**
   *  @name Other functions
   */
  /**
   * @{
   */

  /**
   * Test whether the cell has children
   * (this is the criterion for activity
   * of a cell).
   *
   * See the @ref GlossActive "glossary"
   * for more information.
   */
  bool active () const;

  /**
   * Return whether this cell is owned by the current processor
   * or is owned by another processor. The function always returns
   * true if applied to an object of type dealii::Triangulation,
   * but may yield false if the triangulation is of type
   * parallel::distributed::Triangulation.
   *
   * See the @ref GlossGhostCell
   * "glossary" and the @ref
   * distributed module for more
   * information.
   *
   * @post The returned value is equal to <code>!is_ghost() &&
   * !is_artificial()</code>.
   *
   * @note Whether a cell is a ghost cell, artificial, or is
   * locally owned or is a property that only pertains to cells
   * that are active. Consequently,
   * you can only call this function if the cell it refers to has
   * no children.
   */
  bool is_locally_owned () const;

  /**
   * Return true if either the Triangulation is not distributed or if
   * level_subdomain_id() is equal to the id of the current processor.
   */
  bool is_locally_owned_on_level () const;

  /**
   * Return whether this cell
   * exists in the global mesh but
   * (i) is owned by another
   * processor, i.e. has a
   * subdomain_id different from
   * the one the current processor
   * owns and (ii) is adjacent to a
   * cell owned by the current
   * processor.
   *
   * This function only makes sense
   * if the triangulation used is
   * of kind
   * parallel::distributed::Triangulation. In
   * all other cases, the returned
   * value is always false.
   *
   * See the @ref GlossGhostCell
   * "glossary" and the @ref
   * distributed module for more
   * information.
   *
   * @post The returned value is equal to
   * <code>!is_locally_owned() &&
   * !is_artificial()</code>.
   *
   * @note Whether a cell is a ghost cell, artificial, or is
   * locally owned or is a property that only pertains to cells
   * that are active. Consequently,
   * you can only call this function if the cell it refers to has
   * no children.
   */
  bool is_ghost () const;

  /**
   * Return whether this cell is
   * artificial, i.e. it isn't one
   * of the cells owned by the
   * current processor, and it also
   * doesn't border on one. As a
   * consequence, it exists in the
   * mesh to ensure that each
   * processor has all coarse mesh
   * cells and that the 2:1 ratio
   * of neighboring cells is
   * maintained, but it is not one
   * of the cells we should work on
   * on the current processor. In
   * particular, there is no
   * guarantee that this cell
   * isn't, in fact, further
   * refined on one of the other
   * processors.
   *
   * This function only makes sense
   * if the triangulation used is
   * of kind
   * parallel::distributed::Triangulation. In
   * all other cases, the returned
   * value is always false.
   *
   * See the @ref
   * GlossArtificialCell "glossary"
   * and the @ref distributed
   * module for more information.
   *
   * @post The returned value is equal to
   * <code>!is_ghost() &&
   * !is_artificial()</code>.
   *
   * @note Whether a cell is a ghost cell, artificial, or is
   * locally owned or is a property that only pertains to cells
   * that are active. Consequently,
   * you can only call this function if the cell it refers to has
   * no children.
   */
  bool is_artificial () const;

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
   *
   * In case of codim>0, the point is first projected
   * to the manifold where the cell is embedded and
   * then check if this projection is inside the cell.
   */
  bool point_inside (const Point<spacedim> &p) const;

  /**
   *  Set the neighbor @p i of
   *  this cell to the cell pointed
   *  to by @p pointer.
   *
   * This function shouldn't really be
   * public (but needs to for various
   * reasons in order not to make a long
   * list of functions friends): it
   * modifies internal data structures and
   * may leave things. Do not use it from
   * application codes.
   */
  void set_neighbor (const unsigned int i,
                     const TriaIterator<CellAccessor<dim, spacedim> > &pointer) const;

  /**
   * Return a unique ID for the current cell. This ID is constructed from the
   * path in the hierarchy from the coarse father cell and works correctly in
   * parallel computations using objects of type
   * parallel::distributed::Triangulation. This function is therefore useful
   * in providing a unique identifier for cells (active or not) that also
   * works for parallel triangulations. See the documentation of the CellId
   * class for more information.
   *
   * @note This operation takes O(level) time to compute. In most practicaly
   * cases, the number of levels of a triangulation will depend
   * logarithmically on the number of cells in the triangulation.
   */
  CellId id() const;

  /**
   * @}
   */


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

protected:
  /**
   * This function assumes that the
   * neighbor is not coarser than
   * the current cell. In this case
   * it returns the
   * neighbor_of_neighbor() value.
   * If, however, the neighbor is
   * coarser this function returns
   * an
   * <code>invalid_unsigned_int</code>.
   *
   * This function is not for
   * public use. Use the function
   * neighbor_of_neighbor() instead
   * which throws an exception if
   * called for a coarser
   * neighbor. If neighbor is
   * indeed coarser (you get to
   * know this by e.g. the
   * neighbor_is_coarser()
   * function) then the
   * neighbor_of_coarser_neighbor()
   * function should be call. If
   * you'd like to know only the
   * <code>face_no</code> which is
   * required to get back from the
   * neighbor to the present cell
   * then simply use the
   * neighbor_face_no() function
   * which can be used for coarser
   * as well as noncoarser
   * neighbors.
   */
  unsigned int neighbor_of_neighbor_internal (const unsigned int neighbor) const;

  /**
     As for any codim>0 we can use a similar code
     and c++ does not allow partial templates.
     we use this auxiliary function that is then
     called from point_inside.
  */
  template<int dim_,int spacedim_ >
  bool point_inside_codim(const Point<spacedim_> &p) const;



private:
  /**
   * Set the parent of a cell.
   */
  void set_parent (const unsigned int parent_index);

  /**
   * Set the orientation of this
   * cell.
   *
   * For the meaning of this flag, see
   * @ref GlossDirectionFlag .
   */
  void set_direction_flag (const bool new_direction_flag) const;
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
  void operator = (const CellAccessor<dim, spacedim> &);

  template <int, int> friend class Triangulation;

  friend struct dealii::internal::Triangulation::Implementation;
};



/* -------------- declaration of explicit
   specializations and general templates ------------- */


template <int structdim, int dim, int spacedim>
template <typename OtherAccessor>
InvalidAccessor<structdim, dim, spacedim>::
InvalidAccessor (const OtherAccessor &)
{
  Assert (false,
          ExcMessage ("You are attempting an illegal conversion between "
                      "iterator/accessor types. The constructor you call "
                      "only exists to make certain template constructs "
                      "easier to write as dimension independent code but "
                      "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
TriaAccessor<structdim,dim,spacedim>::
TriaAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &)
{
  Assert (false,
          ExcMessage ("You are attempting an illegal conversion between "
                      "iterator/accessor types. The constructor you call "
                      "only exists to make certain template constructs "
                      "easier to write as dimension independent code but "
                      "the conversion is not valid in the current context."));
}



template <int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
CellAccessor<dim,spacedim>::
CellAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &)
{
  Assert (false,
          ExcMessage ("You are attempting an illegal conversion between "
                      "iterator/accessor types. The constructor you call "
                      "only exists to make certain template constructs "
                      "easier to write as dimension independent code but "
                      "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
TriaAccessor<structdim,dim,spacedim>::
TriaAccessor (const TriaAccessor<structdim2,dim2,spacedim2> &)
{
  Assert (false,
          ExcMessage ("You are attempting an illegal conversion between "
                      "iterator/accessor types. The constructor you call "
                      "only exists to make certain template constructs "
                      "easier to write as dimension independent code but "
                      "the conversion is not valid in the current context."));
}



template <int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
CellAccessor<dim,spacedim>::
CellAccessor (const TriaAccessor<structdim2,dim2,spacedim2> &)
{
  Assert (false,
          ExcMessage ("You are attempting an illegal conversion between "
                      "iterator/accessor types. The constructor you call "
                      "only exists to make certain template constructs "
                      "easier to write as dimension independent code but "
                      "the conversion is not valid in the current context."));
}

template <int dim, int spacedim>
CellId
CellAccessor<dim,spacedim>::id() const
{
  std::vector<unsigned char> id(this->level(), -1);
  unsigned int coarse_index;

  CellAccessor<dim,spacedim> ptr = *this;
  while (ptr.level()>0)
    {
      // find the 'v'st child of our parent we are
      unsigned char v=-1;
      for (unsigned int c=0; c<ptr.parent()->n_children(); ++c)
        {
          if (ptr.parent()->child_index(c)==ptr.index())
            {
              v = c;
              break;
            }
        }

      Assert(v != (unsigned char)-1, ExcInternalError());
      id[ptr.level()-1] = v;

      ptr.copy_from( *(ptr.parent()));
    }

  Assert(ptr.level()==0, ExcInternalError());
  coarse_index = ptr.index();

  return CellId(coarse_index, id);
}


#ifndef DOXYGEN

template <> bool CellAccessor<1,1>::point_inside (const Point<1> &) const;
template <> bool CellAccessor<2,2>::point_inside (const Point<2> &) const;
template <> bool CellAccessor<3,3>::point_inside (const Point<3> &) const;
template <> bool CellAccessor<1,2>::point_inside (const Point<2> &) const;
template <> bool CellAccessor<1,3>::point_inside (const Point<3> &) const;
template <> bool CellAccessor<2,3>::point_inside (const Point<3> &) const;
// -------------------------------------------------------------------

template <> void TriaAccessor<3,3,3>::set_all_manifold_ids (const types::manifold_id) const;

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

// include more templates in debug and optimized mode
#  include "tria_accessor.templates.h"

#endif
