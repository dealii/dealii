// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tria_accessor_h
#define dealii_tria_accessor_h


#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria_faces.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_iterator_base.h>
#include <deal.II/grid/tria_iterator_selector.h>
#include <deal.II/grid/tria_levels.h>

#include <boost/container/small_vector.hpp>

#include <cmath>
#include <limits>
#include <utility>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class Triangulation;
template <typename Accessor>
class TriaRawIterator;
template <typename Accessor>
class TriaIterator;
template <typename Accessor>
class TriaActiveIterator;

namespace parallel
{
  template <int dim, int spacedim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  class TriangulationBase;
}

template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;
template <int dim, int spacedim, bool lda>
class DoFCellAccessor;


template <int dim, int spacedim>
class Manifold;

template <int dim, int spacedim>
class Mapping;
#endif

namespace internal
{
  namespace TriangulationImplementation
  {
    class TriaObjects;
    struct Implementation;
    struct ImplementationMixedMesh;
  } // namespace TriangulationImplementation

  namespace TriaAccessorImplementation
  {
    struct Implementation;

    /**
     * Implementation of a type with which to store the level of an accessor
     * object. We only need it for the case that <tt>structdim == dim</tt>.
     * Otherwise, an empty object is sufficient.
     */
    template <int structdim, int dim>
    struct PresentLevelType
    {
      struct type
      {
        /**
         * Default constructor.
         */
        type() = default;

        /**
         * Dummy constructor. Only level zero is allowed.
         */
        type(const int level)
        {
          Assert(level == 0, ExcInternalError());
          (void)level; // removes -Wunused-parameter warning in optimized mode
        }

        /**
         * Dummy conversion operator. Returns level zero.
         */
        operator int() const
        {
          return 0;
        }

        void
        operator++() const
        {
          DEAL_II_ASSERT_UNREACHABLE();
        }

        void
        operator--() const
        {
          DEAL_II_ASSERT_UNREACHABLE();
        }
      };
    };


    /**
     * Implementation of a type with which to store the level of an accessor
     * object. We only need it for the case that <tt>structdim == dim</tt>.
     * Otherwise, an empty object is sufficient.
     */
    template <int dim>
    struct PresentLevelType<dim, dim>
    {
      using type = int;
    };
  } // namespace TriaAccessorImplementation
} // namespace internal
template <int structdim, int dim, int spacedim>
class TriaAccessor;
template <int dim, int spacedim>
class TriaAccessor<0, dim, spacedim>;
template <int spacedim>
class TriaAccessor<0, 1, spacedim>;

/**
 * A namespace that contains exception classes used by the accessor classes.
 */
namespace TriaAccessorExceptions
{
  /**
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcCellNotUsed,
                   "The operation you are attempting can only be performed for "
                   "(cell, face, or edge) iterators that point to valid "
                   "objects. These objects need not necessarily be active, "
                   "i.e., have no children, but they need to be part of a "
                   "triangulation. (The objects pointed to by an iterator "
                   "may -- after coarsening -- also be objects that used "
                   "to be part of a triangulation, but are now no longer "
                   "used. Their memory location may have been retained "
                   "for re-use upon the next mesh refinement, but is "
                   "currently unused.)");
  /**
   * The cell is not an
   * @ref GlossActive "active"
   * cell, but it already has children. Some operations, like setting
   * refinement flags or accessing degrees of freedom are only possible on
   * active cells.
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcCellNotActive,
                   "The operation you are attempting can only be performed for "
                   "(cell, face, or edge) iterators that point to 'active' "
                   "objects. 'Active' objects are those that do not have "
                   "children (in the case of cells), or that are part of "
                   "an active cell (in the case of faces or edges). However, "
                   "the object on which you are trying the current "
                   "operation is not 'active' in this sense.");
  /**
   * Trying to access the children of a cell which is in fact active.
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcCellHasNoChildren,
                   "The operation you are attempting can only be performed for "
                   "(cell, face, or edge) iterators that have children, "
                   "but the object on which you are trying the current "
                   "operation does not have any.");
  /**
   * Trying to access the parent of a cell which is in the coarsest level of
   * the triangulation.
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcCellHasNoParent,
                   "The operation you are attempting can only be performed for "
                   "(cell, face, or edge) iterators that have a parent object, "
                   "but the object on which you are trying the current "
                   "operation does not have one -- i.e., it is on the "
                   "coarsest level of the triangulation.");
  /**
   * @ingroup Exceptions
   */
  DeclException1(ExcCantSetChildren,
                 int,
                 << "You can only set the child index if the cell does not "
                 << "currently have children registered; or you can clear it. "
                 << "The given index was " << arg1
                 << " (-1 means: clear children).");
  /**
   * @ingroup Exceptions
   */
  template <typename AccessorType>
  DeclException1(ExcDereferenceInvalidObject,
                 AccessorType,
                 << "You tried to dereference an iterator for which this "
                 << "is not possible. More information on this iterator: "
                 << "index=" << arg1.index() << ", state="
                 << (arg1.state() == IteratorState::valid ?
                       "valid" :
                       (arg1.state() == IteratorState::past_the_end ?
                          "past_the_end" :
                          "invalid")));
  /**
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcCantCompareIterators,
                   "Iterators can only be compared if they point to the same "
                   "triangulation, or if neither of them are associated "
                   "with a triangulation.");
  // TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException0(ExcNeighborIsCoarser);
  // TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException0(ExcNeighborIsNotCoarser);
  /**
   * You are trying to access the level of a face, but faces have no inherent
   * level. The level of a face can only be determined by the level of an
   * adjacent face, which in turn implies that a face can have several levels.
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcFacesHaveNoLevel);
  /**
   * You are trying to get the periodic neighbor for a face, which does not
   * have a periodic neighbor. For more information on this, refer to
   * @ref GlossPeriodicConstraints "entry for periodic boundaries".
   * @ingroup Exceptions
   */
  DeclException0(ExcNoPeriodicNeighbor);
  // TODO: Write documentation!
  /**
   * @ingroup Exceptions
   */
  DeclException1(
    ExcSetOnlyEvenChildren,
    int,
    << "You can only set the child index of an even numbered child."
    << "The number of the child given was " << arg1 << '.');
} // namespace TriaAccessorExceptions


/**
 * A base class for the accessor classes used by TriaRawIterator and derived
 * classes.
 *
 * This class offers only the basic functionality required by the iterators
 * (stores the necessary data members, offers comparison operators and the
 * like), but has no functionality to actually dereference data. This is done
 * in the derived classes.
 *
 * In the implementation, the behavior of this class differs between the cases
 * where <tt>structdim==dim</tt> (cells of a mesh) and
 * <tt>structdim&lt;dim</tt> (faces and edges). For the latter, #present_level
 * is always equal to zero and the constructors may not receive a positive
 * value there. For cells, any level is possible, but only those within the
 * range of the levels of the Triangulation are reasonable. Furthermore, the
 * function objects() returns either the container with all cells on the same
 * level or the container with all objects of this dimension
 * (<tt>structdim&lt;dim</tt>).
 *
 * Some internals of this class are discussed in
 * @ref IteratorAccessorInternals.
 *
 * @ingroup grid
 * @ingroup Accessors
 */
template <int structdim, int dim, int spacedim = dim>
class TriaAccessorBase
{
public:
  /**
   * Dimension of the space the object represented by this accessor lives in.
   * For example, if this accessor represents a quadrilateral that is part of
   * a two-dimensional surface in four-dimensional space, then this value is
   * four.
   */
  static constexpr unsigned int space_dimension = spacedim;

  /**
   * Dimensionality of the object that the thing represented by this accessor
   * is part of. For example, if this accessor represents a line that is part
   * of a hexahedron, then this value will be three.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Dimensionality of the current object represented by this accessor. For
   * example, if it is line (irrespective of whether it is part of a 2d or 3d
   * subobject), then this value equals 1.
   */
  static const unsigned int structure_dimension = structdim;

  /**
   * Copy operator. These operators are usually used in a context like
   * <tt>iterator a,b; *a=*b;</tt>. Presumably, the intent here is to copy the
   * object pointed to
   * by @p b to the object pointed to by @p a. However, the result of
   * dereferencing an iterator is not an object but an accessor; consequently,
   * this operation is not useful for iterators on triangulations.
   * Consequently, this operator is declared as deleted and can not be used.
   */
  void
  operator=(const TriaAccessorBase *) = delete;

protected:
  /**
   * Declare the data type that this accessor class expects to get passed from
   * the iterator classes. Since the pure triangulation iterators need no
   * additional data, this data type is @p void.
   */
  using AccessorData = void;

  /**
   * Constructor. Protected, thus only callable from friend classes.
   */
  TriaAccessorBase(const Triangulation<dim, spacedim> *parent = nullptr,
                   const int                           level  = -1,
                   const int                           index  = -1,
                   const AccessorData                       * = nullptr);

  /**
   * Copy constructor. Creates an object with exactly the same data.
   */
  TriaAccessorBase(const TriaAccessorBase &);

  /**
   * Copy operator. Since this is only called from iterators, do not return
   * anything, since the iterator will return itself.
   *
   * This method is protected, since it is only to be called from the iterator
   * class.
   */
  void
  copy_from(const TriaAccessorBase &);

  /**
   * Copy operator. Creates an object with exactly the same data.
   */
  TriaAccessorBase &
  operator=(const TriaAccessorBase &);

  /**
   * Comparison operator for accessors. This operator is used when comparing
   * iterators into objects of a triangulation, for example when putting
   * them into a `std::map`.
   *
   * If #structure_dimension is less than #dimension, we simply compare the
   * index of such an object because faces and edges do not have levels. If
   * #structure_dimension equals #dimension, we compare the level first, and
   * the index only if levels are equal.
   */
  bool
  operator<(const TriaAccessorBase &other) const;

protected:
  /**
   * Compare for equality.
   */
  bool
  operator==(const TriaAccessorBase &) const;

  /**
   * Compare for inequality.
   */
  bool
  operator!=(const TriaAccessorBase &) const;

  /**
   * @name Advancement of iterators
   */
  /**
   * @{
   */
  /**
   * This operator advances the iterator to the next element.
   *
   * For @p dim=1 only: The next element is next on this level if there are
   * more. If the present element is the last on this level, the first on the
   * next level is accessed.
   */
  void
  operator++();

  /**
   * This operator moves the iterator to the previous element.
   *
   * For @p dim=1 only: The previous element is previous on this level if
   * <tt>index>0</tt>. If the present element is the first on this level, the
   * last on the previous level is accessed.
   */
  void
  operator--();
  /**
   * @}
   */

  /**
   * Access to the other objects of a Triangulation with same dimension.
   */
  dealii::internal::TriangulationImplementation::TriaObjects &
  objects() const;

public:
  /**
   * Data type to be used for passing parameters from iterators to the
   * accessor classes in a unified way, no matter what the type of number of
   * these parameters is.
   */
  using LocalData = void *;

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
   * refinement level and the latter is the index of the cell within this
   * refinement level (the former being what this function returns).
   * Consequently, there may be multiple cells on different refinement levels
   * but with the same index within their level. Contrary to this, if the
   * current object corresponds to a face or edge, then the object is uniquely
   * identified solely by its index as faces and edges do not have a
   * refinement level. For these objects, the current function always returns
   * zero as the level.
   */
  int
  level() const;

  /**
   * Return the index of the element presently pointed to on the present
   * level.
   *
   * Within a Triangulation object, cells are uniquely identified by a pair
   * <code>(level, index)</code> where the former is the cell's refinement
   * level and the latter is the index of the cell within this refinement
   * level (the latter being what this function returns). Consequently, there
   * may be multiple cells on different refinement levels but with the same
   * index within their level. Contrary to this, if the current object
   * corresponds to a face or edge, then the object is uniquely identified
   * solely by its index as faces and edges do not have a refinement level.
   *
   * @note The indices objects returned by this function are not a contiguous
   * set of numbers on each level: going from cell to cell, some of the
   * indices in a level may be unused.
   *
   * @note If the triangulation is actually of type
   * parallel::distributed::Triangulation then the indices are relatively only
   * to that part of the distributed triangulation that is stored on the
   * current processor. In other words, cells living in the partitions of the
   * triangulation stored on different processors may have the same index even
   * if they refer to the same cell, and the may have different indices even
   * if they do refer to the same cell (e.g., if a cell is owned by one
   * processor but is a ghost cell on another).
   */
  int
  index() const;

  /**
   * Return the state of the iterator.  For the different states an accessor
   * can be in, refer to the TriaRawIterator documentation.
   */
  IteratorState::IteratorStates
  state() const;

  /**
   * Return a reference to the triangulation which the object pointed to by this
   * class belongs to.
   */
  const Triangulation<dim, spacedim> &
  get_triangulation() const;

  /**
   * @}
   */
protected:
  /**
   * The level if this is a cell (<tt>structdim==dim</tt>). Else, contains
   * zero.
   */
  typename dealii::internal::TriaAccessorImplementation::
    PresentLevelType<structdim, dim>::type present_level;

  /**
   * Used to store the index of the element presently pointed to on the level
   * presently used.
   */
  int present_index;

  /**
   * Pointer to the triangulation which we act on.
   */
  const Triangulation<dim, spacedim> *tria;

private:
  template <typename Accessor>
  friend class TriaRawIterator;
  template <typename Accessor>
  friend class TriaIterator;
  template <typename Accessor>
  friend class TriaActiveIterator;
};



/**
 * A class that represents accessor objects to iterators that don't make sense
 * such as quad iterators in on 1d meshes.  This class can not be used to
 * create objects (it will in fact throw an exception if this should ever be
 * attempted but it sometimes allows code to be written in a simpler way in a
 * dimension independent way. For example, it allows to write code that works
 * on quad iterators that is dimension independent -- i.e., also compiles
 * in 1d -- because quad iterators
 * (via the current class) exist and are syntactically correct. You can not
 * expect, however, to ever create an actual object of one of these iterators
 * in 1d, meaning you need to expect to wrap the code block in which you use
 * quad iterators into something like <code>if (dim@>1)</code> -- which makes
 * eminent sense anyway.
 *
 * This class provides the minimal interface necessary for Accessor classes to
 * interact with Iterator classes. However, this is only for syntactic
 * correctness, none of the functions do anything but generate errors.
 *
 * @ingroup Accessors
 */
template <int structdim, int dim, int spacedim = dim>
class InvalidAccessor
{
public:
  /**
   * Dimension of the space the object represented by this accessor lives in.
   * For example, if this accessor represents a quadrilateral that is part of
   * a two-dimensional surface in four-dimensional space, then this value is
   * four.
   */
  static constexpr unsigned int space_dimension = spacedim;

  /**
   * Dimensionality of the object that the thing represented by this accessor
   * is part of. For example, if this accessor represents a line that is part
   * of a hexahedron, then this value will be three.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Dimensionality of the current object represented by this accessor. For
   * example, if it is line (irrespective of whether it is part of a 2d or 3d
   * subobject), then this value equals 1.
   */
  static const unsigned int structure_dimension = structdim;

  /**
   * Declare the data type that this accessor class expects to get passed from
   * the iterator classes. Since the pure triangulation iterators need no
   * additional data, this data type is @p void.
   */
  using AccessorData = void;

  /**
   * Constructor.  This class is used for iterators that do not make
   * sense in a given dimension, for example quads for 1d meshes. Consequently,
   * while the creation of such objects is syntactically valid, they make no
   * semantic sense, and we generate an exception when such an object is
   * actually generated.
   */
  InvalidAccessor(const void         *parent     = nullptr,
                  const int           level      = -1,
                  const int           index      = -1,
                  const AccessorData *local_data = nullptr);

  /**
   * Copy constructor.  This class is used for iterators that do not make
   * sense in a given dimension, for example quads for 1d meshes. Consequently,
   * while the creation of such objects is syntactically valid, they make no
   * semantic sense, and we generate an exception when such an object is
   * actually generated.
   */
  InvalidAccessor(const InvalidAccessor &);

  /**
   * Conversion from other accessors to the current invalid one. This of
   * course also leads to a run-time error.
   */
  template <typename OtherAccessor>
  InvalidAccessor(const OtherAccessor &);

  /**
   * Dummy copy operation.
   */
  void
  copy_from(const InvalidAccessor &);

  /**
   * Dummy comparison operators.
   */
  bool
  operator==(const InvalidAccessor &) const;
  bool
  operator!=(const InvalidAccessor &) const;

  /**
   * Dummy operators to make things compile. Does nothing.
   */
  void
  operator++() const;
  void
  operator--() const;

  /**
   * Return the state of the iterator.  For the different states an accessor
   * can be in, refer to the TriaRawIterator documentation.
   */
  static IteratorState::IteratorStates
  state();


  /**
   * Level of this object. Vertices have no level, so this function always
   * returns zero.
   */
  static int
  level();

  /**
   * Index of this object. Returns the global index of the vertex this object
   * points to.
   */
  static int
  index();

  /**
   * Dummy function representing whether the accessor points to a used or an
   * unused object.
   */
  bool
  used() const;

  /**
   * Dummy function representing whether the accessor points to an object that
   * has children.
   */
  bool
  has_children() const;

  /**
   * Dummy function that always returns numbers::flat_manifold_id.
   */
  types::manifold_id
  manifold_id() const;

  /**
   * Dummy function that always returns numbers::invalid_unsigned_int.
   */
  unsigned int
  user_index() const;

  /**
   * Dummy function that always throws.
   */
  void
  set_user_index(const unsigned int p) const;

  /**
   * Dummy function that always throws.
   */
  void
  set_manifold_id(const types::manifold_id) const;

  /**
   * Dummy function to extract vertices. Returns the origin.
   */
  Point<spacedim> &
  vertex(const unsigned int i) const;

  /**
   * Dummy function to extract lines. Returns a default-constructed line
   * iterator.
   */
  void *
  line(const unsigned int i) const;

  /**
   * Dummy function to extract quads. Returns a default-constructed quad
   * iterator.
   */
  void *
  quad(const unsigned int i) const;
};



/**
 * A class that provides access to objects in a triangulation such as its
 * vertices, sub-objects, children, geometric information, etc. This class
 * represents objects of dimension <code>structdim</code> (i.e. 1 for lines, 2
 * for quads, 3 for hexes) in a triangulation of dimensionality
 * <code>dim</code> (i.e. 1 for a triangulation of lines, 2 for a
 * triangulation of quads, and 3 for a triangulation of hexes) that is
 * embedded in a space of dimensionality <code>spacedim</code> (for
 * <code>spacedim==dim</code> the triangulation represents a domain in
 * $R^{dim}$, for <code>spacedim@>dim</code> the triangulation is of a
 * manifold embedded in a higher dimensional space).
 *
 * There is a specialization of this class for the case where
 * @p structdim equals zero, i.e., for vertices of a triangulation.
 *
 * @ingroup Accessors
 */
template <int structdim, int dim, int spacedim>
class TriaAccessor : public TriaAccessorBase<structdim, dim, spacedim>
{
public:
  /**
   * Propagate alias from base class to this class.
   */
  using AccessorData =
    typename TriaAccessorBase<structdim, dim, spacedim>::AccessorData;

  /**
   * Constructor.
   */
  TriaAccessor(const Triangulation<dim, spacedim> *parent     = nullptr,
               const int                           level      = -1,
               const int                           index      = -1,
               const AccessorData                 *local_data = nullptr);

  /**
   * The copy constructor is not deleted but copied constructed elements should
   * not be modified, also the comments to the copy assignment operator.
   */
  TriaAccessor(const TriaAccessor &) = default;

  /**
   * Move constructor.
   */
  TriaAccessor(TriaAccessor &&) = default; // NOLINT

  /**
   * Conversion constructor. This constructor exists to make certain
   * constructs simpler to write in dimension independent code. For example,
   * it allows assigning a face iterator to a line iterator, an operation that
   * is useful in 2d but doesn't make any sense in 3d. The constructor here
   * exists for the purpose of making the code conform to C++ but it will
   * unconditionally abort; in other words, assigning a face iterator to a
   * line iterator is better put into an if-statement that checks that the
   * dimension is two, and assign to a quad iterator in 3d (an operator that,
   * without this constructor would be illegal if we happen to compile for
   * 2d).
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * Another conversion operator between objects that don't make sense, just
   * like the previous one.
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor(const TriaAccessor<structdim2, dim2, spacedim2> &);

  /**
   * Copy operator. These operators are usually used in a context like
   * <tt>iterator a,b; *a=*b;</tt>. Presumably, the intent here is to copy the
   * object pointed to
   * by @p b to the object pointed to by @p a. However, the result of
   * dereferencing an iterator is not an object but an accessor; consequently,
   * this operation is not useful for iterators on triangulations.
   * Consequently, this operator is declared as deleted and can not be used.
   */
  TriaAccessor &
  operator=(const TriaAccessor &) = delete;

  /**
   * Move assignment operator. Moving is allowed.
   */
  TriaAccessor &
  operator=(TriaAccessor &&) = default; // NOLINT

  /**
   * Defaulted destructor.
   */
  ~TriaAccessor() = default;

  /**
   * Test for the element being used or not.  The return value is @p true for
   * all iterators that are either normal iterators or active iterators, only
   * raw iterators can return @p false. Since raw iterators are only used in
   * the interiors of the library, you will not usually need this function.
   */
  bool
  used() const;

  /**
   * @name Accessing sub-objects
   */
  /**
   * @{
   */

  /**
   * Pointer to the @p ith vertex bounding this object. Throw an exception if
   * <code>dim=1</code>.
   */
  TriaIterator<TriaAccessor<0, dim, spacedim>>
  vertex_iterator(const unsigned int i) const;

  /**
   * Return the global index of i-th vertex of the current object. The
   * convention regarding the numbering of vertices is laid down in the
   * documentation of the GeometryInfo class.
   *
   * Note that the returned value is only the index of the geometrical vertex.
   * It has nothing to do with possible degrees of freedom associated with it.
   * For this, see the @p DoFAccessor::vertex_dof_index functions.
   *
   * @note Despite the name, the index returned here is only global in the
   * sense that it is specific to a particular Triangulation object or, in the
   * case the triangulation is actually of type
   * parallel::distributed::Triangulation, specific to that part of the
   * distributed triangulation stored on the current processor.
   */
  unsigned int
  vertex_index(const unsigned int i) const;

  /**
   * Return a reference to the @p ith vertex. The reference is not const,
   * i.e., it is possible to call this function on the left hand side of an
   * assignment, thereby moving the vertex of a cell within the triangulation.
   * Of course, doing so requires that you ensure that the new location of the
   * vertex remains useful -- for example, avoiding inverted or otherwise
   * distorted (see also
   * @ref GlossDistorted "this glossary entry").
   *
   * @note When a cell is refined, its children inherit the position of the
   * vertex positions of those vertices they share with the parent cell (plus
   * the locations of the new vertices on edges, faces, and cell interiors
   * that are created for the new child cells). If the vertex of a cell is
   * moved, this implies that its children will also use these new locations.
   * On the other hand, imagine a 2d situation where you have one cell that is
   * refined (with four children) and then you move the central vertex
   * connecting all four children. If you coarsen these four children again to
   * the parent cell, then the location of the moved vertex is lost and if, in
   * a later step, you refine the parent cell again, the then again new vertex
   * will be placed again at the same position as the first time around --
   * i.e., not at the location you had previously moved it to.
   *
   * @note The behavior described above is relevant if you have a
   * parallel::distributed::Triangulation object. There, refining a mesh
   * always involves a re-partitioning. In other words, vertices of locally
   * owned cells (see
   * @ref GlossLocallyOwnedCell "this glossary entry")
   * that you may have moved to a different location on one processor may be
   * moved to a different processor upon mesh refinement (even if these
   * particular cells were not refined) which will re-create their position
   * based on the position of the coarse cells they previously had, not based
   * on the position these vertices had on the processor that previously owned
   * them. In other words, in parallel computations, you will probably have to
   * move nodes explicitly after every mesh refinement because vertex
   * positions may or may not be preserved across the re-partitioning that
   * accompanies mesh refinement.
   */
  Point<spacedim> &
  vertex(const unsigned int i) const;

  /**
   * Pointer to the @p ith line bounding this object.
   */
  typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::line_iterator
    line(const unsigned int i) const;

  /**
   * Line index of the @p ith line bounding this object.
   *
   * Implemented only for <tt>structdim>1</tt>, otherwise an exception
   * generated.
   */
  unsigned int
  line_index(const unsigned int i) const;

  /**
   * Pointer to the @p ith quad bounding this object.
   */
  typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::quad_iterator
    quad(const unsigned int i) const;

  /**
   * Quad index of the @p ith quad bounding this object.
   *
   * Implemented only for <tt>structdim>2</tt>, otherwise an exception
   * generated.
   */
  unsigned int
  quad_index(const unsigned int i) const;
  /**
   * @}
   */

  /**
   * @name Orientation of sub-objects
   */
  /**
   * @{
   */

  /**
   * Return an integer representation that uniquely encodes the orientation,
   * flip, and rotation of a @p face.
   *
   * @ingroup reordering
   */
  types::geometric_orientation
  combined_face_orientation(const unsigned int face) const;

  /**
   * Return whether the face with index @p face has its normal pointing in the
   * standard direction (@p true) or whether it is the opposite (@p false).
   * Which is the standard direction is documented with the GeometryInfo
   * class. In 1d and 2d, this is always @p true, but in 3d it may be
   * different, see the respective discussion in the documentation of the
   * GeometryInfo class.
   *
   * This function is really only for internal use in the library unless you
   * absolutely know what this is all about.
   */
  bool
  face_orientation(const unsigned int face) const;

  /**
   * Return whether the face with index @p face is rotated by 180 degrees (@p
   * true) or not (@p false). In 1d and 2d, this is always @p false, but in
   * 3d it may be different, see the respective discussion in the
   * documentation of the GeometryInfo class.
   *
   * This function is really only for internal use in the library unless you
   * absolutely know what this is all about.
   */
  bool
  face_flip(const unsigned int face) const;

  /**
   * Return whether the face with index @p face is rotated by 90 degrees (@p
   * true) or not (@p false). In 1d and 2d, this is always @p false, but in
   * 3d it may be different, see the respective discussion in the
   * documentation of the GeometryInfo class.
   *
   * This function is really only for internal use in the library unless you
   * absolutely know what this is all about.
   */
  bool
  face_rotation(const unsigned int face) const;

  /**
   * Return whether the line with index @p line is oriented in standard
   * direction.
   *
   * @warning This function is really only for internal use in the library
   * unless you absolutely know what this is all about.
   *
   * @note This function queries ReferenceCell::face_to_cell_line_orientation().
   */
  types::geometric_orientation
  line_orientation(const unsigned int line) const;
  /**
   * @}
   */

  /**
   * @name Accessing children
   */
  /**
   * @{
   */

  /**
   * Test whether the object has children.
   */
  bool
  has_children() const;

  /**
   * Return the number of immediate children of this object. The number of
   * children of an unrefined cell is zero.
   */
  unsigned int
  n_children() const;

  /**
   * Compute and return the number of active descendants of this objects. For
   * example, if all of the eight children of a hex are further refined
   * isotropically exactly once, the returned number will be 64, not 80.
   *
   * If the present cell is not refined, one is returned.
   *
   * If one considers a triangulation as a forest where the root of each tree
   * are the coarse mesh cells and nodes have descendants (the children of a
   * cell), then this function returns the number of terminal nodes in the
   * sub-tree originating from the current object; consequently, if the
   * current object is not further refined, the answer is one.
   */
  unsigned int
  n_active_descendants() const;

  /**
   * Return the number of times that this object is refined. Note that not all
   * its children are refined that often (which is why we prepend @p max_),
   * the returned number is rather the maximum number of refinement in any
   * branch of children of this object.
   *
   * For example, if this object is refined, and one of its children is
   * refined exactly one more time, then <tt>max_refinement_depth</tt> should
   * return 2.
   *
   * If this object is not refined (i.e. it is active), then the return value
   * is zero.
   */
  unsigned int
  max_refinement_depth() const;

  /**
   * Return an iterator to the @p ith child.
   */
  TriaIterator<TriaAccessor<structdim, dim, spacedim>>
  child(const unsigned int i) const;

  /**
   * Return the child number of @p child on the current cell. This is the
   * inverse function of TriaAccessor::child().
   */
  unsigned int
  child_iterator_to_index(
    const TriaIterator<TriaAccessor<structdim, dim, spacedim>> &child) const;

  /**
   * Return an iterator to that object that is identical to the ith child for
   * isotropic refinement. If the current object is refined isotropically,
   * then the returned object is the ith child. If the current object is
   * refined anisotropically, the returned child may in fact be a grandchild
   * of the object, or may not exist at all (in which case an exception is
   * generated).
   */
  TriaIterator<TriaAccessor<structdim, dim, spacedim>>
  isotropic_child(const unsigned int i) const;

  /**
   * Return the RefinementCase of this cell.
   */
  RefinementCase<structdim>
  refinement_case() const;

  /**
   * Index of the @p ith child. The level of the child is one higher than that
   * of the present cell, if the children of a cell are accessed. The children
   * of faces have no level. If the child does not exist, -1 is returned.
   */
  int
  child_index(const unsigned int i) const;

  /**
   * Index of the @p ith isotropic child. See the isotropic_child() function
   * for a definition of this concept.  If the child does not exist, -1 is
   * returned.
   */
  int
  isotropic_child_index(const unsigned int i) const;
  /**
   * @}
   */

  /**
   * @name Dealing with boundary indicators
   */
  /**
   * @{
   */

  /**
   * Return the boundary indicator of this object.
   *
   * If the return value is the special value
   * numbers::internal_face_boundary_id, then this object is in the interior
   * of the domain.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  types::boundary_id
  boundary_id() const;

  /**
   * Set the boundary indicator of the current object. The same applies as for
   * the boundary_id() function.
   *
   * This function only sets the boundary object of the current object itself,
   * not the indicators of the ones that bound it. For example, in 3d, if this
   * function is called on a face, then the boundary indicator of the 4 edges
   * that bound the face remain unchanged. If you want to set the boundary
   * indicators of face and edges at the same time, use the
   * set_all_boundary_ids() function. You can see the result of not using the
   * correct function in the results section of step-49.
   *
   * @warning You should never set the boundary indicator of an interior face
   * (a face not at the boundary of the domain), or set the boundary
   * indicator of an exterior face to numbers::internal_face_boundary_id (this
   * value is reserved for another purpose). Algorithms may not work or
   * produce very confusing results if boundary cells have a boundary
   * indicator of numbers::internal_face_boundary_id or if interior cells have
   * boundary indicators other than numbers::internal_face_boundary_id.
   * Unfortunately, the current object has no means of finding out whether it
   * really is at the boundary of the domain and so cannot determine whether
   * the value you are trying to set makes sense under the current
   * circumstances.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  void
  set_boundary_id(const types::boundary_id) const;

  /**
   * Do as set_boundary_id() but also set the boundary indicators of the
   * objects that bound the current object. For example, in 3d, if
   * set_boundary_id() is called on a face, then the boundary indicator of the
   * 4 edges that bound the face remain unchanged. In contrast, if you call
   * the current function, the boundary indicators of face and edges are all
   * set to the given value.
   *
   * This function is useful if you set boundary indicators of faces in 3d (in
   * 2d, the function does the same as set_boundary_id()) and you do so
   * because you want a curved boundary object to represent the part of the
   * boundary that corresponds to the current face. In that case, the
   * Triangulation class needs to figure out where to put new vertices upon
   * mesh refinement, and higher order Mapping objects also need to figure out
   * where new interpolation points for a curved boundary approximation should
   * be. In either case, the two classes first determine where interpolation
   * points on the edges of a boundary face should be, asking the boundary
   * object, before asking the boundary object for the interpolation points
   * corresponding to the interior of the boundary face. For this to work
   * properly, it is not sufficient to have set the boundary indicator for the
   * face alone, but you also need to set the boundary indicators of the edges
   * that bound the face. This function does all of this at once. You can see
   * the result of not using the correct function in the results section of
   * step-49.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  void
  set_all_boundary_ids(const types::boundary_id) const;

  /**
   * Return whether this object is at the boundary. Obviously, the use of this
   * function is only possible for <tt>dim@>structdim</tt>; however, for
   * <tt>dim==structdim</tt>, an object is a cell and the CellAccessor class
   * offers another possibility to determine whether a cell is at the boundary
   * or not.
   */
  bool
  at_boundary() const;

  /**
   * Return a constant reference to the manifold object used for this object.
   *
   * As explained in the
   * @ref manifold
   * topic, the process involved in finding the appropriate manifold
   * description involves querying both the manifold or boundary
   * indicators. See there for more information.
   */
  const Manifold<dim, spacedim> &
  get_manifold() const;

  /**
   * @}
   */

  /**
   * @name Dealing with manifold indicators
   */
  /**
   * @{
   */

  /**
   * Return the manifold indicator of this object.
   *
   * If the return value is the special value numbers::flat_manifold_id, then
   * this object is associated with a standard Cartesian Manifold Description.
   *
   * @see
   * @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  types::manifold_id
  manifold_id() const;

  /**
   * Set the manifold indicator.  The same applies as for the
   * <tt>manifold_id()</tt> function.
   *
   * Note that it only sets the manifold object of the current object itself,
   * not the indicators of the ones that bound it, nor of its children. For
   * example, in 3d, if this function is called on a face, then the manifold
   * indicator of the 4 edges that bound the face remain unchanged. If you
   * want to set the manifold indicators of face, edges and all children at
   * the same time, use the set_all_manifold_ids() function.
   *
   *
   * @ingroup manifold
   *
   * @see
   * @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  void
  set_manifold_id(const types::manifold_id) const;

  /**
   * Do as set_manifold_id() but also set the manifold indicators of the
   * objects that bound the current object. For example, in 3d, if
   * set_manifold_id() is called on a face, then the manifold indicator of the
   * 4 edges that bound the face remain unchanged. On the other hand, the
   * manifold indicators of face and edges are all set at the same time using
   * the current function.
   *
   * @ingroup manifold
   *
   * @see
   * @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  void
  set_all_manifold_ids(const types::manifold_id) const;

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
   * Read the user flag. See
   * @ref GlossUserFlags
   * for more information.
   */
  bool
  user_flag_set() const;

  /**
   * Set the user flag. See
   * @ref GlossUserFlags
   * for more information.
   */
  void
  set_user_flag() const;

  /**
   * Clear the user flag. See
   * @ref GlossUserFlags
   * for more information.
   */
  void
  clear_user_flag() const;

  /**
   * Set the user flag for this and all descendants. See
   * @ref GlossUserFlags
   * for more information.
   */
  void
  recursively_set_user_flag() const;

  /**
   * Clear the user flag for this and all descendants. See
   * @ref GlossUserFlags
   * for more information.
   */
  void
  recursively_clear_user_flag() const;

  /**
   * Reset the user data to zero, independent if pointer or index. See
   * @ref GlossUserData
   * for more information.
   */
  void
  clear_user_data() const;

  /**
   * Set the user pointer to @p p.
   *
   * @note User pointers and user indices are mutually exclusive. Therefore,
   * you can only use one of them, unless you call
   * Triangulation::clear_user_data() in between.
   *
   * See
   * @ref GlossUserData
   * for more information.
   */
  void
  set_user_pointer(void *p) const;

  /**
   * Reset the user pointer to `nullptr`. See
   * @ref GlossUserData
   * for more information.
   */
  void
  clear_user_pointer() const;

  /**
   * Access the value of the user pointer. It is in the responsibility of the
   * user to make sure that the pointer points to something useful and always
   * requires casting to a known type, e.g.,
   *
   * @code
   * auto *a = static_cast<A*>(cell->user_pointer());
   * @endcode
   *
   * @note User pointers and user indices are mutually exclusive. Therefore, you
   * can only use one of them, unless you call Triangulation::clear_user_data()
   * in between.
   *
   * See
   * @ref GlossUserData
   * for more information.
   */
  void *
  user_pointer() const;

  /**
   * Set the user pointer of this object and all its children to the given
   * value. This is useful for example if all cells of a certain subdomain, or
   * all faces of a certain part of the boundary should have user pointers
   * pointing to objects describing this part of the domain or boundary.
   *
   * Note that the user pointer is not inherited under mesh refinement, so
   * after mesh refinement there might be cells or faces that don't have user
   * pointers pointing to the describing object. In this case, simply loop
   * over all the elements of the coarsest level that has this information,
   * and use this function to recursively set the user pointer of all finer
   * levels of the triangulation.
   *
   * @note User pointers and user indices are mutually exclusive. Therefore,
   * you can only use one of them, unless you call
   * Triangulation::clear_user_data() in between.
   *
   * See
   * @ref GlossUserData
   * for more information.
   */
  void
  recursively_set_user_pointer(void *p) const;

  /**
   * Clear the user pointer of this object and all of its descendants. The
   * same holds as said for the recursively_set_user_pointer() function. See
   * @ref GlossUserData
   * for more information.
   */
  void
  recursively_clear_user_pointer() const;

  /**
   * Set the user index to @p p.
   *
   * @note User pointers and user indices are mutually exclusive. Therefore,
   * you can only use one of them, unless you call
   * Triangulation::clear_user_data() in between. See
   * @ref GlossUserData
   * for more information.
   */
  void
  set_user_index(const unsigned int p) const;

  /**
   * Reset the user index to 0. See
   * @ref GlossUserData
   * for more information.
   */
  void
  clear_user_index() const;

  /**
   * Access the value of the user index.
   *
   * @note User pointers and user indices are mutually exclusive. Therefore,
   * you can only use one of them, unless you call
   * Triangulation::clear_user_data() in between.
   *
   * See
   * @ref GlossUserData
   * for more information.
   */
  unsigned int
  user_index() const;

  /**
   * Set the user index of this object and all its children.
   *
   * Note that the user index is not inherited under mesh refinement, so after
   * mesh refinement there might be cells or faces that don't have the
   * expected user indices. In this case, simply loop over all the elements of
   * the coarsest level that has this information, and use this function to
   * recursively set the user index of all finer levels of the triangulation.
   *
   * @note User pointers and user indices are mutually exclusive. Therefore,
   * you can only use one of them, unless you call
   * Triangulation::clear_user_data() in between.
   *
   * See
   * @ref GlossUserData
   * for more information.
   */
  void
  recursively_set_user_index(const unsigned int p) const;

  /**
   * Clear the user index of this object and all of its descendants. The same
   * holds as said for the recursively_set_user_index() function.
   *
   * See
   * @ref GlossUserData
   * for more information.
   */
  void
  recursively_clear_user_index() const;
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
   * The diameter of an object is computed to be the largest diagonal of the
   * current object. If this object is a quadrilateral, then there are two
   * such diagonal, and if it is a hexahedron, then there are four diagonals
   * that connect "opposite" points. For triangles and tetrahedra, the function
   * simply returns the length of the longest edge.
   *
   * The situation is more difficult for wedges and pyramids: For wedges, we
   * return the length of the longest diagonal of the three quadrilateral faces
   * or the longest edge length of the two triangular faces. For pyramids,
   * the same principle is applied.
   *
   * In all of these cases, this definition of "diameter" is
   * not necessarily the true diameter in the sense of the largest distance
   * between points inside the object. Indeed, one can often construct objects
   * for which it is not, though these are generally quite deformed compared to
   * the reference shape. Furthermore, for objects that may use higher order
   * mappings, one may have bulging faces that also create trouble for
   * computing an exact representation of the diameter of the object. That said,
   * the definition used above is completely sufficient for most computations.
   */
  double
  diameter() const;

  /**
   * Return a pair of Point and double corresponding to the center and
   * the radius of a reasonably small enclosing ball of the object.
   *
   * The function implements Ritter's O(n) algorithm to get a reasonably
   * small enclosing ball around the vertices of the object.
   * The initial guess for the enclosing ball is taken to be the ball
   * which contains the largest diagonal of the object as its diameter.
   * Starting from such an initial guess, the algorithm tests whether all
   * the vertices of the object (except the vertices of the largest diagonal)
   * are geometrically within the ball.
   * If any vertex (v) is found to be geometrically outside the ball,
   * a new iterate of the ball is constructed by shifting its center and
   * increasing the radius so as to geometrically enclose both the previous
   * ball and the vertex (v). The algorithm terminates when all the vertices
   * are geometrically inside the ball.
   *
   * If a vertex (v) is geometrically inside a particular iterate of the ball,
   * then it will continue to be so in the subsequent iterates of
   * the ball (this is true \a by \a construction).
   *
   * @note This function assumes d-linear mapping from the reference cell.
   *
   * <a href="http://geomalgorithms.com/a08-_containers.html">see this</a> and
   * [Ritter 1990]
   */
  std::pair<Point<spacedim>, double>
  enclosing_ball() const;

  /**
   * Return the smallest bounding box that encloses the object.
   *
   * Notice that this method is not aware of any mapping you may be using to
   * do your computations. If you are using a mapping object that modifies the
   * position of the vertices, like MappingQEulerian, or MappingFEField, then
   * you should call the function Mapping::get_bounding_box() instead.
   */
  BoundingBox<spacedim>
  bounding_box() const;

  /**
   * Length of an object in the direction of the given axis, specified in the
   * local coordinate system. See the documentation of GeometryInfo for the
   * meaning and enumeration of the local axes.
   *
   * Note that the "length" of an object can be interpreted in a variety of
   * ways. Here, we choose it as the maximal length of any of the edges of the
   * object that are parallel to the chosen axis on the reference cell.
   */
  double
  extent_in_direction(const unsigned int axis) const;

  /**
   * Return the minimal distance between any two vertices.
   */
  double
  minimum_vertex_distance() const;

  /**
   * Return a point belonging to the Manifold<dim,spacedim> where this object
   * lives, given its parametric coordinates on the reference @p structdim
   * cell. This function queries the underlying manifold object, and can be
   * used to obtain the exact geometrical location of arbitrary points on this
   * object.
   *
   * Notice that the argument @p coordinates are the coordinates on the
   * <em>reference cell</em>, given in reference coordinates. In other words,
   * the argument provides a weighting between the different vertices. For
   * example, for lines, calling this function with argument Point<1>(.5), is
   * equivalent to asking the line for its center.
   */
  Point<spacedim>
  intermediate_point(const Point<structdim> &coordinates) const;

  /**
   * This function computes a fast approximate transformation from the real to
   * the unit cell by inversion of an affine approximation of the $d$-linear
   * function from the reference $d$-dimensional cell.
   *
   * The affine approximation of the unit to real cell mapping is found by a
   * least squares fit of an affine function to the $2^d$ vertices of the
   * present object. For any valid mesh cell whose geometry is not degenerate,
   * this operation results in a unique affine mapping. Thus, this function
   * will return a finite result for all given input points, even in cases
   * where the actual transformation by an actual bi-/trilinear or higher
   * order mapping might be singular. Besides only approximating the mapping
   * from the vertex points, this function also ignores the attached manifold
   * descriptions. The result is only exact in case the transformation from
   * the unit to the real cell is indeed affine, such as in one dimension or
   * for Cartesian and affine (parallelogram) meshes in 2d/3d.
   *
   * For exact transformations to the unit cell, use
   * Mapping::transform_real_to_unit_cell().
   *
   * @note If dim<spacedim we first project p onto the plane.
   */
  Point<structdim>
  real_to_unit_cell_affine_approximation(const Point<spacedim> &point) const;

  /**
   * Center of the object. The center of an object is defined to be the
   * average of the locations of the vertices, which is also where a $Q_1$
   * mapping would map the center of the reference cell. However, you can also
   * ask this function to instead return the average of the vertices as
   * computed by the underlying Manifold object associated with the current
   * object, by setting to true the optional parameter @p respect_manifold.
   * Manifolds would then typically pull back the coordinates of the vertices
   * to a reference domain (not necessarily the reference cell), compute the
   * average there, and then push forward the coordinates of the averaged
   * point to the physical space again; the resulting point is guaranteed to
   * lie within the manifold, even if the manifold is curved.
   *
   * When the object uses a different manifold description as its surrounding,
   * like when part of the bounding objects of this TriaAccessor use a
   * non-flat manifold description but the object itself is flat, the result
   * given by the TriaAccessor::center() function may not be accurate enough,
   * even when parameter @p respect_manifold is set to true. If you find this
   * to be case, than you can further refine the computation of the center by
   * setting to true the second additional parameter @p
   * interpolate_from_surrounding. This computes the location of the center by
   * a so-called transfinite interpolation from the center of all the bounding
   * objects. For a 2d object, it puts a weight of <code>1/2</code> on each of
   * the four surrounding lines and a weight <code>-1/4</code> on the four
   * vertices. This corresponds to a linear interpolation between the
   * descriptions of the four faces, subtracting the contribution of the
   * vertices that is added twice when coming through both lines adjacent to
   * the vertex. In 3d, the weights for faces are <code>1/2</code>, the
   * weights for lines are <code>-1/4</code>, and the weights for vertices are
   * <code>1/8</code>. For further information, also confer to the
   * TransfiniteInterpolationManifold class that is able to not only apply
   * this beneficial description to a single cell but all children of a coarse
   * cell.
   */
  Point<spacedim>
  center(const bool respect_manifold             = false,
         const bool interpolate_from_surrounding = false) const;

  /**
   * Return the barycenter (also called centroid)
   * of the object. The barycenter for an object $K$
   * of dimension $d$ in $D$ space dimensions is given by the $D$-dimensional
   * vector $\mathbf x_K$ defined by
   * @f[
   *   \mathbf x_K = \frac{1}{|K|} \int_K \mathbf x \; \textrm{d}x
   * @f]
   * where the measure of the object is given by
   * @f[
   *   |K| = \int_K \mathbf 1 \; \textrm{d}x.
   * @f]
   * This function assumes that $K$ is mapped by a $d$-linear function from
   * the reference $d$-dimensional cell. Then the integrals above can be
   * pulled back to the reference cell and evaluated exactly (if through
   * lengthy and, compared to the center() function, expensive computations).
   */
  Point<spacedim>
  barycenter() const;

  /**
   * Compute the `structdim`-dimensional measure of the object. For a
   * `dim`-dimensional cell in `dim`-dimensional space, this equals its volume.
   * On the other hand, for a 2d cell in 3d space, or if the current object
   * pointed to is a 2d face of a 3d cell in 3d space, then the function
   * computes the area the object occupies. For a one-dimensional (i.e.,
   * `structdim = 1`) object, regardless of `dim` and `spacedim`, return its
   * length. Similarly, the measure of any vertex (i.e., `structdim = 0`
   * objects) is 1.
   *
   * The function only computes the measure of cells, faces or edges assumed
   * to be represented by (bi-/tri-)linear mappings. In other words, it only
   * takes into account the locations of the vertices that bound the current
   * object but not how the interior of the object may actually be mapped. In
   * most simple cases, this is exactly what you want. However, for objects
   * that are not "straight", e.g. 2d cells embedded in 3d space as part of a
   * triangulation of a curved domain, two-dimensional faces of 3d cells that
   * are not just parallelograms, or for faces that are at the boundary of a
   * domain that is not just bounded by straight line segments or planes, this
   * function only computes the dim-dimensional measure of a (bi-/tri-)linear
   * interpolation of the "real" object as defined by the manifold or boundary
   * object describing the real geometry of the object in question. If you
   * want to consider the "real" geometry, you will need to compute the
   * measure by integrating a function equal to one over the object, which
   * after applying quadrature equals the summing the JxW values returned by
   * the FEValues or FEFaceValues object you will want to use for the
   * integral.
   *
   * @note There is no analytic formula for the area of a bilinear face (i.e.,
   * something with a quadrilateral reference cell) in 3D. This function uses
   * the quadrature defined by QGauss<2>(4) to approximate the measure in this
   * case.
   */
  double
  measure() const;

  /**
   * Return true if the current object is a translation of the given argument.
   *
   * @note For the purpose of a triangulation, cells, faces, etc are only
   * characterized by their vertices. The current function therefore only
   * compares the locations of vertices. For many practical applications,
   * however, it is not only the vertices that determine whether one cell is a
   * translation of another, but also how the cell is mapped from the
   * reference cell to its location in real space. For example, if we are
   * using higher order mappings, then not only do the vertices have to be
   * translations of each other, but also the points along edges. In these
   * questions, therefore, it would be appropriate to ask the mapping, not the
   * current function, whether two objects are translations of each other.
   */
  bool
  is_translation_of(
    const TriaIterator<TriaAccessor<structdim, dim, spacedim>> &o) const;

  /**
   * Reference cell type of the current object.
   */
  ReferenceCell
  reference_cell() const;

  /**
   * Number of vertices.
   */
  unsigned int
  n_vertices() const;

  /**
   * Number of lines.
   */
  unsigned int
  n_lines() const;

  /**
   * Return the number of faces for a cell. This function is only
   * implemented for cells (i.e., `structdim==dim`) to avoid the question
   * of what exactly is meant in a construct such as
   * `cell->face(f)->n_faces()`. If you want to ask how many bounding
   * lines a face of a 3d cell has, use `cell->face(f)->n_lines()`; if
   * you want to ask about the number of vertices of a face of a 2d cell,
   * use `cell->face(f)->n_vertices()`.
   */
  unsigned int
  n_faces() const;

  /**
   * Return an object that can be thought of as an array containing all indices
   * from zero to n_vertices().
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  vertex_indices() const;

  /**
   * Return an object that can be thought of as an array containing all indices
   * from zero to n_lines().
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  line_indices() const;

  /**
   * Return an object that can be thought of as an array containing all indices
   * from zero to n_faces().
   *
   * @note Only implemented for cells (structdim==dim).
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  face_indices() const;

  /**
   * @}
   */

private:
  /**
   * Like set_boundary_id but without checking for internal faces or invalid
   * ids.
   */
  void
  set_boundary_id_internal(const types::boundary_id id) const;

  /**
   * Set the indices of those objects that bound the current
   * object. For example, if the current object represents a cell,
   * then the argument denotes the indices of the faces that bound the
   * cell. If the current object represents a line, the argument
   * denotes the indices of the vertices that bound it. And so on.
   */
  void
  set_bounding_object_indices(
    const std::initializer_list<int> &new_indices) const;

  /**
   * The same as above but for `unsigned int`.
   */
  void
  set_bounding_object_indices(
    const std::initializer_list<unsigned int> &new_indices) const;

  /**
   * Set the flag indicating, what <code>line_orientation()</code> will
   * return.
   *
   * It is only possible to set the line_orientation of faces in 3d (i.e.
   * <code>structdim==2 && dim==3</code>).
   */
  void
  set_line_orientation(const unsigned int                 line,
                       const types::geometric_orientation orientation) const;

  /**
   * Set the combined face orientation (i.e., the integer that uniquely encodes
   * the orientation, flip, and rotation). This function is only implemented for
   * objects which have faces, i.e., for structdim == dim.
   * For more information see the
   * @ref GlossCombinedOrientation "combined orientation glossary entry".
   *
   * @ingroup reordering
   */
  void
  set_combined_face_orientation(
    const unsigned int                 face,
    const types::geometric_orientation combined_orientation) const;

  /**
   * Set the @p used flag. Only for internal use in the library.
   */
  void
  set_used_flag() const;

  /**
   * Clear the @p used flag. Only for internal use in the library.
   */
  void
  clear_used_flag() const;

  /**
   * Set the @p RefinementCase<dim> this TriaObject is refined with. Not
   * defined for <tt>structdim=1</tt> as lines are always refined resulting in
   * 2 children lines (isotropic refinement).
   *
   * You should know quite exactly what you are doing if you touch this
   * function. It is exclusively for internal use in the library.
   */
  void
  set_refinement_case(const RefinementCase<structdim> &ref_case) const;

  /**
   * Clear the RefinementCase<dim> of this TriaObject, i.e. reset it to
   * RefinementCase<dim>::no_refinement.
   *
   * You should know quite exactly what you are doing if you touch this
   * function. It is exclusively for internal use in the library.
   */
  void
  clear_refinement_case() const;

  /**
   * Set the index of the ith child. Since the children come at least in
   * pairs, we need to store the index of only every second child, i.e. of the
   * even numbered children. Make sure, that the index of child i=0 is set
   * first. Calling this function for odd numbered children is not allowed.
   */
  void
  set_children(const unsigned int i, const int index) const;

  /**
   * Clear the child field, i.e. set it to a value which indicates that this
   * cell has no children.
   */
  void
  clear_children() const;

private:
  friend class Triangulation<dim, spacedim>;

  friend struct dealii::internal::TriangulationImplementation::Implementation;
  friend struct dealii::internal::TriangulationImplementation::
    ImplementationMixedMesh;
  friend struct dealii::internal::TriaAccessorImplementation::Implementation;
};



/**
 * This class is a specialization of <code>TriaAccessor<structdim, dim,
 * spacedim></code>
 * for the case that @p structdim is zero. This
 * class represents vertices in a triangulation of dimensionality
 * <code>dim</code> (i.e. 1 for a triangulation of lines, 2 for a
 * triangulation of quads, and 3 for a triangulation of hexes) that is
 * embedded in a space of dimensionality <code>spacedim</code> (for
 * <code>spacedim==dim</code> the triangulation represents a domain in
 * ${\mathbb R}^\text{dim}$, for <code>spacedim@>dim</code> the triangulation
 * is of a manifold embedded in a higher dimensional space).
 *
 * There is a further specialization of this class for the case that
 * @p dim equals one, i.e., for vertices of a one-dimensional triangulation,
 * since in that case vertices are also faces.
 *
 * @ingroup Accessors
 */
template <int dim, int spacedim>
class TriaAccessor<0, dim, spacedim>
{
public:
  /**
   * Dimension of the space the object represented by this accessor lives in.
   * For example, if this accessor represents a quad that is part of a
   * two-dimensional surface in four-dimensional space, then this value is four.
   */
  static constexpr unsigned int space_dimension = spacedim;

  /**
   * Dimensionality of the object that the thing represented by this accessor
   * is part of. For example, if this accessor represents a line that is part
   * of a hexahedron, then this value will be three.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Dimensionality of the current object represented by this accessor. For
   * example, if it is line (irrespective of whether it is part of a quad or
   * hex, and what dimension we are in), then this value equals 1.
   */
  static const unsigned int structure_dimension = 0;

  /**
   * Pointer to internal data.
   */
  using AccessorData = void;

  /**
   * Constructor. The second argument is the global index of the vertex we
   * point to.
   */
  TriaAccessor(const Triangulation<dim, spacedim> *tria,
               const unsigned int                  vertex_index);

  /**
   * Constructor. This constructor exists in order to maintain interface
   * compatibility with the other accessor classes. @p index can be used to
   * set the global index of the vertex we point to.
   */
  TriaAccessor(const Triangulation<dim, spacedim> *tria  = nullptr,
               const int                           level = 0,
               const int                           index = 0,
               const AccessorData                      * = nullptr);

  /**
   * Constructor. Should never be called and thus produces an error.
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor(const TriaAccessor<structdim2, dim2, spacedim2> &);

  /**
   * Constructor. Should never be called and thus produces an error.
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * Return the state of the iterator.
   */
  IteratorState::IteratorStates
  state() const;

  /**
   * Level of this object. Vertices have no level, so this function always
   * returns zero.
   */
  static int
  level();

  /**
   * Index of this object. Returns the global index of the vertex this object
   * points to.
   */
  int
  index() const;

  /**
   * Return a reference to the triangulation which the object pointed to by this
   * class belongs to.
   */
  const Triangulation<dim, spacedim> &
  get_triangulation() const;

  /**
   * @name Advancement of iterators
   */
  /**
   * @{
   */
  /**
   * This operator advances the iterator to the next element.
   */
  void
  operator++();

  /**
   * This operator moves the iterator to the previous element.
   */
  void
  operator--();
  /**
   * Compare for equality.
   */
  bool
  operator==(const TriaAccessor &) const;

  /**
   * Compare for inequality.
   */
  bool
  operator!=(const TriaAccessor &) const;

  /**
   * @}
   */


  /**
   * @name Accessing sub-objects
   */
  /**
   * @{
   */

  /**
   * Return the global index of i-th vertex of the current object. If @p i is
   * zero, this returns the index of the current point to which this object
   * refers. Otherwise, it throws an exception.
   *
   * Note that the returned value is only the index of the geometrical vertex.
   * It has nothing to do with possible degrees of freedom associated with it.
   * For this, see the @p DoFAccessor::vertex_dof_index functions.
   *
   * @note Despite the name, the index returned here is only global in the
   * sense that it is specific to a particular Triangulation object or, in the
   * case the triangulation is actually of type
   * parallel::distributed::Triangulation, specific to that part of the
   * distributed triangulation stored on the current processor.
   */
  unsigned int
  vertex_index(const unsigned int i = 0) const;

  /**
   * Return a reference to the @p ith vertex. If i is zero, this returns a
   * reference to the current point to which this object refers. Otherwise, it
   * throws an exception.
   */
  Point<spacedim> &
  vertex(const unsigned int i = 0) const;

  /**
   * Pointer to the @p ith line bounding this object. Will point to an invalid
   * object.
   */
  typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::line_iterator static line(const unsigned int);

  /**
   * Line index of the @p ith line bounding this object. Throws an exception.
   */
  static unsigned int
  line_index(const unsigned int i);

  /**
   * Pointer to the @p ith quad bounding this object.
   */
  static typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::quad_iterator
    quad(const unsigned int i);

  /**
   * Quad index of the @p ith quad bounding this object. Throws an exception.
   */
  static unsigned int
  quad_index(const unsigned int i);

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
   * Diameter of the object. This function always returns zero.
   */
  double
  diameter() const;

  /**
   * Length of an object in the direction of the given axis, specified in the
   * local coordinate system. See the documentation of GeometryInfo for the
   * meaning and enumeration of the local axes.
   *
   * This function always returns zero.
   */
  double
  extent_in_direction(const unsigned int axis) const;

  /**
   * Return the center of this object, which of course coincides with the
   * location of the vertex this object refers to. The parameters @p
   * respect_manifold and @p interpolate_from_surrounding are not used. They
   * are there to provide the same interface as
   * <code>TriaAccessor<structdim,dim,spacedim></code>.
   */
  Point<spacedim>
  center(const bool respect_manifold             = false,
         const bool interpolate_from_surrounding = false) const;

  /**
   * Compute the `structdim`-dimensional measure of the present object. Since,
   * in this context, `structdim = 0`, this function always returns 1.
   *
   * @note This is consistent with what ReferenceCells::Vertex::volume()
   * returns.
   */
  double
  measure() const;
  /**
   * @}
   */

  /**
   * @name Orientation of sub-objects
   */
  /**
   * @{
   */

  /**
   * @brief Always return 0
   *
   * @ingroup reordering
   */
  static types::geometric_orientation
  combined_face_orientation(const unsigned int face);

  /**
   * @brief Always return false
   */
  static bool
  face_orientation(const unsigned int face);

  /**
   * @brief Always return false
   */
  static bool
  face_flip(const unsigned int face);

  /**
   * @brief Always return false
   */
  static bool
  face_rotation(const unsigned int face);

  /**
   * @brief Always return numbers::reverse_line_orientation
   */
  static types::geometric_orientation
  line_orientation(const unsigned int line);

  /**
   * @}
   */

  /**
   * @name Accessing children
   */
  /**
   * @{
   */

  /**
   * Test whether the object has children. Always false.
   */
  static bool
  has_children();

  /**
   * Return the number of immediate children of this object. This is always
   * zero.
   */
  static unsigned int
  n_children();

  /**
   * Compute and return the number of active descendants of this objects.
   * Always zero.
   */
  static unsigned int
  n_active_descendants();

  /**
   * Return the number of times that this object is refined. Always 0.
   */
  static unsigned int
  max_refinement_depth();

  /**
   * @brief Return an invalid unsigned integer.
   */
  static unsigned int
  child_iterator_to_index(const TriaIterator<TriaAccessor<0, dim, spacedim>> &);

  /**
   * @brief Return an invalid object.
   */
  static TriaIterator<TriaAccessor<0, dim, spacedim>>
  child(const unsigned int);

  /**
   * @brief Return an invalid object.
   */
  static TriaIterator<TriaAccessor<0, dim, spacedim>>
  isotropic_child(const unsigned int);

  /**
   * Always return no refinement.
   */
  static RefinementCase<0>
  refinement_case();

  /**
   * @brief Returns -1
   */
  static int
  child_index(const unsigned int i);

  /**
   * @brief Returns -1
   */
  static int
  isotropic_child_index(const unsigned int i);
  /**
   * @}
   */

  /**
   * Return whether the vertex pointed to here is used.
   */
  bool
  used() const;

protected:
  /**
   * Copy operator. Since this is only called from iterators, do not return
   * anything, since the iterator will return itself.
   *
   * This method is protected, since it is only to be called from the iterator
   * class.
   */
  void
  copy_from(const TriaAccessor &);

  /**
   * Comparison operator for accessors. This operator is used when comparing
   * iterators into objects of a triangulation, for example when putting
   * them into a `std::map`.
   *
   * This operator simply compares the global index of the vertex the
   * current object points to.
   */
  bool
  operator<(const TriaAccessor &other) const;

  /**
   * Pointer to the triangulation we operate on.
   */
  const Triangulation<dim, spacedim> *tria;

  /**
   * The global vertex index of the vertex this object corresponds to.
   */
  unsigned int global_vertex_index;

private:
  template <typename Accessor>
  friend class TriaRawIterator;
  template <typename Accessor>
  friend class TriaIterator;
  template <typename Accessor>
  friend class TriaActiveIterator;
};



/**
 * This class is a specialization of <code>TriaAccessor<structdim, dim,
 * spacedim></code>
 * for the case that @p structdim is zero and @p dim is one. This
 * class represents vertices in a one-dimensional triangulation that is
 * embedded in a space of dimensionality <code>spacedim</code> (for
 * <code>spacedim==dim==1</code> the triangulation represents a domain in
 * ${\mathbb R}^\text{dim}$, for <code>spacedim@>dim==1</code> the triangulation
 * is of a manifold embedded in a higher dimensional space).
 *
 * The current specialization of the TriaAccessor<0,dim,spacedim> class
 * for vertices of a one-dimensional triangulation exists
 * since in the @p dim == 1 case vertices are also faces.
 *
 * @ingroup Accessors
 */
template <int spacedim>
class TriaAccessor<0, 1, spacedim>
{
public:
  /**
   * Dimension of the space the object represented by this accessor lives in.
   * For example, if this accessor represents a quad that is part of a
   * two-dimensional surface in four-dimensional space, then this value is four.
   */
  static constexpr unsigned int space_dimension = spacedim;

  /**
   * Dimensionality of the object that the thing represented by this accessor
   * is part of. For example, if this accessor represents a line that is part
   * of a hexahedron, then this value will be three.
   */
  static constexpr unsigned int dimension = 1;

  /**
   * Dimensionality of the current object represented by this accessor. For
   * example, if it is line (irrespective of whether it is part of a 2d or 3d
   * subobject), then this value equals 1.
   */
  static const unsigned int structure_dimension = 0;

  /**
   * Pointer to internal data.
   */
  using AccessorData = void;

  /**
   * Whether the vertex represented here is at the left end of the domain, the
   * right end, or in the interior.
   */
  enum VertexKind
  {
    /**
     * Left vertex.
     */
    left_vertex,
    /**
     * Interior vertex.
     */
    interior_vertex,
    /**
     * Right vertex.
     */
    right_vertex
  };

  /**
   * Constructor.
   *
   * Since there is no mapping from vertices to cells, an accessor object for
   * a point has no way to figure out whether it is at the boundary of the
   * domain or not. Consequently, the second argument must be passed by the
   * object that generates this accessor -- e.g. a 1d cell that can figure out
   * whether its left or right vertex are at the boundary.
   *
   * The third argument is the global index of the vertex we point to.
   */
  TriaAccessor(const Triangulation<1, spacedim> *tria,
               const VertexKind                  vertex_kind,
               const unsigned int                vertex_index);

  /**
   * Constructor. This constructor exists in order to maintain interface
   * compatibility with the other accessor classes. However, it doesn't do
   * anything useful here and so may not actually be called.
   */
  TriaAccessor(const Triangulation<1, spacedim> *tria = nullptr,
               const int                              = 0,
               const int                              = 0,
               const AccessorData                   * = nullptr);

  /**
   * Constructor. Should never be called and thus produces an error.
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor(const TriaAccessor<structdim2, dim2, spacedim2> &);

  /**
   * Constructor. Should never be called and thus produces an error.
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * Copy operator. Since this is only called from iterators, do not return
   * anything, since the iterator will return itself.
   */
  void
  copy_from(const TriaAccessor &);

  /**
   * Copy operator. We need this function to support generic
   * programming, but it just throws an exception because it cannot do
   * the required operations.
   */
  void
  copy_from(const TriaAccessorBase<0, 1, spacedim> &);

  /**
   * Return the state of the iterator. Since an iterator to points can not be
   * incremented or decremented, its state remains constant, and in particular
   * equal to IteratorState::valid.
   */
  static IteratorState::IteratorStates
  state();

  /**
   * Level of this object. Vertices have no level, so this function always
   * returns zero.
   */
  static int
  level();

  /**
   * Index of this object. Returns the global index of the vertex this object
   * points to.
   */
  int
  index() const;

  /**
   * Return a reference to the triangulation which the object pointed to by this
   * class belongs to.
   */
  const Triangulation<1, spacedim> &
  get_triangulation() const;

  /**
   * @name Advancement of iterators
   */
  /**
   * @{
   */
  /**
   * This operator advances the iterator to the next element. For points, this
   * operation is not defined, so you can't iterate over point iterators.
   */
  void
  operator++() const;

  /**
   * This operator moves the iterator to the previous element. For points,
   * this operation is not defined, so you can't iterate over point iterators.
   */
  void
  operator--() const;
  /**
   * Compare for equality.
   */
  bool
  operator==(const TriaAccessor &) const;

  /**
   * Compare for inequality.
   */
  bool
  operator!=(const TriaAccessor &) const;

  /**
   * Comparison operator for accessors. This operator is used when comparing
   * iterators into objects of a triangulation, for example when putting
   * them into a `std::map`.
   *
   * This operator simply compares the global index of the vertex the
   * current object points to.
   */
  bool
  operator<(const TriaAccessor &other) const;

  /**
   * @}
   */

  /**
   * @name Accessing sub-objects
   */
  /**
   * @{
   */

  /**
   * Return the global index of i-th vertex of the current object. If i is
   * zero, this returns the index of the current point to which this object
   * refers. Otherwise, it throws an exception.
   *
   * Note that the returned value is only the index of the geometrical vertex.
   * It has nothing to do with possible degrees of freedom associated with it.
   * For this, see the @p DoFAccessor::vertex_dof_index functions.
   *
   * @note Despite the name, the index returned here is only global in the
   * sense that it is specific to a particular Triangulation object or, in the
   * case the triangulation is actually of type
   * parallel::distributed::Triangulation, specific to that part of the
   * distributed triangulation stored on the current processor.
   */
  unsigned int
  vertex_index(const unsigned int i = 0) const;

  /**
   * Return a reference to the @p ith vertex. If i is zero, this returns a
   * reference to the current point to which this object refers. Otherwise, it
   * throws an exception.
   */
  Point<spacedim> &
  vertex(const unsigned int i = 0) const;

  /**
   * Return the center of this object, which of course coincides with the
   * location of the vertex this object refers to.
   */
  Point<spacedim>
  center() const;

  /**
   * Pointer to the @p ith line bounding this object. Will point to an invalid
   * object.
   */
  typename dealii::internal::TriangulationImplementation::
    Iterators<1, spacedim>::line_iterator static line(const unsigned int);

  /**
   * Line index of the @p ith line bounding this object.
   *
   * Implemented only for <tt>structdim>1</tt>, otherwise an exception
   * generated.
   */
  static unsigned int
  line_index(const unsigned int i);

  /**
   * Pointer to the @p ith quad bounding this object.
   */
  static typename dealii::internal::TriangulationImplementation::
    Iterators<1, spacedim>::quad_iterator
    quad(const unsigned int i);

  /**
   * Quad index of the @p ith quad bounding this object.
   *
   * Implemented only for <tt>structdim>2</tt>, otherwise an exception
   * generated.
   */
  static unsigned int
  quad_index(const unsigned int i);

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
   * Return 1.
   *
   * @note This is consistent with what ReferenceCells::Vertex::volume()
   * returns.
   */
  static double
  measure();

  /**
   * @}
   */

  /**
   * Return whether this point is at the boundary of the one-dimensional
   * triangulation we deal with here.
   */
  bool
  at_boundary() const;

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
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  types::boundary_id
  boundary_id() const;

  /**
   * Return a constant reference to the manifold object used for this object.
   */
  const Manifold<1, spacedim> &
  get_manifold() const;

  /**
   * Return the manifold indicator of this object.
   *
   * @see
   * @ref GlossManifoldIndicator "Glossary entry on manifold indicators".
   */
  types::manifold_id
  manifold_id() const;


  /**
   * @name User data
   */
  /**
   * @{
   */
  /**
   * Read the user flag. See
   * @ref GlossUserFlags
   * for more information.
   */
  bool
  user_flag_set() const;

  /**
   * Set the user flag. See
   * @ref GlossUserFlags
   * for more information.
   */
  void
  set_user_flag() const;

  /**
   * Clear the user flag. See
   * @ref GlossUserFlags
   * for more information.
   */
  void
  clear_user_flag() const;

  /**
   * Set the user flag for this and all descendants. See
   * @ref GlossUserFlags
   * for more information.
   */
  void
  recursively_set_user_flag() const;

  /**
   * Clear the user flag for this and all descendants. See
   * @ref GlossUserFlags
   * for more information.
   */
  void
  recursively_clear_user_flag() const;

  /**
   * Reset the user data to zero, independent if pointer or index. See
   * @ref GlossUserData
   * for more information.
   */
  void
  clear_user_data() const;

  /**
   * Set the user pointer to @p p.
   *
   * @note User pointers and user indices are mutually exclusive. Therefore,
   * you can only use one of them, unless you call
   * Triangulation::clear_user_data() in between.
   *
   * See
   * @ref GlossUserData
   * for more information.
   */
  void
  set_user_pointer(void *p) const;

  /**
   * Reset the user pointer to `nullptr`. See
   * @ref GlossUserData
   * for more information.
   */
  void
  clear_user_pointer() const;

  /**
   * Access the value of the user pointer. It is in the responsibility of the
   * user to make sure that the pointer points to something useful. You should
   * use the new style cast operator to maintain a minimum of type safety,
   * e.g.
   *
   * @note User pointers and user indices are mutually exclusive. Therefore,
   * you can only use one of them, unless you call
   * Triangulation::clear_user_data() in between. <tt>A
   * *a=static_cast<A*>(cell->user_pointer());</tt>.
   *
   * See
   * @ref GlossUserData
   * for more information.
   */
  void *
  user_pointer() const;

  /**
   * Set the user pointer of this object and all its children to the given
   * value. This is useful for example if all cells of a certain subdomain, or
   * all faces of a certain part of the boundary should have user pointers
   * pointing to objects describing this part of the domain or boundary.
   *
   * Note that the user pointer is not inherited under mesh refinement, so
   * after mesh refinement there might be cells or faces that don't have user
   * pointers pointing to the describing object. In this case, simply loop
   * over all the elements of the coarsest level that has this information,
   * and use this function to recursively set the user pointer of all finer
   * levels of the triangulation.
   *
   * @note User pointers and user indices are mutually exclusive. Therefore,
   * you can only use one of them, unless you call
   * Triangulation::clear_user_data() in between.
   *
   * See
   * @ref GlossUserData
   * for more information.
   */
  void
  recursively_set_user_pointer(void *p) const;

  /**
   * Clear the user pointer of this object and all of its descendants. The
   * same holds as said for the recursively_set_user_pointer() function. See
   * @ref GlossUserData
   * for more information.
   */
  void
  recursively_clear_user_pointer() const;

  /**
   * Set the user index to @p p.
   *
   * @note User pointers and user indices are mutually exclusive. Therefore,
   * you can only use one of them, unless you call
   * Triangulation::clear_user_data() in between. See
   * @ref GlossUserData
   * for more information.
   */
  void
  set_user_index(const unsigned int p) const;

  /**
   * Reset the user index to 0. See
   * @ref GlossUserData
   * for more information.
   */
  void
  clear_user_index() const;

  /**
   * Access the value of the user index.
   *
   * @note User pointers and user indices are mutually exclusive. Therefore,
   * you can only use one of them, unless you call
   * Triangulation::clear_user_data() in between.
   *
   * See
   * @ref GlossUserData
   * for more information.
   */
  unsigned int
  user_index() const;

  /**
   * Set the user index of this object and all its children.
   *
   * Note that the user index is not inherited under mesh refinement, so after
   * mesh refinement there might be cells or faces that don't have the
   * expected user indices. In this case, simply loop over all the elements of
   * the coarsest level that has this information, and use this function to
   * recursively set the user index of all finer levels of the triangulation.
   *
   * @note User pointers and user indices are mutually exclusive. Therefore,
   * you can only use one of them, unless you call
   * Triangulation::clear_user_data() in between.
   *
   * See
   * @ref GlossUserData
   * for more information.
   */
  void
  recursively_set_user_index(const unsigned int p) const;

  /**
   * Clear the user index of this object and all of its descendants. The same
   * holds as said for the recursively_set_user_index() function.
   *
   * See
   * @ref GlossUserData
   * for more information.
   */
  void
  recursively_clear_user_index() const;
  /**
   * @}
   */

  /**
   * @name Orientation of sub-objects
   */
  /**
   * @{
   */

  /**
   * @brief Always return 0
   *
   * @ingroup reordering
   */
  static types::geometric_orientation
  combined_face_orientation(const unsigned int face);

  /**
   * @brief Always return false
   */
  static bool
  face_orientation(const unsigned int face);

  /**
   * @brief Always return false
   */
  static bool
  face_flip(const unsigned int face);

  /**
   * @brief Always return false
   */
  static bool
  face_rotation(const unsigned int face);

  /**
   * @brief Always return numbers::reverse_line_orientation
   */
  static types::geometric_orientation
  line_orientation(const unsigned int line);

  /**
   * @}
   */

  /**
   * @name Accessing children
   */
  /**
   * @{
   */

  /**
   * Test whether the object has children. Always false.
   */
  static bool
  has_children();

  /**
   * Return the number of immediate children of this object.This is always
   * zero in dimension 0.
   */
  static unsigned int
  n_children();

  /**
   * Compute and return the number of active descendants of this objects.
   * Always zero.
   */
  static unsigned int
  n_active_descendants();

  /**
   * Return the number of times that this object is refined. Always 0.
   */
  static unsigned int
  max_refinement_depth();

  /**
   * @brief Return an invalid unsigned integer.
   */
  static unsigned int
  child_iterator_to_index(const TriaIterator<TriaAccessor<0, 1, spacedim>> &);

  /**
   * @brief Return an invalid object
   */
  static TriaIterator<TriaAccessor<0, 1, spacedim>>
  child(const unsigned int);

  /**
   * @brief Return an invalid object
   */
  static TriaIterator<TriaAccessor<0, 1, spacedim>>
  isotropic_child(const unsigned int);

  /**
   * Always return no refinement.
   */
  static RefinementCase<0>
  refinement_case();

  /**
   * @brief Returns -1
   */
  static int
  child_index(const unsigned int i);

  /**
   * @brief Returns -1
   */
  static int
  isotropic_child_index(const unsigned int i);
  /**
   * @}
   */

  /**
   * @name Dealing with boundary indicators
   */
  /**
   * @{
   */

  /**
   * Set the boundary indicator. The same applies as for the
   * <tt>boundary_id()</tt> function.
   *
   * @warning You should never set the boundary indicator of an interior face
   * (a face not at the boundary of the domain), or set the boundary
   * indicator of an exterior face to numbers::internal_face_boundary_id (this
   * value is reserved for another purpose). Algorithms may not work or
   * produce very confusing results if boundary cells have a boundary
   * indicator of numbers::internal_face_boundary_id or if interior cells have
   * boundary indicators other than numbers::internal_face_boundary_id.
   * Unfortunately, the current object has no means of finding out whether it
   * really is at the boundary of the domain and so cannot determine whether
   * the value you are trying to set makes sense under the current
   * circumstances.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  void
  set_boundary_id(const types::boundary_id) const;

  /**
   * Set the manifold indicator of this vertex. This does nothing so far since
   * manifolds are only used to refine and map objects, but vertices are not
   * refined and the mapping is trivial. This function is here only to allow
   * dimension independent programming.
   */
  void
  set_manifold_id(const types::manifold_id);

  /**
   * Set the boundary indicator of this object and all of its
   * lower-dimensional sub-objects.  Since this object only represents a single
   * vertex, there are no lower-dimensional object and this function is
   * equivalent to calling set_boundary_id() with the same argument.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  void
  set_all_boundary_ids(const types::boundary_id) const;

  /**
   * Set the manifold indicator of this object and all of its
   * lower-dimensional sub-objects.  Since this object only represents a single
   * vertex, there are no lower-dimensional object and this function is
   * equivalent to calling set_manifold_id() with the same argument.
   *
   * @ingroup manifold
   *
   * @see
   * @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  void
  set_all_manifold_ids(const types::manifold_id);
  /**
   * @}
   */

  /**
   * Return whether the vertex pointed to here is used.
   */
  bool
  used() const;

  /**
   * Reference cell type of the current object.
   */
  ReferenceCell
  reference_cell() const;

  /**
   * Number of vertices.
   */
  unsigned int
  n_vertices() const;

  /**
   * Number of lines.
   */
  unsigned int
  n_lines() const;

  /**
   * Return an object that can be thought of as an array containing all indices
   * from zero to n_vertices().
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  vertex_indices() const;

  /**
   * Return an object that can be thought of as an array containing all indices
   * from zero to n_lines().
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  line_indices() const;

protected:
  /**
   * Pointer to the triangulation we operate on.
   */
  const Triangulation<1, spacedim> *tria;

  /**
   * Whether this is a left end, right end, or interior vertex. This
   * information is provided by the cell at the time of creation.
   */
  VertexKind vertex_kind;

  /**
   * The global vertex index of the vertex this object corresponds to.
   */
  unsigned int global_vertex_index;
};



/**
 * This class allows access to a cell: a line in one dimension, a quad in two
 * dimension, etc.
 *
 * The following refers to any dimension:
 *
 * This class allows access to a <tt>cell</tt>, which is a line in 1d and a
 * quad in 2d. Cells have more functionality than lines or quads by
 * themselves, for example they can be flagged for refinement, they have
 * neighbors, they have the possibility to check whether they are at the
 * boundary etc. This class offers access to all this data.
 *
 * @ingroup grid
 * @ingroup Accessors
 */
template <int dim, int spacedim = dim>
class CellAccessor : public TriaAccessor<dim, dim, spacedim>
{
public:
  /**
   * Propagate the AccessorData type into the present class.
   */
  using AccessorData = typename TriaAccessor<dim, dim, spacedim>::AccessorData;

  /**
   * Define the type of the container this is part of.
   */
  using Container = Triangulation<dim, spacedim>;

  /**
   * @name Constructors
   */
  /**
   * @{
   */

  /**
   * Constructor.
   */
  CellAccessor(const Triangulation<dim, spacedim> *parent     = nullptr,
               const int                           level      = -1,
               const int                           index      = -1,
               const AccessorData                 *local_data = nullptr);

  /**
   * Copy constructor.
   */
  CellAccessor(const TriaAccessor<dim, dim, spacedim> &cell_accessor);

  /**
   * Conversion constructor. This constructor exists to make certain
   * constructs simpler to write in dimension independent code. For example,
   * it allows assigning a face iterator to a line iterator, an operation that
   * is useful in 2d but doesn't make any sense in 3d. The constructor here
   * exists for the purpose of making the code conform to C++ but it will
   * unconditionally abort; in other words, assigning a face iterator to a
   * line iterator is better put into an if-statement that checks that the
   * dimension is two, and assign to a quad iterator in 3d (an operator that,
   * without this constructor would be illegal if we happen to compile for
   * 2d).
   */
  template <int structdim2, int dim2, int spacedim2>
  CellAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * Another conversion operator between objects that don't make sense, just
   * like the previous one.
   */
  template <int structdim2, int dim2, int spacedim2>
  CellAccessor(const TriaAccessor<structdim2, dim2, spacedim2> &);

  /**
   * Copy constructor.
   */
  CellAccessor(const CellAccessor<dim, spacedim> &) = default;

  /**
   * Move constructor.
   */
  // NOLINTNEXTLINE OSX does not compile with noexcept
  CellAccessor(CellAccessor<dim, spacedim> &&) = default;

  /**
   * Destructor.
   */
  ~CellAccessor() = default;

  /**
   * Copy operator. These operators are usually used in a context like
   * <tt>iterator a,b; *a=*b;</tt>. Presumably, the intent here is to copy the
   * object pointed to
   * by @p b to the object pointed to by @p a. However, the result of
   * dereferencing an iterator is not an object but an accessor; consequently,
   * this operation is not useful for iterators on triangulations.
   * Consequently, this operator is declared as deleted and can not be used.
   */
  CellAccessor<dim, spacedim> &
  operator=(const CellAccessor<dim, spacedim> &) = delete;

  /**
   * Move assignment operator.
   */
  // NOLINTNEXTLINE OSX does not compile with noexcept
  CellAccessor<dim, spacedim> &
  operator=(CellAccessor<dim, spacedim> &&) = default; // NOLINT

  /**
   * @}
   */

  /**
   * @name Converting iterators
   */
  /**
   * @{
   */

  /**
   * A function that converts a Triangulation active cell iterator to a
   * DoFHandler active cell iterator, or a DoFHandler active cell iterator
   * to an active cell iterator of another DoFHandler. The @p iterator must be
   * associated with the triangulation of the @p dof_handler.
   *
   * @param dof_handler The DoFHandler for the output active cell iterator.
   * @return An active cell iterator for the @p dof_handler, matching the cell
   *         referenced by the input @p iterator. The type of the
   *         returned object is a DoFHandler::active_cell_iterator.
   */
  TriaActiveIterator<DoFCellAccessor<dim, spacedim, false>>
  as_dof_handler_iterator(const DoFHandler<dim, spacedim> &dof_handler) const;

  /**
   * A function similar to as_dof_handler_iterator(). It converts
   * a Triangulation active/level cell iterator to a
   * DoFHandler active level cell iterator, or a DoFHandler level active/cell
   * iterator to a level cell iterator of another DoFHandler. The
   * @p iterator must be associated with the triangulation of the
   * @p dof_handler.
   */
  TriaIterator<DoFCellAccessor<dim, spacedim, true>>
  as_dof_handler_level_iterator(
    const DoFHandler<dim, spacedim> &dof_handler) const;


  /**
   * @}
   */

  /**
   * @name Accessing sub-objects and neighbors
   */
  /**
   * @{
   */

  /**
   * Return a pointer to the @p ith child. Overloaded version which returns a
   * more reasonable iterator class.
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  child(const unsigned int i) const;

  /**
   * Return an array of iterators to all children of this cell.
   */
  boost::container::small_vector<TriaIterator<CellAccessor<dim, spacedim>>,
                                 GeometryInfo<dim>::max_children_per_cell>
  child_iterators() const;

  /**
   * Return an iterator to the @p ith face of this cell.
   */
  TriaIterator<TriaAccessor<dim - 1, dim, spacedim>>
  face(const unsigned int i) const;

  /**
   * Return the face number of @p face on the current cell. This is the
   * inverse function of TriaAccessor::face().
   */
  unsigned int
  face_iterator_to_index(
    const TriaIterator<TriaAccessor<dim - 1, dim, spacedim>> &face) const;

  /**
   * Return an array of iterators to all faces of this cell.
   */
  boost::container::small_vector<
    TriaIterator<TriaAccessor<dim - 1, dim, spacedim>>,
#ifndef _MSC_VER // MSVC prior to 2022 cannot use a constexpr function this way
    ReferenceCells::max_n_faces<dim>()
#else
    GeometryInfo<dim>::faces_per_cell
#endif
    >
  face_iterators() const;

  /**
   * Return the (global) index of the @p ith face of this cell.
   *
   * @note Despite the name, the index returned here is only global in the
   * sense that it is specific to a particular Triangulation object or, in the
   * case the triangulation is actually of type
   * parallel::distributed::Triangulation, specific to that part of the
   * distributed triangulation stored on the current processor.
   */
  unsigned int
  face_index(const unsigned int i) const;

  /**
   * Return an iterator to that cell that neighbors the present cell on the
   * given face and subface number.
   *
   * To succeed, the present cell must not be further refined, and the
   * neighbor on the given face must be further refined exactly once; the
   * returned cell is then a child of that neighbor.
   *
   * The function may not be called in 1d, since there we have no subfaces.
   * The implementation of this function is rather straightforward in 2d, by
   * first determining which face of the neighbor cell the present cell is
   * bordering on (this is what the @p neighbor_of_neighbor function does),
   * and then asking @p GeometryInfo::child_cell_on_subface for the index of
   * the child.
   *
   * However, the situation is more complicated in 3d, since there faces may
   * have more than one orientation, and we have to use @p face_orientation,
   * @p face_flip and @p face_rotation for both this and the neighbor cell to
   * figure out which cell we want to have.
   *
   * This can lead to surprising results: if we are sitting on a cell and are
   * asking for a cell behind subface <tt>sf</tt>, then this means that we are
   * considering the subface for the face in the natural direction for the
   * present cell. However, if the face as seen from this cell has
   * <tt>face_orientation()==false</tt>, then the child of the face that
   * separates the present cell from the neighboring cell's child is not
   * necessarily the @p sf-th child of the face of this cell. This is so
   * because the @p subface_no on a cell corresponds to the subface with
   * respect to the intrinsic ordering of the present cell, whereas children
   * of face iterators are computed with respect to the intrinsic ordering of
   * faces; these two orderings are only identical if the face orientation is
   * @p true, and reversed otherwise.
   *
   * Similarly, effects of <tt>face_flip()==true</tt> and
   * <tt>face_rotation()==true()</tt>, both of which indicate a non-standard
   * face have to be considered.
   *
   * Fortunately, this is only very rarely of concern, since usually one
   * simply wishes to loop over all finer neighbors at a given face of an
   * active cell. Only in the process of refinement of a Triangulation we want
   * to set neighbor information for both our child cells and the neighbor's
   * children. Since we can respect orientation of faces from our current cell
   * in that case, we do NOT respect face_orientation, face_flip and
   * face_rotation of the present cell within this function, i.e. the returned
   * neighbor's child is behind subface @p subface concerning the intrinsic
   * ordering of the given face.
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  neighbor_child_on_subface(const unsigned int face_no,
                            const unsigned int subface_no) const;

  /**
   * Return an iterator to the neighboring cell on the other side of the face
   * with number @p face_no. If the neighbor does not exist,
   * i.e., if the face with number @p face_no of the current object is at the boundary, then
   * an invalid iterator is returned. In detail, the smallest cell `neighbor`
   * for which `cell->face(face_no)` is a subset of
   * `neighbor->face(opposite_face_no)`, where `opposite_face_no` is the face
   * number opposite to `face_no`.
   *
   * Consequently, the index @p face_no must be less than n_faces().
   *
   * For example, consider the following situation:
   * @image html limit_level_difference_at_vertices_anisotropic.png ""
   *
   * Here, if you are on cell `1.3` and ask for its left neighbor (which is,
   * according to the conventions spelled out in the GeometryInfo class, its
   * <i>zeroth</i> neighbor), then you will get the parent cell of `3.5`, since
   * this is the smallest cell for which we have `(1.3)->face(0) ==
   * (3.5)->parent()->face(1)`. Note, that you will not obtain the parent cell
   * of `2.8`.
   *
   * Further, if you ask for the right (i.e. the <i>first</i>) neighbor of cell
   * `4.1`, then you will get cell `1.3`. Consequently, there are two
   * neighboring cells that differ by three in their levels. In fact, using
   * anisotropic refinement it is possible to generate arbitrary large
   * differences in the level of neighboring cells. Perform e.g. arbitrarily
   * many `y`-refinements of cell `4.1` and its children. While the second and
   * third neighbors are being refined as well, due to the avoidance of multiple
   * hanging nodes, cell `1.3` always remains as neighbor of the resulting
   * right-most child.
   *
   * On the other hand, if you were at cell `3.3` and ask for its third
   * neighbor, cell `4.1` will be returned, since it is the smallest cell that
   * fulfills the described property. This shows, that the neighbor <i>can</i>
   * have a higher level than the cell itself.
   *
   * However, using only isotropic refinement, the neighbor will have <i>at
   * most</i> the level as the cell itself. This can be verified in the bottom
   * half of the image, where only isotropic refinement was done: The first
   * neighbor of `3.3` is given by `2.6`. Further refinement of `2.6` will
   * generate a new first neighbor with level 3, but any further refinements of
   * that child will not affect the neighbor of cell `3.3`. Due to the avoidance
   * of multiple hanging nodes on a mesh, it is also impossible to obtain a
   * coarser cell than `2.6` as the first neighbor of `3.3`. Consequently, the
   * neighbor of a fully isotropic refined mesh has either the same level as the
   * cell itself, or is exactly one level coarser.
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  neighbor(const unsigned int face_no) const;

  /**
   * Return the cell index of the neighboring cell on the other side of the face
   * with index @p face_no. If the neighbor does not exist, this function returns -1.
   *
   * This function is equivalent to <tt>cell->neighbor(face_no)->index()</tt>.
   * See neighbor() for more details.
   */
  int
  neighbor_index(const unsigned int face_no) const;

  /**
   * Return the level of the neighboring cell on the other side of the face with
   * number @p face_no. If the neighbor does not exist, this function returns -1.
   *
   * This function is equivalent to `cell->neighbor(face_no)->level()`.
   * See neighbor() for more details.
   */
  int
  neighbor_level(const unsigned int face_no) const;

  /**
   * Return the how-many'th neighbor this cell is of
   * <tt>cell->neighbor(face_no)</tt>, i.e. return @p other_face_no such that
   * <tt>cell->neighbor(face_no)->neighbor(other_face_no)==cell</tt>. This
   * function is the right one if you want to know how to get back from a
   * neighbor to the present cell.
   *
   * Note that this operation is only useful if the neighbor is not coarser
   * than the present cell. If the neighbor is coarser this function throws an
   * exception. Use the @p neighbor_of_coarser_neighbor function in that case.
   */
  unsigned int
  neighbor_of_neighbor(const unsigned int face_no) const;

  /**
   * Return, whether the neighbor is coarser then the present cell. This is
   * important in case of anisotropic refinement where this information does
   * not depend on the levels of the cells.
   *
   * Note, that in an anisotropic setting, a cell can only be coarser than
   * another one at a given face, not on a general basis. The face of the
   * finer cell is contained in the corresponding face of the coarser cell,
   * the finer face is either a child or a grandchild of the coarser face.
   */
  bool
  neighbor_is_coarser(const unsigned int face_no) const;

  /**
   * This function is a generalization of the @p neighbor_of_neighbor function
   * for the case of a coarser neighbor. It returns a pair of numbers, face_no
   * and subface_no, with the following property, if the neighbor is not
   * refined: <tt>cell->neighbor(neighbor)->neighbor_child_on_subface(face_no,
   * subface_no)==cell</tt>. In 3d, a coarser neighbor can still be refined.
   * In that case subface_no denotes the child index of the neighbors face
   * that relates to our face:
   * <tt>cell->neighbor(neighbor)->face(face_no)->child(subface_no)==cell->face(neighbor)</tt>.
   * This case in 3d and how it can happen is discussed in the introduction of
   * the step-30 tutorial program.
   *
   * This function is impossible for <tt>dim==1</tt>.
   */
  std::pair<unsigned int, unsigned int>
  neighbor_of_coarser_neighbor(const unsigned int neighbor) const;

  /**
   * This function is a generalization of the @p neighbor_of_neighbor and the
   * @p neighbor_of_coarser_neighbor functions. It checks whether the neighbor
   * is coarser or not and calls the respective function. In both cases, only
   * the face_no is returned.
   */
  unsigned int
  neighbor_face_no(const unsigned int neighbor) const;

  /**
   * Compatibility interface with DoFCellAccessor. Always returns @p false.
   */
  static bool
  is_level_cell();

  /**
   * @}
   */
  /**
   * @name Dealing with periodic neighbors
   */
  /**
   * @{
   */
  /**
   * If the cell has a periodic neighbor at its @c ith face, this function
   * returns true, otherwise, the returned value is false.
   */
  bool
  has_periodic_neighbor(const unsigned int i) const;

  /**
   * For a cell with its @c ith face at a periodic boundary,
   * see
   * @ref GlossPeriodicConstraints "the entry for periodic boundaries",
   * this function returns an iterator to the cell on the other side
   * of the periodic boundary. If there is no periodic boundary at the @c ith
   * face, an exception will be thrown.
   * In order to avoid running into an exception, check the result of
   * has_periodic_neighbor() for the @c ith face prior to using this function.
   * The behavior of periodic_neighbor() is similar to neighbor(), in
   * the sense that the returned cell has at most the same level of refinement
   * as the current cell. On distributed meshes, by calling
   * Triangulation::add_periodicity(),
   * we can make sure that the element on the other side of the periodic
   * boundary exists in this rank as a ghost cell or a locally owned cell.
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  periodic_neighbor(const unsigned int i) const;

  /**
   * For a cell whose @c ith face is not at a boundary, this function returns
   * the same result as neighbor(). If the @c ith face is at a periodic boundary
   * this function returns the same result as periodic_neighbor(). If neither of
   * the aforementioned conditions are met, i.e. the @c ith face is on a
   * nonperiodic boundary, an exception will be thrown.
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  neighbor_or_periodic_neighbor(const unsigned int i) const;

  /**
   * Return an iterator to the periodic neighbor of the cell at a given
   * face and subface number. The general guidelines for using this function
   * is similar to the function neighbor_child_on_subface(). The
   * implementation of this function is consistent with
   * periodic_neighbor_of_coarser_periodic_neighbor(). For instance,
   * assume that we are sitting on a cell named @c cell1, whose neighbor behind
   * the @c ith face is one level coarser. Let us name this coarser neighbor
   * @c cell2. Then, by calling
   * periodic_neighbor_of_coarser_periodic_neighbor(), from @c cell1, we get
   * a @c face_num and a @c subface_num. Now, if we call
   * periodic_neighbor_child_on_subface() from cell2, with the above face_num
   * and subface_num, we get an iterator to @c cell1.
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  periodic_neighbor_child_on_subface(const unsigned int face_no,
                                     const unsigned int subface_no) const;

  /**
   * This function is a generalization of
   * periodic_neighbor_of_periodic_neighbor()
   * for those cells which have a coarser periodic neighbor. The returned
   * pair of numbers can be used in periodic_neighbor_child_on_subface()
   * to get back to the current cell. In other words, the following
   * assertion should be true, for a cell with coarser periodic neighbor:
   * cell->periodic_neighbor(i)->periodic_neighbor_child_on_subface(face_no,
   * subface_no)==cell
   */
  std::pair<unsigned int, unsigned int>
  periodic_neighbor_of_coarser_periodic_neighbor(const unsigned face_no) const;

  /**
   * This function returns the index of the periodic neighbor at the
   * @c ith face of the current cell. If there is no periodic neighbor
   * at the given face, the returned value is -1.
   */
  int
  periodic_neighbor_index(const unsigned int i) const;

  /**
   * This function returns the level of the periodic neighbor at the
   * @c ith face of the current cell. If there is no periodic neighbor
   * at the given face, the returned value is -1.
   */
  int
  periodic_neighbor_level(const unsigned int i) const;

  /**
   * For a cell with a periodic neighbor at its @c ith face, this function
   * returns the face number of that periodic neighbor such that the
   * current cell is the periodic neighbor of that neighbor. In other words
   * the following assertion holds for those cells which have a periodic
   * neighbor with the same or a higher level of refinement as the current
   * cell:
   * @c {cell->periodic_neighbor(i)->
   *     periodic_neighbor(cell->periodic_neighbor_of_periodic_neighbor(i))==cell}
   * For the cells with a coarser periodic neighbor, one should use
   * periodic_neighbor_of_coarser_periodic_neighbor() and
   * periodic_neighbor_child_on_subface()
   * to get back to the current cell.
   */
  unsigned int
  periodic_neighbor_of_periodic_neighbor(const unsigned int i) const;

  /**
   * If a cell has a periodic neighbor at its @c ith face, this function
   * returns the face number of the periodic neighbor, which is connected
   * to the @c ith face of this cell.
   */
  unsigned int
  periodic_neighbor_face_no(const unsigned int i) const;

  /**
   * This function returns true if the element on the other side of the
   * periodic boundary is coarser and returns false otherwise. The
   * implementation allows this function to work in the case of
   * anisotropic refinement.
   */
  bool
  periodic_neighbor_is_coarser(const unsigned int i) const;

  /**
   * @}
   */

  /**
   * @name Dealing with boundary indicators
   */
  /**
   * @{
   */

  /**
   * Return whether the @p ith vertex or face (depending on the dimension) is
   * part of the boundary. This is true, if the @p ith neighbor does not
   * exist.
   */
  bool
  at_boundary(const unsigned int i) const;

  /**
   * Return whether the cell is at the boundary. Being at the boundary is
   * defined by one face being on the boundary. Note that this does not catch
   * cases where only one vertex of a 2d or 3d subobject is at the boundary,
   * or where only one line of a hex is at the boundary while the interiors of
   * all faces are in the interior of the domain. For the latter case, the @p
   * has_boundary_lines function is the right one to ask.
   */
  bool
  at_boundary() const;

  /**
   * This is a slight variation to the @p at_boundary function: for 1 and 2
   * dimensions, it is equivalent, for three dimensions it returns whether at
   * least one of the 12 lines of the hexahedron is at a boundary. This, of
   * course, includes the case where a whole face is at the boundary, but also
   * some other cases.
   */
  bool
  has_boundary_lines() const;
  /**
   * @}
   */

  /**
   * @name Dealing with refinement indicators
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
  RefinementCase<dim>
  refine_flag_set() const;

  /**
   * Flag the cell pointed to for refinement. This function is only allowed
   * for active cells. Keeping the default value for @p ref_case will mark
   * this cell for isotropic refinement.
   *
   * If you choose anisotropic refinement, for example by passing as argument
   * one of the flags RefinementCase::cut_x, RefinementCase::cut_y,
   * RefinementCase::cut_z, or a combination of these, then keep in mind
   * that refining in x-, y-, or z-direction happens with regard to the
   * <em>local</em> coordinate system of the cell. In other words, these
   * flags determine which edges and faces of the cell will be cut into new
   * edges and faces. On the other hand, this process is independent of
   * how the cell is oriented within the <em>global</em> coordinate system,
   * and you should not assume any particular orientation of the cell's
   * local coordinate system within the global coordinate system of the
   * space it lives in.
   */
  void
  set_refine_flag(const RefinementCase<dim> ref_case =
                    RefinementCase<dim>::isotropic_refinement) const;

  /**
   * Clear the refinement flag.
   */
  void
  clear_refine_flag() const;

  /**
   * Return the @p IsotropicRefinementChoices this cell was flagged to be refined
   * with.
   */
  std::uint8_t
  refine_choice() const;

  /**
   * Set the @p IsotropicRefinementChoices this cell is flagged to be refined
   * with.
   */
  void
  set_refine_choice(const std::uint8_t refinement_choice = static_cast<char>(
                      IsotropicRefinementChoice::isotropic_refinement)) const;

  /**
   * Clear the @p IsotropicRefinementChoices flag.
   */
  void
  clear_refine_choice() const;

  /**
   * Modify the refinement flag of the cell to ensure (at least) the given
   * refinement case @p face_refinement_case at face <tt>face_no</tt>, taking
   * into account orientation, flip and rotation of the face. Return, whether
   * the refinement flag had to be modified. This function is only allowed for
   * active cells.
   */
  bool
  flag_for_face_refinement(
    const unsigned int             face_no,
    const RefinementCase<dim - 1> &face_refinement_case =
      RefinementCase<dim - 1>::isotropic_refinement) const;

  /**
   * Modify the refinement flag of the cell to ensure that line
   * <tt>face_no</tt> will be refined. Return, whether the refinement flag had
   * to be modified. This function is only allowed for active cells.
   */
  bool
  flag_for_line_refinement(const unsigned int line_no) const;

  /**
   * Return the SubfaceCase of face <tt>face_no</tt>. Note that this is not
   * identical to asking <tt>cell->face(face_no)->refinement_case()</tt> since
   * the latter returns a RefinementCase<dim-1> and thus only considers one
   * (anisotropic) refinement, whereas this function considers the complete
   * refinement situation including possible refinement of the face's
   * children. This function may only be called for active cells in 2d and 3d.
   */
  internal::SubfaceCase<dim>
  subface_case(const unsigned int face_no) const;

  /**
   * Return whether the coarsen flag is set or not.
   */
  bool
  coarsen_flag_set() const;

  /**
   * Flag the cell pointed to for coarsening. This function is only allowed
   * for active cells.
   */
  void
  set_coarsen_flag() const;

  /**
   * Clear the coarsen flag.
   */
  void
  clear_coarsen_flag() const;
  /**
   * @}
   */

  /**
   * @name Dealing with material indicators
   */
  /**
   * @{
   */

  /**
   * Return the material id of this cell.
   *
   * For a typical use of this function, see the
   * @ref step_28 "step-28"
   * tutorial program.
   *
   * See the
   * @ref GlossMaterialId "glossary"
   * for more information.
   */
  types::material_id
  material_id() const;

  /**
   * Set the material id of this cell.
   *
   * For a typical use of this function, see the
   * @ref step_28 "step-28"
   * tutorial program.
   *
   * See the
   * @ref GlossMaterialId "glossary"
   * for more information.
   */
  void
  set_material_id(const types::material_id new_material_id) const;

  /**
   * Set the material id of this cell and all its children (and
   * grand-children, and so on) to the given value.
   *
   * See the
   * @ref GlossMaterialId "glossary"
   * for more information.
   */
  void
  recursively_set_material_id(const types::material_id new_material_id) const;
  /**
   * @}
   */

  /**
   * @name Dealing with subdomain indicators
   */
  /**
   * @{
   */

  /**
   * Return the subdomain id of this cell.
   *
   * See the
   * @ref GlossSubdomainId "glossary"
   * for more information.
   *
   * @note The subdomain of a cell is a property only defined for active
   * cells, i.e., cells that are not further refined. Consequently, you can
   * only call this function if the cell it refers to has no children. For
   * multigrid methods in parallel, it is also important to know which
   * processor owns non-active cells, and for this you can call
   * level_subdomain_id().
   */
  types::subdomain_id
  subdomain_id() const;

  /**
   * Set the subdomain id of this cell.
   *
   * See the
   * @ref GlossSubdomainId "glossary"
   * for more information. This function should not be called if you use a
   * parallel::distributed::Triangulation object.
   *
   * @note The subdomain of a cell is a property only defined for active
   * cells, i.e., cells that are not further refined. Consequently, you can
   * only call this function if the cell it refers to has no children. For
   * multigrid methods in parallel, it is also important to know which
   * processor owns non-active cells, and for this you can call
   * level_subdomain_id().
   */
  void
  set_subdomain_id(const types::subdomain_id new_subdomain_id) const;

  /**
   * Get the level subdomain id of this cell. This is used for parallel
   * multigrid where not only the global mesh (consisting of the active cells)
   * is partitioned among processors, but also the individual levels of the
   * hierarchy of recursively refined cells that make up the mesh. In
   * other words, the level subdomain id is a property that is also defined
   * for non-active cells if a multigrid hierarchy is used.
   */
  types::subdomain_id
  level_subdomain_id() const;

  /**
   * Set the level subdomain id of this cell. This is used for parallel
   * multigrid.
   */
  void
  set_level_subdomain_id(
    const types::subdomain_id new_level_subdomain_id) const;


  /**
   * Set the subdomain id of this cell (if it is active) or all its terminal
   * children (and grand-children, and so on, as long as they have no children
   * of their own) to the given value. Since the subdomain id is a concept
   * that is only defined for cells that are active (i.e., have no children of
   * their own), this function only sets the subdomain ids for all children
   * and grand children of this cell that are actually active, skipping
   * intermediate child cells.
   *
   * See the
   * @ref GlossSubdomainId "glossary"
   * for more information. This function should not be called if you use a
   * parallel::distributed::Triangulation object since there the subdomain id
   * is implicitly defined by which processor you're on.
   */
  void
  recursively_set_subdomain_id(
    const types::subdomain_id new_subdomain_id) const;
  /**
   * @}
   */

  /**
   * Return a globally unique cell index for the current cell,
   * assuming it is not artificial. The value is identical to
   * active_cell_index() if the cell is part of a serial
   * triangulation.
   *
   * In the context of parallel triangulations, locally-owned cells
   * are enumerated contiguously within each subdomain of the
   * mesh. This ensures that the index returned by this function can
   * be used as the index into vectors with a total of
   * Triangulation::n_globally_active_cells() entries, and for which
   * every process stores a contiguous part.  If such a cell-data
   * vector has been set up with
   * parallel::TriangulationBase::global_active_cell_index_partitioner(),
   * the index returned by this function can then be used to access
   * the correct vector entry.
   */
  types::global_cell_index
  global_active_cell_index() const;

  /**
   * Return a globally unique index for a non-artificial level cell.
   *
   * @note Similar to global_active_cell_index(), with the difference
   * that the cell-data vector has been set up with
   * parallel::TriangulationBase::global_level_cell_index_partitioner().
   */
  types::global_cell_index
  global_level_cell_index() const;

  /**
   * @name Dealing with codim 1 cell orientation
   */
  /**
   * @{
   */

  /**
   * Return the orientation of this cell. This function always returns
   * `true` if `dim==spacedim`. It can return `true` or `false` if
   * `dim==spacedim-1`. The function cannot be called (and will abort
   * with an error) if called for `dim<spacedim-1`.
   *
   * For the meaning of this flag, see
   * @ref GlossDirectionFlag.
   */
  bool
  direction_flag() const;

  /**
   * Return the how many-th active cell the current cell is (assuming the
   * current cell is indeed active). This is useful, for example, if you are
   * accessing the elements of a vector with as many entries as there are
   * active cells. Such vectors are used for estimating the error on each cell
   * of a triangulation, for specifying refinement criteria passed to the
   * functions in GridRefinement, and for generating cell-wise output.
   *
   * The function throws an exception if the current cell is not active.
   *
   * @note If the triangulation this function is called on is of type
   * parallel::distributed::Triangulation, then active cells may be locally
   * owned, ghost cells, or artificial (see
   * @ref GlossLocallyOwnedCell,
   * @ref GlossGhostCell,
   * and
   * @ref GlossArtificialCell).
   * This function counts over all of them, including ghost and artificial
   * active cells. This implies that the index returned by this function
   * uniquely identifies a cell within the triangulation on a single
   * processor, but does not uniquely identify the cell among the (parts of
   * the) triangulation that is shared among processors. If you would like to
   * identify active cells across processors, you need to consider the CellId
   * of a cell returned by CellAccessor::id().
   */
  unsigned int
  active_cell_index() const;

  /**
   * Return the index of the parent of this cell within the level of the
   * triangulation to which the parent cell belongs. The level of the parent
   * is of course one lower than that of the present cell. If the parent does
   * not exist (i.e., if the object is at the coarsest level of the mesh
   * hierarchy), an exception is generated.
   */
  int
  parent_index() const;

  /**
   * Return an iterator to the parent. If the parent does not exist (i.e., if
   * the object is at the coarsest level of the mesh hierarchy), an exception
   * is generated.
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  parent() const;

  /**
   * @}
   */

  /**
   * @name Other functions
   */
  /**
   * @{
   */

  /**
   * Test that the cell has no children (this is the criterion for whether a
   * cell is called "active").
   *
   * See the
   * @ref GlossActive "glossary"
   * for more information.
   */
  bool
  is_active() const;

  /**
   * Return whether this cell is owned by the current processor or is owned by
   * another processor. The function always returns true if applied to an
   * object of type dealii::Triangulation, but may yield false if the
   * triangulation is of type parallel::distributed::Triangulation.
   *
   * See the
   * @ref GlossGhostCell "glossary"
   * and the
   * @ref distributed
   * topic for more information.
   *
   * @post The returned value is equal to <code>!is_ghost() &&
   * !is_artificial()</code>.
   *
   * @note Whether a cell is a ghost cell, artificial, or is locally owned or
   * is a property that only pertains to cells that are active. Consequently,
   * you can only call this function if the cell it refers to has no children.
   */
  bool
  is_locally_owned() const;

  /**
   * Return true if either the Triangulation is not distributed or if
   * level_subdomain_id() is equal to the id of the current processor.
   */
  bool
  is_locally_owned_on_level() const;

  /**
   * Return true if:
   * <ol>
   * <li>This cell exists in the global mesh (i.e., it is not artificial),
   * and</li>
   * <li>This cell is owned by another processor (i.e., has a subdomain_id
   * different from Triangulation::locally_owned_subdomain())</li>
   * </ol>
   *
   * In all other cases the returned value is false. In particular, only
   * parallel Triangulations (i.e., Triangulations inheriting from
   * parallel::TriangulationBase) can have ghost cells, so for a serial
   * Triangulation the returned value is false.
   *
   * See the
   * @ref GlossGhostCell "glossary"
   * and the
   * @ref distributed
   * topic for more information.
   *
   * @post The returned value is equal to <code>!is_locally_owned() &&
   * !is_artificial()</code>.
   *
   * @note For parallel::distributed::Triangulation and
   * parallel::fullydistributed::Triangulation, ghost cells are always adjacent
   * to locally owned cells. For parallel::shared::Triangulation they may not
   * be, dependent on whether or not the triangulation uses artificial cells -
   * see parallel::shared::Triangulation::Triangulation() for more information.
   *
   * @note Whether a cell is a ghost cell, artificial, or is locally owned
   * is a property that only pertains to cells that are active. Consequently,
   * you can only call this function if the cell it refers to has no children.
   */
  bool
  is_ghost() const;

  /**
   * Return true if either the Triangulation is not distributed or if the
   * cell is not artificial and the level_subdomain_id() is not equal to the id
   * of the current processor.
   */
  bool
  is_ghost_on_level() const;

  /**
   * Return whether this cell is artificial, i.e. it isn't one of the cells
   * owned by the current processor, and it also doesn't border on one. As a
   * consequence, it exists in the mesh to ensure that each processor has all
   * coarse mesh cells and that the 2:1 ratio of neighboring cells is
   * maintained, but it is not one of the cells we should work on on the
   * current processor. In particular, there is no guarantee that this cell
   * isn't, in fact, further refined on one of the other processors.
   *
   * This function only makes sense if the triangulation used is of kind
   * parallel::distributed::Triangulation. In all other cases, the returned
   * value is always false.
   *
   * See the
   * @ref GlossArtificialCell "glossary"
   * and the
   * @ref distributed
   * topic for more information.
   *
   * @post The returned value is equal to <code>!is_ghost() &&
   * !is_locally_owned()</code>.
   *
   * @note Whether a cell is a ghost cell, artificial, or is locally owned is
   * a property that only pertains to cells that are active. Consequently, you
   * can only call this function if the cell it refers to has no children.
   */
  bool
  is_artificial() const;

  /**
   * Similar to is_artificial() but checking the conditions on the levels.
   *
   * @post The returned value is equal to <code>!is_ghost_on_level() &&
   * !is_locally_owned_on_level()</code>.
   */
  bool
  is_artificial_on_level() const;

  /**
   * Test whether the point @p p is inside this cell. Points on the boundary
   * are counted as being inside the cell.
   *
   * Note that this function assumes that the mapping between unit cell and
   * real cell is (bi-, tri-)linear, i.e. that faces in 2d and edges in 3d are
   * straight lines. If you have higher order transformations, results may be
   * different as to whether a point is in- or outside the cell in real space.
   *
   * In case of codim>0, the point is first projected to the manifold where
   * the cell is embedded and then check if this projection is inside the
   * cell.
   */
  bool
  point_inside(const Point<spacedim> &p) const;

  /**
   * Set the neighbor @p i of this cell to the cell pointed to by @p pointer.
   *
   * This function shouldn't really be public (but needs to for various
   * reasons in order not to make a long list of functions friends): it
   * modifies internal data structures and may leave things. Do not use it
   * from application codes.
   */
  void
  set_neighbor(const unsigned int                               i,
               const TriaIterator<CellAccessor<dim, spacedim>> &pointer) const;

  /**
   * Return a unique ID for the current cell. This ID is constructed from the
   * path in the hierarchy from the coarse parent cell and works correctly in
   * parallel computations using objects of type
   * parallel::distributed::Triangulation. This function is therefore useful
   * in providing a unique identifier for cells (active or not) that also
   * works for parallel triangulations. See the documentation of the CellId
   * class for more information.
   *
   * @note This operation takes O(level) time to compute. In most practical
   * cases, the number of levels of a triangulation will depend
   * logarithmically on the number of cells in the triangulation.
   */
  CellId
  id() const;

  using TriaAccessor<dim, dim, spacedim>::diameter;

  /**
   * The same as TriaAccessor::diameter() but also taking a Mapping class.
   */
  double
  diameter(const Mapping<dim, spacedim> &mapping) const;

  /**
   * @}
   */


  /**
   * @ingroup Exceptions
   */
  DeclException0(ExcRefineCellNotActive);
  /**
   * @ingroup Exceptions
   */
  DeclException0(ExcCellFlaggedForRefinement);
  /**
   * @ingroup Exceptions
   */
  DeclException0(ExcCellFlaggedForCoarsening);

protected:
  /**
   * This function assumes that the neighbor is not coarser than the current
   * cell. In this case it returns the neighbor_of_neighbor() value. If,
   * however, the neighbor is coarser this function returns an
   * <code>invalid_unsigned_int</code>.
   *
   * This function is not for public use. Use the function
   * neighbor_of_neighbor() instead which throws an exception if called for a
   * coarser neighbor. If neighbor is indeed coarser (you get to know this by
   * e.g. the neighbor_is_coarser() function) then the
   * neighbor_of_coarser_neighbor() function should be call. If you'd like to
   * know only the <code>face_no</code> which is required to get back from the
   * neighbor to the present cell then simply use the neighbor_face_no()
   * function which can be used for coarser as well as non-coarser neighbors.
   */
  unsigned int
  neighbor_of_neighbor_internal(const unsigned int neighbor) const;

  /**
   * As for any codim>0 we can use a similar code and c++ does not allow
   * partial templates. we use this auxiliary function that is then called
   * from point_inside.
   */
  template <int dim_, int spacedim_>
  bool
  point_inside_codim(const Point<spacedim_> &p) const;



private:
  /**
   * Set the active cell index of a cell. This is done at the end of
   * refinement.
   */
  void
  set_active_cell_index(const unsigned int active_cell_index) const;

  /**
   * Set global active cell index for a cell.
   */
  void
  set_global_active_cell_index(const types::global_cell_index index) const;

  /**
   * Set global level cell index for a level cell.
   */
  void
  set_global_level_cell_index(const types::global_cell_index index) const;

  /**
   * Set the parent of a cell.
   */
  void
  set_parent(const unsigned int parent_index);

  /**
   * Set the orientation of this cell. This function can only be
   * called if the argument is `true` if `dim==spacedim`. It can be
   * called with either `true` or `false` if `dim==spacedim-1`. The
   * function cannot be called (and will abort with an error) if called
   * for `dim<spacedim-1`.
   *
   * For the meaning of this flag, see
   * @ref GlossDirectionFlag.
   */
  void
  set_direction_flag(const bool new_direction_flag) const;

  friend class Triangulation<dim, spacedim>;

  friend class parallel::TriangulationBase<dim, spacedim>;

  friend struct dealii::internal::TriangulationImplementation::Implementation;
  friend struct dealii::internal::TriangulationImplementation::
    ImplementationMixedMesh;
};



/* ----- declaration of explicit specializations and general templates ----- */


template <int structdim, int dim, int spacedim>
template <typename OtherAccessor>
InvalidAccessor<structdim, dim, spacedim>::InvalidAccessor(
  const OtherAccessor &)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
TriaAccessor<structdim, dim, spacedim>::TriaAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
CellAccessor<dim, spacedim>::CellAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
TriaAccessor<structdim, dim, spacedim>::TriaAccessor(
  const TriaAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
CellAccessor<dim, spacedim>::CellAccessor(
  const TriaAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}


#ifndef DOXYGEN

template <>
bool
CellAccessor<1, 1>::point_inside(const Point<1> &) const;
template <>
bool
CellAccessor<2, 2>::point_inside(const Point<2> &) const;
template <>
bool
CellAccessor<3, 3>::point_inside(const Point<3> &) const;
template <>
bool
CellAccessor<1, 2>::point_inside(const Point<2> &) const;
template <>
bool
CellAccessor<1, 3>::point_inside(const Point<3> &) const;
template <>
bool
CellAccessor<2, 3>::point_inside(const Point<3> &) const;
// -------------------------------------------------------------------

template <>
void
TriaAccessor<3, 3, 3>::set_all_manifold_ids(const types::manifold_id) const;



namespace internal
{
  namespace TriaAccessorImplementation
  {
    /**
     * Compute the diameter for a given set of vertices. The vertices are
     * interpreted, depending on their count, as the vertices of a particular
     * reference cell.
     */
    template <int dim, int spacedim>
    inline double
    diameter(
      const boost::container::small_vector<Point<spacedim>,
                                           GeometryInfo<dim>::vertices_per_cell>
        vertices)
    {
      const ReferenceCell reference_cell =
        ReferenceCell::n_vertices_to_type(dim, vertices.size());

      if (reference_cell == ReferenceCells::Line)
        // Return the distance between the two vertices
        return (vertices[1] - vertices[0]).norm();
      else if (reference_cell == ReferenceCells::Triangle)
        // Return the longest of the three edges
        return std::max({(vertices[1] - vertices[0]).norm(),
                         (vertices[2] - vertices[1]).norm(),
                         (vertices[2] - vertices[0]).norm()});
      else if (reference_cell == ReferenceCells::Quadrilateral)
        // Return the longer one of the two diagonals of the quadrilateral
        return std::max({(vertices[3] - vertices[0]).norm(),
                         (vertices[2] - vertices[1]).norm()});
      else if (reference_cell == ReferenceCells::Tetrahedron)
        // Return the longest of the six edges of the tetrahedron
        return std::max({(vertices[1] - vertices[0]).norm(),
                         (vertices[2] - vertices[0]).norm(),
                         (vertices[2] - vertices[1]).norm(),
                         (vertices[3] - vertices[0]).norm(),
                         (vertices[3] - vertices[1]).norm(),
                         (vertices[3] - vertices[2]).norm()});
      else if (reference_cell == ReferenceCells::Pyramid)
        // Return ...
        return std::max({// the longest diagonal of the quadrilateral base
                         // of the pyramid or ...
                         (vertices[3] - vertices[0]).norm(),
                         (vertices[2] - vertices[1]).norm(),
                         // the longest edge connected with the apex of the
                         // pyramid
                         (vertices[4] - vertices[0]).norm(),
                         (vertices[4] - vertices[1]).norm(),
                         (vertices[4] - vertices[2]).norm(),
                         (vertices[4] - vertices[3]).norm()});
      else if (reference_cell == ReferenceCells::Wedge)
        // Return ...
        return std::max({// the longest of the 2*3=6 diagonals of the three
                         // quadrilateral sides of the wedge or ...
                         (vertices[4] - vertices[0]).norm(),
                         (vertices[3] - vertices[1]).norm(),
                         (vertices[5] - vertices[1]).norm(),
                         (vertices[4] - vertices[2]).norm(),
                         (vertices[5] - vertices[0]).norm(),
                         (vertices[3] - vertices[2]).norm(),
                         // the longest of the 3*2=6 edges of the two
                         // triangular faces of the wedge
                         (vertices[1] - vertices[0]).norm(),
                         (vertices[2] - vertices[1]).norm(),
                         (vertices[2] - vertices[0]).norm(),
                         (vertices[4] - vertices[3]).norm(),
                         (vertices[5] - vertices[4]).norm(),
                         (vertices[5] - vertices[3]).norm()});
      else if (reference_cell == ReferenceCells::Hexahedron)
        // Return the longest of the four diagonals of the hexahedron
        return std::max({(vertices[7] - vertices[0]).norm(),
                         (vertices[6] - vertices[1]).norm(),
                         (vertices[2] - vertices[5]).norm(),
                         (vertices[3] - vertices[4]).norm()});

      DEAL_II_NOT_IMPLEMENTED();
      return -1e10;
    }
  } // namespace TriaAccessorImplementation
} // namespace internal


/*--------------------- Functions: TriaAccessorBase -------------------------*/

template <int structdim, int dim, int spacedim>
inline TriaAccessorBase<structdim, dim, spacedim>::TriaAccessorBase(
  const Triangulation<dim, spacedim> *tria,
  const int                           level,
  const int                           index,
  const AccessorData *)
  : present_level((structdim == dim) ? level : 0)
  , present_index(index)
  , tria(tria)
{
  // non-cells have no level, so a 0
  // should have been passed, or a -1
  // for an end-iterator, or -2 for
  // an invalid (default constructed)
  // iterator
  if (structdim != dim)
    {
      Assert((level == 0) || (level == -1) || (level == -2),
             ExcInternalError());
    }
}


template <int structdim, int dim, int spacedim>
inline TriaAccessorBase<structdim, dim, spacedim>::TriaAccessorBase(
  const TriaAccessorBase<structdim, dim, spacedim> &a)
  : present_level(a.present_level)
  , present_index(a.present_index)
  , tria(a.tria)
{}


template <int structdim, int dim, int spacedim>
inline void
TriaAccessorBase<structdim, dim, spacedim>::copy_from(
  const TriaAccessorBase<structdim, dim, spacedim> &a)
{
  present_level = a.present_level;
  present_index = a.present_index;
  tria          = a.tria;

  if (structdim != dim)
    {
      Assert((present_level == 0) || (present_level == -1) ||
               (present_level == -2),
             ExcInternalError());
    }
}



template <int structdim, int dim, int spacedim>
inline TriaAccessorBase<structdim, dim, spacedim> &
TriaAccessorBase<structdim, dim, spacedim>::operator=(
  const TriaAccessorBase<structdim, dim, spacedim> &a)
{
  present_level = a.present_level;
  present_index = a.present_index;
  tria          = a.tria;

  if (structdim != dim)
    {
      Assert((present_level == 0) || (present_level == -1) ||
               (present_level == -2),
             ExcInternalError());
    }
  return *this;
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessorBase<structdim, dim, spacedim>::operator==(
  const TriaAccessorBase<structdim, dim, spacedim> &a) const
{
  Assert(tria == a.tria || tria == nullptr || a.tria == nullptr,
         TriaAccessorExceptions::ExcCantCompareIterators());
  return ((tria == a.tria) && (present_level == a.present_level) &&
          (present_index == a.present_index));
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessorBase<structdim, dim, spacedim>::operator!=(
  const TriaAccessorBase<structdim, dim, spacedim> &a) const
{
  Assert(tria == a.tria || tria == nullptr || a.tria == nullptr,
         TriaAccessorExceptions::ExcCantCompareIterators());
  return ((tria != a.tria) || (present_level != a.present_level) ||
          (present_index != a.present_index));
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessorBase<structdim, dim, spacedim>::operator<(
  const TriaAccessorBase<structdim, dim, spacedim> &other) const
{
  Assert(tria == other.tria, TriaAccessorExceptions::ExcCantCompareIterators());

  if (present_level != other.present_level)
    return (present_level < other.present_level);

  return (present_index < other.present_index);
}



template <int structdim, int dim, int spacedim>
inline int
TriaAccessorBase<structdim, dim, spacedim>::level() const
{
  // This is always zero or invalid
  // if the object is not a cell
  return present_level;
}



template <int structdim, int dim, int spacedim>
inline int
TriaAccessorBase<structdim, dim, spacedim>::index() const
{
  return present_index;
}



template <int structdim, int dim, int spacedim>
inline IteratorState::IteratorStates
TriaAccessorBase<structdim, dim, spacedim>::state() const
{
  if ((present_level >= 0) && (present_index >= 0))
    return IteratorState::valid;
  else if (present_index == -1)
    return IteratorState::past_the_end;
  else
    return IteratorState::invalid;
}



template <int structdim, int dim, int spacedim>
inline const Triangulation<dim, spacedim> &
TriaAccessorBase<structdim, dim, spacedim>::get_triangulation() const
{
  return *tria;
}



template <int structdim, int dim, int spacedim>
inline void
TriaAccessorBase<structdim, dim, spacedim>::operator++()
{
  // this iterator is used for
  // objects without level
  ++this->present_index;

  if (structdim != dim)
    {
      // is index still in the range of
      // the vector? (note that we don't
      // have to set the level, since
      // dim!=1 and the object therefore
      // has no level)
      if (this->present_index >= static_cast<int>(objects().n_objects()))
        this->present_index = -1;
    }
  else
    {
      while (this->present_index >=
             static_cast<int>(
               this->tria->levels[this->present_level]->cells.n_objects()))
        {
          // no -> go one level up until we find
          // one with more than zero cells
          ++this->present_level;
          this->present_index = 0;
          // highest level reached?
          if (this->present_level >=
              static_cast<int>(this->tria->levels.size()))
            {
              // return with past the end pointer
              this->present_level = this->present_index = -1;
              return;
            }
        }
    }
}


template <int structdim, int dim, int spacedim>
inline void
TriaAccessorBase<structdim, dim, spacedim>::operator--()
{
  // same as operator++
  --this->present_index;

  if (structdim != dim)
    {
      if (this->present_index < 0)
        this->present_index = -1;
    }
  else
    {
      while (this->present_index < 0)
        {
          // no -> go one level down
          --this->present_level;
          // lowest level reached?
          if (this->present_level == -1)
            {
              // return with past the end pointer
              this->present_level = this->present_index = -1;
              return;
            }
          // else
          this->present_index =
            this->tria->levels[this->present_level]->cells.n_objects() - 1;
        }
    }
}



template <int structdim, int dim, int spacedim>
inline dealii::internal::TriangulationImplementation::TriaObjects &
TriaAccessorBase<structdim, dim, spacedim>::objects() const
{
  if (structdim == dim)
    return this->tria->levels[this->present_level]->cells;

  if (structdim == 1 && dim > 1)
    return this->tria->faces->lines;

  if (structdim == 2 && dim > 2)
    return this->tria->faces->quads;

  DEAL_II_ASSERT_UNREACHABLE();

  return this->tria->levels[this->present_level]->cells;
}



/*---------------------- Functions: InvalidAccessor -------------------------*/

template <int structdim, int dim, int spacedim>
InvalidAccessor<structdim, dim, spacedim>::InvalidAccessor(const void *,
                                                           const int,
                                                           const int,
                                                           const AccessorData *)
{
  Assert(false,
         ExcMessage("You are attempting an invalid conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
InvalidAccessor<structdim, dim, spacedim>::InvalidAccessor(
  const InvalidAccessor &)
{
  Assert(false,
         ExcMessage("You are attempting an invalid conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
void
InvalidAccessor<structdim, dim, spacedim>::copy_from(const InvalidAccessor &)
{
  // nothing to do here. we could
  // throw an exception but we can't
  // get here without first creating
  // an object which would have
  // already thrown
}



template <int structdim, int dim, int spacedim>
bool
InvalidAccessor<structdim, dim, spacedim>::operator==(
  const InvalidAccessor &) const
{
  // nothing to do here. we could
  // throw an exception but we can't
  // get here without first creating
  // an object which would have
  // already thrown
  return false;
}



template <int structdim, int dim, int spacedim>
bool
InvalidAccessor<structdim, dim, spacedim>::operator!=(
  const InvalidAccessor &) const
{
  // nothing to do here. we could
  // throw an exception but we can't
  // get here without first creating
  // an object which would have
  // already thrown
  return true;
}



template <int structdim, int dim, int spacedim>
bool
InvalidAccessor<structdim, dim, spacedim>::used() const
{
  // nothing to do here. we could
  // throw an exception but we can't
  // get here without first creating
  // an object which would have
  // already thrown
  return false;
}



template <int structdim, int dim, int spacedim>
bool
InvalidAccessor<structdim, dim, spacedim>::has_children() const
{
  // nothing to do here. we could
  // throw an exception but we can't
  // get here without first creating
  // an object which would have
  // already thrown
  return false;
}



template <int structdim, int dim, int spacedim>
void
InvalidAccessor<structdim, dim, spacedim>::operator++() const
{}



template <int structdim, int dim, int spacedim>
void
InvalidAccessor<structdim, dim, spacedim>::operator--() const
{}



template <int structdim, int dim, int spacedim>
types::manifold_id
InvalidAccessor<structdim, dim, spacedim>::manifold_id() const
{
  return numbers::flat_manifold_id;
}



template <int structdim, int dim, int spacedim>
unsigned int
InvalidAccessor<structdim, dim, spacedim>::user_index() const
{
  return numbers::invalid_unsigned_int;
}



template <int structdim, int dim, int spacedim>
void
InvalidAccessor<structdim, dim, spacedim>::set_user_index(
  const unsigned int) const
{
  Assert(false,
         ExcMessage("You are trying to set the user index of an "
                    "invalid object."));
}



template <int structdim, int dim, int spacedim>
void
InvalidAccessor<structdim, dim, spacedim>::set_manifold_id(
  const types::manifold_id) const
{
  Assert(false,
         ExcMessage("You are trying to set the manifold id of an "
                    "invalid object."));
}



template <int structdim, int dim, int spacedim>
inline Point<spacedim> &
InvalidAccessor<structdim, dim, spacedim>::vertex(const unsigned int) const
{
  // nothing to do here. we could throw an exception but we can't get here
  // without first creating an object which would have already thrown
  static Point<spacedim> invalid_vertex;
  return invalid_vertex;
}


template <int structdim, int dim, int spacedim>
inline void *
InvalidAccessor<structdim, dim, spacedim>::line(const unsigned int) const
{
  // nothing to do here. we could throw an exception but we can't get here
  // without first creating an object which would have already thrown
  return nullptr;
}



template <int structdim, int dim, int spacedim>
inline void *
InvalidAccessor<structdim, dim, spacedim>::quad(const unsigned int) const
{
  // nothing to do here. we could throw an exception but we can't get here
  // without first creating an object which would have already thrown
  return nullptr;
}


/*------------------------ Functions: TriaAccessor ---------------------------*/


namespace internal
{
  namespace TriaAccessorImplementation
  {
    // make sure that if in the following we
    // write TriaAccessor
    // we mean the *class*
    // dealii::TriaAccessor, not the
    // enclosing namespace
    // dealii::internal::TriaAccessor
    using dealii::TriaAccessor;

    /**
     * A class with the same purpose as the similarly named class of the
     * Triangulation class. See there for more information.
     */
    struct Implementation
    {
      /**
       * Implementation of the function of some name in the parent class.
       */
      template <int structdim, int dim, int spacedim>
      inline static void
      set_combined_face_orientation(
        const TriaAccessor<structdim, dim, spacedim> &accessor,
        const unsigned int                            face_no,
        const types::geometric_orientation            combined_orientation)
      {
        Assert(structdim == dim,
               ExcMessage("This function can only be used on objects that are "
                          "cells and not on objects which bound cells."));
        AssertIndexRange(face_no, accessor.n_faces());
        AssertIndexRange(combined_orientation,
                         accessor.reference_cell().n_face_orientations(
                           face_no));

        // face_orientations is not set up in 1d
        if (dim != 1)
          accessor.tria->levels[accessor.present_level]
            ->face_orientations.set_combined_orientation(
              accessor.present_index * ReferenceCells::max_n_faces<dim>() +
                face_no,
              combined_orientation);
      }



      template <int dim, int spacedim>
      static std::array<unsigned int, 1>
      get_line_indices_of_cell(const TriaAccessor<1, dim, spacedim> &)
      {
        DEAL_II_ASSERT_UNREACHABLE();
        return {};
      }



      template <int structdim, int dim, int spacedim>
      static std::array<unsigned int, 4>
      get_line_indices_of_cell(const TriaAccessor<2, dim, spacedim> &cell)
      {
        // For 2d cells the access cell->line_orientation() is already
        // efficient
        std::array<unsigned int, 4> line_indices = {};
        for (const unsigned int line : cell.line_indices())
          line_indices[line] = cell.line_index(line);
        return line_indices;
      }

      /**
       * A helper function to provide faster access to cell->line_index() in
       * 3d
       */
      template <int structdim, int dim, int spacedim>
      static std::array<unsigned int, 12>
      get_line_indices_of_cell(
        const TriaAccessor<structdim, dim, spacedim> &cell)
      {
        std::array<unsigned int, 12> line_indices = {};

        // For hexahedra, the classical access via quads -> lines is too
        // inefficient. Unroll this code here to allow the compiler to inline
        // the necessary functions.
        const auto ref_cell = cell.reference_cell();
        if (ref_cell == ReferenceCells::Hexahedron)
          {
            for (unsigned int f = 4; f < 6; ++f)
              {
                const auto orientation =
                  cell.get_triangulation()
                    .levels[cell.level()]
                    ->face_orientations.get_combined_orientation(
                      cell.index() * ReferenceCells::max_n_faces<dim>() + f);

                // It might seem superfluous to spell out the four indices
                // that get later consumed by a for loop over these four
                // elements; however, for the compiler it is easier to inline
                // the statement of standard_to_real_face_line() when next to
                // each other, as opposed to be interleaved with a
                // line_index() call.
                const std::array<unsigned int, 4> my_indices{
                  {ref_cell.standard_to_real_face_line(0, f, orientation),
                   ref_cell.standard_to_real_face_line(1, f, orientation),
                   ref_cell.standard_to_real_face_line(2, f, orientation),
                   ref_cell.standard_to_real_face_line(3, f, orientation)}};
                const auto quad = cell.quad(f);
                for (unsigned int l = 0; l < 4; ++l)
                  line_indices[4 * (f - 4) + l] =
                    quad->line_index(my_indices[l]);
              }
            for (unsigned int f = 0; f < 2; ++f)
              {
                const auto orientation =
                  cell.get_triangulation()
                    .levels[cell.level()]
                    ->face_orientations.get_combined_orientation(
                      cell.index() * ReferenceCells::max_n_faces<dim>() + f);
                const std::array<unsigned int, 2> my_indices{
                  {ref_cell.standard_to_real_face_line(0, f, orientation),
                   ref_cell.standard_to_real_face_line(1, f, orientation)}};
                const auto quad      = cell.quad(f);
                line_indices[8 + f]  = quad->line_index(my_indices[0]);
                line_indices[10 + f] = quad->line_index(my_indices[1]);
              }
          }
        else if (ref_cell == ReferenceCells::Tetrahedron)
          {
            std::array<unsigned int, 3> orientations{
              {cell.combined_face_orientation(0),
               cell.combined_face_orientation(1),
               cell.combined_face_orientation(2)}};
            const std::array<unsigned int, 6> my_indices{
              {ref_cell.standard_to_real_face_line(0, 0, orientations[0]),
               ref_cell.standard_to_real_face_line(1, 0, orientations[0]),
               ref_cell.standard_to_real_face_line(2, 0, orientations[0]),
               ref_cell.standard_to_real_face_line(1, 1, orientations[1]),
               ref_cell.standard_to_real_face_line(2, 1, orientations[1]),
               ref_cell.standard_to_real_face_line(1, 2, orientations[2])}};
            line_indices[0] = cell.quad(0)->line_index(my_indices[0]);
            line_indices[1] = cell.quad(0)->line_index(my_indices[1]);
            line_indices[2] = cell.quad(0)->line_index(my_indices[2]);
            line_indices[3] = cell.quad(1)->line_index(my_indices[3]);
            line_indices[4] = cell.quad(1)->line_index(my_indices[4]);
            line_indices[5] = cell.quad(2)->line_index(my_indices[5]);
          }
        else
          // For other shapes (wedges, pyramids), we do not currently
          // implement an optimized function.
          for (unsigned int l = 0; l < std::min(12U, cell.n_lines()); ++l)
            line_indices[l] = cell.line_index(l);

        return line_indices;
      }



      /**
       * A helper function to provide faster access to
       * cell->line_orientation(), 1d specialization
       */
      template <int dim, int spacedim>
      static std::array<types::geometric_orientation, 1>
      get_line_orientations_of_cell(const TriaAccessor<1, dim, spacedim> &)
      {
        DEAL_II_ASSERT_UNREACHABLE();
        return {};
      }



      /**
       * A helper function to provide faster access to
       * cell->line_orientation(), 2d specialization
       */
      template <int dim, int spacedim>
      static std::array<types::geometric_orientation, 4>
      get_line_orientations_of_cell(const TriaAccessor<2, dim, spacedim> &cell)
      {
        // For 2d cells the access cell->line_orientation() is already
        // efficient
        std::array<types::geometric_orientation, 4> line_orientations = {};
        for (const unsigned int line : cell.line_indices())
          line_orientations[line] = cell.line_orientation(line);
        return line_orientations;
      }



      /**
       * A helper function to provide faster access to
       * cell->line_orientation(), 3d specialization
       */
      template <int dim, int spacedim>
      static std::array<types::geometric_orientation, 12>
      get_line_orientations_of_cell(const TriaAccessor<3, dim, spacedim> &cell)
      {
        std::array<types::geometric_orientation, 12> line_orientations = {};

        // For hexahedra, the classical access via quads -> lines is too
        // inefficient. Unroll this code here to allow the compiler to inline
        // the necessary functions.
        const auto ref_cell = cell.reference_cell();
        if (ref_cell == ReferenceCells::Hexahedron)
          {
            for (unsigned int f = 4; f < 6; ++f)
              {
                const auto orientation =
                  cell.get_triangulation()
                    .levels[cell.level()]
                    ->face_orientations.get_combined_orientation(
                      cell.index() * ReferenceCells::max_n_faces<dim>() + f);

                // It might seem superfluous to spell out the four indices and
                // orientations that get later consumed by a for loop over
                // these four elements; however, for the compiler it is easier
                // to inline the statement of standard_to_real_face_line()
                // when next to each other, as opposed to be interleaved with
                // a line_index() call.
                const std::array<unsigned int, 4> my_indices{
                  {ref_cell.standard_to_real_face_line(0, f, orientation),
                   ref_cell.standard_to_real_face_line(1, f, orientation),
                   ref_cell.standard_to_real_face_line(2, f, orientation),
                   ref_cell.standard_to_real_face_line(3, f, orientation)}};
                const auto quad = cell.quad(f);
                const std::array<types::geometric_orientation, 4>
                  my_orientations{{ref_cell.face_to_cell_line_orientation(
                                     0,
                                     f,
                                     orientation,
                                     quad->line_orientation(my_indices[0])),
                                   ref_cell.face_to_cell_line_orientation(
                                     1,
                                     f,
                                     orientation,
                                     quad->line_orientation(my_indices[1])),
                                   ref_cell.face_to_cell_line_orientation(
                                     2,
                                     f,
                                     orientation,
                                     quad->line_orientation(my_indices[2])),
                                   ref_cell.face_to_cell_line_orientation(
                                     3,
                                     f,
                                     orientation,
                                     quad->line_orientation(my_indices[3]))}};
                for (unsigned int l = 0; l < 4; ++l)
                  line_orientations[4 * (f - 4) + l] = my_orientations[l];
              }
            for (unsigned int f = 0; f < 2; ++f)
              {
                const auto orientation =
                  cell.get_triangulation()
                    .levels[cell.level()]
                    ->face_orientations.get_combined_orientation(
                      cell.index() * ReferenceCells::max_n_faces<3>() + f);
                const std::array<unsigned int, 2> my_indices{
                  {ref_cell.standard_to_real_face_line(0, f, orientation),
                   ref_cell.standard_to_real_face_line(1, f, orientation)}};
                const auto quad = cell.quad(f);
                const std::array<types::geometric_orientation, 2>
                  my_orientations{{ref_cell.face_to_cell_line_orientation(
                                     0,
                                     f,
                                     orientation,
                                     quad->line_orientation(my_indices[0])),
                                   ref_cell.face_to_cell_line_orientation(
                                     1,
                                     f,
                                     orientation,
                                     quad->line_orientation(my_indices[1]))}};
                line_orientations[8 + f]  = my_orientations[0];
                line_orientations[10 + f] = my_orientations[1];
              }
          }
        else if (ref_cell == ReferenceCells::Tetrahedron)
          {
            std::array<unsigned int, 3> orientations{
              {cell.combined_face_orientation(0),
               cell.combined_face_orientation(1),
               cell.combined_face_orientation(2)}};
            const std::array<unsigned int, 6> my_indices{
              {ref_cell.standard_to_real_face_line(0, 0, orientations[0]),
               ref_cell.standard_to_real_face_line(1, 0, orientations[0]),
               ref_cell.standard_to_real_face_line(2, 0, orientations[0]),
               ref_cell.standard_to_real_face_line(1, 1, orientations[1]),
               ref_cell.standard_to_real_face_line(2, 1, orientations[1]),
               ref_cell.standard_to_real_face_line(1, 2, orientations[2])}};
            line_orientations[0] = ref_cell.face_to_cell_line_orientation(
              0,
              0,
              orientations[0],
              cell.quad(0)->line_orientation(my_indices[0]));
            line_orientations[1] = ref_cell.face_to_cell_line_orientation(
              1,
              0,
              orientations[0],
              cell.quad(0)->line_orientation(my_indices[1]));
            line_orientations[2] = ref_cell.face_to_cell_line_orientation(
              2,
              0,
              orientations[0],
              cell.quad(0)->line_orientation(my_indices[2]));
            line_orientations[3] = ref_cell.face_to_cell_line_orientation(
              1,
              1,
              orientations[1],
              cell.quad(1)->line_orientation(my_indices[3]));
            line_orientations[4] = ref_cell.face_to_cell_line_orientation(
              2,
              1,
              orientations[1],
              cell.quad(1)->line_orientation(my_indices[4]));
            line_orientations[5] = ref_cell.face_to_cell_line_orientation(
              1,
              2,
              orientations[2],
              cell.quad(2)->line_orientation(my_indices[5]));
          }
        else
          // For other shapes (wedges, pyramids), we do not currently
          // implement an optimized function
          for (unsigned int l = 0; l < std::min(12U, cell.n_lines()); ++l)
            line_orientations[l] = cell.line_orientation(l);

        return line_orientations;
      }
    };
  } // namespace TriaAccessorImplementation
} // namespace internal



template <int structdim, int dim, int spacedim>
inline TriaAccessor<structdim, dim, spacedim>::TriaAccessor(
  const Triangulation<dim, spacedim> *parent,
  const int                           level,
  const int                           index,
  const AccessorData                 *local_data)
  : TriaAccessorBase<structdim, dim, spacedim>(parent, level, index, local_data)
{}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::used() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));
  return this->objects().used[this->present_index];
}



template <int structdim, int dim, int spacedim>
inline TriaIterator<TriaAccessor<0, dim, spacedim>>
TriaAccessor<structdim, dim, spacedim>::vertex_iterator(
  const unsigned int i) const
{
  return TriaIterator<TriaAccessor<0, dim, spacedim>>(this->tria,
                                                      0,
                                                      vertex_index(i));
}



template <int structdim, int dim, int spacedim>
inline ReferenceCell
TriaAccessor<structdim, dim, spacedim>::reference_cell() const
{
  if (structdim == 0)
    return ReferenceCells::Vertex;
  else if (structdim == 1)
    return ReferenceCells::Line;
  else if (structdim == dim)
    return this->tria->levels[this->present_level]
      ->reference_cell[this->present_index];
  else
    return this->tria->faces->get_quad_type(this->present_index);
}



template <int structdim, int dim, int spacedim>
inline unsigned int
TriaAccessor<structdim, dim, spacedim>::vertex_index(
  const unsigned int corner) const
{
  AssertIndexRange(corner, this->n_vertices());

  if constexpr (structdim == 1)
    {
      // This branch needs to be first (and not combined with the structdim ==
      // dim branch) so that we can get line vertex indices when setting up the
      // cell vertex index cache
      return this->objects()
        .cells[this->present_index * ReferenceCells::max_n_faces<1>() + corner];
    }
  else if constexpr (structdim == dim)
    {
      // This branch should only be used after the cell vertex index cache is
      // set up
      const auto my_index = static_cast<std::size_t>(this->present_index) *
                            ReferenceCells::max_n_vertices<dim>();
      AssertIndexRange(my_index + corner,
                       this->tria->levels[this->present_level]
                         ->cell_vertex_indices_cache.size());
      const unsigned int vertex_index =
        this->tria->levels[this->present_level]
          ->cell_vertex_indices_cache[my_index + corner];
      Assert(vertex_index != numbers::invalid_unsigned_int, ExcInternalError());
      return vertex_index;
    }
  else if constexpr (structdim == 2)
    {
      const auto [line_index, vertex_index] =
        this->reference_cell().standard_vertex_to_face_and_vertex_index(corner);
      const auto vertex_within_line_index =
        this->reference_cell().standard_to_real_face_vertex(
          vertex_index, line_index, this->line_orientation(line_index));

      return this->line(line_index)->vertex_index(vertex_within_line_index);
    }
  else
    {
      DEAL_II_ASSERT_UNREACHABLE();
      return numbers::invalid_unsigned_int;
    }
}



template <int structdim, int dim, int spacedim>
inline Point<spacedim> &
TriaAccessor<structdim, dim, spacedim>::vertex(const unsigned int i) const
{
  return const_cast<Point<spacedim> &>(this->tria->vertices[vertex_index(i)]);
}



template <int structdim, int dim, int spacedim>
inline typename dealii::internal::TriangulationImplementation::
  Iterators<dim, spacedim>::line_iterator
  TriaAccessor<structdim, dim, spacedim>::line(const unsigned int i) const
{
  // checks happen in line_index
  return typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::line_iterator(this->tria, 0, line_index(i));
}



template <int structdim, int dim, int spacedim>
inline unsigned int
TriaAccessor<structdim, dim, spacedim>::line_index(const unsigned int i) const
{
  (void)i;
  AssertIndexRange(i, this->n_lines());
  Assert(structdim != 1,
         ExcMessage("You can't ask for the index of a line bounding a "
                    "one-dimensional cell because it is not bounded by "
                    "lines."));

  if constexpr (structdim == 2)
    {
      return this->objects()
        .cells[this->present_index * ReferenceCells::max_n_faces<2>() + i];
    }
  else if constexpr (structdim == 3)
    {
      const auto [face_index, line_index] =
        this->reference_cell().standard_line_to_face_and_line_index(i);
      const auto line_within_face_index =
        this->reference_cell().standard_to_real_face_line(
          line_index, face_index, this->combined_face_orientation(face_index));

      return this->quad(face_index)->line_index(line_within_face_index);
    }

  DEAL_II_ASSERT_UNREACHABLE();
  return numbers::invalid_unsigned_int;
}



template <int structdim, int dim, int spacedim>
inline typename dealii::internal::TriangulationImplementation::
  Iterators<dim, spacedim>::quad_iterator
  TriaAccessor<structdim, dim, spacedim>::quad(const unsigned int i) const
{
  // checks happen in quad_index
  return typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::quad_iterator(this->tria, 0, quad_index(i));
}



template <int structdim, int dim, int spacedim>
inline unsigned int
TriaAccessor<structdim, dim, spacedim>::quad_index(const unsigned int i) const
{
  Assert(structdim == 3,
         ExcMessage("You can't ask for the index of a quad bounding "
                    "a one- or two-dimensional cell because it is not "
                    "bounded by quads."));
  // work around a bogus GCC-9 warning which considers i unused except in 3d
  (void)i;
  if constexpr (structdim == 3)
    {
      AssertIndexRange(i, n_faces());
      return this->tria->levels[this->present_level]
        ->cells
        .cells[this->present_index * ReferenceCells::max_n_faces<3>() + i];
    }
  else
    return numbers::invalid_unsigned_int;
}



template <int structdim, int dim, int spacedim>
inline types::geometric_orientation
TriaAccessor<structdim, dim, spacedim>::combined_face_orientation(
  const unsigned int face) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());
  AssertIndexRange(face, n_faces());
  Assert(structdim == dim,
         ExcMessage("This function can only be used on objects "
                    "that are cells, but not on faces or edges "
                    "that bound cells."));
  // work around a bogus GCC-9 warning which considers face unused except in 3d
  (void)face;

  if constexpr (structdim == 1)
    return numbers::default_geometric_orientation;
  else if constexpr (structdim == 2)
    {
      // if all elements are quads (or if we have a very special consistently
      // oriented triangular mesh) then we do not store this array
      if (this->tria->levels[this->present_level]
            ->face_orientations.n_objects() == 0)
        return numbers::default_geometric_orientation;
      else
        return this->tria->levels[this->present_level]
          ->face_orientations.get_combined_orientation(
            this->present_index * ReferenceCells::max_n_faces<dim>() + face);
    }
  else
    return this->tria->levels[this->present_level]
      ->face_orientations.get_combined_orientation(
        this->present_index * ReferenceCells::max_n_faces<dim>() + face);
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::face_orientation(
  const unsigned int face) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());
  AssertIndexRange(face, n_faces());
  Assert(structdim == dim,
         ExcMessage("This function can only be used on objects "
                    "that are cells, but not on faces or edges "
                    "that bound cells."));
  // work around a bogus GCC-9 warning which considers face unused in 1d
  (void)face;

  if constexpr (structdim == 1)
    // in 1d 'faces' are vertices and those are always consistently oriented
    return true;
  else if constexpr (structdim == 2)
    return this->line_orientation(face) ==
           numbers::default_geometric_orientation;
  else
    return this->tria->levels[this->present_level]
      ->face_orientations.get_orientation(
        this->present_index * ReferenceCells::max_n_faces<structdim>() + face);
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::face_flip(const unsigned int face) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert(structdim == dim,
         ExcMessage("This function can only be used on objects "
                    "that are cells, but not on faces or edges "
                    "that bound cells."));
  AssertIndexRange(face, n_faces());
  // work around a bogus GCC-9 warning which considers face unused except in 3d
  (void)face;

  if constexpr (structdim == 3)
    return this->tria->levels[this->present_level]->face_orientations.get_flip(
      this->present_index * ReferenceCells::max_n_faces<structdim>() + face);
  else
    // In 1d and 2d, face_flip is always false as faces can only be
    // 'flipped' in 3d.
    return false;
}


template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::face_rotation(
  const unsigned int face) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert(structdim == dim,
         ExcMessage("This function can only be used on objects "
                    "that are cells, but not on faces or edges "
                    "that bound cells."));
  AssertIndexRange(face, n_faces());
  // work around a bogus GCC-9 warning which considers face unused except in 3d
  (void)face;

  if constexpr (structdim == 3)
    return this->tria->levels[this->present_level]
      ->face_orientations.get_rotation(
        this->present_index * ReferenceCells::max_n_faces<structdim>() + face);
  else
    // In 1d and 2d, face_rotation is always false as faces can only be
    // 'rotated' in 3d.
    return false;
}



template <int structdim, int dim, int spacedim>
inline types::geometric_orientation
TriaAccessor<structdim, dim, spacedim>::line_orientation(
  const unsigned int line) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());
  AssertIndexRange(line, this->n_lines());
  // work around a bogus GCC-9 warning which considers line unused in 1d
  (void)line;

  if constexpr (structdim == 1)
    return numbers::default_geometric_orientation;
  else if constexpr (structdim == 2 && dim == 2)
    // lines in 2d are faces
    {
      const auto combined_orientation = combined_face_orientation(line);
      Assert(combined_orientation == numbers::default_geometric_orientation ||
               combined_orientation == numbers::reverse_line_orientation,
             ExcInternalError());
      return combined_orientation;
    }
  else if constexpr (structdim == 2 && dim == 3)
    {
      // line orientations in 3d are stored in their own array as bools: here
      // 'true' is the default orientation and 'false' is the reversed one
      // (which matches set_line_orientation())
      const auto index =
        this->present_index * ReferenceCells::max_n_lines<2>() + line;
      Assert(index < this->tria->faces->quads_line_orientations.size(),
             ExcInternalError());
      return this->tria->faces->quads_line_orientations[index] ?
               numbers::default_geometric_orientation :
               numbers::reverse_line_orientation;
    }
  else if constexpr (structdim == 3 && dim == 3)
    {
      const auto reference_cell = this->reference_cell();
      // First pick a face on which this line is a part of, and the
      // index of the line within.
      const auto [face_index, line_index] =
        reference_cell.standard_line_to_face_and_line_index(line);
      const auto line_within_face_index =
        reference_cell.standard_to_real_face_line(
          line_index, face_index, this->combined_face_orientation(face_index));

      // Then query how that line is oriented within that face:
      return reference_cell.face_to_cell_line_orientation(
        line_index,
        face_index,
        this->combined_face_orientation(face_index),
        this->quad(face_index)->line_orientation(line_within_face_index));
    }
  else
    {
      DEAL_II_ASSERT_UNREACHABLE();
      return false;
    }
}



template <int structdim, int dim, int spacedim>
inline void
TriaAccessor<structdim, dim, spacedim>::set_line_orientation(
  const unsigned int                 line,
  const types::geometric_orientation value) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());
  AssertIndexRange(line, this->n_lines());
  Assert(dim != 1,
         ExcMessage("In 1d lines are cells and thus do not need to have their "
                    "orientations set."));
  Assert(dim != 2,
         ExcMessage("In 2d lines are faces, and, for compatibility with other "
                    "dimensions, their orientations should be set via "
                    "set_combined_face_orientation()."));
  // work around a bogus GCC-9 warning which considers line and value unused
  // except in 3d
  (void)line;
  (void)value;

  if constexpr (dim == 3)
    {
      // We set line orientations per face, not per cell, so this only works for
      // faces in 3d.
      Assert(structdim == 2, ExcNotImplemented());
      const auto index =
        this->present_index * ReferenceCells::max_n_lines<2>() + line;
      Assert(index < this->tria->faces->quads_line_orientations.size(),
             ExcInternalError());
      this->tria->faces->quads_line_orientations[index] =
        value == numbers::default_geometric_orientation;
    }
}



template <int structdim, int dim, int spacedim>
inline void
TriaAccessor<structdim, dim, spacedim>::set_combined_face_orientation(
  const unsigned int                 face,
  const types::geometric_orientation combined_orientation) const
{
  Assert(used(), TriaAccessorExceptions::ExcCellNotUsed());
  AssertIndexRange(face, this->n_faces());

  dealii::internal::TriaAccessorImplementation::Implementation::
    set_combined_face_orientation(*this, face, combined_orientation);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_used_flag() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));
  this->objects().used[this->present_index] = true;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_used_flag() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));
  this->objects().used[this->present_index] = false;
}


template <int structdim, int dim, int spacedim>
int
TriaAccessor<structdim, dim, spacedim>::child_index(const unsigned int i) const
{
  Assert(has_children(), TriaAccessorExceptions::ExcCellHasNoChildren());
  AssertIndexRange(i, n_children());

  // each set of two children are stored
  // consecutively, so we only have to find
  // the location of the set of children
  const unsigned int n_sets_of_two =
    GeometryInfo<structdim>::max_children_per_cell / 2;
  return this->objects().children[n_sets_of_two * this->present_index + i / 2] +
         i % 2;
}



template <int structdim, int dim, int spacedim>
int
TriaAccessor<structdim, dim, spacedim>::isotropic_child_index(
  const unsigned int i) const
{
  AssertIndexRange(i, GeometryInfo<structdim>::max_children_per_cell);

  switch (structdim)
    {
      case 1:
        return child_index(i);
      case 2:
        {
          const RefinementCase<2> this_refinement_case(
            static_cast<std::uint8_t>(refinement_case()));

          Assert(this_refinement_case != RefinementCase<2>::no_refinement,
                 TriaAccessorExceptions::ExcCellHasNoChildren());

          if (this_refinement_case == RefinementCase<2>::cut_xy)
            return child_index(i);
          else if ((this_refinement_case == RefinementCase<2>::cut_x) &&
                   (child(i % 2)->refinement_case() ==
                    RefinementCase<2>::cut_y))
            return child(i % 2)->child_index(i / 2);
          else if ((this_refinement_case == RefinementCase<2>::cut_y) &&
                   (child(i / 2)->refinement_case() ==
                    RefinementCase<2>::cut_x))
            return child(i / 2)->child_index(i % 2);
          else
            Assert(
              false,
              ExcMessage(
                "This cell has no grandchildren equivalent to isotropic refinement"));
          break;
        }

      case 3:
        DEAL_II_NOT_IMPLEMENTED();
    }
  return -1;
}



template <int structdim, int dim, int spacedim>
RefinementCase<structdim>
TriaAccessor<structdim, dim, spacedim>::refinement_case() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));

  switch (structdim)
    {
      case 1:
        return (RefinementCase<structdim>(
          this->objects().children[this->present_index] != -1 ?
            // cast the branches
            // here first to uchar
            // and then (above) to
            // RefinementCase<structdim>
            // so that the
            // conversion is valid
            // even for the case
            // structdim>1 (for
            // which this part of
            // the code is dead
            // anyway)
            static_cast<std::uint8_t>(RefinementCase<1>::cut_x) :
            static_cast<std::uint8_t>(RefinementCase<1>::no_refinement)));

      default:
        Assert(static_cast<unsigned int>(this->present_index) <
                 this->objects().refinement_cases.size(),
               ExcIndexRange(this->present_index,
                             0,
                             this->objects().refinement_cases.size()));

        return (static_cast<RefinementCase<structdim>>(
          this->objects().refinement_cases[this->present_index]));
    }
}



template <int structdim, int dim, int spacedim>
inline TriaIterator<TriaAccessor<structdim, dim, spacedim>>
TriaAccessor<structdim, dim, spacedim>::child(const unsigned int i) const

{
  // checking of 'i' happens in child_index
  const TriaIterator<TriaAccessor<structdim, dim, spacedim>> q(
    this->tria, (dim == structdim ? this->level() + 1 : 0), child_index(i));

  Assert((q.state() == IteratorState::past_the_end) || q->used(),
         ExcInternalError());

  return q;
}



template <int structdim, int dim, int spacedim>
inline unsigned int
TriaAccessor<structdim, dim, spacedim>::child_iterator_to_index(
  const TriaIterator<TriaAccessor<structdim, dim, spacedim>> &child) const
{
  const auto n_children = this->n_children();
  for (unsigned int child_n = 0; child_n < n_children; ++child_n)
    if (this->child(child_n) == child)
      return child_n;

  Assert(false,
         ExcMessage("The given child is not a child of the current object."));
  return numbers::invalid_unsigned_int;
}



template <int structdim, int dim, int spacedim>
inline TriaIterator<TriaAccessor<structdim, dim, spacedim>>
TriaAccessor<structdim, dim, spacedim>::isotropic_child(
  const unsigned int i) const
{
  // checking of 'i' happens in child() or
  // child_index() called below
  switch (structdim)
    {
      case 1:
        // no anisotropic refinement in 1d
        return child(i);

      case 2:
        {
          const RefinementCase<2> this_refinement_case(
            static_cast<std::uint8_t>(refinement_case()));

          Assert(this_refinement_case != RefinementCase<2>::no_refinement,
                 TriaAccessorExceptions::ExcCellHasNoChildren());

          if (this_refinement_case == RefinementCase<2>::cut_xy)
            return child(i);
          else if ((this_refinement_case == RefinementCase<2>::cut_x) &&
                   (child(i % 2)->refinement_case() ==
                    RefinementCase<2>::cut_y))
            return child(i % 2)->child(i / 2);
          else if ((this_refinement_case == RefinementCase<2>::cut_y) &&
                   (child(i / 2)->refinement_case() ==
                    RefinementCase<2>::cut_x))
            return child(i / 2)->child(i % 2);
          else
            Assert(
              false,
              ExcMessage(
                "This cell has no grandchildren equivalent to isotropic refinement"));
          break;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
  // we don't get here but have to return
  // something...
  return child(0);
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::has_children() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));

  // each set of two children are stored
  // consecutively, so we only have to find
  // the location of the set of children
  const unsigned int n_sets_of_two =
    GeometryInfo<structdim>::max_children_per_cell / 2;
  return (this->objects().children[n_sets_of_two * this->present_index] != -1);
}



template <int structdim, int dim, int spacedim>
inline unsigned int
TriaAccessor<structdim, dim, spacedim>::n_children() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));
  if (reference_cell() == ReferenceCells::Tetrahedron)
    return GeometryInfo<structdim>::n_children(
      RefinementCase<structdim>::isotropic_refinement);
  else
    return GeometryInfo<structdim>::n_children(refinement_case());
}



template <int structdim, int dim, int spacedim>
inline void
TriaAccessor<structdim, dim, spacedim>::set_refinement_case(
  const RefinementCase<structdim> &refinement_case) const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));
  Assert(static_cast<unsigned int>(this->present_index) <
           this->objects().refinement_cases.size(),
         ExcIndexRange(this->present_index,
                       0,
                       this->objects().refinement_cases.size()));

  this->objects().refinement_cases[this->present_index] = refinement_case;
}


template <int structdim, int dim, int spacedim>
inline void
TriaAccessor<structdim, dim, spacedim>::clear_refinement_case() const
{
  Assert(this->state() == IteratorState::valid,
         TriaAccessorExceptions::ExcDereferenceInvalidObject<TriaAccessor>(
           *this));
  Assert(static_cast<unsigned int>(this->present_index) <
           this->objects().refinement_cases.size(),
         ExcIndexRange(this->present_index,
                       0,
                       this->objects().refinement_cases.size()));

  this->objects().refinement_cases[this->present_index] =
    RefinementCase<structdim>::no_refinement;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_children(const unsigned int i,
                                                     const int index) const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert(i % 2 == 0, TriaAccessorExceptions::ExcSetOnlyEvenChildren(i));

  // each set of two children are stored
  // consecutively, so we only have to find
  // the location of the set of children
  const unsigned int n_sets_of_two =
    GeometryInfo<structdim>::max_children_per_cell / 2;

  Assert(
    // clearing the child index for a cell
    (index == -1) ||
      // if setting the child index for the i'th child (with i==0),
      // then the index must be a non-negative number
      (i == 0 && !this->has_children() && (index >= 0)) ||
      // if setting the child index for the i'th child (with i>0),
      // then the previously stored index must be the invalid
      // index
      (i > 0 && this->has_children() && (index >= 0) &&
       this->objects().children[n_sets_of_two * this->present_index + i / 2] ==
         -1),
    TriaAccessorExceptions::ExcCantSetChildren(index));

  this->objects().children[n_sets_of_two * this->present_index + i / 2] = index;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_children() const
{
  // each set of two children are stored
  // consecutively, so we only have to find
  // the location of the set of children
  const unsigned int n_sets_of_two =
    GeometryInfo<structdim>::max_children_per_cell / 2;

  for (unsigned int i = 0; i < n_sets_of_two; ++i)
    set_children(2 * i, -1);
}



template <int structdim, int dim, int spacedim>
inline bool
TriaAccessor<structdim, dim, spacedim>::user_flag_set() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->objects().user_flags[this->present_index];
}



template <int structdim, int dim, int spacedim>
inline void
TriaAccessor<structdim, dim, spacedim>::set_user_flag() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_flags[this->present_index] = true;
}



template <int structdim, int dim, int spacedim>
inline void
TriaAccessor<structdim, dim, spacedim>::clear_user_flag() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_flags[this->present_index] = false;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::recursively_set_user_flag() const
{
  set_user_flag();

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->recursively_set_user_flag();
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::recursively_clear_user_flag() const
{
  clear_user_flag();

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->recursively_clear_user_flag();
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_user_data() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().clear_user_data(this->present_index);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_user_pointer(void *p) const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_pointer(this->present_index) = p;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_user_pointer() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_pointer(this->present_index) = nullptr;
}



template <int structdim, int dim, int spacedim>
void *
TriaAccessor<structdim, dim, spacedim>::user_pointer() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->objects().user_pointer(this->present_index);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::recursively_set_user_pointer(
  void *p) const
{
  set_user_pointer(p);

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->recursively_set_user_pointer(p);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::recursively_clear_user_pointer() const
{
  clear_user_pointer();

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->recursively_clear_user_pointer();
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_user_index(
  const unsigned int p) const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_index(this->present_index) = p;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::clear_user_index() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->objects().user_index(this->present_index) = 0;
}



template <int structdim, int dim, int spacedim>
unsigned int
TriaAccessor<structdim, dim, spacedim>::user_index() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->objects().user_index(this->present_index);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::recursively_set_user_index(
  const unsigned int p) const
{
  set_user_index(p);

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->recursively_set_user_index(p);
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::recursively_clear_user_index() const
{
  clear_user_index();

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->recursively_clear_user_index();
}



template <int structdim, int dim, int spacedim>
inline unsigned int
TriaAccessor<structdim, dim, spacedim>::max_refinement_depth() const
{
  if (!this->has_children())
    return 0;

  unsigned int max_depth = 1;
  for (unsigned int c = 0; c < n_children(); ++c)
    max_depth = std::max(max_depth, child(c)->max_refinement_depth() + 1);
  return max_depth;
}



template <int structdim, int dim, int spacedim>
unsigned int
TriaAccessor<structdim, dim, spacedim>::n_active_descendants() const
{
  if (!this->has_children())
    return 1;
  else
    {
      unsigned int sum = 0;
      for (unsigned int c = 0; c < n_children(); ++c)
        sum += this->child(c)->n_active_descendants();
      return sum;
    }
}



template <int structdim, int dim, int spacedim>
types::boundary_id
TriaAccessor<structdim, dim, spacedim>::boundary_id() const
{
  Assert(structdim < dim, ExcImpossibleInDim(dim));
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  return this->objects()
    .boundary_or_material_id[this->present_index]
    .boundary_id;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_boundary_id(
  const types::boundary_id boundary_ind) const
{
  Assert(structdim < dim, ExcImpossibleInDim(dim));
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert(boundary_ind != numbers::internal_face_boundary_id,
         ExcMessage("You are trying to set the boundary_id to an invalid "
                    "value (numbers::internal_face_boundary_id is reserved)."));
  Assert(this->at_boundary(),
         ExcMessage("You are trying to set the boundary_id of an "
                    "internal object, which is not allowed!"));

  this->objects().boundary_or_material_id[this->present_index].boundary_id =
    boundary_ind;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_boundary_id_internal(
  const types::boundary_id boundary_ind) const
{
  Assert(structdim < dim, ExcImpossibleInDim(dim));
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  this->objects().boundary_or_material_id[this->present_index].boundary_id =
    boundary_ind;
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_all_boundary_ids(
  const types::boundary_id boundary_ind) const
{
  set_boundary_id(boundary_ind);

  switch (structdim)
    {
      case 1:
        // 1d objects have no sub-objects
        // where we have to do anything
        break;

      case 2:
        // for boundary quads also set
        // boundary_id of bounding lines
        for (unsigned int i = 0; i < this->n_lines(); ++i)
          this->line(i)->set_boundary_id(boundary_ind);
        break;

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int structdim, int dim, int spacedim>
bool
TriaAccessor<structdim, dim, spacedim>::at_boundary() const
{
  // error checking is done
  // in boundary_id()
  return (boundary_id() != numbers::internal_face_boundary_id);
}



template <int structdim, int dim, int spacedim>
const Manifold<dim, spacedim> &
TriaAccessor<structdim, dim, spacedim>::get_manifold() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->get_manifold(this->manifold_id());
}


template <int structdim, int dim, int spacedim>
types::manifold_id
TriaAccessor<structdim, dim, spacedim>::manifold_id() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  return this->objects().manifold_id[this->present_index];
}



template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_manifold_id(
  const types::manifold_id manifold_ind) const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());

  this->objects().manifold_id[this->present_index] = manifold_ind;
}


template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::set_all_manifold_ids(
  const types::manifold_id manifold_ind) const
{
  set_manifold_id(manifold_ind);

  if (this->has_children())
    for (unsigned int c = 0; c < this->n_children(); ++c)
      this->child(c)->set_all_manifold_ids(manifold_ind);

  switch (structdim)
    {
      case 1:
        if (dim == 1)
          {
            (*this->tria->vertex_to_manifold_id_map_1d)[vertex_index(0)] =
              manifold_ind;
            (*this->tria->vertex_to_manifold_id_map_1d)[vertex_index(1)] =
              manifold_ind;
          }
        break;

      case 2:
        // for quads/simplices also set manifold_id of bounding lines
        for (unsigned int i = 0; i < this->n_lines(); ++i)
          this->line(i)->set_manifold_id(manifold_ind);
        break;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int structdim, int dim, int spacedim>
double
TriaAccessor<structdim, dim, spacedim>::diameter() const
{
  boost::container::small_vector<Point<spacedim>,
#  ifndef _MSC_VER
                                 ReferenceCells::max_n_vertices<structdim>()
#  else
                                 GeometryInfo<structdim>::vertices_per_cell
#  endif
                                 >
    vertices(this->n_vertices());

  for (unsigned int v = 0; v < vertices.size(); ++v)
    vertices[v] = this->vertex(v);

  return internal::TriaAccessorImplementation::diameter<structdim, spacedim>(
    vertices);
}



template <int dim, int spacedim>
double
CellAccessor<dim, spacedim>::diameter(
  const Mapping<dim, spacedim> &mapping) const
{
  return internal::TriaAccessorImplementation::diameter<dim, spacedim>(
    mapping.get_vertices(typename Triangulation<dim, spacedim>::cell_iterator(
      this->tria, this->level(), this->index())));
}



template <int structdim, int dim, int spacedim>
std::pair<Point<spacedim>, double>
TriaAccessor<structdim, dim, spacedim>::enclosing_ball() const
{
  // If the object is one dimensional,
  // the enclosing ball is the initial iterate
  // i.e., the ball's center and diameter are
  // the center and the diameter of the object.
  if (structdim == 1)
    return std::make_pair((this->vertex(1) + this->vertex(0)) * 0.5,
                          (this->vertex(1) - this->vertex(0)).norm() * 0.5);

  // The list is_initial_guess_vertex contains bool values and has
  // the same size as the number of vertices per object.
  // The entries of is_initial_guess_vertex are set true only for those
  // two vertices corresponding to the largest diagonal which is being used
  // to construct the initial ball.
  // We employ this mask to skip these two vertices while enlarging the ball.
  std::vector<bool> is_initial_guess_vertex(this->n_vertices());

  // First let all the vertices be outside
  std::fill(is_initial_guess_vertex.begin(),
            is_initial_guess_vertex.end(),
            false);

  // Get an initial guess by looking at the largest diagonal
  Point<spacedim> center;
  double          radius = 0;

  switch (structdim)
    {
      case 2:
        {
          const Point<spacedim> p30(this->vertex(3) - this->vertex(0));
          const Point<spacedim> p21(this->vertex(2) - this->vertex(1));
          if (p30.norm() > p21.norm())
            {
              center                     = this->vertex(0) + 0.5 * p30;
              radius                     = p30.norm() / 2.;
              is_initial_guess_vertex[3] = true;
              is_initial_guess_vertex[0] = true;
            }
          else
            {
              center                     = this->vertex(1) + 0.5 * p21;
              radius                     = p21.norm() / 2.;
              is_initial_guess_vertex[2] = true;
              is_initial_guess_vertex[1] = true;
            }
          break;
        }
      case 3:
        {
          const Point<spacedim>     p70(this->vertex(7) - this->vertex(0));
          const Point<spacedim>     p61(this->vertex(6) - this->vertex(1));
          const Point<spacedim>     p25(this->vertex(2) - this->vertex(5));
          const Point<spacedim>     p34(this->vertex(3) - this->vertex(4));
          const std::vector<double> diagonals = {p70.norm(),
                                                 p61.norm(),
                                                 p25.norm(),
                                                 p34.norm()};
          const std::vector<double>::const_iterator it =
            std::max_element(diagonals.begin(), diagonals.end());
          if (it == diagonals.begin())
            {
              center                     = this->vertex(0) + 0.5 * p70;
              is_initial_guess_vertex[7] = true;
              is_initial_guess_vertex[0] = true;
            }
          else if (it == diagonals.begin() + 1)
            {
              center                     = this->vertex(1) + 0.5 * p61;
              is_initial_guess_vertex[6] = true;
              is_initial_guess_vertex[1] = true;
            }
          else if (it == diagonals.begin() + 2)
            {
              center                     = this->vertex(5) + 0.5 * p25;
              is_initial_guess_vertex[2] = true;
              is_initial_guess_vertex[5] = true;
            }
          else
            {
              center                     = this->vertex(4) + 0.5 * p34;
              is_initial_guess_vertex[3] = true;
              is_initial_guess_vertex[4] = true;
            }
          radius = *it * 0.5;
          break;
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
        return std::pair<Point<spacedim>, double>();
    }

  // For each vertex that is found to be geometrically outside the ball
  // enlarge the ball  so that the new ball contains both the previous ball
  // and the given vertex.
  for (const unsigned int v : this->vertex_indices())
    if (!is_initial_guess_vertex[v])
      {
        const double distance = center.distance(this->vertex(v));
        if (distance > radius)
          {
            // we found a vertex which is outside of the ball
            // extend it (move center and change radius)
            const Point<spacedim> pCV(center - this->vertex(v));
            radius = (distance + radius) * 0.5;
            center = this->vertex(v) + pCV * (radius / distance);

            // Now the new ball constructed in this block
            // encloses the vertex (v) that was found to be geometrically
            // outside the old ball.
          }
      }
  if constexpr (running_in_debug_mode())
    {
      bool all_vertices_within_ball = true;

      // Set all_vertices_within_ball false if any of the vertices of the object
      // are geometrically outside the ball
      for (const unsigned int v : this->vertex_indices())
        if (center.distance(this->vertex(v)) >
            radius + 100. * std::numeric_limits<double>::epsilon())
          {
            all_vertices_within_ball = false;
            break;
          }
      // If all the vertices are not within the ball throw error
      Assert(all_vertices_within_ball, ExcInternalError());
    }
  return std::make_pair(center, radius);
}


template <int structdim, int dim, int spacedim>
double
TriaAccessor<structdim, dim, spacedim>::minimum_vertex_distance() const
{
  switch (structdim)
    {
      case 1:
        return (this->vertex(1) - this->vertex(0)).norm();
      case 2:
      case 3:
        {
          double min = std::numeric_limits<double>::max();
          for (const unsigned int i : this->vertex_indices())
            for (unsigned int j = i + 1; j < this->n_vertices(); ++j)
              min = std::min(min,
                             (this->vertex(i) - this->vertex(j)) *
                               (this->vertex(i) - this->vertex(j)));
          return std::sqrt(min);
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
        return -1e10;
    }
}


template <int structdim, int dim, int spacedim>
bool
TriaAccessor<structdim, dim, spacedim>::is_translation_of(
  const TriaIterator<TriaAccessor<structdim, dim, spacedim>> &o) const
{
  // go through the vertices and check... The
  // cell is a translation of the previous
  // one in case the distance between the
  // individual vertices in the two cell is
  // the same for all the vertices. So do the
  // check by first getting the distance on
  // the first vertex, and then checking
  // whether all others have the same down to
  // rounding errors (we have to be careful
  // here because the calculation of the
  // displacement between one cell and the
  // next can already result in the loss of
  // one or two digits), so we choose 1e-12
  // times the distance between the zeroth
  // vertices here.
  bool                      is_translation = true;
  const Tensor<1, spacedim> dist           = o->vertex(0) - this->vertex(0);
  const double              tol_square     = 1e-24 * dist.norm_square();
  for (unsigned int i = 1; i < this->n_vertices(); ++i)
    {
      const Tensor<1, spacedim> dist_new =
        (o->vertex(i) - this->vertex(i)) - dist;
      if (dist_new.norm_square() > tol_square)
        {
          is_translation = false;
          break;
        }
    }
  return is_translation;
}



template <int structdim, int dim, int spacedim>
unsigned int
TriaAccessor<structdim, dim, spacedim>::n_vertices() const
{
  return this->reference_cell().n_vertices();
}



template <int structdim, int dim, int spacedim>
unsigned int
TriaAccessor<structdim, dim, spacedim>::n_lines() const
{
  return this->reference_cell().n_lines();
}



template <int structdim, int dim, int spacedim>
unsigned int
TriaAccessor<structdim, dim, spacedim>::n_faces() const
{
  Assert(structdim == dim,
         ExcMessage("This function can only be used on objects "
                    "that are cells, but not on faces or edges "
                    "that bound cells."));

  return this->reference_cell().n_faces();
}



template <int structdim, int dim, int spacedim>
std_cxx20::ranges::iota_view<unsigned int, unsigned int>
TriaAccessor<structdim, dim, spacedim>::vertex_indices() const
{
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(0U,
                                                                  n_vertices());
}



template <int structdim, int dim, int spacedim>
std_cxx20::ranges::iota_view<unsigned int, unsigned int>
TriaAccessor<structdim, dim, spacedim>::line_indices() const
{
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(0U,
                                                                  n_lines());
}



template <int structdim, int dim, int spacedim>
std_cxx20::ranges::iota_view<unsigned int, unsigned int>
TriaAccessor<structdim, dim, spacedim>::face_indices() const
{
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(0U,
                                                                  n_faces());
}



/*----------------- Functions: TriaAccessor<0,dim,spacedim> -----------------*/

template <int dim, int spacedim>
inline TriaAccessor<0, dim, spacedim>::TriaAccessor(
  const Triangulation<dim, spacedim> *tria,
  const unsigned int                  vertex_index)
  : tria(tria)
  , global_vertex_index(vertex_index)
{}



template <int dim, int spacedim>
inline TriaAccessor<0, dim, spacedim>::TriaAccessor(
  const Triangulation<dim, spacedim> *tria,
  const int /*level*/,
  const int index,
  const AccessorData *)
  : tria(tria)
  , global_vertex_index(index)
{}



template <int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
inline TriaAccessor<0, dim, spacedim>::TriaAccessor(
  const TriaAccessor<structdim2, dim2, spacedim2> &)
  : tria(nullptr)
  , global_vertex_index(numbers::invalid_unsigned_int)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
inline TriaAccessor<0, dim, spacedim>::TriaAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
  : tria(nullptr)
  , global_vertex_index(numbers::invalid_unsigned_int)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int dim, int spacedim>
inline void
TriaAccessor<0, dim, spacedim>::copy_from(const TriaAccessor &t)
{
  tria                = t.tria;
  global_vertex_index = t.global_vertex_index;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::operator<(
  const TriaAccessor<0, dim, spacedim> &other) const
{
  Assert(tria == other.tria, TriaAccessorExceptions::ExcCantCompareIterators());

  return (global_vertex_index < other.global_vertex_index);
}



template <int dim, int spacedim>
inline IteratorState::IteratorStates
TriaAccessor<0, dim, spacedim>::state() const
{
  if (global_vertex_index != numbers::invalid_unsigned_int)
    return IteratorState::valid;
  else
    return IteratorState::past_the_end;
}



template <int dim, int spacedim>
inline int
TriaAccessor<0, dim, spacedim>::level()
{
  return 0;
}



template <int dim, int spacedim>
inline int
TriaAccessor<0, dim, spacedim>::index() const
{
  return global_vertex_index;
}



template <int dim, int spacedim>
inline const Triangulation<dim, spacedim> &
TriaAccessor<0, dim, spacedim>::get_triangulation() const
{
  return *tria;
}



template <int dim, int spacedim>
inline void
TriaAccessor<0, dim, spacedim>::operator++()
{
  ++global_vertex_index;
  if (global_vertex_index >= tria->n_vertices())
    global_vertex_index = numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline void
TriaAccessor<0, dim, spacedim>::operator--()
{
  if (global_vertex_index != numbers::invalid_unsigned_int)
    {
      if (global_vertex_index != 0)
        --global_vertex_index;
      else
        global_vertex_index = numbers::invalid_unsigned_int;
    }
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::operator==(const TriaAccessor &t) const
{
  const bool result =
    ((tria == t.tria) && (global_vertex_index == t.global_vertex_index));

  return result;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::operator!=(const TriaAccessor &t) const
{
  return !(*this == t);
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::vertex_index(const unsigned int) const
{
  return global_vertex_index;
}



template <int dim, int spacedim>
inline Point<spacedim> &
TriaAccessor<0, dim, spacedim>::vertex(const unsigned int) const
{
  return const_cast<Point<spacedim> &>(
    this->tria->vertices[global_vertex_index]);
}



template <int dim, int spacedim>
inline typename dealii::internal::TriangulationImplementation::
  Iterators<dim, spacedim>::line_iterator
  TriaAccessor<0, dim, spacedim>::line(const unsigned int)
{
  return typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::line_iterator();
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::line_index(const unsigned int)
{
  Assert(false, ExcImpossibleInDim(0));
  return numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline typename dealii::internal::TriangulationImplementation::
  Iterators<dim, spacedim>::quad_iterator
  TriaAccessor<0, dim, spacedim>::quad(const unsigned int)
{
  return typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::quad_iterator();
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::quad_index(const unsigned int)
{
  Assert(false, ExcImpossibleInDim(0));
  return numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline double
TriaAccessor<0, dim, spacedim>::diameter() const
{
  return 0.;
}



template <int dim, int spacedim>
inline double
TriaAccessor<0, dim, spacedim>::extent_in_direction(const unsigned int) const
{
  return 0.;
}



template <int dim, int spacedim>
inline Point<spacedim>
TriaAccessor<0, dim, spacedim>::center(const bool, const bool) const
{
  return this->tria->vertices[global_vertex_index];
}



template <int dim, int spacedim>
inline double
TriaAccessor<0, dim, spacedim>::measure() const
{
  return 1.0;
}



template <int dim, int spacedim>
inline types::geometric_orientation
TriaAccessor<0, dim, spacedim>::combined_face_orientation(
  const unsigned int /*face*/)
{
  return numbers::default_geometric_orientation;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::face_orientation(const unsigned int /*face*/)
{
  return false;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::face_flip(const unsigned int /*face*/)
{
  return false;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::face_rotation(const unsigned int /*face*/)
{
  return false;
}



template <int dim, int spacedim>
inline types::geometric_orientation
TriaAccessor<0, dim, spacedim>::line_orientation(const unsigned int /*line*/)
{
  return numbers::reverse_line_orientation;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::has_children()
{
  return false;
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::n_children()
{
  return 0;
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::n_active_descendants()
{
  return 0;
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::max_refinement_depth()
{
  return 0;
}



template <int dim, int spacedim>
inline unsigned int
TriaAccessor<0, dim, spacedim>::child_iterator_to_index(
  const TriaIterator<TriaAccessor<0, dim, spacedim>> &)
{
  return numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline TriaIterator<TriaAccessor<0, dim, spacedim>>
TriaAccessor<0, dim, spacedim>::child(const unsigned int)
{
  return TriaIterator<TriaAccessor<0, dim, spacedim>>();
}



template <int dim, int spacedim>
inline TriaIterator<TriaAccessor<0, dim, spacedim>>
TriaAccessor<0, dim, spacedim>::isotropic_child(const unsigned int)
{
  return TriaIterator<TriaAccessor<0, dim, spacedim>>();
}



template <int dim, int spacedim>
inline RefinementCase<0>
TriaAccessor<0, dim, spacedim>::refinement_case()
{
  return {RefinementPossibilities<0>::no_refinement};
}



template <int dim, int spacedim>
inline int
TriaAccessor<0, dim, spacedim>::child_index(const unsigned int)
{
  return -1;
}



template <int dim, int spacedim>
inline int
TriaAccessor<0, dim, spacedim>::isotropic_child_index(const unsigned int)
{
  return -1;
}



template <int dim, int spacedim>
inline bool
TriaAccessor<0, dim, spacedim>::used() const
{
  return tria->vertex_used(global_vertex_index);
}



/*------------------- Functions: TriaAccessor<0,1,spacedim> -----------------*/

template <int spacedim>
inline TriaAccessor<0, 1, spacedim>::TriaAccessor(
  const Triangulation<1, spacedim> *tria,
  const VertexKind                  vertex_kind,
  const unsigned int                vertex_index)
  : tria(tria)
  , vertex_kind(vertex_kind)
  , global_vertex_index(vertex_index)
{}



template <int spacedim>
inline TriaAccessor<0, 1, spacedim>::TriaAccessor(
  const Triangulation<1, spacedim> *tria,
  const int                         level,
  const int                         index,
  const AccessorData *)
  : tria(tria)
  , vertex_kind(interior_vertex)
  , global_vertex_index(numbers::invalid_unsigned_int)
{
  // in general, calling this constructor should yield an error -- users should
  // instead call the one immediately above. however, if you create something
  // like Triangulation<1>::face_iterator() then this calls the default
  // constructor of the iterator which calls the accessor with argument list
  // (0,-2,-2,0), so in this particular case accept this call and create an
  // object that corresponds to the default constructed (invalid) vertex
  // accessor
  (void)level;
  (void)index;
  Assert((level == -2) && (index == -2),
         ExcMessage(
           "This constructor can not be called for face iterators in 1d, "
           "except to default-construct iterator objects."));
}



template <int spacedim>
template <int structdim2, int dim2, int spacedim2>
inline TriaAccessor<0, 1, spacedim>::TriaAccessor(
  const TriaAccessor<structdim2, dim2, spacedim2> &)
  : tria(nullptr)
  , vertex_kind(interior_vertex)
  , global_vertex_index(numbers::invalid_unsigned_int)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int spacedim>
template <int structdim2, int dim2, int spacedim2>
inline TriaAccessor<0, 1, spacedim>::TriaAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
  : tria(nullptr)
  , vertex_kind(interior_vertex)
  , global_vertex_index(numbers::invalid_unsigned_int)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::copy_from(const TriaAccessor &t)
{
  tria                = t.tria;
  vertex_kind         = t.vertex_kind;
  global_vertex_index = t.global_vertex_index;
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::copy_from(
  const TriaAccessorBase<0, 1, spacedim> &)
{
  // We cannot convert from TriaAccessorBase to
  // TriaAccessor<0,1,spacedim> because the latter is not derived from
  // the former. We should never get here.
  DEAL_II_ASSERT_UNREACHABLE();
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::operator<(
  const TriaAccessor<0, 1, spacedim> &other) const
{
  Assert(tria == other.tria, TriaAccessorExceptions::ExcCantCompareIterators());

  return (global_vertex_index < other.global_vertex_index);
}



template <int spacedim>
inline IteratorState::IteratorStates
TriaAccessor<0, 1, spacedim>::state()
{
  return IteratorState::valid;
}


template <int spacedim>
inline int
TriaAccessor<0, 1, spacedim>::level()
{
  return 0;
}



template <int spacedim>
inline int
TriaAccessor<0, 1, spacedim>::index() const
{
  return global_vertex_index;
}



template <int spacedim>
inline const Triangulation<1, spacedim> &
TriaAccessor<0, 1, spacedim>::get_triangulation() const
{
  return *tria;
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::operator++() const
{
  DEAL_II_NOT_IMPLEMENTED();
}


template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::operator--() const
{
  DEAL_II_NOT_IMPLEMENTED();
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::operator==(const TriaAccessor &t) const
{
  const bool result =
    ((tria == t.tria) && (global_vertex_index == t.global_vertex_index));
  // if we point to the same vertex,
  // make sure we know the same about
  // it
  if (result == true)
    Assert(vertex_kind == t.vertex_kind, ExcInternalError());

  return result;
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::operator!=(const TriaAccessor &t) const
{
  return !(*this == t);
}



template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::vertex_index(const unsigned int i) const
{
  AssertIndexRange(i, 1);
  (void)i;
  return global_vertex_index;
}



template <int spacedim>
inline Point<spacedim> &
TriaAccessor<0, 1, spacedim>::vertex(const unsigned int i) const
{
  AssertIndexRange(i, 1);
  (void)i;
  return const_cast<Point<spacedim> &>(
    this->tria->vertices[global_vertex_index]);
}



template <int spacedim>
inline Point<spacedim>
TriaAccessor<0, 1, spacedim>::center() const
{
  return this->tria->vertices[global_vertex_index];
}



template <int spacedim>
inline typename dealii::internal::TriangulationImplementation::
  Iterators<1, spacedim>::line_iterator
  TriaAccessor<0, 1, spacedim>::line(const unsigned int)
{
  return {};
}


template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::line_index(const unsigned int)
{
  Assert(false, ExcImpossibleInDim(0));
  return numbers::invalid_unsigned_int;
}


template <int spacedim>
inline typename dealii::internal::TriangulationImplementation::
  Iterators<1, spacedim>::quad_iterator
  TriaAccessor<0, 1, spacedim>::quad(const unsigned int)
{
  return {};
}



template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::quad_index(const unsigned int)
{
  Assert(false, ExcImpossibleInDim(0));
  return numbers::invalid_unsigned_int;
}



template <int spacedim>
inline double
TriaAccessor<0, 1, spacedim>::measure()
{
  return 1.0;
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::at_boundary() const
{
  return vertex_kind != interior_vertex;
}


template <int spacedim>
inline types::boundary_id
TriaAccessor<0, 1, spacedim>::boundary_id() const
{
  switch (vertex_kind)
    {
      case left_vertex:
      case right_vertex:
        {
          Assert(tria->vertex_to_boundary_id_map_1d->find(
                   this->vertex_index()) !=
                   tria->vertex_to_boundary_id_map_1d->end(),
                 ExcInternalError());

          return (*tria->vertex_to_boundary_id_map_1d)[this->vertex_index()];
        }

      default:
        return numbers::internal_face_boundary_id;
    }
}



template <int spacedim>
inline const Manifold<1, spacedim> &
TriaAccessor<0, 1, spacedim>::get_manifold() const
{
  return this->tria->get_manifold(this->manifold_id());
}



template <int spacedim>
inline types::manifold_id
TriaAccessor<0, 1, spacedim>::manifold_id() const
{
  if (tria->vertex_to_manifold_id_map_1d->find(this->vertex_index()) !=
      tria->vertex_to_manifold_id_map_1d->end())
    return (*tria->vertex_to_manifold_id_map_1d)[this->vertex_index()];
  else
    return numbers::flat_manifold_id;
}


template <int spacedim>
inline types::geometric_orientation
TriaAccessor<0, 1, spacedim>::combined_face_orientation(
  const unsigned int /*face*/)
{
  return numbers::reverse_line_orientation;
}


template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::face_orientation(const unsigned int /*face*/)
{
  return false;
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::face_flip(const unsigned int /*face*/)
{
  return false;
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::face_rotation(const unsigned int /*face*/)
{
  return false;
}



template <int spacedim>
inline types::geometric_orientation
TriaAccessor<0, 1, spacedim>::line_orientation(const unsigned int /*line*/)
{
  return numbers::reverse_line_orientation;
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::has_children()
{
  return false;
}



template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::n_children()
{
  return 0;
}



template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::n_active_descendants()
{
  return 0;
}



template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::max_refinement_depth()
{
  return 0;
}



template <int spacedim>
inline unsigned int
TriaAccessor<0, 1, spacedim>::child_iterator_to_index(
  const TriaIterator<TriaAccessor<0, 1, spacedim>> &)
{
  return numbers::invalid_unsigned_int;
}



template <int spacedim>
inline TriaIterator<TriaAccessor<0, 1, spacedim>>
TriaAccessor<0, 1, spacedim>::child(const unsigned int)
{
  return TriaIterator<TriaAccessor<0, 1, spacedim>>();
}


template <int spacedim>
inline TriaIterator<TriaAccessor<0, 1, spacedim>>
TriaAccessor<0, 1, spacedim>::isotropic_child(const unsigned int)
{
  return TriaIterator<TriaAccessor<0, 1, spacedim>>();
}


template <int spacedim>
inline RefinementCase<0>
TriaAccessor<0, 1, spacedim>::refinement_case()
{
  return {RefinementPossibilities<0>::no_refinement};
}

template <int spacedim>
inline int
TriaAccessor<0, 1, spacedim>::child_index(const unsigned int)
{
  return -1;
}


template <int spacedim>
inline int
TriaAccessor<0, 1, spacedim>::isotropic_child_index(const unsigned int)
{
  return -1;
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::set_boundary_id(const types::boundary_id b) const
{
  Assert(tria->vertex_to_boundary_id_map_1d->find(this->vertex_index()) !=
           tria->vertex_to_boundary_id_map_1d->end(),
         ExcMessage("You can't set the boundary_id of a face of a cell that is "
                    "not actually at the boundary."));

  (*tria->vertex_to_boundary_id_map_1d)[this->vertex_index()] = b;
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::set_manifold_id(const types::manifold_id b)
{
  (*tria->vertex_to_manifold_id_map_1d)[this->vertex_index()] = b;
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::set_all_boundary_ids(
  const types::boundary_id b) const
{
  set_boundary_id(b);
}



template <int spacedim>
inline void
TriaAccessor<0, 1, spacedim>::set_all_manifold_ids(const types::manifold_id b)
{
  set_manifold_id(b);
}



template <int spacedim>
inline bool
TriaAccessor<0, 1, spacedim>::used() const
{
  return tria->vertex_used(global_vertex_index);
}



template <int spacedim>
inline ReferenceCell
TriaAccessor<0, 1, spacedim>::reference_cell() const
{
  return ReferenceCells::Vertex;
}



template <int spacedim>
unsigned int
TriaAccessor<0, 1, spacedim>::n_vertices() const
{
  return 1;
}



template <int spacedim>
unsigned int
TriaAccessor<0, 1, spacedim>::n_lines() const
{
  return 0;
}



template <int spacedim>
std_cxx20::ranges::iota_view<unsigned int, unsigned int>
TriaAccessor<0, 1, spacedim>::vertex_indices() const
{
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(0U,
                                                                  n_vertices());
}



template <int spacedim>
std_cxx20::ranges::iota_view<unsigned int, unsigned int>
TriaAccessor<0, 1, spacedim>::line_indices() const
{
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(0U,
                                                                  n_lines());
}

/*------------------ Functions: CellAccessor<dim,spacedim> ------------------*/


template <int dim, int spacedim>
inline CellAccessor<dim, spacedim>::CellAccessor(
  const Triangulation<dim, spacedim> *parent,
  const int                           level,
  const int                           index,
  const AccessorData                 *local_data)
  : TriaAccessor<dim, dim, spacedim>(parent, level, index, local_data)
{}



template <int dim, int spacedim>
inline CellAccessor<dim, spacedim>::CellAccessor(
  const TriaAccessor<dim, dim, spacedim> &cell_accessor)
  : TriaAccessor<dim, dim, spacedim>(
      static_cast<const TriaAccessor<dim, dim, spacedim> &>(cell_accessor))
{}



template <int dim, int spacedim>
inline TriaIterator<CellAccessor<dim, spacedim>>
CellAccessor<dim, spacedim>::child(const unsigned int i) const
{
  TriaIterator<CellAccessor<dim, spacedim>> q(this->tria,
                                              this->present_level + 1,
                                              this->child_index(i));

  Assert((q.state() == IteratorState::past_the_end) || q->used(),
         ExcInternalError());

  return q;
}



template <int dim, int spacedim>
inline boost::container::small_vector<TriaIterator<CellAccessor<dim, spacedim>>,
                                      GeometryInfo<dim>::max_children_per_cell>
CellAccessor<dim, spacedim>::child_iterators() const
{
  boost::container::small_vector<TriaIterator<CellAccessor<dim, spacedim>>,
                                 GeometryInfo<dim>::max_children_per_cell>
    child_iterators(this->n_children());

  for (unsigned int i = 0; i < this->n_children(); ++i)
    child_iterators[i] = this->child(i);

  return child_iterators;
}



template <int dim, int spacedim>
inline TriaIterator<TriaAccessor<dim - 1, dim, spacedim>>
CellAccessor<dim, spacedim>::face(const unsigned int i) const
{
  AssertIndexRange(i, this->n_faces());
  if constexpr (dim == 1)
    {
      using VertexKind = typename TriaAccessor<0, 1, spacedim>::VertexKind;
      VertexKind vertex_kind = VertexKind::interior_vertex;
      if (i == 0 && at_boundary(0))
        vertex_kind = VertexKind::left_vertex;
      if (i == 1 && at_boundary(1))
        vertex_kind = VertexKind::right_vertex;
      TriaAccessor<0, 1, spacedim> a(&this->get_triangulation(),
                                     vertex_kind,
                                     this->vertex_index(i));
      return TriaIterator<TriaAccessor<0, 1, spacedim>>(a);
    }
  else if constexpr (dim == 2)
    return this->line(i);
  else if constexpr (dim == 3)
    return this->quad(i);
  else
    {
      Assert(false, ExcNotImplemented());
      return {};
    }
}



template <int dim, int spacedim>
inline unsigned int
CellAccessor<dim, spacedim>::face_iterator_to_index(
  const TriaIterator<TriaAccessor<dim - 1, dim, spacedim>> &face) const
{
  for (const unsigned int face_n : this->face_indices())
    if (this->face(face_n) == face)
      return face_n;

  Assert(false,
         ExcMessage("The given face is not a face of the current cell."));
  return numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
inline boost::container::small_vector<
  TriaIterator<TriaAccessor<dim - 1, dim, spacedim>>,
#  ifndef _MSC_VER
  ReferenceCells::max_n_faces<dim>()
#  else
  GeometryInfo<dim>::faces_per_cell
#  endif
  >
CellAccessor<dim, spacedim>::face_iterators() const
{
  boost::container::small_vector<
    TriaIterator<TriaAccessor<dim - 1, dim, spacedim>>,
#  ifndef _MSC_VER
    ReferenceCells::max_n_faces<dim>()
#  else
    GeometryInfo<dim>::faces_per_cell
#  endif
    >
    face_iterators(this->n_faces());

  for (const unsigned int i : this->face_indices())
    face_iterators[i] = this->face(i);

  return face_iterators;
}



template <int dim, int spacedim>
inline unsigned int
CellAccessor<dim, spacedim>::face_index(const unsigned int i) const
{
  switch (dim)
    {
      case 1:
        return this->vertex_index(i);

      case 2:
        return this->line_index(i);

      case 3:
        return this->quad_index(i);

      default:
        return numbers::invalid_unsigned_int;
    }
}



template <int dim, int spacedim>
inline int
CellAccessor<dim, spacedim>::neighbor_index(const unsigned int face_no) const
{
  AssertIndexRange(face_no, this->n_faces());
  return this->tria->levels[this->present_level]
    ->neighbors[this->present_index * ReferenceCells::max_n_faces<dim>() +
                face_no]
    .second;
}



template <int dim, int spacedim>
inline int
CellAccessor<dim, spacedim>::neighbor_level(const unsigned int face_no) const
{
  AssertIndexRange(face_no, this->n_faces());
  return this->tria->levels[this->present_level]
    ->neighbors[this->present_index * ReferenceCells::max_n_faces<dim>() +
                face_no]
    .first;
}



template <int dim, int spacedim>
inline RefinementCase<dim>
CellAccessor<dim, spacedim>::refine_flag_set() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  // cells flagged for refinement must be active
  // (the @p set_refine_flag function checks this,
  // but activity may change when refinement is
  // executed and for some reason the refine
  // flag is not cleared).
  Assert(this->is_active() || !this->tria->levels[this->present_level]
                                 ->refine_flags[this->present_index],
         ExcRefineCellNotActive());
  return RefinementCase<dim>(
    this->tria->levels[this->present_level]->refine_flags[this->present_index]);
}



template <int dim, int spacedim>
inline void
CellAccessor<dim, spacedim>::set_refine_flag(
  const RefinementCase<dim> refinement_case) const
{
  Assert(this->used() && this->is_active(), ExcRefineCellNotActive());
  Assert(!coarsen_flag_set(), ExcCellFlaggedForCoarsening());

  this->tria->levels[this->present_level]->refine_flags[this->present_index] =
    refinement_case;
}



template <int dim, int spacedim>
inline void
CellAccessor<dim, spacedim>::clear_refine_flag() const
{
  Assert(this->used() && this->is_active(), ExcRefineCellNotActive());
  this->tria->levels[this->present_level]->refine_flags[this->present_index] =
    RefinementCase<dim>::no_refinement;
}


template <int dim, int spacedim>
inline std::uint8_t
CellAccessor<dim, spacedim>::refine_choice() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  if (this->tria->levels[this->present_level]->refine_choice.size() == 0)
    return 0U;
  return this->tria->levels[this->present_level]
    ->refine_choice[this->present_index];
}


template <int dim, int spacedim>
inline void
CellAccessor<dim, spacedim>::set_refine_choice(
  const std::uint8_t refinement_choice) const
{
  Assert(this->used() && this->is_active(), ExcRefineCellNotActive());
  if (this->tria->levels[this->present_level]->refine_choice.size() != 0)
    this->tria->levels[this->present_level]
      ->refine_choice[this->present_index] = refinement_choice;
}


template <int dim, int spacedim>
inline void
CellAccessor<dim, spacedim>::clear_refine_choice() const
{
  Assert(this->used() && this->is_active(), ExcRefineCellNotActive());
  if (this->tria->levels[this->present_level]->refine_choice.size() != 0)
    this->tria->levels[this->present_level]
      ->refine_choice[this->present_index] =
      static_cast<char>(IsotropicRefinementChoice::isotropic_refinement);
}


template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::flag_for_face_refinement(
  const unsigned int             face_no,
  const RefinementCase<dim - 1> &face_refinement_case) const
{
  Assert(dim > 1, ExcImpossibleInDim(dim));
  AssertIndexRange(face_no, this->n_faces());
  AssertIndexRange(face_refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);

  // the new refinement case is a combination
  // of the minimum required one for the given
  // face refinement and the already existing
  // flagged refinement case
  RefinementCase<dim> old_ref_case = refine_flag_set();
  RefinementCase<dim> new_ref_case =
    (old_ref_case |
     GeometryInfo<dim>::min_cell_refinement_case_for_face_refinement(
       face_refinement_case,
       face_no,
       this->face_orientation(face_no),
       this->face_flip(face_no),
       this->face_rotation(face_no)));
  set_refine_flag(new_ref_case);
  // return, whether we had to change the
  // refinement flag
  return new_ref_case != old_ref_case;
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::flag_for_line_refinement(
  const unsigned int line_no) const
{
  Assert(dim > 1, ExcImpossibleInDim(dim));
  AssertIndexRange(line_no, this->n_lines());

  // the new refinement case is a combination
  // of the minimum required one for the given
  // line refinement and the already existing
  // flagged refinement case
  RefinementCase<dim>
    old_ref_case = refine_flag_set(),
    new_ref_case =
      old_ref_case |
      GeometryInfo<dim>::min_cell_refinement_case_for_line_refinement(line_no);
  set_refine_flag(new_ref_case);
  // return, whether we had to change the
  // refinement flag
  return new_ref_case != old_ref_case;
}



template <int dim, int spacedim>
inline internal::SubfaceCase<dim>
CellAccessor<dim, spacedim>::subface_case(const unsigned int face_no) const
{
  Assert(is_active(), TriaAccessorExceptions::ExcCellNotActive());
  AssertIndexRange(face_no, this->n_faces());

  if constexpr (dim == 1)
    return internal::SubfaceCase<1>::case_none;
  else if constexpr (dim == 2)
    return ((face(face_no)->has_children()) ?
              internal::SubfaceCase<2>::case_x :
              internal::SubfaceCase<2>::case_none);
  else if constexpr (dim == 3)
    {
      switch (static_cast<std::uint8_t>(face(face_no)->refinement_case()))
        {
          case RefinementCase<3>::no_refinement:
            return internal::SubfaceCase<3>::case_none;
          case RefinementCase<3>::cut_x:
            if (face(face_no)->child(0)->has_children())
              {
                Assert(face(face_no)->child(0)->refinement_case() ==
                         RefinementCase<2>::cut_y,
                       ExcInternalError());
                if (face(face_no)->child(1)->has_children())
                  {
                    Assert(face(face_no)->child(1)->refinement_case() ==
                             RefinementCase<2>::cut_y,
                           ExcInternalError());
                    return internal::SubfaceCase<3>::case_x1y2y;
                  }
                else
                  return internal::SubfaceCase<3>::case_x1y;
              }
            else
              {
                if (face(face_no)->child(1)->has_children())
                  {
                    Assert(face(face_no)->child(1)->refinement_case() ==
                             RefinementCase<2>::cut_y,
                           ExcInternalError());
                    return internal::SubfaceCase<3>::case_x2y;
                  }
                else
                  return internal::SubfaceCase<3>::case_x;
              }
          case RefinementCase<3>::cut_y:
            if (face(face_no)->child(0)->has_children())
              {
                Assert(face(face_no)->child(0)->refinement_case() ==
                         RefinementCase<2>::cut_x,
                       ExcInternalError());
                if (face(face_no)->child(1)->has_children())
                  {
                    Assert(face(face_no)->child(1)->refinement_case() ==
                             RefinementCase<2>::cut_x,
                           ExcInternalError());
                    return internal::SubfaceCase<3>::case_y1x2x;
                  }
                else
                  return internal::SubfaceCase<3>::case_y1x;
              }
            else
              {
                if (face(face_no)->child(1)->has_children())
                  {
                    Assert(face(face_no)->child(1)->refinement_case() ==
                             RefinementCase<2>::cut_x,
                           ExcInternalError());
                    return internal::SubfaceCase<3>::case_y2x;
                  }
                else
                  return internal::SubfaceCase<3>::case_y;
              }
          case RefinementCase<3>::cut_xy:
            return internal::SubfaceCase<3>::case_xy;
          default:
            DEAL_II_ASSERT_UNREACHABLE();
        }
    }

  // we should never get here
  DEAL_II_ASSERT_UNREACHABLE();
  return internal::SubfaceCase<dim>::case_none;
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::coarsen_flag_set() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  // cells flagged for coarsening must be active
  // (the @p set_refine_flag function checks this,
  // but activity may change when refinement is
  // executed and for some reason the refine
  // flag is not cleared).
  Assert(this->is_active() || !this->tria->levels[this->present_level]
                                 ->coarsen_flags[this->present_index],
         ExcRefineCellNotActive());
  return this->tria->levels[this->present_level]
    ->coarsen_flags[this->present_index];
}



template <int dim, int spacedim>
inline void
CellAccessor<dim, spacedim>::set_coarsen_flag() const
{
  Assert(this->used() && this->is_active(), ExcRefineCellNotActive());
  Assert(!refine_flag_set(), ExcCellFlaggedForRefinement());

  this->tria->levels[this->present_level]->coarsen_flags[this->present_index] =
    true;
}



template <int dim, int spacedim>
inline void
CellAccessor<dim, spacedim>::clear_coarsen_flag() const
{
  Assert(this->used() && this->is_active(), ExcRefineCellNotActive());
  this->tria->levels[this->present_level]->coarsen_flags[this->present_index] =
    false;
}



template <int dim, int spacedim>
inline TriaIterator<CellAccessor<dim, spacedim>>
CellAccessor<dim, spacedim>::neighbor(const unsigned int face_no) const
{
  TriaIterator<CellAccessor<dim, spacedim>> q(this->tria,
                                              neighbor_level(face_no),
                                              neighbor_index(face_no));

  Assert((q.state() == IteratorState::past_the_end) || q->used(),
         ExcInternalError());

  return q;
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_active() const
{
  return !this->has_children();
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_locally_owned() const
{
  Assert(this->is_active(),
         ExcMessage("is_locally_owned() can only be called on active cells!"));
#  ifndef DEAL_II_WITH_MPI
  return true;
#  else

  // Serial triangulations report invalid_subdomain_id as their locally owned
  // subdomain, so the first condition checks whether we have a serial
  // triangulation, in which case all cells are locally owned. The second
  // condition compares the subdomain id in the parallel case.
  const types::subdomain_id locally_owned_subdomain =
    this->tria->locally_owned_subdomain();
  return (locally_owned_subdomain == numbers::invalid_subdomain_id ||
          this->subdomain_id() == locally_owned_subdomain);

#  endif
}


template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_locally_owned_on_level() const
{
#  ifndef DEAL_II_WITH_MPI
  return true;
#  else

  // Serial triangulations report invalid_subdomain_id as their locally owned
  // subdomain, so the first condition checks whether we have a serial
  // triangulation, in which case all cells are locally owned. The second
  // condition compares the subdomain id in the parallel case.
  const types::subdomain_id locally_owned_subdomain =
    this->tria->locally_owned_subdomain();
  return (locally_owned_subdomain == numbers::invalid_subdomain_id ||
          this->level_subdomain_id() == locally_owned_subdomain);

#  endif
}


template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_ghost() const
{
  Assert(this->is_active(),
         ExcMessage("is_ghost() can only be called on active cells!"));
  if (this->has_children())
    return false;

#  ifndef DEAL_II_WITH_MPI
  return false;
#  else

  // Serial triangulations report invalid_subdomain_id as their locally owned
  // subdomain, so the first condition rules out that case as all cells to a
  // serial triangulation are locally owned and none is ghosted. The second
  // and third conditions check whether the cell's subdomain is not the
  // locally owned one and not artificial.
  const types::subdomain_id locally_owned_subdomain =
    this->tria->locally_owned_subdomain();
  const types::subdomain_id subdomain_id = this->subdomain_id();
  return (locally_owned_subdomain != numbers::invalid_subdomain_id &&
          subdomain_id != locally_owned_subdomain &&
          subdomain_id != numbers::artificial_subdomain_id);

#  endif
}


template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_ghost_on_level() const
{
#  ifndef DEAL_II_WITH_MPI
  return false;
#  else

  // Serial triangulations report invalid_subdomain_id as their locally owned
  // subdomain, so the first condition checks whether we have a serial
  // triangulation, in which case all cells are locally owned. The second
  // condition compares the subdomain id in the parallel case.
  const types::subdomain_id locally_owned_subdomain =
    this->tria->locally_owned_subdomain();
  const types::subdomain_id subdomain_id = this->level_subdomain_id();
  return (locally_owned_subdomain != numbers::invalid_subdomain_id &&
          subdomain_id != locally_owned_subdomain &&
          subdomain_id != numbers::artificial_subdomain_id);

#  endif
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_artificial() const
{
  Assert(this->is_active(),
         ExcMessage("is_artificial() can only be called on active cells!"));
#  ifndef DEAL_II_WITH_MPI
  return false;
#  else

  // Serial triangulations report invalid_subdomain_id as their locally owned
  // subdomain, so the first condition rules out that case as all cells to a
  // serial triangulation are locally owned and none is artificial.
  return (this->tria->locally_owned_subdomain() !=
            numbers::invalid_subdomain_id &&
          this->subdomain_id() == numbers::artificial_subdomain_id);

#  endif
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_artificial_on_level() const
{
#  ifndef DEAL_II_WITH_MPI
  return false;
#  else
  return (this->tria->locally_owned_subdomain() !=
            numbers::invalid_subdomain_id &&
          this->level_subdomain_id() == numbers::artificial_subdomain_id);
#  endif
}



template <int dim, int spacedim>
inline types::subdomain_id
CellAccessor<dim, spacedim>::subdomain_id() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert(this->is_active(),
         ExcMessage("subdomain_id() can only be called on active cells!"));
  return this->tria->levels[this->present_level]
    ->subdomain_ids[this->present_index];
}



template <int dim, int spacedim>
inline types::subdomain_id
CellAccessor<dim, spacedim>::level_subdomain_id() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->levels[this->present_level]
    ->level_subdomain_ids[this->present_index];
}



template <int dim, int spacedim>
inline unsigned int
CellAccessor<dim, spacedim>::neighbor_face_no(const unsigned int neighbor) const
{
  const unsigned int n2 = neighbor_of_neighbor_internal(neighbor);
  if (n2 != numbers::invalid_unsigned_int)
    // return this value as the
    // neighbor is not coarser
    return n2;
  else
    // the neighbor is coarser
    return neighbor_of_coarser_neighbor(neighbor).first;
}



template <int dim, int spacedim>
inline bool
CellAccessor<dim, spacedim>::is_level_cell()
{
  return false;
}



template <int dim, int spacedim>
inline unsigned int
CellAccessor<dim, spacedim>::active_cell_index() const
{
  Assert(this->is_active(), TriaAccessorExceptions::ExcCellNotActive());
  return this->tria->levels[this->present_level]
    ->active_cell_indices[this->present_index];
}



template <int dim, int spacedim>
inline types::global_cell_index
CellAccessor<dim, spacedim>::global_active_cell_index() const
{
  Assert(this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert(this->is_active(),
         ExcMessage(
           "global_active_cell_index() can only be called on active cells!"));

  return this->tria->levels[this->present_level]
    ->global_active_cell_indices[this->present_index];
}



template <int dim, int spacedim>
inline types::global_cell_index
CellAccessor<dim, spacedim>::global_level_cell_index() const
{
  return this->tria->levels[this->present_level]
    ->global_level_cell_indices[this->present_index];
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE


#endif
