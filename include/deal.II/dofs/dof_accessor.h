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

#ifndef dealii_dof_accessor_h
#define dealii_dof_accessor_h


#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_faces.h>
#include <deal.II/dofs/dof_iterator_selector.h>
#include <deal.II/dofs/dof_levels.h>

#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/read_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include <boost/container/small_vector.hpp>

#include <limits>
#include <set>
#include <type_traits>
#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename number>
class FullMatrix;
template <typename number>
class Vector;
template <typename number>
class AffineConstraints;

template <typename Accessor>
class TriaRawIterator;

template <int, int>
class FiniteElement;

namespace internal
{
  namespace DoFCellAccessorImplementation
  {
    struct Implementation;
  }

  namespace DoFHandlerImplementation
  {
    struct Implementation;
    namespace Policy
    {
      struct Implementation;
    }
  } // namespace DoFHandlerImplementation

  namespace hp
  {
    namespace DoFHandlerImplementation
    {
      struct Implementation;
    }
  } // namespace hp
} // namespace internal
#endif


namespace internal
{
  namespace DoFAccessorImplementation
  {
    /**
     * This is a switch class which only declares an @p alias. It is meant to
     * determine which class a DoFAccessor class is to be derived from. By
     * default, <tt>DoFAccessor@<structdim,dim,spacedim@></tt> derives from
     * the alias in the general
     * <tt>Inheritance@<structdim,dim,spacedim@></tt> class, which is
     * <tt>TriaAccessor@<structdim,dim,spacedim@></tt>, but if
     * <tt>structdim==dim</tt>, then the specialization
     * <tt>Inheritance@<dim,dim,spacedim@></tt> is used which declares its
     * local type to be <tt>CellAccessor@<dim,spacedim@></tt>. Therefore, the
     * inheritance is automatically chosen to be from CellAccessor if the
     * object under consideration has full dimension, i.e. constitutes a cell.
     *
     * @ingroup dofs
     * @ingroup Accessors
     */
    template <int structdim, int dim, int spacedim>
    struct Inheritance
    {
      /**
       * Declaration of the @p alias. See the full documentation for more
       * information.
       */
      using BaseClass = dealii::TriaAccessor<structdim, dim, spacedim>;
    };


    /**
     * This is the specialization of the general template used for the case
     * where an object has full dimension, i.e. is a cell. See the general
     * template for more details.
     */
    template <int dim, int spacedim>
    struct Inheritance<dim, dim, spacedim>
    {
      /**
       * Declaration of the @p alias. See the full documentation for more
       * information.
       */
      using BaseClass = dealii::CellAccessor<dim, spacedim>;
    };

    struct Implementation;
  } // namespace DoFAccessorImplementation
} // namespace internal


/* -------------------------------------------------------------------------- */



/**
 * A class that gives access to the degrees of freedom stored in a DoFHandler
 * object. Accessors are used to access the data that pertains to edges,
 * faces, and cells of a triangulation. The concept is explained in more
 * detail in connection to
 * @ref Iterators.
 *
 * This class follows mainly the route laid out by the accessor library
 * declared in the triangulation library (TriaAccessor). It enables the user
 * to access the degrees of freedom on lines, quads, or hexes. The first
 * template argument of this class determines the dimensionality of the object
 * under consideration: 1 for lines, 2 for quads, and 3 for hexes. The second
 * argument denotes the type of DoFHandler we should work on. From the second
 * template argument we also deduce the dimension of the Triangulation this
 * object refers to as well as the dimension of the space into which it is
 * embedded. Finally, the template argument <code>level_dof_access</code>
 * governs the behavior of the function get_active_or_mg_dof_indices(). See
 * the section on Generic loops below.
 *
 * <h3>Alias</h3>
 *
 * Usage is best to happen through the alias to the various kinds of iterators
 * provided by the DoFHandler class, since they are more secure to changes in
 * the class naming and template interface as well as providing easier typing
 * (much less complicated names!).
 *
 * <h3>Generic loops and the third template argument</h3>
 *
 * Many loops look very similar, whether they operate on the active dofs of
 * the active cells of the Triangulation or on the level dofs of a single
 * level or the whole grid hierarchy. In order to use polymorphism in such
 * loops, they access degrees of freedom through the function
 * get_active_or_mg_dof_indices(), which changes behavior according to the
 * third template argument.  If the argument is false, then the active dofs of
 * active cells are accessed. If it is true, the level dofs are used.
 * DoFHandler has functions, for instance begin() and begin_mg(), which return
 * either type or the other. Additionally, they can be cast into each other,
 * in case this is needed, since they access the same data.
 *
 * It is recommended to use the function get_active_or_mg_dof_indices()
 * in generic loops in lieu of get_dof_indices() or get_mg_dof_indices().
 *
 * <h3>Inheritance</h3>
 *
 * If the structural dimension given by the first template argument equals the
 * dimension of the DoFHandler (given as the second template argument), then
 * we are obviously dealing with cells, rather than lower-dimensional objects.
 * In that case, inheritance is from CellAccessor, to provide access to all
 * the cell specific information afforded by that class. Otherwise, i.e. for
 * lower-dimensional objects, inheritance is from TriaAccessor.
 *
 * There is a DoFCellAccessor class that provides the equivalent to the
 * CellAccessor class.
 *
 * @tparam structdim The dimensionality of the objects the accessor
 *   represents. For example, points have @p structdim equal to zero,
 *   edges have @p structdim equal to one, etc.
 * @tparam dim Dimension of the underlying DoFHandler.
 * @tparam spacedim Space dimension of the underlying DoFHandler.
 * @tparam level_dof_access If @p false, then the accessor simply represents
 *   a cell, face, or edge in a DoFHandler for which degrees of freedom only
 *   exist on the finest level. Some operations are not allowed in this case,
 *   such as asking for DoF indices on non-active cells. On the other hand,
 *   if this template argument is @p true, then the accessor represents an
 *   object in a multilevel hierarchy of degrees of freedom. In this case,
 *   accessing DoF indices of <i>any</i> cell is possible, and will return
 *   the <i>level</i> indices (which, for active cells, may be different
 *   from the <i>global</i> indices).
 *
 * @ingroup dofs
 * @ingroup Accessors
 */
template <int structdim, int dim, int spacedim, bool level_dof_access>
class DoFAccessor : public dealii::internal::DoFAccessorImplementation::
                      Inheritance<structdim, dim, spacedim>::BaseClass
{
public:
  /**
   * A static variable that allows users of this class to discover the value
   * of the second template argument.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * A static variable that allows users of this class to discover the value
   * of the third template argument.
   */
  static constexpr unsigned int space_dimension = spacedim;

  /**
   * Declare an alias to the base class to make accessing some of the
   * exception classes simpler.
   */
  using BaseClass = typename dealii::internal::DoFAccessorImplementation::
    Inheritance<structdim, dimension, space_dimension>::BaseClass;

  /**
   * Data type passed by the iterator class.
   */
  using AccessorData = DoFHandler<dimension, space_dimension>;

  /**
   * @name Constructors
   */
  /**
   * @{
   */

  /**
   * Default constructor. Provides an accessor that can't be used.
   */
  DoFAccessor();

  /**
   * Constructor that generates an access that points to a particular cell or
   * face or edge in a DoFHandler.
   *
   * @param tria The triangulation into which this accessor points.
   * @param level The level within the mesh hierarchy of the object pointed
   *   to. For example, coarse mesh cells will have level zero, their children
   *   level one, and so on. This argument is ignored for faces and edges
   *   which do not have a level.
   * @param index The index of the object pointed to within the specified
   *   refinement level.
   * @param dof_handler A pointer to the DoFHandler object to which the
   *   accessor shall refer. This DoFHandler object must of course be built on
   *   the same triangulation as the one specified in the first argument.
   */
  DoFAccessor(const Triangulation<dim, spacedim> *tria,
              const int                           level,
              const int                           index,
              const DoFHandler<dim, spacedim>    *dof_handler);

  /**
   * Copy constructor.
   */
  DoFAccessor(const DoFAccessor<structdim, dim, spacedim, level_dof_access> &) =
    default;

  /**
   * Move constructor.
   */
  DoFAccessor(                                                    // NOLINT
    DoFAccessor<structdim, dim, spacedim, level_dof_access> &&) = // NOLINT
    default;                                                      // NOLINT

  /**
   * Destructor.
   */
  ~DoFAccessor() = default;

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
  DoFAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * Another conversion operator between objects that don't make sense, just
   * like the previous one.
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  DoFAccessor(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &);

  /**
   * Copy constructor allowing to switch level access and active access.
   */
  template <bool level_dof_access2>
  DoFAccessor(const DoFAccessor<structdim, dim, spacedim, level_dof_access2> &);

  /**
   * Copy operator. These operators are usually used in a context like
   * <tt>iterator a,b; *a=*b;</tt>. Presumably, the intent here is to copy the
   * object pointed to
   * by @p b to the object pointed to by @p a. However, the result of
   * dereferencing an iterator is not an object but an accessor; consequently,
   * this operation is not useful for iterators on DoF handler objects.
   * Consequently, this operator is declared as deleted and can not be used.
   */
  DoFAccessor<structdim, dim, spacedim, level_dof_access> &
  operator=(const DoFAccessor<structdim, dim, spacedim, level_dof_access> &da) =
    delete;

  /**
   * Move assignment operator.
   */
  DoFAccessor<structdim, dim, spacedim, level_dof_access> &       // NOLINT
  operator=(                                                      // NOLINT
    DoFAccessor<structdim, dim, spacedim, level_dof_access> &&) = // NOLINT
    default;                                                      // NOLINT

  /**
   * @}
   */

  /**
   * Return a handle on the DoFHandler object which we are using.
   */
  const DoFHandler<dim, spacedim> &
  get_dof_handler() const;

  /**
   * Implement the copy operator needed for the iterator classes.
   */
  template <bool level_dof_access2>
  void
  copy_from(const DoFAccessor<structdim, dim, spacedim, level_dof_access2> &a);

  /**
   * Copy operator used by the iterator class. Keeps the previously set dof
   * handler, but sets the object coordinates of the TriaAccessor.
   */
  void
  copy_from(const TriaAccessorBase<structdim, dim, spacedim> &da);

  /**
   * Tell the caller whether get_active_or_mg_dof_indices() accesses active or
   * level dofs.
   */
  static bool
  is_level_cell();

  /**
   * @name Accessing sub-objects
   */
  /**
   * @{
   */

  /**
   * Return an iterator pointing to the @p c-th child.
   */
  TriaIterator<DoFAccessor<structdim, dim, spacedim, level_dof_access>>
  child(const unsigned int c) const;

  /**
   * Pointer to the @p ith line bounding this object. If the current object is
   * a line itself, then the only valid index is @p i equals to zero, and the
   * function returns an iterator to itself.
   */
  typename dealii::internal::DoFHandlerImplementation::
    Iterators<dim, spacedim, level_dof_access>::line_iterator
    line(const unsigned int i) const;

  /**
   * Pointer to the @p ith quad bounding this object. If the current object is
   * a quad itself, then the only valid index is @p i equals to zero, and the
   * function returns an iterator to itself.
   */
  typename dealii::internal::DoFHandlerImplementation::
    Iterators<dim, spacedim, level_dof_access>::quad_iterator
    quad(const unsigned int i) const;

  /**
   * @}
   */

  /**
   * @name Accessing the DoF indices of this object
   */
  /**
   * @{
   */

  /**
   * Return the <i>global</i> indices of the degrees of freedom located on
   * this object in the standard ordering defined by the finite element (i.e.,
   * dofs on vertex 0, dofs on vertex 1, etc, dofs on line 0, dofs on line 1,
   * etc, dofs on quad 0, etc.) This function is only available on
   * <i>active</i> objects (see
   * @ref GlossActive "this glossary entry").
   *
   * The cells needs to be an active cell (and not artificial in a parallel
   * distributed computation).
   *
   * The vector has to have the right size before being passed to this
   * function.
   *
   * The last argument denotes the finite element index. For the standard
   * ::DoFHandler class, this value must be equal to its default value since
   * that class only supports the same finite element on all cells anyway.
   *
   * However, when the relevant DoFHandler object has hp-capabilities enabled,
   * different finite element objects may be used on different cells. On faces
   * between two cells, as well as vertices, there may therefore be two sets
   * of degrees of freedom, one for each of the finite elements used on the
   * adjacent cells. In order to specify which set of degrees of freedom to
   * work on, the last argument is used to disambiguate. Finally, if this
   * function is called for a cell object, there can only be a single set of
   * degrees of freedom, and fe_index has to match the result of
   * active_fe_index().
   *
   * For cells, there is only a single possible finite element index (namely
   * the one for that cell, returned by <code>cell-@>active_fe_index</code>.
   * Consequently, the derived DoFCellAccessor class has an overloaded version
   * of this function that calls the present function with
   * <code>cell-@>active_fe_index</code> as last argument.
   */
  void
  get_dof_indices(
    std::vector<types::global_dof_index> &dof_indices,
    const types::fe_index fe_index = numbers::invalid_fe_index) const;

  /**
   * Return the global multilevel indices of the degrees of freedom that live
   * on the current object with respect to the given level within the
   * multigrid hierarchy. The indices refer to the local numbering for the
   * level this line lives on.
   */
  void
  get_mg_dof_indices(
    const int                             level,
    std::vector<types::global_dof_index> &dof_indices,
    const types::fe_index fe_index = numbers::invalid_fe_index) const;

  /**
   * Set the level DoF indices that are returned by get_mg_dof_indices.
   */
  void
  set_mg_dof_indices(
    const int                                   level,
    const std::vector<types::global_dof_index> &dof_indices,
    const types::fe_index fe_index = numbers::invalid_fe_index);

  /**
   * Global DoF index of the <i>i</i> degree associated with the @p vertexth
   * vertex of the present cell.
   *
   * The last argument denotes the finite element index. For the standard
   * ::DoFHandler class, this value must be equal to its default value since
   * that class only supports the same finite element on all cells anyway.
   *
   * However, when hp-capabilities are enabled, different finite
   * element objects may be used on different cells. On faces between
   * two cells, as well as vertices, there may therefore be two sets
   * of degrees of freedom, one for each of the finite elements used
   * on the adjacent cells.  In order to specify which set of degrees
   * of freedom to work on, the last argument is used to
   * disambiguate. Finally, if this function is called for a cell
   * object, there can only be a single set of degrees of freedom, and
   * `fe_index` has to match the result of
   * `cell->active_fe_index()`. Alternatively, if `fe_index` is left
   * to its default value when this function is called on a cell, then
   * this is interpreted as equal to `cell->active_fe_index()`.
   */
  types::global_dof_index
  vertex_dof_index(
    const unsigned int    vertex,
    const unsigned int    i,
    const types::fe_index fe_index = numbers::invalid_fe_index) const;

  /**
   * Return the global DoF index of the <code>i</code>th degree of freedom
   * associated with the <code>vertex</code>th vertex on level @p level. Also
   * see vertex_dof_index().
   */
  types::global_dof_index
  mg_vertex_dof_index(
    const int             level,
    const unsigned int    vertex,
    const unsigned int    i,
    const types::fe_index fe_index = numbers::invalid_fe_index) const;

  /**
   * Index of the <i>i</i>th degree of freedom of this object.
   *
   * The last argument denotes the finite element index. For the standard
   * ::DoFHandler class, this value must be equal to its default value since
   * that class only supports the same finite element on all cells anyway.
   *
   * However, when hp-capabilities are enabled, different finite element
   * objects may be used on different cells. On faces between two cells, as
   * well as vertices, there may therefore be two sets of degrees of freedom,
   * one for each of the finite elements used on the adjacent cells.  In order
   * to specify which set of degrees of freedom to work on, the last argument
   * is used to disambiguate. Finally, if this function is called for a cell
   * object, there can only be a single set of degrees of freedom, and
   * fe_index has to match the result of active_fe_index().
   *
   * @note While the get_dof_indices() function returns an array that contains
   * the indices of all degrees of freedom that somehow live on this object
   * (i.e. on the vertices, edges or interior of this object), the current
   * dof_index() function only considers the DoFs that really belong to this
   * particular object's interior. In other words, as an example, if the
   * current object refers to a quad (a cell in 2d, a face in 3d) and the
   * finite element associated with it is a bilinear one, then the
   * get_dof_indices() will return an array of size 4 while dof_index() will
   * produce an exception because no degrees are defined in the interior of
   * the face.
   */
  types::global_dof_index
  dof_index(const unsigned int    i,
            const types::fe_index fe_index = numbers::invalid_fe_index) const;

  /**
   * Return the dof_index on the given level. Also see dof_index.
   */
  types::global_dof_index
  mg_dof_index(const int level, const unsigned int i) const;

  /**
   * @}
   */

  /**
   * @name Accessing the finite element associated with this object
   */
  /**
   * @{
   */

  /**
   * Return the number of finite elements that are active on a given object.
   *
   * When hp-capabilities are disabled the answer is, of course, always one.
   * However, when hp-capabilities are enabled, this isn't the case: If this
   * is a cell, the answer is of course one. If it is a face, the answer may
   * be one or two, depending on whether the two adjacent cells use the same
   * finite element or not. If it is an edge in 3d, the possible return value
   * may be one or any other value larger than that.
   */
  unsigned int
  n_active_fe_indices() const;

  /**
   * Return the @p n-th active FE index on this object. For cells and all
   * non-hp-objects, there is only a single active FE index, so the argument
   * must be equal to zero. For lower-dimensional hp-objects, there are
   * n_active_fe_indices() active finite elements, and this function can be
   * queried for their indices.
   */
  types::fe_index
  nth_active_fe_index(const unsigned int n) const;

  /**
   * Returns all active FE indices on this object.
   *
   * The size of the returned set equals the number of finite elements that
   * are active on this object.
   */
  std::set<types::fe_index>
  get_active_fe_indices() const;

  /**
   * Return true if the finite element with given index is active on the
   * present object. When the current DoFHandler does not have
   * hp-capabilities, this is of course the case only if @p fe_index equals
   * zero. For cells, it is the case if @p fe_index equals active_fe_index()
   * of this cell. For faces and other lower- dimensional objects, there may
   * be more than one @p fe_index that are active on any given object (see
   * n_active_fe_indices()).
   */
  bool
  fe_index_is_active(const types::fe_index fe_index) const;

  /**
   * Return a reference to the finite element used on this object with the
   * given @p fe_index. @p fe_index must be used on this object, i.e.
   * <code>fe_index_is_active(fe_index)</code> must return true.
   */
  const FiniteElement<dim, spacedim> &
  get_fe(const types::fe_index fe_index) const;

  /**
   * @}
   */

  /**
   * Exceptions for child classes
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcInvalidObject,
                   "This accessor object has not been "
                   "associated with any DoFHandler object.");
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcVectorNotEmpty);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcVectorDoesNotMatch);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcMatrixDoesNotMatch);
  /**
   * A function has been called for a cell which should be
   * @ref GlossActive "active",
   * but is refined.
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcNotActive);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcCantCompareIterators);

protected:
  /**
   * Store the address of the DoFHandler object to be accessed.
   */
  DoFHandler<dim, spacedim> *dof_handler;

public:
  /**
   * Compare for equality. Return <tt>true</tt> if the two accessors refer to
   * the same object.
   *
   * The template parameters of this function allow for a comparison of very
   * different objects. Therefore, some of them are disabled. Namely, if the
   * dimension, or the dof handler of the two objects differ, an exception is
   * generated. It can be expected that this is an unwanted comparison.
   *
   * The template parameter <tt>level_dof_access2</tt> is ignored, such that
   * an iterator with level access can be equal to one with access to the
   * active degrees of freedom.
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  bool
  operator==(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &) const;

  /**
   * Compare for inequality. The boolean not of operator==().
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  bool
  operator!=(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &) const;

protected:
  /**
   * Reset the DoF handler pointer.
   */
  void
  set_dof_handler(DoFHandler<dim, spacedim> *dh);

  /**
   * Set the index of the <i>i</i>th degree of freedom of this object to @p
   * index.
   *
   * The last argument denotes the finite element index. For the standard
   * ::DoFHandler class, this value must be equal to its default value since
   * that class only supports the same finite element on all cells anyway.
   *
   * However, when the relevant DoFHandler has hp-capabilities, different
   * finite element objects may be used on different cells. On faces between
   * two cells, as well as vertices, there may therefore be two sets of
   * degrees of freedom, one for each of the finite elements used on the
   * adjacent cells.  In order to specify which set of degrees of freedom to
   * work on, the last argument is used to disambiguate. Finally, if this
   * function is called for a cell object, there can only be a single set of
   * degrees of freedom, and fe_index has to match the result of
   * active_fe_index().
   */
  void
  set_dof_index(
    const unsigned int            i,
    const types::global_dof_index index,
    const types::fe_index         fe_index = numbers::invalid_fe_index) const;

  void
  set_mg_dof_index(const int                     level,
                   const unsigned int            i,
                   const types::global_dof_index index) const;

  void
  set_mg_vertex_dof_index(
    const int                     level,
    const unsigned int            vertex,
    const unsigned int            i,
    const types::global_dof_index index,
    const types::fe_index         fe_index = numbers::invalid_fe_index) const;

  // Iterator classes need to be friends because they need to access
  // operator== and operator!=.
  template <typename>
  friend class TriaRawIterator;
  template <int, int, int, bool>
  friend class DoFAccessor;

private:
  // Make the DoFHandler class a friend so that it can call the set_xxx()
  // functions.
  friend class DoFHandler<dim, spacedim>;

  friend struct dealii::internal::DoFHandlerImplementation::Policy::
    Implementation;
  friend struct dealii::internal::DoFHandlerImplementation::Implementation;
  friend struct dealii::internal::hp::DoFHandlerImplementation::Implementation;
  friend struct dealii::internal::DoFCellAccessorImplementation::Implementation;
  friend struct dealii::internal::DoFAccessorImplementation::Implementation;
};



/**
 * Specialization of the general DoFAccessor class template for the case of
 * zero-dimensional objects (a vertex) that are the face of a one-dimensional
 * cell in spacedim space dimensions. Since vertices function differently than
 * general faces, this class does a few things differently than the general
 * template, but the interface should look the same.
 */
template <int spacedim, bool level_dof_access>
class DoFAccessor<0, 1, spacedim, level_dof_access>
  : public TriaAccessor<0, 1, spacedim>
{
public:
  /**
   * A static variable that allows users of this class to discover the value
   * of the second template argument.
   */
  static constexpr unsigned int dimension = 1;

  /**
   * A static variable that allows users of this class to discover the value
   * of the third template argument.
   */
  static constexpr unsigned int space_dimension = spacedim;

  /**
   * Declare an alias to the base class to make accessing some of the
   * exception classes simpler.
   */
  using BaseClass = TriaAccessor<0, 1, spacedim>;

  /**
   * Data type passed by the iterator class.
   */
  using AccessorData = DoFHandler<1, spacedim>;

  /**
   * @name Constructors
   */
  /**
   * @{
   */

  /**
   * Default constructor. Provides an accessor that can't be used.
   */
  DoFAccessor();

  /**
   * Constructor to be used if the object here refers to a vertex of a
   * one-dimensional triangulation, i.e. a face of the triangulation.
   *
   * Since there is no mapping from vertices to cells, an accessor object for
   * a point has no way to figure out whether it is at the boundary of the
   * domain or not. Consequently, the second argument must be passed by the
   * object that generates this accessor -- e.g. a 1d cell that can figure out
   * whether its left or right vertex are at the boundary.
   *
   * The third argument is the global index of the vertex we point to.
   *
   * The fourth argument is a pointer to the DoFHandler object.
   *
   * This iterator can only be called for one-dimensional triangulations.
   */
  DoFAccessor(
    const Triangulation<1, spacedim>                       *tria,
    const typename TriaAccessor<0, 1, spacedim>::VertexKind vertex_kind,
    const unsigned int                                      vertex_index,
    const DoFHandler<1, spacedim>                          *dof_handler);

  /**
   * Constructor. This constructor exists in order to maintain interface
   * compatibility with the other accessor classes. However, it doesn't do
   * anything useful here and so may not actually be called except to
   * default-construct iterator objects.
   */
  DoFAccessor(const Triangulation<1, spacedim> *,
              const int                                  = 0,
              const int                                  = 0,
              const DoFHandler<1, spacedim> *dof_handler = 0);

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
  DoFAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * Another conversion operator between objects that don't make sense, just
   * like the previous one.
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  DoFAccessor(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &);

  /**
   * Copy constructor.
   */
  DoFAccessor(const DoFAccessor<0, 1, spacedim, level_dof_access> &) = default;

  /**
   * Move constructor.
   */
  // NOLINTNEXTLINE OSX does not compile with noexcept
  DoFAccessor(DoFAccessor<0, 1, spacedim, level_dof_access> &&) = default;

  /**
   * Destructor.
   */
  ~DoFAccessor() = default;

  /**
   * Copy operator. These operators are usually used in a context like
   * <tt>iterator a,b; *a=*b;</tt>. Presumably, the intent here is to copy the
   * object pointed to
   * by @p b to the object pointed to by @p a. However, the result of
   * dereferencing an iterator is not an object but an accessor; consequently,
   * this operation is not useful for iterators on DoF handler objects.
   * Consequently, this operator is declared as deleted and can not be used.
   */
  DoFAccessor<0, 1, spacedim, level_dof_access> &
  operator=(const DoFAccessor<0, 1, spacedim, level_dof_access> &da) = delete;

  /**
   * Move assignment operator.
   */
  DoFAccessor<0, 1, spacedim, level_dof_access> &
  operator=(DoFAccessor<0, 1, spacedim, level_dof_access> &&) noexcept =
    default;

  /**
   * @}
   */

  /**
   * Return a handle on the DoFHandler object which we are using.
   */
  const DoFHandler<1, spacedim> &
  get_dof_handler() const;

  /**
   * Implement the copy operator needed for the iterator classes.
   */
  template <bool level_dof_access2>
  void
  copy_from(const DoFAccessor<0, 1, spacedim, level_dof_access2> &a);

  /**
   * Copy operator used by the iterator class. Keeps the previously set dof
   * handler, but sets the object coordinates of the TriaAccessor.
   */
  void
  copy_from(const TriaAccessorBase<0, 1, spacedim> &da);

  /**
   * @name Accessing sub-objects
   */
  /**
   * @{
   */

  /**
   * Return an invalid iterator of a type that represents pointing to a child
   * of the current object. The object is invalid because points (as
   * represented by the current class) do not have children.
   */
  TriaIterator<DoFAccessor<0, 1, spacedim, level_dof_access>>
  child(const unsigned int c) const;

  /**
   * Pointer to the @p ith line bounding this object.
   *
   * Since meshes with dimension 1 do not have quads this method just throws
   * an exception.
   */
  typename dealii::internal::DoFHandlerImplementation::
    Iterators<1, spacedim, level_dof_access>::line_iterator
    line(const unsigned int i) const;

  /**
   * Pointer to the @p ith quad bounding this object.
   *
   * Since meshes with dimension 1 do not have quads this method just throws
   * an exception.
   */
  typename dealii::internal::DoFHandlerImplementation::
    Iterators<1, spacedim, level_dof_access>::quad_iterator
    quad(const unsigned int i) const;

  /**
   * @}
   */

  /**
   * @name Accessing the DoF indices of this object
   */
  /**
   * @{
   */

  /**
   * Return the <i>global</i> indices of the degrees of freedom located on
   * this object in the standard ordering defined by the finite element. This
   * function is only available on <i>active</i> objects (see
   * @ref GlossActive "this glossary entry").
   *
   * The present vertex must belong to an active cell (and not artificial in a
   * parallel distributed computation).
   *
   * The vector has to have the right size before being passed to this
   * function.
   *
   * The last argument denotes the finite element index. For the standard
   * ::DoFHandler class, this value must be equal to its default value since
   * that class only supports the same finite element on all cells anyway.
   *
   * However, when the relevant DoFHandler has hp-capabilities, different
   * finite element objects may be used on different cells. On faces between
   * two cells, as well as vertices, there may therefore be two sets of
   * degrees of freedom, one for each of the finite elements used on the
   * adjacent cells.  In order to specify which set of degrees of freedom to
   * work on, the last argument is used to disambiguate. Finally, if this
   * function is called for a cell object, there can only be a single set of
   * degrees of freedom, and fe_index has to match the result of
   * active_fe_index().
   *
   * For cells, there is only a single possible finite element index (namely
   * the one for that cell, returned by <code>cell-@>active_fe_index</code>.
   * Consequently, the derived DoFCellAccessor class has an overloaded version
   * of this function that calls the present function with
   * <code>cell-@>active_fe_index</code> as last argument.
   */
  void
  get_dof_indices(
    std::vector<types::global_dof_index> &dof_indices,
    const types::fe_index fe_index = numbers::invalid_fe_index) const;

  /**
   * Return the global multilevel indices of the degrees of freedom that live
   * on the current object with respect to the given level within the
   * multigrid hierarchy. The indices refer to the local numbering for the
   * level this line lives on.
   */
  void
  get_mg_dof_indices(
    const int                             level,
    std::vector<types::global_dof_index> &dof_indices,
    const types::fe_index fe_index = numbers::invalid_fe_index) const;

  /**
   * Global DoF index of the <i>i</i> degree associated with the @p vertexth
   * vertex of the present cell.
   *
   * The last argument denotes the finite element index. For the standard
   * ::DoFHandler class, this value must be equal to its default value since
   * that class only supports the same finite element on all cells anyway.
   *
   * However, when the relevant DoFHandler has hp-capabilities, different
   * finite element objects may be used on different cells. On faces between
   * two cells, as well as vertices, there may therefore be two sets of
   * degrees of freedom, one for each of the finite elements used on the
   * adjacent cells.  In order to specify which set of degrees of freedom to
   * work on, the last argument is used to disambiguate. Finally, if this
   * function is called for a cell object, there can only be a single set of
   * degrees of freedom, and fe_index has to match the result of
   * active_fe_index().
   */
  types::global_dof_index
  vertex_dof_index(
    const unsigned int    vertex,
    const unsigned int    i,
    const types::fe_index fe_index = numbers::invalid_fe_index) const;

  /**
   * Index of the <i>i</i>th degree of freedom of this object.
   *
   * The last argument denotes the finite element index. For the standard
   * ::DoFHandler class, this value must be equal to its default value since
   * that class only supports the same finite element on all cells anyway.
   *
   * However, when the relevant DoFHandler has hp-capabilities, different
   * finite element objects may be used on different cells. On faces between
   * two cells, as well as vertices, there may therefore be two sets of
   * degrees of freedom, one for each of the finite elements used on the
   * adjacent cells.  In order to specify which set of degrees of freedom to
   * work on, the last argument is used to disambiguate. Finally, if this
   * function is called for a cell object, there can only be a single set of
   * degrees of freedom, and fe_index has to match the result of
   * active_fe_index().
   */
  types::global_dof_index
  dof_index(const unsigned int    i,
            const types::fe_index fe_index = numbers::invalid_fe_index) const;

  /**
   * @}
   */

  /**
   * @name Accessing the finite element associated with this object
   */
  /**
   * @{
   */

  /**
   * Return the number of finite elements that are active on a given object.
   *
   * Since vertices do not store the information necessary for this to be
   * calculated, this method just raises an exception and only exists to
   * enable dimension-independent programming.
   */
  unsigned int
  n_active_fe_indices() const;

  /**
   * Return the @p n-th active FE index on this object.
   *
   * Since vertices do not store the information necessary for this to be
   * calculated, this method just raises an exception and only exists to
   * enable dimension-independent programming.
   */
  types::fe_index
  nth_active_fe_index(const unsigned int n) const;

  /**
   * Return true if the finite element with given index is active on the
   * present object.
   *
   * Since vertices do not store the information necessary for this to be
   * calculated, this method just raises an exception and only exists to
   * enable dimension-independent programming.
   */
  bool
  fe_index_is_active(const types::fe_index fe_index) const;

  /**
   * Return a reference to the finite element used on this object with the
   * given @p fe_index. @p fe_index must be used on this object, i.e.
   * <code>fe_index_is_active(fe_index)</code> must return true.
   */
  const FiniteElement<1, spacedim> &
  get_fe(const types::fe_index fe_index) const;

  /**
   * @}
   */

  /**
   * Exceptions for child classes
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcInvalidObject,
                   "This accessor object has not been "
                   "associated with any DoFHandler object.");
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcVectorNotEmpty);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcVectorDoesNotMatch);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcMatrixDoesNotMatch);
  /**
   * A function has been called for a cell which should be
   * @ref GlossActive "active",
   * but is refined.
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcNotActive);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcCantCompareIterators);

protected:
  /**
   * Store the address of the DoFHandler object to be accessed.
   */
  DoFHandler<1, spacedim> *dof_handler;

  /**
   * Compare for equality.
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  bool
  operator==(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &) const;

  /**
   * Compare for inequality.
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  bool
  operator!=(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &) const;

  /**
   * Reset the DoF handler pointer.
   */
  void
  set_dof_handler(DoFHandler<1, spacedim> *dh);

  /**
   * Set the index of the <i>i</i>th degree of freedom of this object to @p
   * index.
   *
   * The last argument denotes the finite element index. For the standard
   * ::DoFHandler class, this value must be equal to its default value since
   * that class only supports the same finite element on all cells anyway.
   *
   * However, when the relevant DoFHandler has hp-capabilities, different
   * finite element objects may be used on different cells. On faces between
   * two cells, as well as vertices, there may therefore be two sets of
   * degrees of freedom, one for each of the finite elements used on the
   * adjacent cells.  In order to specify which set of degrees of freedom to
   * work on, the last argument is used to disambiguate. Finally, if this
   * function is called for a cell object, there can only be a single set of
   * degrees of freedom, and fe_index has to match the result of
   * active_fe_index().
   */
  void
  set_dof_index(
    const unsigned int            i,
    const types::global_dof_index index,
    const types::fe_index         fe_index = numbers::invalid_fe_index) const;

  // Iterator classes need to be friends because they need to access
  // operator== and operator!=.
  template <typename>
  friend class TriaRawIterator;


  // Make the DoFHandler class a friend so that it can call the set_xxx()
  // functions.
  friend class DoFHandler<1, spacedim>;

  friend struct dealii::internal::DoFHandlerImplementation::Policy::
    Implementation;
  friend struct dealii::internal::DoFHandlerImplementation::Implementation;
  friend struct dealii::internal::hp::DoFHandlerImplementation::Implementation;
  friend struct dealii::internal::DoFCellAccessorImplementation::Implementation;
};



/* -------------------------------------------------------------------------- */


/**
 * A class that represents DoF accessor objects to iterators that don't make
 * sense such as quad iterators in on 1d meshes.  This class can not be used to
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
class DoFInvalidAccessor : public InvalidAccessor<structdim, dim, spacedim>
{
public:
  /**
   * Propagate alias from base class to this class.
   */
  using AccessorData =
    typename InvalidAccessor<structdim, dim, spacedim>::AccessorData;

  /**
   * Constructor.  This class is used for iterators that do not make
   * sense in a given dimension, for example quads for 1d meshes. Consequently,
   * while the creation of such objects is syntactically valid, they make no
   * semantic sense, and we generate an exception when such an object is
   * actually generated.
   */
  DoFInvalidAccessor(const void         *parent     = nullptr,
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
  DoFInvalidAccessor(const DoFInvalidAccessor<structdim, dim, spacedim> &);

  /**
   * Conversion from other accessors to the current invalid one. This of
   * course also leads to a run-time error.
   */
  template <typename OtherAccessor>
  DoFInvalidAccessor(const OtherAccessor &);

  /**
   * Return the index of the <i>i</i>th degree of freedom of this object to
   * @p index. Since the current object doesn't point to anything useful, like
   * all other functions in this class this function only throws an exception.
   */
  types::global_dof_index
  dof_index(const unsigned int    i,
            const types::fe_index fe_index =
              DoFHandler<dim, spacedim>::default_fe_index) const;

  /**
   * Set the index of the <i>i</i>th degree of freedom of this object to @p
   * index. Since the current object doesn't point to anything useful, like
   * all other functions in this class this function only throws an exception.
   */
  void
  set_dof_index(
    const unsigned int            i,
    const types::global_dof_index index,
    const types::fe_index         fe_index = numbers::invalid_fe_index) const;
};



/* -------------------------------------------------------------------------- */


/**
 * Grant access to the degrees of freedom on a cell.
 *
 * Note that since for the class we derive from, i.e.
 * <tt>DoFAccessor<dim></tt>, the two template parameters are equal, the base
 * class is actually derived from CellAccessor, which makes the functions of
 * this class available to the DoFCellAccessor class as well.
 *
 * @ingroup dofs
 * @ingroup Accessors
 */
template <int dimension_, int space_dimension_, bool level_dof_access>
class DoFCellAccessor : public DoFAccessor<dimension_,
                                           dimension_,
                                           space_dimension_,
                                           level_dof_access>
{
public:
  /**
   * Extract dimension from DoFHandler.
   */
  static const unsigned int dim = dimension_;

  /**
   * Extract space dimension from DoFHandler.
   */
  static const unsigned int spacedim = space_dimension_;


  /**
   * Data type passed by the iterator class.
   */
  using AccessorData = DoFHandler<dimension_, space_dimension_>;

  /**
   * Declare an alias to the base class to make accessing some of the
   * exception classes simpler.
   */
  using BaseClass =
    DoFAccessor<dimension_, dimension_, space_dimension_, level_dof_access>;

  /**
   * Define the type of the container this is part of.
   */
  using Container = DoFHandler<dimension_, space_dimension_>;

  /**
   * A type for an iterator over the faces of a cell. This is what the face()
   * function returns.
   */
  using face_iterator = TriaIterator<DoFAccessor<dimension_ - 1,
                                                 dimension_,
                                                 space_dimension_,
                                                 level_dof_access>>;

  /**
   * @name Constructors and initialization
   */
  /**
   * @{
   */

  /**
   * Constructor
   */
  DoFCellAccessor(const Triangulation<dimension_, space_dimension_> *tria,
                  const int                                          level,
                  const int                                          index,
                  const AccessorData *local_data);

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
  DoFCellAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * Another conversion operator between objects that don't make sense, just
   * like the previous one.
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  explicit DoFCellAccessor(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &);

  /**
   * Copy constructor.
   */
  DoFCellAccessor(
    const DoFCellAccessor<dimension_, space_dimension_, level_dof_access> &) =
    default;

  /**
   * Move constructor.
   */
  DoFCellAccessor(                                                  // NOLINT
    DoFCellAccessor<dimension_, space_dimension_, level_dof_access> // NOLINT
      &&) = default;                                                // NOLINT

  /**
   * Destructor
   */
  ~DoFCellAccessor() = default;

  /**
   * Copy operator. These operators are usually used in a context like
   * <tt>iterator a,b; *a=*b;</tt>. Presumably, the intent here is to copy the
   * object pointed to
   * by @p b to the object pointed to by @p a. However, the result of
   * dereferencing an iterator is not an object but an accessor; consequently,
   * this operation is not useful for iterators on DoF handler objects.
   * Consequently, this operator is declared as deleted and can not be used.
   */
  DoFCellAccessor<dimension_, space_dimension_, level_dof_access> &
  operator=(
    const DoFCellAccessor<dimension_, space_dimension_, level_dof_access> &da) =
    delete;

  /**
   * Move assignment operator.
   */
  DoFCellAccessor<dimension_, space_dimension_, level_dof_access> & // NOLINT
  operator=(                                                        // NOLINT
    DoFCellAccessor<dimension_, space_dimension_, level_dof_access> // NOLINT
      &&) = default;                                                // NOLINT

  /**
   * @}
   */

  /**
   * Return the parent of this cell as a DoF cell iterator. If the parent does
   * not exist (i.e., if the object is at the coarsest level of the mesh
   * hierarchy), an exception is generated.
   *
   * This function is needed since the parent function of the base class
   * CellAccessor returns a triangulation cell accessor without access to the
   * DoF data.
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  parent() const;

  /**
   * @name Accessing sub-objects and neighbors
   */
  /**
   * @{
   */

  /**
   * Return the @p ith neighbor as a DoF cell iterator. This function is
   * needed since the neighbor function of the base class returns a cell
   * accessor without access to the DoF data.
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  neighbor(const unsigned int i) const;

  /**
   * Return the @p ith periodic neighbor as a DoF cell iterator. This function
   * is needed since the neighbor function of the base class returns a cell
   * accessor without access to the DoF data.
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  periodic_neighbor(const unsigned int i) const;

  /**
   * Return the @p ith neighbor or periodic neighbor as a DoF cell iterator.
   * This function is needed since the neighbor function of the base class
   * returns a cell accessor without access to the DoF data.
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  neighbor_or_periodic_neighbor(const unsigned int i) const;

  /**
   * Return the @p ith child as a DoF cell iterator. This function is needed
   * since the child function of the base class returns a cell accessor
   * without access to the DoF data.
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  child(const unsigned int i) const;

  /**
   * Return an array of iterators to all children of this cell.
   */
  boost::container::small_vector<
    TriaIterator<
      DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>,
    GeometryInfo<dimension_>::max_children_per_cell>
  child_iterators() const;

  /**
   * Return an iterator to the @p ith face of this cell.
   *
   * This function returns a DoFAccessor with <code>structdim == 0</code> in
   * 1d, a DoFAccessor::line in 2d, and a DoFAccessor::quad in 3d.
   */
  face_iterator
  face(const unsigned int i) const;

  /**
   * Return an array of iterators to all faces of this cell.
   */
  boost::container::small_vector<face_iterator,
                                 GeometryInfo<dimension_>::faces_per_cell>
  face_iterators() const;

  /**
   * Return the result of the @p neighbor_child_on_subface function of the
   * base class, but convert it so that one can also access the DoF data (the
   * function in the base class only returns an iterator with access to the
   * triangulation data).
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  neighbor_child_on_subface(const unsigned int face_no,
                            const unsigned int subface_no) const;

  /**
   * Return the result of the @p periodic_neighbor_child_on_subface function
   * of the base class, but convert it so that one can also access the DoF
   * data (the function in the base class only returns an iterator with access
   * to the triangulation data).
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  periodic_neighbor_child_on_subface(const unsigned int face_no,
                                     const unsigned int subface_no) const;

  /**
   * @}
   */

  /**
   * @name Extracting values from global vectors
   */
  /**
   * @{
   */

  /**
   * Collect the values of the given vector restricted to the dofs of this cell
   * in the standard ordering: dofs on vertex 0, dofs on vertex 1, etc, dofs
   * on line 0, dofs on line 1, etc, dofs on quad 0, etc. In other
   * words, this function implements a
   * [gather
   * operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing)).
   *
   * The vector has to have the right size before being passed to this
   * function. This function is only callable for active cells.
   *
   * The input vector may be either a <tt>Vector<float></tt>, Vector<double>,
   * or a BlockVector<double>, or a PETSc or Trilinos vector if deal.II is
   * compiled to support these libraries. It is in the responsibility of the
   * caller to assure that the types of the numbers stored in input and output
   * vectors are compatible and with similar accuracy.
   */
  template <class InputVector, typename number>
  void
  get_dof_values(const InputVector &values, Vector<number> &local_values) const;

  /**
   * Collect the values of the given vector restricted to the dofs of this cell
   * in the standard ordering: dofs on vertex 0, dofs on vertex 1, etc, dofs
   * on line 0, dofs on line 1, etc, dofs on quad 0, etc. In other
   * words, this function implements a
   * [gather
   * operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing)).
   *
   * The vector has to have the right size before being passed to this
   * function. This function is only callable for active cells.
   *
   * The input vector may be either a <tt>Vector<float></tt>, Vector<double>,
   * or a BlockVector<double>, or a PETSc or Trilinos vector if deal.II is
   * compiled to support these libraries. It is in the responsibility of the
   * caller to assure that the types of the numbers stored in input and output
   * vectors are compatible and with similar accuracy.
   */
  template <typename Number, typename ForwardIterator>
  void
  get_dof_values(const ReadVector<Number> &values,
                 ForwardIterator           local_values_begin,
                 ForwardIterator           local_values_end) const;

  /**
   * Collect the values of the given vector restricted to the dofs of this cell
   * in the standard ordering: dofs on vertex 0, dofs on vertex 1, etc, dofs
   * on line 0, dofs on line 1, etc, dofs on quad 0, etc. In other
   * words, this function implements a
   * [gather
   * operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing)).
   *
   * The vector has to have the right size before being passed to this
   * function. This function is only callable for active cells.
   *
   * The input vector may be either a <tt>Vector<float></tt>, Vector<double>,
   * or a BlockVector<double>, or a PETSc or Trilinos vector if deal.II is
   * compiled to support these libraries. It is in the responsibility of the
   * caller to assure that the types of the numbers stored in input and output
   * vectors are compatible and with similar accuracy. The
   * AffineConstraints object passed as an argument to this function makes
   * sure that constraints are correctly distributed when the dof values
   * are calculated.
   */
  template <class InputVector, typename ForwardIterator>
  void
  get_dof_values(
    const AffineConstraints<typename InputVector::value_type> &constraints,
    const InputVector                                         &values,
    ForwardIterator local_values_begin,
    ForwardIterator local_values_end) const;

  /**
   * This function is the counterpart to get_dof_values(): it takes a vector
   * of values for the degrees of freedom of the cell pointed to by this
   * iterator and writes these values into the global data vector @p values.
   * In other words, this function implements a
   * [scatter
   * operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing)).
   * This function is only callable for active cells.
   *
   * Note that for continuous finite elements, calling this function affects
   * the dof values on neighboring cells as well. It may also violate
   * continuity requirements for hanging nodes, if neighboring cells are less
   * refined than the present one. These requirements are not taken care of
   * and must be enforced by the user afterwards.
   *
   * The vector has to have the right size before being passed to this
   * function.
   *
   * The output vector may be either a Vector<float>, Vector<double>, or a
   * BlockVector<double>, or a PETSc vector if deal.II is compiled to support
   * these libraries. It is in the responsibility of the caller to assure that
   * the types of the numbers stored in input and output vectors are
   * compatible and with similar accuracy.
   */
  template <class OutputVector, typename number>
  void
  set_dof_values(const Vector<number> &local_values,
                 OutputVector         &values) const;

  /**
   * Return the interpolation of the given finite element function to the
   * present cell. In the simplest case, the cell is a terminal one, i.e., it
   * has no children; then, the returned value is the vector of nodal values
   * on that cell. You could as well get the desired values through the @p
   * get_dof_values function. In the other case, when the cell has children,
   * we use the restriction matrices provided by the finite element class to
   * compute the interpolation from the children to the present cell.
   *
   * If the cell is part of a DoFHandler with hp-capabilities, cells only have
   * an associated finite element space if they are active. However, this
   * function is supposed to also provide information on inactive cells with
   * children. Consequently, it carries a third argument that can be used in
   * the hp-context that denotes the finite element space we are supposed to
   * interpolate onto. If the cell is active, this function then obtains the
   * finite element function from the <code>values</code> vector on this cell
   * and interpolates it onto the space described by the
   * <code>fe_index</code>th element of the hp::FECollection associated with
   * the DoFHandler of which this cell is a part of. If the cell is not
   * active, then we first perform this interpolation on all of its terminal
   * children and then interpolate this function down to the cell requested
   * keeping the function space the same.
   *
   * It is assumed that both input vectors already have the right size
   * beforehand.
   *
   * @note Unlike the get_dof_values() function, this function is only
   * available on cells, rather than on lines, quads, and hexes, since
   * interpolation is presently only provided for cells by the finite element
   * classes.
   */
  template <typename Number>
  void
  get_interpolated_dof_values(
    const ReadVector<Number> &values,
    Vector<Number>           &interpolated_values,
    const types::fe_index     fe_index = numbers::invalid_fe_index) const;

  /**
   * This function is the counterpart to get_interpolated_dof_values(): you
   * specify the dof values on a cell and these are interpolated to the
   * children of the present cell and set on the terminal cells.
   *
   * In principle, it works as follows: if the cell pointed to by this object
   * is terminal (i.e., has no children), then the dof values are set in the
   * global data vector by calling the set_dof_values() function; otherwise,
   * the values are prolonged to each of the children and this function is
   * called for each of them.
   *
   * Using the get_interpolated_dof_values() and this function, you can
   * compute the interpolation of a finite element function to a coarser grid
   * by first getting the interpolated solution on a cell of the coarse grid
   * and afterwards redistributing it using this function.
   *
   * Note that for continuous finite elements, calling this function affects
   * the dof values on neighboring cells as well. It may also violate
   * continuity requirements for hanging nodes, if neighboring cells are less
   * refined than the present one, or if their children are less refined than
   * the children of this cell. These requirements are not taken care of and
   * must be enforced by the user afterward.
   *
   * If the cell is part of a DoFHandler with hp-capabilities, cells only have
   * an associated finite element space if they are active. However, this
   * function is supposed to also work on inactive cells with children.
   * Consequently, it carries a third argument that can be used in the
   * hp-context that denotes the finite element space we are supposed to
   * interpret the input vector of this function in. If the cell is active,
   * this function then interpolates the input vector interpreted as an
   * element of the space described by the <code>fe_index</code>th element of
   * the hp::FECollection associated with the DoFHandler of which this
   * cell is a part of, and interpolates it into the space that is associated
   * with this cell. On the other hand, if the cell is not active, then we
   * first perform this interpolation from this cell to its children using the
   * given <code>fe_index</code> until we end up on an active cell, at which
   * point we follow the procedure outlined at the beginning of the paragraph.
   *
   * It is assumed that both vectors already have the right size beforehand.
   * This function relies on the existence of a natural interpolation property
   * of finite element spaces of a cell to its children, denoted by the
   * prolongation matrices of finite element classes. For some elements, the
   * spaces on coarse and fine grids are not nested, in which case the
   * interpolation to a child is not the identity; refer to the documentation
   * of the respective finite element class for a description of what the
   * prolongation matrices represent in this case.
   *
   * @note Cells set the values of DoFs independently and might overwrite
   * previously set values in the global vector, for example when calling the
   * same function earlier from a different cell. By setting @p perform_check,
   * you can enable a check that the previous value and the one to be set here
   * are at least roughly the same. In practice, they might be slightly
   * different because they are computed in a way that theoretically ensures
   * that they are the same, but in practice they are only equal up to
   * round-off.
   *
   * @note Unlike the get_dof_values() function, this function is only
   * available on cells, rather than on lines, quads, and hexes, since
   * interpolation is presently only provided for cells by the finite element
   * classes.
   */
  template <class OutputVector, typename number>
  void
  set_dof_values_by_interpolation(
    const Vector<number> &local_values,
    OutputVector         &values,
    const types::fe_index fe_index      = numbers::invalid_fe_index,
    const bool            perform_check = false) const;

  /**
   * Similar to set_dof_values_by_interpolation() with the difference that
   * values are added into the vector.
   *
   * @note In parallel::distributed::SolutionTransfer, this function is used
   *   to accumulate the contributions of all cells to a DoF; with a
   *   subsequent multiplication with the inverse of the valence, finally,
   *   the average value is obtained.
   */
  template <class OutputVector, typename number>
  void
  distribute_local_to_global_by_interpolation(
    const Vector<number> &local_values,
    OutputVector         &values,
    const types::fe_index fe_index = numbers::invalid_fe_index) const;

  /**
   * Distribute a local (cell based) vector to a global one by mapping the
   * local numbering of the degrees of freedom to the global one and entering
   * the local values into the global vector. In other words, this function
   * implements a
   * [scatter
   * operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing)).
   *
   * The elements are <em>added</em> to the existing elements in the global
   * vector, rather than just set, since this is usually what one wants. You may
   * also want to take a look at the
   * AffineConstraints::distribute_local_to_global() function if you need to
   * deal with constraints.
   */
  template <typename number, typename OutputVector>
  void
  distribute_local_to_global(const Vector<number> &local_source,
                             OutputVector         &global_destination) const;

  /**
   * Distribute a local (cell based) vector in iterator format to a global one
   * by mapping the local numbering of the degrees of freedom to the global
   * one and entering the local values into the global vector.
   * In other words, this function implements a
   * [scatter
   * operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing)).
   *
   * The elements are <em>added</em> to the existing elements in the global
   * vector, rather than just set, since this is usually what one wants. You may
   * also want to take a look at the
   * AffineConstraints::distribute_local_to_global() function if you need to
   * deal with constraints.
   */
  template <typename ForwardIterator, typename OutputVector>
  void
  distribute_local_to_global(ForwardIterator local_source_begin,
                             ForwardIterator local_source_end,
                             OutputVector   &global_destination) const;

  /**
   * Distribute a local (cell based) vector in iterator format to a global one
   * by mapping the local numbering of the degrees of freedom to the global
   * one and entering the local values into the global vector.
   * In other words, this function implements a
   * [scatter
   * operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing)).
   *
   * The elements are <em>added</em> up to the elements in the global vector,
   * rather than just set, since this is usually what one wants. Moreover, the
   * AffineConstraints object passed to this function makes sure that also
   * constraints are eliminated in this process.
   */
  template <typename ForwardIterator, typename OutputVector>
  void
  distribute_local_to_global(
    const AffineConstraints<typename OutputVector::value_type> &constraints,
    ForwardIterator local_source_begin,
    ForwardIterator local_source_end,
    OutputVector   &global_destination) const;

  /**
   * This function does much the same as the
   * <tt>distribute_local_to_global(Vector,Vector)</tt> function, but operates
   * on matrices instead of vectors. If the matrix type is a sparse matrix
   * then it is supposed to have non-zero entry slots where required.
   */
  template <typename number, typename OutputMatrix>
  void
  distribute_local_to_global(const FullMatrix<number> &local_source,
                             OutputMatrix &global_destination) const;

  /**
   * This function does what the two <tt>distribute_local_to_global</tt>
   * functions with vector and matrix argument do, but all at once.
   */
  template <typename number, typename OutputMatrix, typename OutputVector>
  void
  distribute_local_to_global(const FullMatrix<number> &local_matrix,
                             const Vector<number>     &local_vector,
                             OutputMatrix             &global_matrix,
                             OutputVector             &global_vector) const;

  /**
   * @}
   */

  /**
   * @name Accessing the DoF indices of this object
   */

  /**
   * @{
   */

  /**
   * Obtain the global indices of the local degrees of freedom on this cell.
   *
   * If this object accesses a level cell (indicated by the third template
   * argument or #is_level_cell), then return the result of
   * get_mg_dof_indices(), else return get_dof_indices().
   *
   * You will get a level_cell_iterator when calling begin_mg() and a normal
   * one otherwise.
   *
   * Examples for this use are in the implementation of DoFRenumbering.
   */
  void
  get_active_or_mg_dof_indices(
    std::vector<types::global_dof_index> &dof_indices) const;

  /**
   * Return the <i>global</i> indices of the degrees of freedom located on
   * this object in the standard ordering defined by the finite element (i.e.,
   * dofs on vertex 0, dofs on vertex 1, etc, dofs on line 0, dofs on line 1,
   * etc, dofs on quad 0, etc.) This function is only available on
   * <i>active</i> objects (see
   * @ref GlossActive "this glossary entry").
   *
   * @param[out] dof_indices The vector into which the indices will be
   * written. It has to have the right size (namely,
   * <code>fe.n_dofs_per_cell()</code>, <code>fe.dofs_per_face</code>, or
   * <code>fe.dofs_per_line</code>, depending on which kind of object this
   * function is called) before being passed to this function.
   *
   * This function reimplements the same function in the base class. In
   * contrast to the function in the base class, we do not need the
   * <code>fe_index</code> here because there is always a unique finite
   * element index on cells.
   *
   * This is a function which requires that the cell is active.
   *
   * Also see get_active_or_mg_dof_indices().
   *
   * @note In many places in the tutorial and elsewhere in the library, the
   * argument to this function is called <code>local_dof_indices</code> by
   * convention. The name is not meant to indicate the <i>local</i> numbers of
   * degrees of freedom (which are always between zero and
   * <code>fe.n_dofs_per_cell()</code>) but instead that the returned values are
   * the <i>global</i> indices of those degrees of freedom that are located
   * locally on the current cell.
   */
  void
  get_dof_indices(std::vector<types::global_dof_index> &dof_indices) const;

  /**
   * Retrieve the global indices of the degrees of freedom on this cell in the
   * level vector associated to the level of the cell.
   */
  void
  get_mg_dof_indices(std::vector<types::global_dof_index> &dof_indices) const;

  /**
   * @}
   */

  /**
   * @name Accessing the finite element associated with this object
   */
  /**
   * @{
   */

  /**
   * Return the finite element that is used on the cell pointed to by this
   * iterator. For DoFHandler objects without hp-capabilities, this is of
   * course always the same element, independent of the cell we are presently
   * on, but for hp-DoFHandler objects this may change from cell to cell.
   *
   * @note Since degrees of freedom only exist on active cells for DoFHandler
   * objects with hp-capabilities (i.e., there is currently no implementation
   * of multilevel such objects), it does not make sense to query the finite
   * element on non-active cells since they do not have finite element spaces
   * associated with them without having any degrees of freedom. Consequently,
   * this function will produce an exception when called on non-active cells.
   */
  const FiniteElement<dimension_, space_dimension_> &
  get_fe() const;

  /**
   * Return the index inside the hp::FECollection of the FiniteElement used
   * for this cell. This function is only useful if the DoFHandler object
   * associated with the current cell has hp-capabilities enabled.
   *
   * @note Since degrees of freedom only exist on active cells for DoFHandler
   * objects with hp-capabilities (i.e., there is currently no implementation
   * of multilevel such objects), it does not make sense to query the finite
   * element on non-active cells since they do not have finite element spaces
   * associated with them without having any degrees of freedom. Consequently,
   * this function will produce an exception when called on non-active cells.
   *
   * @note When using parallel meshes, either through the
   * parallel::shared::Triangulation or parallel::distributed::Triangulation
   * classes, it is only allowed to call this function on locally
   * owned or ghost cells. No information is available on artificial cells.
   * Furthermore, @p active_fe_index information is only exchanged from locally
   * owned cells on one processor to other processors where they may be
   * ghost cells, during the call to DoFHandler::set_fe() and
   * DoFHandler::distribute_dofs(). Be aware that if you call
   * set_active_fe_index() on a cell after calling one of these functions, then
   * this information will not be propagated to other processors who may have
   * this cell as a ghost cell. See the documentation of DoFHandler for more
   * information.
   */
  types::fe_index
  active_fe_index() const;

  /**
   * Set the index of the FiniteElement used for this cell. This determines
   * which element in an hp::FECollection to use. This function is only useful
   * if the DoF handler object associated with the current cell has
   * hp-capabilities enabled.
   *
   * @note Since degrees of freedom only exist on active cells for DoFHandler
   * objects with hp-capabilities (i.e., there is currently no implementation
   * of multilevel such objects), it does not make sense to query the finite
   * element on non-active cells since they do not have finite element spaces
   * associated with them without having any degrees of freedom. Consequently,
   * this function will produce an exception when called on non-active cells.
   *
   * @note When using parallel meshes, either through the
   * parallel::shared::Triangulation or parallel::distributed::Triangulation
   * classes, it is only allowed to call this function on locally
   * owned or ghost cells. No information is available on artificial cells.
   * Furthermore, @p active_fe_index information is only exchanged from locally
   * owned cells on one processor to other processors where they may be
   * ghost cells, during the call to DoFHandler::set_fe() and
   * DoFHandler::distribute_dofs(). Be aware that if you call
   * set_active_fe_index() on a cell after calling one of these functions, then
   * this information will not be propagated to other processors who may have
   * this cell as a ghost cell. See the documentation of DoFHandler for more
   * information.
   */
  void
  set_active_fe_index(const types::fe_index i) const;
  /**
   * @}
   */

  /**
   * Set the DoF indices of this cell to the given values. This function
   * bypasses the DoF cache, if one exists for the given DoF handler class.
   */
  void
  set_dof_indices(const std::vector<types::global_dof_index> &dof_indices);

  /**
   * Set the Level DoF indices of this cell to the given values.
   */
  void
  set_mg_dof_indices(const std::vector<types::global_dof_index> &dof_indices);

  /**
   * @name Dealing with refinement indicators
   */
  /**
   * @{
   */

  /**
   * Return the finite element that will be assigned to this cell next time the
   * triangulation gets refined and coarsened. If no future finite element has
   * been specified for this cell via the set_future_fe_index() function, the
   * active one will remain unchanged, in which case the active finite element
   * will be returned.
   *
   * For DoFHandlers without hp-capabilities enabled, this is of course always
   * the same element, independent of the cell we are presently on, but for
   * hp-DoFHandler objects this may change from cell to cell.
   *
   * @note Since degrees of freedom only exist on active cells for DoFHandler
   * objects with hp-capabilities (i.e., there is currently no implementation
   * of multilevel such objects), it does not make sense to query the finite
   * element on non-active cells since they do not have finite element spaces
   * associated with them without having any degrees of freedom. Consequently,
   * this function will produce an exception when called on non-active cells.
   */
  const FiniteElement<dimension_, space_dimension_> &
  get_future_fe() const;

  /**
   * Return the fe_index of the finite element that will be assigned to this
   * cell next time the triangulation gets refined and coarsened. If no future
   * finite element has been specified for this cell via the
   * set_future_fe_index() function, the active one will remain unchanged, in
   * which case the fe_index of the active finite element will be returned.
   *
   * @note Since degrees of freedom only exist on active cells for DoFHandler
   * objects with hp-capabilities (i.e., there is currently no implementation
   * of multilevel such objects), it does not make sense to query the finite
   * element on non-active cells since they do not have finite element spaces
   * associated with them without having any degrees of freedom. Consequently,
   * this function will produce an exception when called on non-active cells.
   *
   * @note When using parallel meshes, either through the
   * parallel::shared::Triangulation or parallel::distributed::Triangulation
   * classes, it is only allowed to call this function on locally owned cells.
   */
  types::fe_index
  future_fe_index() const;

  /**
   * Set the fe_index of the finite element that will be assigned to this
   * cell next time the triangulation gets refined and coarsened. A previously
   * assigned future finite element will be overwritten.
   *
   * See notes of future_fe_index() for information about restrictions on this
   * functionality.
   */
  void
  set_future_fe_index(const types::fe_index i) const;

  /**
   * Return whether a future finite element has been set.
   *
   * See notes of future_fe_index() for information about restrictions on this
   * functionality.
   */
  bool
  future_fe_index_set() const;

  /**
   * Revoke the future finite element assigned. Thus, the active finite element
   * will remain unchanged next time the triangulation gets refined and
   * coarsened.
   *
   * See notes on future_fe_index() for information about restrictions on this
   * functionality.
   */
  void
  clear_future_fe_index() const;
  /**
   * @}
   */

private:
  friend struct dealii::internal::DoFCellAccessorImplementation::Implementation;
};


template <int structdim, int dim, int spacedim, bool level_dof_access>
inline bool
DoFAccessor<structdim, dim, spacedim, level_dof_access>::is_level_cell()
{
  return level_dof_access;
}



template <int structdim, int dim, int spacedim>
template <typename OtherAccessor>
DoFInvalidAccessor<structdim, dim, spacedim>::DoFInvalidAccessor(
  const OtherAccessor &)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



/*------------------------- Functions: DoFAccessor ---------------------------*/


template <int structdim, int dim, int spacedim, bool level_dof_access>
inline DoFAccessor<structdim, dim, spacedim, level_dof_access>::DoFAccessor()
{
  Assert(false, ExcInvalidObject());
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline DoFAccessor<structdim, dim, spacedim, level_dof_access>::DoFAccessor(
  const Triangulation<dim, spacedim> *tria,
  const int                           level,
  const int                           index,
  const DoFHandler<dim, spacedim>    *dof_handler)
  : dealii::internal::DoFAccessorImplementation::
      Inheritance<structdim, dim, spacedim>::BaseClass(tria, level, index)
  , dof_handler(const_cast<DoFHandler<dim, spacedim> *>(dof_handler))
{
  Assert(
    tria == nullptr || &dof_handler->get_triangulation() == tria,
    ExcMessage(
      "You can't create a DoF accessor in which the DoFHandler object "
      "uses a different triangulation than the one you pass as argument."));
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2>
DoFAccessor<structdim, dim, spacedim, level_dof_access>::DoFAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false, ExcInvalidObject());
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline DoFAccessor<structdim, dim, spacedim, level_dof_access>::DoFAccessor(
  const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &other)
  : BaseClass(other)
  , dof_handler(nullptr)
{
  Assert(false,
         ExcMessage(
           "You are trying to assign iterators that are incompatible. "
           "The reason for incompatibility is that they refer to objects of "
           "different dimensionality (e.g., assigning a line iterator "
           "to a quad iterator)."));
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
template <bool level_dof_access2>
inline DoFAccessor<structdim, dim, spacedim, level_dof_access>::DoFAccessor(
  const DoFAccessor<structdim, dim, spacedim, level_dof_access2> &other)
  : BaseClass(other)
  , dof_handler(const_cast<DoFHandler<dim, spacedim> *>(other.dof_handler))
{}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::set_dof_handler(
  DoFHandler<dim, spacedim> *dh)
{
  Assert(dh != nullptr, ExcInvalidObject());
  this->dof_handler = dh;
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline const DoFHandler<dim, spacedim> &
DoFAccessor<structdim, dim, spacedim, level_dof_access>::get_dof_handler() const
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  return *this->dof_handler;
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::copy_from(
  const TriaAccessorBase<structdim, dim, spacedim> &da)
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  BaseClass::copy_from(da);
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
template <bool level_dof_access2>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::copy_from(
  const DoFAccessor<structdim, dim, spacedim, level_dof_access2> &a)
{
  BaseClass::copy_from(a);
  this->dof_handler = a.dof_handler;
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline bool
DoFAccessor<structdim, dim, spacedim, level_dof_access>::operator==(
  const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &a) const
{
  Assert(structdim == structdim2, ExcCantCompareIterators());
  Assert(this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator==(a));
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline bool
DoFAccessor<structdim, dim, spacedim, level_dof_access>::operator!=(
  const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &a) const
{
  Assert(structdim == structdim2, ExcCantCompareIterators());
  Assert(this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator!=(a));
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline TriaIterator<DoFAccessor<structdim, dim, spacedim, level_dof_access>>
DoFAccessor<structdim, dim, spacedim, level_dof_access>::child(
  const unsigned int i) const
{
  Assert(static_cast<unsigned int>(this->present_level) <
           this->dof_handler->object_dof_indices.size(),
         ExcMessage("DoFHandler not initialized"));

  TriaIterator<TriaAccessor<structdim, dim, spacedim>> t =
    TriaAccessor<structdim, dim, spacedim>::child(i);

  TriaIterator<DoFAccessor<structdim, dim, spacedim, level_dof_access>> q(
    *t, this->dof_handler);
  return q;
}


namespace internal
{
  namespace DoFAccessorImplementation
  {
    /**
     * Convert an FE index that might contain the right value but also
     * numbers::invalid_fe_index to a right value if needed/possible.
     */
    template <int structdim, int dim, int spacedim, bool level_dof_access>
    types::fe_index
    get_fe_index_or_default(
      const DoFAccessor<structdim, dim, spacedim, level_dof_access> &cell,
      const types::fe_index                                          fe_index)
    {
      if (cell.get_dof_handler().has_hp_capabilities() == false)
        {
          // No hp enabled, and the argument is at its default value -> we
          // can translate to the default active fe index
          Assert(
            (fe_index == numbers::invalid_fe_index) ||
              (fe_index == DoFHandler<dim, spacedim>::default_fe_index),
            ExcMessage(
              "It is not possible to specify a FE index if no hp support is used!"));

          return DoFHandler<dim, spacedim>::default_fe_index;
        }
      else
        {
          // Otherwise: If anything other than the default is provided by
          // the caller, then we should take just that. As an exception, if
          // we are on a cell (rather than a face/edge/vertex), then we know
          // that there is only one active fe index on this cell and we can
          // use that:
          if ((dim == structdim) && (fe_index == numbers::invalid_fe_index))
            {
              AssertDimension(cell.n_active_fe_indices(), 1);

              return cell.nth_active_fe_index(0);
            }

          Assert((fe_index != numbers::invalid_fe_index),
                 ExcMessage(
                   "You need to specify a FE index if hp support is used!"));

          return fe_index;
        }
    }

    /**
     * A class like the one with same name in tria.cc. See there for more
     * information.
     */
    struct Implementation
    {
      /**
       * In several applications of DoFAccessor::get_dof_values(), we want to
       * extract some indices without having to allocate memory. We do this by
       * setting a boost small_vector with 27 elements on the stack, and only
       * allocate when we exceed 27. The number 27 is heuristic and allows up
       * to quadratic shape functions on scalar problems in 3d, or linear
       * shape functions on systems (elasticity).
       */
      using dof_index_vector_type =
        boost::container::small_vector<dealii::types::global_dof_index, 27>;

      /**
       * Process the @p local_index-th degree of freedom corresponding to the
       * finite element specified by @p fe_index on the vertex with global
       * number @p vertex_index to @p global_index.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim,
                int spacedim,
                int structdim,
                typename GlobalIndexType,
                typename DoFPProcessor>
      static void
      process_dof_index(const DoFHandler<dim, spacedim> &dof_handler,
                        const unsigned int               obj_level,
                        const unsigned int               obj_index,
                        const types::fe_index            fe_index,
                        const unsigned int               local_index,
                        const std::integral_constant<int, structdim> &,
                        GlobalIndexType     &global_index,
                        const DoFPProcessor &process)
      {
        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // 1) no hp used -> fe_index == 0
        if (dof_handler.hp_capability_enabled == false)
          {
            AssertDimension(fe_index,
                            (DoFHandler<dim, spacedim>::default_fe_index));

            process(
              dof_handler.object_dof_indices
                [obj_level][structdim]
                [dof_handler.object_dof_ptr[obj_level][structdim][obj_index] +
                 local_index],
              global_index);

            return;
          }

        // 2) cell and hp is used -> there is only one fe_index
        if (structdim == dim)
          {
            process(
              dof_handler.object_dof_indices
                [obj_level][structdim]
                [dof_handler.object_dof_ptr[obj_level][structdim][obj_index] +
                 local_index],
              global_index);
            return;
          }

        // 3) general entity and hp is used
        AssertIndexRange(obj_level, dof_handler.object_dof_indices.size());
        AssertIndexRange(structdim,
                         dof_handler.object_dof_indices[obj_level].size());

        Assert(dof_handler.hp_capability_enabled, ExcInternalError());

        AssertIndexRange(structdim, dof_handler.hp_object_fe_ptr.size());
        AssertIndexRange(obj_index,
                         dof_handler.hp_object_fe_ptr[structdim].size());

        const auto ptr =
          std::find(dof_handler.hp_object_fe_indices[structdim].begin() +
                      dof_handler.hp_object_fe_ptr[structdim][obj_index],
                    dof_handler.hp_object_fe_indices[structdim].begin() +
                      dof_handler.hp_object_fe_ptr[structdim][obj_index + 1],
                    fe_index);

        Assert(ptr != dof_handler.hp_object_fe_indices[structdim].begin() +
                        dof_handler.hp_object_fe_ptr[structdim][obj_index + 1],
               ExcMessage(
                 "You are requesting an active FE index that is not assigned "
                 "to any of the cells connected to this entity."));

        const types::fe_index fe_index_ =
          std::distance(dof_handler.hp_object_fe_indices[structdim].begin() +
                          dof_handler.hp_object_fe_ptr[structdim][obj_index],
                        ptr);

        AssertIndexRange(
          dof_handler.hp_capability_enabled ?
            (dof_handler.hp_object_fe_ptr[structdim][obj_index] + fe_index_) :
            obj_index,
          dof_handler.object_dof_ptr[obj_level][structdim].size());

        AssertIndexRange(
          dof_handler.object_dof_ptr
              [obj_level][structdim]
              [dof_handler.hp_capability_enabled ?
                 (dof_handler.hp_object_fe_ptr[structdim][obj_index] +
                  fe_index_) :
                 obj_index] +
            local_index,
          dof_handler.object_dof_indices[obj_level][structdim].size());

        process(dof_handler.object_dof_indices
                  [obj_level][structdim]
                  [dof_handler.object_dof_ptr
                     [obj_level][structdim]
                     [dof_handler.hp_capability_enabled ?
                        (dof_handler.hp_object_fe_ptr[structdim][obj_index] +
                         fe_index_) :
                        obj_index] +
                   local_index],
                global_index);
      }

      /**
       * Determine start index and number of dofs of object in global data
       * structure.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim, int spacedim, int structdim>
      static std::pair<unsigned int, unsigned int>
      process_object_range(const DoFHandler<dim, spacedim> &dof_handler,
                           const unsigned int               obj_level,
                           const unsigned int               obj_index,
                           const types::fe_index            fe_index,
                           const std::integral_constant<int, structdim> &)
      {
        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // determine range of dofs in global data structure
        // 1) cell
        if (structdim == dim)
          {
            const unsigned int ptr_0 =
              dof_handler.object_dof_ptr[obj_level][structdim][obj_index];
            const unsigned int length =
              dof_handler.get_fe(fe_index).template n_dofs_per_object<dim>(0);

            return {ptr_0, length};
          }

        // 2) hp is not used -> fe_index == 0
        if (dof_handler.hp_capability_enabled == false)
          {
            AssertDimension(fe_index,
                            (DoFHandler<dim, spacedim>::default_fe_index));

            const unsigned int ptr_0 =
              dof_handler.object_dof_ptr[obj_level][structdim][obj_index];
            const unsigned int length =
              dof_handler.object_dof_ptr[obj_level][structdim][obj_index + 1] -
              ptr_0;

            return {ptr_0, length};
          }

        // 3) hp is used
        AssertIndexRange(obj_level, dof_handler.object_dof_indices.size());
        AssertIndexRange(structdim,
                         dof_handler.object_dof_indices[obj_level].size());

        AssertIndexRange(structdim, dof_handler.hp_object_fe_ptr.size());
        AssertIndexRange(obj_index,
                         dof_handler.hp_object_fe_ptr[structdim].size());

        const auto fe_index_local_ptr =
          std::find(dof_handler.hp_object_fe_indices[structdim].begin() +
                      dof_handler.hp_object_fe_ptr[structdim][obj_index],
                    dof_handler.hp_object_fe_indices[structdim].begin() +
                      dof_handler.hp_object_fe_ptr[structdim][obj_index + 1],
                    fe_index);

        Assert(fe_index_local_ptr !=
                 dof_handler.hp_object_fe_indices[structdim].begin() +
                   dof_handler.hp_object_fe_ptr[structdim][obj_index + 1],
               ExcMessage(
                 "You tried to call a function accessing DoF indices, but "
                 "they appear not be available (yet) or inconsistent. "
                 "Did you call distribute_dofs() first? Alternatively, if "
                 "you are using different elements on different cells (i.e., "
                 "you are using the hp capabilities of deal.II), did you "
                 "change the active_fe_index of a cell since you last "
                 "called distribute_dofs()?"));

        const types::fe_index fe_index_local =
          std::distance(dof_handler.hp_object_fe_indices[structdim].begin() +
                          dof_handler.hp_object_fe_ptr[structdim][obj_index],
                        fe_index_local_ptr);

        AssertIndexRange(
          dof_handler.hp_object_fe_ptr[structdim][obj_index] + fe_index_local,
          dof_handler.object_dof_ptr[obj_level][structdim].size());

        const unsigned int ptr_0 =
          dof_handler
            .object_dof_ptr[obj_level][structdim]
                           [dof_handler.hp_object_fe_ptr[structdim][obj_index] +
                            fe_index_local];
        const unsigned int ptr_1 =
          dof_handler
            .object_dof_ptr[obj_level][structdim]
                           [dof_handler.hp_object_fe_ptr[structdim][obj_index] +
                            fe_index_local + 1];

        return {ptr_0, ptr_1 - ptr_0};
      }

      template <int dim, int spacedim, int structdim, bool level_dof_access>
      static std::pair<unsigned int, unsigned int>
      process_object_range(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                              accessor,
        const types::fe_index fe_index)
      {
        return process_object_range(accessor.get_dof_handler(),
                                    accessor.level(),
                                    accessor.index(),
                                    fe_index,
                                    std::integral_constant<int, structdim>());
      }

      template <int dim, int spacedim, int structdim>
      static std::pair<unsigned int, unsigned int>
      process_object_range(dealii::DoFInvalidAccessor<structdim, dim, spacedim>,
                           const unsigned int)
      {
        DEAL_II_ASSERT_UNREACHABLE();

        return {0, 0};
      }



      /**
       * Process all dofs of an object.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim,
                int spacedim,
                int structdim,
                typename DoFProcessor,
                typename DoFMapping>
      static DEAL_II_ALWAYS_INLINE void
      process_object(const DoFHandler<dim, spacedim>              &dof_handler,
                     const unsigned int                            obj_level,
                     const unsigned int                            obj_index,
                     const types::fe_index                         fe_index,
                     const DoFMapping                             &mapping,
                     const std::integral_constant<int, structdim> &dd,
                     types::global_dof_index *&dof_indices_ptr,
                     const DoFProcessor       &process)
      {
        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // determine range of dofs in global data structure
        const auto range =
          process_object_range(dof_handler, obj_level, obj_index, fe_index, dd);
        if (range.second == 0)
          return;

        std::vector<types::global_dof_index> &object_dof_indices =
          dof_handler
            .object_dof_indices[structdim < dim ? 0 : obj_level][structdim];
        AssertIndexRange(range.first, object_dof_indices.size());
        types::global_dof_index *DEAL_II_RESTRICT stored_indices =
          object_dof_indices.data() + range.first;

        // process dofs
        for (unsigned int i = 0; i < range.second; ++i)
          {
            process(
              stored_indices[(structdim == 0 || structdim == dim) ? i :
                                                                    mapping(i)],
              dof_indices_ptr);
            if (dof_indices_ptr != nullptr)
              ++dof_indices_ptr;
          }
      }



      /**
       * Set the @p local_index-th degree of freedom corresponding to the
       * finite element specified by @p fe_index on the vertex with global
       * number @p vertex_index to @p global_index.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim, int spacedim, int structdim>
      static void
      set_dof_index(const DoFHandler<dim, spacedim>              &dof_handler,
                    const unsigned int                            obj_level,
                    const unsigned int                            obj_index,
                    const types::fe_index                         fe_index,
                    const unsigned int                            local_index,
                    const std::integral_constant<int, structdim> &dd,
                    const types::global_dof_index                 global_index)
      {
        process_dof_index(dof_handler,
                          obj_level,
                          obj_index,
                          fe_index,
                          local_index,
                          dd,
                          global_index,
                          [](auto &ptr, const auto &value) { ptr = value; });
      }


      /**
       * Get the @p local_index-th degree of freedom corresponding to the
       * finite element specified by @p fe_index on the vertex with global
       * number @p vertex_index to @p global_index.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim, int spacedim, int structdim>
      static types::global_dof_index
      get_dof_index(const DoFHandler<dim, spacedim>              &dof_handler,
                    const unsigned int                            obj_level,
                    const unsigned int                            obj_index,
                    const types::fe_index                         fe_index,
                    const unsigned int                            local_index,
                    const std::integral_constant<int, structdim> &dd)
      {
        types::global_dof_index global_index;
        process_dof_index(dof_handler,
                          obj_level,
                          obj_index,
                          fe_index,
                          local_index,
                          dd,
                          global_index,
                          [](const auto &ptr, auto &value) { value = ptr; });
        return global_index;
      }


      template <int dim, int spacedim>
      static types::global_dof_index &
      mg_vertex_dof_index(DoFHandler<dim, spacedim> &dof_handler,
                          const int                  level,
                          const unsigned int         vertex_index,
                          const unsigned int         i)
      {
        Assert(dof_handler.hp_capability_enabled == false,
               ExcMessage(
                 "DoFHandler in hp-mode does not implement multilevel DoFs."));

        return dof_handler.mg_vertex_dofs[vertex_index].access_index(
          level, i, dof_handler.get_fe().n_dofs_per_vertex());
      }



      /**
       * Return the number of different finite elements that are active on a
       * given object such as a vertex, line, or cell.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim, int spacedim, int structdim>
      static unsigned int
      n_active_fe_indices(const DoFHandler<dim, spacedim> &dof_handler,
                          const unsigned int               obj_level,
                          const unsigned int               obj_index,
                          const std::integral_constant<int, structdim> &)
      {
        (void)obj_level;

        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // 1) no hp used -> fe_index == 0
        if (dof_handler.hp_capability_enabled == false)
          return 1;

        // 2) cell and hp is used -> there is only one fe_index
        if (structdim == dim)
          return 1;

        // 3) general entity and hp is used
        AssertIndexRange(structdim, dof_handler.hp_object_fe_ptr.size());
        AssertIndexRange(obj_index + 1,
                         dof_handler.hp_object_fe_ptr[structdim].size());

        return dof_handler.hp_object_fe_ptr[structdim][obj_index + 1] -
               dof_handler.hp_object_fe_ptr[structdim][obj_index];
      }



      /**
       * Return the FE index of the n-th finite element active on a given
       * object such as a vertex, line, or cell.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim, int spacedim, int structdim>
      static types::fe_index
      nth_active_fe_index(const DoFHandler<dim, spacedim> &dof_handler,
                          const unsigned int               obj_level,
                          const unsigned int               obj_index,
                          const unsigned int               local_index,
                          const std::integral_constant<int, structdim> &)
      {
        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // for cells only one active FE index available
        Assert(((structdim == dim) &&
                (local_index != DoFHandler<dim, spacedim>::default_fe_index)) ==
                 false,
               ExcNotImplemented());

        // 1) no hp used -> fe_index == 0
        if (dof_handler.hp_capability_enabled == false)
          return DoFHandler<dim, spacedim>::default_fe_index;

        // 2) cell and hp is used -> there is only one fe_index
        if (structdim == dim)
          return dof_handler.hp_cell_active_fe_indices[obj_level][obj_index];

        // 3) general entity and hp is used
        AssertIndexRange(structdim, dof_handler.hp_object_fe_indices.size());
        AssertIndexRange(structdim, dof_handler.hp_object_fe_ptr.size());
        AssertIndexRange(obj_index,
                         dof_handler.hp_object_fe_ptr[structdim].size());
        AssertIndexRange(dof_handler.hp_object_fe_ptr[structdim][obj_index] +
                           local_index,
                         dof_handler.hp_object_fe_indices[structdim].size());

        return dof_handler.hp_object_fe_indices
          [structdim]
          [dof_handler.hp_object_fe_ptr[structdim][obj_index] + local_index];
      }



      /**
       * Returns all active FE indices on a given object such as a
       * vertex, line, or cell.
       *
       * The size of the returned set equals the number of finite elements that
       * are active on this vertex.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim, int spacedim, int structdim>
      static std::set<types::fe_index>
      get_active_fe_indices(const DoFHandler<dim, spacedim> &dof_handler,
                            const unsigned int               obj_level,
                            const unsigned int               obj_index,
                            const std::integral_constant<int, structdim> &t)
      {
        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // 1) no hp used -> fe_index == 0
        if (dof_handler.hp_capability_enabled == false)
          return {DoFHandler<dim, spacedim>::default_fe_index};

        // 2) cell and hp is used -> there is only one fe_index
        if (structdim == dim)
          return {dof_handler.hp_cell_active_fe_indices[obj_level][obj_index]};

        // 3) general entity and hp is used
        std::set<types::fe_index> active_fe_indices;
        for (unsigned int i = 0;
             i < n_active_fe_indices(dof_handler, obj_level, obj_index, t);
             ++i)
          active_fe_indices.insert(
            nth_active_fe_index(dof_handler, obj_level, obj_index, i, t));
        return active_fe_indices;
      }



      template <int dim, int spacedim, int structdim>
      static bool
      fe_index_is_active(const DoFHandler<dim, spacedim> &dof_handler,
                         const unsigned int               obj_level,
                         const unsigned int               obj_index,
                         const types::fe_index            fe_index,
                         const std::integral_constant<int, structdim> &)
      {
        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // 1) no hp used -> fe_index == 0
        if (dof_handler.hp_capability_enabled == false)
          return (fe_index == DoFHandler<dim, spacedim>::default_fe_index);

        // 2) cell and hp is used -> there is only one fe_index
        if (structdim == dim)
          return dof_handler.hp_cell_active_fe_indices[obj_level][obj_index] ==
                 fe_index;

        // 3) general entity and hp is used
        return std::find(
                 dof_handler.hp_object_fe_indices[structdim].begin() +
                   dof_handler.hp_object_fe_ptr[structdim][obj_index],
                 dof_handler.hp_object_fe_indices[structdim].begin() +
                   dof_handler.hp_object_fe_ptr[structdim][obj_index + 1],
                 fe_index) !=
               (dof_handler.hp_object_fe_indices[structdim].begin() +
                dof_handler.hp_object_fe_ptr[structdim][obj_index + 1]);
      }



      template <typename InputVector, typename ForwardIterator>
      static void
      extract_subvector_to(const InputVector             &values,
                           const types::global_dof_index *cache,
                           const types::global_dof_index *cache_end,
                           ForwardIterator                local_values_begin)
      {
        values.extract_subvector_to(cache, cache_end, local_values_begin);
      }



#ifdef DEAL_II_WITH_TRILINOS
      static std::vector<unsigned int>
      sort_indices(const types::global_dof_index *v_begin,
                   const types::global_dof_index *v_end)
      {
        // initialize original index locations
        std::vector<unsigned int> idx(v_end - v_begin);
        std::iota(idx.begin(), idx.end(), 0u);

        // sort indices based on comparing values in v
        std::sort(idx.begin(),
                  idx.end(),
                  [&v_begin](unsigned int i1, unsigned int i2) {
                    return *(v_begin + i1) < *(v_begin + i2);
                  });

        return idx;
      }



#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
      template <typename ForwardIterator, typename Number, typename MemorySpace>
      static void
      extract_subvector_to(
        const LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>
                                      &values,
        const types::global_dof_index *cache_begin,
        const types::global_dof_index *cache_end,
        ForwardIterator                local_values_begin)
      {
        std::vector<unsigned int> sorted_indices_pos =
          sort_indices(cache_begin, cache_end);
        const unsigned int cache_size = cache_end - cache_begin;
        std::vector<types::global_dof_index> cache_indices(cache_size);
        for (unsigned int i = 0; i < cache_size; ++i)
          cache_indices[i] = *(cache_begin + sorted_indices_pos[i]);

        IndexSet index_set(cache_indices.back() + 1);
        index_set.add_indices(cache_indices.begin(), cache_indices.end());
        index_set.compress();
        LinearAlgebra::ReadWriteVector<Number> read_write_vector(index_set);
        read_write_vector.import_elements(values, VectorOperation::insert);

        // Copy the elements from read_write_vector and reorder them.
        for (unsigned int i = 0; i < cache_size; ++i, ++local_values_begin)
          *local_values_begin = read_write_vector[sorted_indices_pos[i]];
      }
#  endif



      template <typename ForwardIterator>
      static void
      extract_subvector_to(const LinearAlgebra::EpetraWrappers::Vector &values,
                           const types::global_dof_index *cache_begin,
                           const types::global_dof_index *cache_end,
                           ForwardIterator                local_values_begin)
      {
        std::vector<unsigned int> sorted_indices_pos =
          sort_indices(cache_begin, cache_end);
        const unsigned int cache_size = cache_end - cache_begin;
        std::vector<types::global_dof_index> cache_indices(cache_size);
        for (unsigned int i = 0; i < cache_size; ++i)
          cache_indices[i] = *(cache_begin + sorted_indices_pos[i]);

        IndexSet index_set(cache_indices.back() + 1);
        index_set.add_indices(cache_indices.begin(), cache_indices.end());
        index_set.compress();
        LinearAlgebra::ReadWriteVector<double> read_write_vector(index_set);
        read_write_vector.import_elements(values, VectorOperation::insert);

        // Copy the elements from read_write_vector and reorder them.
        for (unsigned int i = 0; i < cache_size; ++i, ++local_values_begin)
          *local_values_begin = read_write_vector[sorted_indices_pos[i]];
      }
#endif

      /**
       * Loop over all degrees of freedom of the object described by the
       * provided @p accessor and @p fe_index and count them.
       */
      template <int dim, int spacedim, bool level_dof_access, int structdim>
      static unsigned int
      n_dof_indices(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                             &accessor,
        const types::fe_index fe_index_,
        const bool            count_level_dofs)
      {
        // note: we cannot rely on the template parameter level_dof_access here,
        // since the function get_mg_dof_indices()/set_mg_dof_indices() can be
        // called even if level_dof_access==false.
        if (count_level_dofs)
          {
            const auto &fe = accessor.get_fe(fe_index_);

            const unsigned int                                   //
              dofs_per_vertex = fe.n_dofs_per_vertex(),          //
              dofs_per_line   = fe.n_dofs_per_line(),            //
              dofs_per_quad   = fe.n_dofs_per_quad(0 /*dummy*/), //
              dofs_per_hex    = fe.n_dofs_per_hex();             //

            unsigned int index = 0;

            // 1) VERTEX dofs
            index += dofs_per_vertex * accessor.n_vertices();

            // 2) LINE dofs
            if (structdim == 2 || structdim == 3)
              index += dofs_per_line * accessor.n_lines();

            // 3) FACE dofs
            if (structdim == 3)
              index += dofs_per_quad * accessor.n_faces();

            // 4) INNER dofs
            const unsigned int interior_dofs =
              structdim == 1 ? dofs_per_line :
                               (structdim == 2 ? dofs_per_quad : dofs_per_hex);

            index += interior_dofs;

            return index;
          }
        else
          {
            const auto fe_index =
              internal::DoFAccessorImplementation::get_fe_index_or_default(
                accessor, fe_index_);

            unsigned int index = 0;

            // 1) VERTEX dofs
            for (const auto vertex : accessor.vertex_indices())
              index += process_object_range(accessor.get_dof_handler(),
                                            0,
                                            accessor.vertex_index(vertex),
                                            fe_index,
                                            std::integral_constant<int, 0>())
                         .second;

            // 2) LINE dofs
            if (structdim == 2 || structdim == 3)
              for (const auto line : accessor.line_indices())
                index +=
                  process_object_range(*accessor.line(line), fe_index).second;

            // 3) FACE dofs
            if (structdim == 3)
              for (const auto face : accessor.face_indices())
                index +=
                  process_object_range(*accessor.quad(face), fe_index).second;

            // 4) INNER dofs
            index += process_object_range(accessor, fe_index).second;

            return index;
          }
      }



      // The next few internal helper functions are needed to support various
      // DoFIndicesType kinds, e.g. actual vectors of DoFIndices or empty
      // types that we use when we only want to work on the internally stored
      // DoFs and never extract any number.
      template <typename ArrayType>
      static unsigned int
      get_array_length(const ArrayType &array)
      {
        return array.size();
      }

      static unsigned int
      get_array_length(const std::tuple<> &)
      {
        return 0;
      }

      template <typename ArrayType>
      static types::global_dof_index *
      get_array_ptr(const ArrayType &array)
      {
        return const_cast<types::global_dof_index *>(array.data());
      }

      static types::global_dof_index *
      get_array_ptr(const std::tuple<> &)
      {
        return nullptr;
      }



      /**
       * Loop over all degrees of freedom of the object described by the
       * provided @p accessor and @p fe_index and perform the static functions
       * provided by DoFOperation (set/get) on these.
       */
      template <int  dim,
                int  spacedim,
                bool level_dof_access,
                int  structdim,
                typename DoFIndicesType,
                typename DoFOperation,
                typename DoFProcessor>
      static void
      process_dof_indices(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                             &accessor,
        const DoFIndicesType &const_dof_indices,
        const types::fe_index fe_index_,
        const DoFOperation   &dof_operation,
        const DoFProcessor   &dof_processor,
        const bool            count_level_dofs)
      {
        const types::fe_index fe_index =
          internal::DoFAccessorImplementation::get_fe_index_or_default(
            accessor, fe_index_);

        // we cannot rely on the template parameter level_dof_access here, since
        // the function get_mg_dof_indices()/set_mg_dof_indices() can be called
        // even if level_dof_access==false.
        (void)count_level_dofs;

        const auto &fe = accessor.get_fe(fe_index);

        // we want to pass in rvalue 'std::tuple<>' types as `DoFIndicesType`,
        // but we need non-const references for std::vector<> types, so get in
        // a const reference here and immediately cast the constness away -
        // note that any use of the dereferenced invalid type will result in a
        // segfault
        types::global_dof_index *dof_indices_ptr =
          get_array_ptr(const_dof_indices);
        types::global_dof_index *end_dof_indices =
          dof_indices_ptr + get_array_length(const_dof_indices);

        // 1) VERTEX dofs, only step into the functions if we actually have
        // DoFs on them
        if (fe.n_dofs_per_vertex() > 0)
          for (const auto vertex : accessor.vertex_indices())
            dof_operation.process_vertex_dofs(*accessor.dof_handler,
                                              accessor.vertex_index(vertex),
                                              fe_index,
                                              dof_indices_ptr,
                                              dof_processor);

        // 2) copy dof numbers from the LINE, accounting for the possibility of
        // reversed line orientations.
        if (structdim > 1 && fe.n_dofs_per_line() > 0)
          {
            const auto line_indices = internal::TriaAccessorImplementation::
              Implementation::get_line_indices_of_cell(accessor);
            const auto line_orientations =
              internal::TriaAccessorImplementation::Implementation::
                get_line_orientations_of_cell(accessor);

            for (const auto line : accessor.line_indices())
              {
                const auto line_orientation = line_orientations[line];
                if (line_orientation == numbers::default_geometric_orientation)
                  dof_operation.process_dofs(
                    accessor.get_dof_handler(),
                    0,
                    line_indices[line],
                    fe_index,
                    [](const auto d) { return d; },
                    std::integral_constant<int, 1>(),
                    dof_indices_ptr,
                    dof_processor);
                else
                  {
                    Assert(line_orientation ==
                             numbers::reverse_line_orientation,
                           ExcInternalError());
                    dof_operation.process_dofs(
                      accessor.get_dof_handler(),
                      0,
                      line_indices[line],
                      fe_index,
                      [&fe, line_orientation](const auto d) {
                        return fe.adjust_line_dof_index_for_line_orientation(
                          d, line_orientation);
                      },
                      std::integral_constant<int, 1>(),
                      dof_indices_ptr,
                      dof_processor);
                  }
              }
          }

        // 3) copy dof numbers from the QUAD (i.e., 3d faces). Like lines we
        // only adjust dof indices (which has some cost) if we are not in the
        // default orientation.
        if (structdim == 3 && fe.max_dofs_per_quad() > 0)
          for (const auto face_no : accessor.face_indices())
            {
              const auto combined_orientation =
                accessor.combined_face_orientation(face_no);
              const unsigned int quad_index = accessor.quad_index(face_no);
              if (combined_orientation ==
                  numbers::default_geometric_orientation)
                dof_operation.process_dofs(
                  accessor.get_dof_handler(),
                  0,
                  quad_index,
                  fe_index,
                  [](const auto d) { return d; },
                  std::integral_constant<int, 2>(),
                  dof_indices_ptr,
                  dof_processor);
              else
                dof_operation.process_dofs(
                  accessor.get_dof_handler(),
                  0,
                  quad_index,
                  fe_index,
                  [&](const auto d) {
                    return fe.adjust_quad_dof_index_for_face_orientation(
                      d, face_no, combined_orientation);
                  },
                  std::integral_constant<int, 2>(),
                  dof_indices_ptr,
                  dof_processor);
            }

        // 4) INNER dofs (i.e., line dofs in 1d, quad dofs in 2d, or hex dofs in
        // 3d) - here we need to make sure that the shortcut to not run the
        // function does not miss the faces of wedge and pyramid elements where
        // n_dofs_per_object might not return the largest possible value
        if (((dim == 3 && structdim == 2) ?
               fe.max_dofs_per_quad() :
               fe.template n_dofs_per_object<structdim>()) > 0)
          dof_operation.process_dofs(
            accessor.get_dof_handler(),
            accessor.level(),
            accessor.index(),
            fe_index,
            [&](const auto d) { return d; },
            std::integral_constant<int, structdim>(),
            dof_indices_ptr,
            dof_processor);

        if (dof_indices_ptr != nullptr)
          {
            AssertDimension(n_dof_indices(accessor, fe_index, count_level_dofs),
                            dof_indices_ptr - get_array_ptr(const_dof_indices));
          }

        // PM: This is a part that should not be reached since it indicates that
        // an object (and/or its subobjects) is not active. However,
        // unfortunately this function is called by
        // DoFTools::set_periodicity_constraints() indirectly by
        // get_dof_indices() also for artificial faces to determine if a face
        // is artificial.
        types::global_dof_index invalid_index = numbers::invalid_dof_index;
        for (; dof_indices_ptr < end_dof_indices; ++dof_indices_ptr)
          dof_processor(invalid_index, dof_indices_ptr);
      }



      /**
       * An internal struct encapsulating the task of getting (vertex)
       * DoF indices.
       */
      template <int dim, int spacedim>
      struct DoFIndexProcessor
      {
        /**
         * Return vertex DoF indices.
         */
        template <typename DoFProcessor>
        DEAL_II_ALWAYS_INLINE void
        process_vertex_dofs(DoFHandler<dim, spacedim> &dof_handler,
                            const unsigned int         vertex_index,
                            const types::fe_index      fe_index,
                            types::global_dof_index  *&dof_indices_ptr,
                            const DoFProcessor        &dof_processor) const
        {
          process_object(
            dof_handler,
            0,
            vertex_index,
            fe_index,
            [](const auto d) {
              DEAL_II_ASSERT_UNREACHABLE();
              return d;
            },
            std::integral_constant<int, 0>(),
            dof_indices_ptr,
            dof_processor);
        }

        /**
         * Return DoF indices for lines, quads, and inner degrees of freedom.
         */
        template <int structdim, typename DoFMapping, typename DoFProcessor>
        DEAL_II_ALWAYS_INLINE void
        process_dofs(const DoFHandler<dim, spacedim> &dof_handler,
                     const unsigned int               obj_level,
                     const unsigned int               obj_index,
                     const types::fe_index            fe_index,
                     const DoFMapping                &mapping,
                     const std::integral_constant<int, structdim>,
                     types::global_dof_index *&dof_indices_ptr,
                     const DoFProcessor       &dof_processor) const
        {
          process_object(
            dof_handler,
            obj_level,
            obj_index,
            fe_index,
            mapping,
            std::integral_constant<int, std::min(structdim, dim)>(),
            dof_indices_ptr,
            dof_processor);
        }
      };



      /**
       * An internal struct encapsulating the task of getting level (vertex)
       * DoF indices.
       */
      template <int dim, int spacedim>
      struct MGDoFIndexProcessor
      {
        /**
         * Constructor.
         */
        MGDoFIndexProcessor(const unsigned int level)
          : level(level)
        {}

        /**
         * Return vertex DoF indices.
         */
        template <typename DoFProcessor>
        DEAL_II_ALWAYS_INLINE void
        process_vertex_dofs(DoFHandler<dim, spacedim> &dof_handler,
                            const unsigned int         vertex_index,
                            const types::fe_index,
                            types::global_dof_index *&dof_indices_ptr,
                            const DoFProcessor       &dof_processor) const
        {
          const unsigned int n_indices =
            dof_handler.get_fe(0).template n_dofs_per_object<0>();
          types::global_dof_index *stored_indices =
            &dof_handler.mg_vertex_dofs[vertex_index].access_index(level,
                                                                   0,
                                                                   n_indices);
          for (unsigned int d = 0; d < n_indices; ++d, ++dof_indices_ptr)
            dof_processor(stored_indices[d], dof_indices_ptr);
        }

        /**
         * Return DoF indices for lines, quads, and inner degrees of freedom.
         */
        template <int structdim, typename DoFMapping, typename DoFProcessor>
        DEAL_II_ALWAYS_INLINE void
        process_dofs(const DoFHandler<dim, spacedim> &dof_handler,
                     const unsigned int,
                     const unsigned int    obj_index,
                     const types::fe_index fe_index,
                     const DoFMapping     &mapping,
                     const std::integral_constant<int, structdim>,
                     types::global_dof_index *&dof_indices_ptr,
                     const DoFProcessor       &dof_processor) const
        {
          const unsigned int n_indices =
            dof_handler.get_fe(0).template n_dofs_per_object<structdim>();
          types::global_dof_index *stored_indices = &get_mg_dof_index(
            dof_handler,
            dof_handler.mg_levels[level],
            dof_handler.mg_faces,
            obj_index,
            fe_index,
            0,
            std::integral_constant<int, std::min(structdim, dim)>());
          for (unsigned int d = 0; d < n_indices; ++d, ++dof_indices_ptr)
            dof_processor(stored_indices[structdim < dim ? mapping(d) : d],
                          dof_indices_ptr);
        }

      private:
        const unsigned int level;
      };



      template <int dim, int spacedim, bool level_dof_access, int structdim>
      static void
      get_dof_indices(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                                             &accessor,
        std::vector<types::global_dof_index> &dof_indices,
        const types::fe_index                 fe_index)
      {
        process_dof_indices(
          accessor,
          dof_indices,
          fe_index,
          DoFIndexProcessor<dim, spacedim>(),
          [](auto stored_index, auto dof_ptr) { *dof_ptr = stored_index; },
          false);
      }



      template <int dim, int spacedim, bool level_dof_access, int structdim>
      static void
      set_dof_indices(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                                                   &accessor,
        const std::vector<types::global_dof_index> &dof_indices,
        const types::fe_index                       fe_index)
      {
        // Note: this function is as general as `get_dof_indices()`. This
        // assert is placed here since it is currently only used by the
        // function DoFCellAccessor::set_dof_indices(), which is called by
        // internal::DoFHandlerImplementation::Policy::Implementation::distribute_dofs().
        // In the case of new use cases, this assert can be removed.
        Assert(
          dim == structdim,
          ExcMessage(
            "This function is intended to be used for DoFCellAccessor, i.e., "
            "dimension == structdim."));

        process_dof_indices(
          accessor,
          dof_indices,
          fe_index,
          DoFIndexProcessor<dim, spacedim>(),
          [](auto &stored_index, auto dof_ptr) { stored_index = *dof_ptr; },
          false);
      }



      template <int dim, int spacedim, bool level_dof_access, int structdim>
      static void
      get_mg_dof_indices(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                                             &accessor,
        const int                             level,
        std::vector<types::global_dof_index> &dof_indices,
        const types::fe_index                 fe_index)
      {
        Assert((fe_index == DoFHandler<dim, spacedim>::default_fe_index),
               ExcMessage("MG DoF indices cannot be queried in hp case"));
        process_dof_indices(
          accessor,
          dof_indices,
          fe_index,
          MGDoFIndexProcessor<dim, spacedim>(level),
          [](auto stored_index, auto dof_ptr) { *dof_ptr = stored_index; },
          true);
      }



      template <int dim, int spacedim, bool level_dof_access, int structdim>
      static void
      set_mg_dof_indices(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                                                   &accessor,
        const int                                   level,
        const std::vector<types::global_dof_index> &dof_indices,
        const types::fe_index                       fe_index)
      {
        Assert((fe_index == DoFHandler<dim, spacedim>::default_fe_index),
               ExcMessage("MG DoF indices cannot be queried in hp case"));

        // Note: this function is as general as `get_mg_dof_indices()`. This
        // assert is placed here since it is currently only used by the
        // function DoFCellAccessor::set_mg_dof_indices(), which is called by
        // internal::DoFHandlerImplementation::Policy::Implementation::distribute_mg_dofs().
        // In the case of new use cases, this assert can be removed.
        Assert(dim == structdim,
               ExcMessage("This function is intended to be used for "
                          "DoFCellAccessor, i.e., dimension == structdim."));

        process_dof_indices(
          accessor,
          dof_indices,
          fe_index,
          MGDoFIndexProcessor<dim, spacedim>(level),
          [](auto &stored_index, auto dof_ptr) { stored_index = *dof_ptr; },
          true);
      }



      template <int dim, int spacedim>
      static types::global_dof_index &
      get_mg_dof_index(
        const DoFHandler<dim, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<dim>>
          &mg_level,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<dim>>
          &,
        const unsigned int    obj_index,
        const types::fe_index fe_index,
        const unsigned int    local_index,
        const std::integral_constant<int, dim>)
      {
        Assert(dof_handler.hp_capability_enabled == false,
               (typename DoFHandler<dim, spacedim>::ExcNotImplementedWithHP()));

        return mg_level->dof_object.access_dof_index(
          static_cast<const DoFHandler<dim, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index);
      }



      template <int dim, int spacedim, std::enable_if_t<(dim > 1), int> = 0>
      static types::global_dof_index &
      get_mg_dof_index(
        const DoFHandler<dim, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<dim>>
          &,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<dim>>
                             &mg_faces,
        const unsigned int    obj_index,
        const types::fe_index fe_index,
        const unsigned int    local_index,
        const std::integral_constant<int, 1>)
      {
        return mg_faces->lines.access_dof_index(
          static_cast<const DoFHandler<dim, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index);
      }



      template <int spacedim>
      static types::global_dof_index &
      get_mg_dof_index(
        const DoFHandler<3, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<3>>
          &,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<3>>
                             &mg_faces,
        const unsigned int    obj_index,
        const types::fe_index fe_index,
        const unsigned int    local_index,
        const std::integral_constant<int, 2>)
      {
        Assert(dof_handler.hp_capability_enabled == false,
               (typename DoFHandler<3, spacedim>::ExcNotImplementedWithHP()));
        return mg_faces->quads.access_dof_index(
          static_cast<const DoFHandler<3, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index);
      }
    };



    template <int dim, int spacedim, bool level_dof_access>
    void
    get_cell_dof_indices(
      const dealii::DoFCellAccessor<dim, spacedim, level_dof_access> &accessor,
      Implementation::dof_index_vector_type &dof_indices,
      const unsigned int                     fe_index);
  } // namespace DoFAccessorImplementation
} // namespace internal



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<structdim, dim, spacedim, level_dof_access>::dof_index(
  const unsigned int    i,
  const types::fe_index fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  // access the respective DoF
  return dealii::internal::DoFAccessorImplementation::Implementation::
    get_dof_index(*this->dof_handler,
                  this->level(),
                  this->index(),
                  fe_index,
                  i,
                  std::integral_constant<int, structdim>());
}


template <int structdim, int dim, int spacedim, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<structdim, dim, spacedim, level_dof_access>::mg_dof_index(
  const int          level,
  const unsigned int i) const
{
  return internal::DoFAccessorImplementation::Implementation::get_mg_dof_index(
    *this->dof_handler,
    this->dof_handler->mg_levels[level],
    this->dof_handler->mg_faces,
    this->index(),
    0,
    i,
    std::integral_constant<int, structdim>());
}


template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::set_dof_index(
  const unsigned int            i,
  const types::global_dof_index index,
  const types::fe_index         fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  // access the respective DoF
  dealii::internal::DoFAccessorImplementation::Implementation::set_dof_index(
    *this->dof_handler,
    this->level(),
    this->index(),
    fe_index,
    i,
    std::integral_constant<int, structdim>(),
    index);
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::set_mg_dof_index(
  const int                     level,
  const unsigned int            i,
  const types::global_dof_index index) const
{
  internal::DoFAccessorImplementation::Implementation::get_mg_dof_index(
    *this->dof_handler,
    this->dof_handler->mg_levels[level],
    this->dof_handler->mg_faces,
    this->index(),
    0,
    i,
    std::integral_constant<int, structdim>()) = index;
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline unsigned int
DoFAccessor<structdim, dim, spacedim, level_dof_access>::n_active_fe_indices()
  const
{
  // access the respective DoF
  return dealii::internal::DoFAccessorImplementation::Implementation::
    n_active_fe_indices(*this->dof_handler,
                        this->level(),
                        this->index(),
                        std::integral_constant<int, structdim>());
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline types::fe_index
DoFAccessor<structdim, dim, spacedim, level_dof_access>::nth_active_fe_index(
  const unsigned int n) const
{
  // access the respective DoF
  return dealii::internal::DoFAccessorImplementation::Implementation::
    nth_active_fe_index(*this->dof_handler,
                        this->level(),
                        this->index(),
                        n,
                        std::integral_constant<int, structdim>());
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline std::set<types::fe_index>
DoFAccessor<structdim, dim, spacedim, level_dof_access>::get_active_fe_indices()
  const
{
  std::set<types::fe_index> active_fe_indices;
  for (unsigned int i = 0; i < n_active_fe_indices(); ++i)
    active_fe_indices.insert(nth_active_fe_index(i));
  return active_fe_indices;
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline bool
DoFAccessor<structdim, dim, spacedim, level_dof_access>::fe_index_is_active(
  const types::fe_index fe_index) const
{
  // access the respective DoF
  return dealii::internal::DoFAccessorImplementation::Implementation::
    fe_index_is_active(*this->dof_handler,
                       this->level(),
                       this->index(),
                       fe_index,
                       std::integral_constant<int, structdim>());
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<structdim, dim, spacedim, level_dof_access>::vertex_dof_index(
  const unsigned int    vertex,
  const unsigned int    i,
  const types::fe_index fe_index_) const
{
  const types::fe_index fe_index =
    (((this->dof_handler->hp_capability_enabled == false) &&
      (fe_index_ == numbers::invalid_fe_index)) ?
       // No hp enabled, and the argument is at its default value -> we
       // can translate to the default active fe index
       DoFHandler<dim, spacedim>::default_fe_index :
       // Otherwise: If anything other than the default is provided by
       // the caller, then we should take just that. As an exception, if
       // we are on a cell (rather than a face/edge/vertex), then we know
       // that there is only one active fe index on this cell and we can
       // use that:
       ((dim == structdim) && (fe_index_ == numbers::invalid_fe_index) ?
          this->nth_active_fe_index(0) :
          fe_index_));

  return dealii::internal::DoFAccessorImplementation::Implementation::
    get_dof_index(*this->dof_handler,
                  0,
                  this->vertex_index(vertex),
                  fe_index,
                  i,
                  std::integral_constant<int, 0>());
}


template <int structdim, int dim, int spacedim, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<structdim, dim, spacedim, level_dof_access>::mg_vertex_dof_index(
  const int             level,
  const unsigned int    vertex,
  const unsigned int    i,
  const types::fe_index fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);
  (void)fe_index;
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  Assert(this->dof_handler->mg_vertex_dofs.size() > 0,
         ExcMessage("Multigrid DoF indices can only be accessed after "
                    "DoFHandler::distribute_mg_dofs() has been called!"));
  AssertIndexRange(vertex, this->n_vertices());
  AssertIndexRange(i, this->dof_handler->get_fe(fe_index).n_dofs_per_vertex());

  Assert(dof_handler->hp_capability_enabled == false,
         ExcMessage(
           "DoFHandler in hp-mode does not implement multilevel DoFs."));

  return this->dof_handler->mg_vertex_dofs[this->vertex_index(vertex)]
    .access_index(level, i, this->dof_handler->get_fe().n_dofs_per_vertex());
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::
  set_mg_vertex_dof_index(const int                     level,
                          const unsigned int            vertex,
                          const unsigned int            i,
                          const types::global_dof_index index,
                          const types::fe_index         fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);
  (void)fe_index;
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  AssertIndexRange(vertex, this->n_vertices());
  AssertIndexRange(i, this->dof_handler->get_fe(fe_index).n_dofs_per_vertex());

  Assert(dof_handler->hp_capability_enabled == false,
         ExcMessage(
           "DoFHandler in hp-mode does not implement multilevel DoFs."));

  this->dof_handler->mg_vertex_dofs[this->vertex_index(vertex)].access_index(
    level, i, this->dof_handler->get_fe().n_dofs_per_vertex()) = index;
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline const FiniteElement<dim, spacedim> &
DoFAccessor<structdim, dim, spacedim, level_dof_access>::get_fe(
  const types::fe_index fe_index) const
{
  Assert(fe_index_is_active(fe_index) == true,
         ExcMessage("This function can only be called for active FE indices"));

  return this->dof_handler->get_fe(fe_index);
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline typename dealii::internal::DoFHandlerImplementation::
  Iterators<dim, spacedim, level_dof_access>::line_iterator
  DoFAccessor<structdim, dim, spacedim, level_dof_access>::line(
    const unsigned int i) const
{
  // if we are asking for a particular line and this object refers to
  // a line, then the only valid index is i==0 and we should return
  // *this
  if (structdim == 1)
    {
      Assert(i == 0,
             ExcMessage("You can only ask for line zero if the "
                        "current object is a line itself."));
      return typename dealii::internal::DoFHandlerImplementation::
        Iterators<dim, spacedim, level_dof_access>::cell_iterator(
          &this->get_triangulation(),
          this->level(),
          this->index(),
          &this->get_dof_handler());
    }

  // otherwise we need to be in structdim>=2
  Assert(structdim > 1, ExcImpossibleInDim(structdim));
  Assert(dim > 1, ExcImpossibleInDim(dim));

  // checking of 'i' happens in line_index(i)
  return typename dealii::internal::DoFHandlerImplementation::
    Iterators<dim, spacedim, level_dof_access>::line_iterator(
      this->tria,
      0, // only sub-objects are allowed, which have no level
      this->line_index(i),
      this->dof_handler);
}


template <int structdim, int dim, int spacedim, bool level_dof_access>
inline typename dealii::internal::DoFHandlerImplementation::
  Iterators<dim, spacedim, level_dof_access>::quad_iterator
  DoFAccessor<structdim, dim, spacedim, level_dof_access>::quad(
    const unsigned int i) const
{
  // if we are asking for a
  // particular quad and this object
  // refers to a quad, then the only
  // valid index is i==0 and we
  // should return *this
  if (structdim == 2)
    {
      Assert(i == 0,
             ExcMessage("You can only ask for quad zero if the "
                        "current object is a quad itself."));
      return typename dealii::internal::DoFHandlerImplementation::
        Iterators<dim, spacedim>::cell_iterator(&this->get_triangulation(),
                                                this->level(),
                                                this->index(),
                                                &this->get_dof_handler());
    }

  // otherwise we need to be in structdim>=3
  Assert(structdim > 2, ExcImpossibleInDim(structdim));
  Assert(dim > 2, ExcImpossibleInDim(dim));

  // checking of 'i' happens in quad_index(i)
  return typename dealii::internal::DoFHandlerImplementation::
    Iterators<dim, spacedim, level_dof_access>::quad_iterator(
      this->tria,
      0, // only sub-objects are allowed, which have no level
      this->quad_index(i),
      this->dof_handler);
}


/*----------------- Functions: DoFAccessor<0,1,spacedim> --------------------*/


template <int spacedim, bool level_dof_access>
inline DoFAccessor<0, 1, spacedim, level_dof_access>::DoFAccessor()
{
  Assert(false, ExcInvalidObject());
}



template <int spacedim, bool level_dof_access>
inline DoFAccessor<0, 1, spacedim, level_dof_access>::DoFAccessor(
  const Triangulation<1, spacedim>                       *tria,
  const typename TriaAccessor<0, 1, spacedim>::VertexKind vertex_kind,
  const unsigned int                                      vertex_index,
  const DoFHandler<1, spacedim>                          *dof_handler)
  : BaseClass(tria, vertex_kind, vertex_index)
  , dof_handler(const_cast<DoFHandler<1, spacedim> *>(dof_handler))
{}



template <int spacedim, bool level_dof_access>
inline DoFAccessor<0, 1, spacedim, level_dof_access>::DoFAccessor(
  const Triangulation<1, spacedim> *tria,
  const int                         level,
  const int                         index,
  const DoFHandler<1, spacedim>    *dof_handler)
  // This is the constructor signature for "ordinary" (non-vertex)
  // accessors and we shouldn't be calling it altogether. But it is also
  // the constructor that the default-constructor of TriaRawIterator
  // calls when default-constructing an iterator object. If so, this
  // happens with level==-2 and index==-2, and this is the only case we
  // would like to support. We do this by just forwarding to the
  // other constructor of this class, and then asserting the condition
  // on level and index.
  : DoFAccessor<0, 1, spacedim, level_dof_access>(
      tria,
      TriaAccessor<0, 1, spacedim>::interior_vertex,
      0U,
      dof_handler)
{
  (void)level;
  (void)index;
  Assert((tria == nullptr) && (level == -2) && (index == -2) &&
           (dof_handler == nullptr),
         ExcMessage(
           "This constructor can not be called for face iterators in 1d, "
           "except to default-construct iterator objects."));
}



template <int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2>
DoFAccessor<0, 1, spacedim, level_dof_access>::DoFAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false, ExcInvalidObject());
}



template <int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline DoFAccessor<0, 1, spacedim, level_dof_access>::DoFAccessor(
  const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &)
{
  Assert(false, ExcInvalidObject());
}



template <int spacedim, bool level_dof_access>
inline void
DoFAccessor<0, 1, spacedim, level_dof_access>::set_dof_handler(
  DoFHandler<1, spacedim> *dh)
{
  Assert(dh != nullptr, ExcInvalidObject());
  this->dof_handler = dh;
}



template <int spacedim, bool level_dof_access>
inline void
DoFAccessor<0, 1, spacedim, level_dof_access>::set_dof_index(
  const unsigned int /*i*/,
  const types::global_dof_index /*index*/,
  const types::fe_index /*fe_index*/) const
{
  DEAL_II_NOT_IMPLEMENTED();
}



template <int spacedim, bool level_dof_access>
inline const DoFHandler<1, spacedim> &
DoFAccessor<0, 1, spacedim, level_dof_access>::get_dof_handler() const
{
  return *this->dof_handler;
}



template <int spacedim, bool level_dof_access>
inline void
DoFAccessor<0, 1, spacedim, level_dof_access>::get_dof_indices(
  std::vector<types::global_dof_index> &dof_indices,
  const types::fe_index                 fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  for (unsigned int i = 0; i < dof_indices.size(); ++i)
    dof_indices[i] = dealii::internal::DoFAccessorImplementation::
      Implementation::get_dof_index(*dof_handler,
                                    0,
                                    this->global_vertex_index,
                                    fe_index,
                                    i,
                                    std::integral_constant<int, 0>());
}



template <int spacedim, bool level_dof_access>
inline void
DoFAccessor<0, 1, spacedim, level_dof_access>::get_mg_dof_indices(
  const int                             level,
  std::vector<types::global_dof_index> &dof_indices,
  const types::fe_index                 fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);
  (void)fe_index;
  AssertDimension(fe_index, (DoFHandler<1, spacedim>::default_fe_index));

  for (unsigned int i = 0; i < dof_indices.size(); ++i)
    dof_indices[i] =
      dealii::internal::DoFAccessorImplementation::Implementation::
        mg_vertex_dof_index(*dof_handler, level, this->global_vertex_index, i);
}



template <int spacedim, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<0, 1, spacedim, level_dof_access>::vertex_dof_index(
  const unsigned int    vertex,
  const unsigned int    i,
  const types::fe_index fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  (void)vertex;
  AssertIndexRange(vertex, 1);
  return dealii::internal::DoFAccessorImplementation::Implementation::
    get_dof_index(*dof_handler,
                  0,
                  this->global_vertex_index,
                  fe_index,
                  i,
                  std::integral_constant<int, 0>());
}



template <int spacedim, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<0, 1, spacedim, level_dof_access>::dof_index(
  const unsigned int    i,
  const types::fe_index fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  return dealii::internal::DoFAccessorImplementation::Implementation::
    get_dof_index(*this->dof_handler,
                  0,
                  this->vertex_index(0),
                  fe_index,
                  i,
                  std::integral_constant<int, 0>());
}



template <int spacedim, bool level_dof_access>
inline unsigned int
DoFAccessor<0, 1, spacedim, level_dof_access>::n_active_fe_indices() const
{
  return 1;
}



template <int spacedim, bool level_dof_access>
inline types::fe_index
DoFAccessor<0, 1, spacedim, level_dof_access>::nth_active_fe_index(
  const unsigned int /*n*/) const
{
  return 0;
}



template <int spacedim, bool level_dof_access>
inline bool
DoFAccessor<0, 1, spacedim, level_dof_access>::fe_index_is_active(
  const types::fe_index /*fe_index*/) const
{
  DEAL_II_NOT_IMPLEMENTED();
  return false;
}



template <int spacedim, bool level_dof_access>
inline const FiniteElement<1, spacedim> &
DoFAccessor<0, 1, spacedim, level_dof_access>::get_fe(
  const types::fe_index fe_index) const
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  return dof_handler->get_fe(fe_index);
}



template <int spacedim, bool level_dof_access>
inline void
DoFAccessor<0, 1, spacedim, level_dof_access>::copy_from(
  const TriaAccessorBase<0, 1, spacedim> &da)
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  BaseClass::copy_from(da);
}



template <int spacedim, bool level_dof_access>
template <bool level_dof_access2>
inline void
DoFAccessor<0, 1, spacedim, level_dof_access>::copy_from(
  const DoFAccessor<0, 1, spacedim, level_dof_access2> &a)
{
  BaseClass::copy_from(a);
  set_dof_handler(a.dof_handler);
}



template <int spacedim, bool level_dof_access>
inline TriaIterator<DoFAccessor<0, 1, spacedim, level_dof_access>>
DoFAccessor<0, 1, spacedim, level_dof_access>::child(
  const unsigned int /*i*/) const
{
  return TriaIterator<DoFAccessor<0, 1, spacedim, level_dof_access>>();
}



template <int spacedim, bool level_dof_access>
inline typename dealii::internal::DoFHandlerImplementation::
  Iterators<1, spacedim, level_dof_access>::line_iterator
  DoFAccessor<0, 1, spacedim, level_dof_access>::line(
    const unsigned int /*c*/) const
{
  DEAL_II_NOT_IMPLEMENTED();
  return typename dealii::internal::DoFHandlerImplementation::
    Iterators<1, spacedim, level_dof_access>::line_iterator();
}



template <int spacedim, bool level_dof_access>
inline typename dealii::internal::DoFHandlerImplementation::
  Iterators<1, spacedim, level_dof_access>::quad_iterator
  DoFAccessor<0, 1, spacedim, level_dof_access>::quad(
    const unsigned int /*c*/) const
{
  DEAL_II_NOT_IMPLEMENTED();
  return typename dealii::internal::DoFHandlerImplementation::
    Iterators<1, spacedim, level_dof_access>::quad_iterator();
}



template <int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline bool
DoFAccessor<0, 1, spacedim, level_dof_access>::operator==(
  const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &a) const
{
  Assert(structdim2 == 0, ExcCantCompareIterators());
  Assert(this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator==(a));
}



template <int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline bool
DoFAccessor<0, 1, spacedim, level_dof_access>::operator!=(
  const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &a) const
{
  Assert(structdim2 == 0, ExcCantCompareIterators());
  Assert(this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator!=(a));
}



/*------------------------- Functions: DoFCellAccessor -----------------------*/


namespace internal
{
  namespace DoFCellAccessorImplementation
  {
    /**
     * A class with the same purpose as the similarly named class of the
     * Triangulation class. See there for more information.
     */
    struct Implementation
    {
      /**
       * Do what the active_fe_index function in the parent class is supposed to
       * do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static types::fe_index
      active_fe_index(
        const DoFCellAccessor<dim, spacedim, level_dof_access> &accessor)
      {
        if (accessor.dof_handler->hp_capability_enabled == false)
          return DoFHandler<dim, spacedim>::default_fe_index;

        Assert(accessor.dof_handler != nullptr,
               (typename std::decay_t<decltype(accessor)>::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->hp_cell_future_fe_indices.size(),
               ExcMessage("DoFHandler not initialized"));

        return accessor.dof_handler
          ->hp_cell_active_fe_indices[accessor.level()][accessor.present_index];
      }



      /**
       * Do what the set_active_fe_index function in the parent class is
       * supposed to do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static void
      set_active_fe_index(
        const DoFCellAccessor<dim, spacedim, level_dof_access> &accessor,
        const types::fe_index                                   i)
      {
        if (accessor.dof_handler->hp_capability_enabled == false)
          {
            AssertDimension(i, (DoFHandler<dim, spacedim>::default_fe_index));
            return;
          }

        Assert(accessor.dof_handler != nullptr,
               (typename std::decay_t<decltype(accessor)>::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->hp_cell_future_fe_indices.size(),
               ExcMessage("DoFHandler not initialized"));
        Assert(i != numbers::invalid_fe_index,
               ExcMessage("Invalid finite element index."));

        accessor.dof_handler
          ->hp_cell_active_fe_indices[accessor.level()]
                                     [accessor.present_index] = i;
      }



      /**
       * Do what the future_fe_index function in the parent class is supposed to
       * do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static types::fe_index
      future_fe_index(
        const DoFCellAccessor<dim, spacedim, level_dof_access> &accessor)
      {
        if (accessor.dof_handler->hp_capability_enabled == false)
          return DoFHandler<dim, spacedim>::default_fe_index;

        Assert(accessor.dof_handler != nullptr,
               (typename std::decay_t<decltype(accessor)>::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->hp_cell_future_fe_indices.size(),
               ExcMessage("DoFHandler not initialized"));

        if (future_fe_index_set(accessor))
          return accessor.dof_handler
            ->hp_cell_future_fe_indices[accessor.level()]
                                       [accessor.present_index];
        else
          return accessor.dof_handler
            ->hp_cell_active_fe_indices[accessor.level()]
                                       [accessor.present_index];
      }


      /**
       * Do what the set_future_fe_index function in the parent class is
       * supposed to do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static void
      set_future_fe_index(
        const DoFCellAccessor<dim, spacedim, level_dof_access> &accessor,
        const types::fe_index                                   i)
      {
        if (accessor.dof_handler->hp_capability_enabled == false)
          {
            AssertDimension(i, (DoFHandler<dim, spacedim>::default_fe_index));
            return;
          }

        Assert(accessor.dof_handler != nullptr,
               (typename std::decay_t<decltype(accessor)>::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->hp_cell_future_fe_indices.size(),
               ExcMessage("DoFHandler not initialized"));
        Assert(i != numbers::invalid_fe_index,
               ExcMessage("Invalid finite element index."));

        accessor.dof_handler
          ->hp_cell_future_fe_indices[accessor.level()]
                                     [accessor.present_index] = i;
      }



      /**
       * Do what the future_fe_index_set function in the parent class is
       * supposed to do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static bool
      future_fe_index_set(
        const DoFCellAccessor<dim, spacedim, level_dof_access> &accessor)
      {
        if (accessor.dof_handler->hp_capability_enabled == false)
          return false;

        Assert(accessor.dof_handler != nullptr,
               (typename std::decay_t<decltype(accessor)>::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->hp_cell_future_fe_indices.size(),
               ExcMessage("DoFHandler not initialized"));

        return accessor.dof_handler
                 ->hp_cell_future_fe_indices[accessor.level()]
                                            [accessor.present_index] !=
               numbers::invalid_fe_index;
      }



      /**
       * Do what the clear_fe_index function in the parent class is supposed to
       * do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static void
      clear_future_fe_index(
        const DoFCellAccessor<dim, spacedim, level_dof_access> &accessor)
      {
        if (accessor.dof_handler->hp_capability_enabled == false)
          return;

        Assert(accessor.dof_handler != nullptr,
               (typename std::decay_t<decltype(accessor)>::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->hp_cell_future_fe_indices.size(),
               ExcMessage("DoFHandler not initialized"));

        accessor.dof_handler
          ->hp_cell_future_fe_indices[accessor.level()]
                                     [accessor.present_index] =
          numbers::invalid_fe_index;
      }
    };
  } // namespace DoFCellAccessorImplementation
} // namespace internal



template <int dimension_, int space_dimension_, bool level_dof_access>
inline DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  DoFCellAccessor(const Triangulation<dimension_, space_dimension_> *tria,
                  const int                                          level,
                  const int                                          index,
                  const AccessorData                                *local_data)
  : DoFAccessor<dimension_, dimension_, space_dimension_, level_dof_access>(
      tria,
      level,
      index,
      local_data)
{}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2>
inline DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  DoFCellAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false, typename BaseClass::ExcInvalidObject());
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  DoFCellAccessor(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &other)
  : BaseClass(other)
{}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline TriaIterator<
  DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::neighbor(
  const unsigned int i) const
{
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
    q(this->tria,
      this->neighbor_level(i),
      this->neighbor_index(i),
      this->dof_handler);

  if constexpr (running_in_debug_mode())
    {
      if (q.state() != IteratorState::past_the_end)
        Assert(q->used(), ExcInternalError());
    }
  return q;
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline TriaIterator<
  DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::child(
  const unsigned int i) const
{
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
    q(this->tria, this->level() + 1, this->child_index(i), this->dof_handler);

  if constexpr (running_in_debug_mode())
    {
      if (q.state() != IteratorState::past_the_end)
        Assert(q->used(), ExcInternalError());
    }
  return q;
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline boost::container::small_vector<
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>,
  GeometryInfo<dimension_>::max_children_per_cell>
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  child_iterators() const
{
  boost::container::small_vector<
    TriaIterator<
      DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>,
    GeometryInfo<dimension_>::max_children_per_cell>
    child_iterators(this->n_children());

  for (unsigned int i = 0; i < this->n_children(); ++i)
    child_iterators[i] = this->child(i);

  return child_iterators;
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline TriaIterator<
  DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::parent() const
{
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
    q(this->tria, this->level() - 1, this->parent_index(), this->dof_handler);

  return q;
}



namespace internal
{
  namespace DoFCellAccessorImplementation
  {
    template <int dim, int spacedim, bool level_dof_access>
    inline TriaIterator<
      dealii::DoFAccessor<dim - 1, dim, spacedim, level_dof_access>>
    get_face(
      const dealii::DoFCellAccessor<dim, spacedim, level_dof_access> &cell,
      const unsigned int                                              i,
      const std::integral_constant<int, 1>)
    {
      dealii::DoFAccessor<0, dim, spacedim, level_dof_access> a(
        &cell.get_triangulation(),
        ((i == 0) && cell.at_boundary(0) ?
           dealii::TriaAccessor<0, 1, spacedim>::left_vertex :
           ((i == 1) && cell.at_boundary(1) ?
              dealii::TriaAccessor<0, 1, spacedim>::right_vertex :
              dealii::TriaAccessor<0, 1, spacedim>::interior_vertex)),
        cell.vertex_index(i),
        &cell.get_dof_handler());
      return dealii::TriaIterator<
        dealii::DoFAccessor<0, dim, spacedim, level_dof_access>>(a);
    }


    template <int dim, int spacedim, bool level_dof_access>
    inline TriaIterator<
      dealii::DoFAccessor<dim - 1, dim, spacedim, level_dof_access>>
    get_face(
      const dealii::DoFCellAccessor<dim, spacedim, level_dof_access> &cell,
      const unsigned int                                              i,
      const std::integral_constant<int, 2>)
    {
      return cell.line(i);
    }


    template <int dim, int spacedim, bool level_dof_access>
    inline TriaIterator<
      dealii::DoFAccessor<dim - 1, dim, spacedim, level_dof_access>>
    get_face(
      const dealii::DoFCellAccessor<dim, spacedim, level_dof_access> &cell,
      const unsigned int                                              i,
      const std::integral_constant<int, 3>)
    {
      return cell.quad(i);
    }
  } // namespace DoFCellAccessorImplementation
} // namespace internal



template <int dimension_, int space_dimension_, bool level_dof_access>
inline typename DoFCellAccessor<dimension_,
                                space_dimension_,
                                level_dof_access>::face_iterator
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::face(
  const unsigned int i) const
{
  AssertIndexRange(i, this->n_faces());

  return dealii::internal::DoFCellAccessorImplementation::get_face(
    *this, i, std::integral_constant<int, dimension_>());
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline boost::container::small_vector<
  typename DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
    face_iterator,
  GeometryInfo<dimension_>::faces_per_cell>
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  face_iterators() const
{
  boost::container::small_vector<
    typename DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
      face_iterator,
    GeometryInfo<dimension_>::faces_per_cell>
    face_iterators(this->n_faces());

  for (const unsigned int i : this->face_indices())
    face_iterators[i] =
      dealii::internal::DoFCellAccessorImplementation::get_face(
        *this, i, std::integral_constant<int, dimension_>());

  return face_iterators;
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  get_active_or_mg_dof_indices(
    std::vector<types::global_dof_index> &dof_indices) const
{
  if (level_dof_access)
    get_mg_dof_indices(dof_indices);
  else
    get_dof_indices(dof_indices);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <class InputVector, typename number>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::get_dof_values(
  const InputVector &values,
  Vector<number>    &local_values) const
{
  get_dof_values(values, local_values.begin(), local_values.end());
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <typename Number, typename ForwardIterator>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::get_dof_values(
  const ReadVector<Number> &values,
  ForwardIterator           local_values_begin,
  ForwardIterator           local_values_end) const
{
  (void)local_values_end;
  Assert(this->is_artificial() == false,
         ExcMessage("Can't ask for DoF indices on artificial cells."));
  Assert(this->is_active(), ExcMessage("Cell must be active."));
  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());

  Assert(static_cast<unsigned int>(local_values_end - local_values_begin) ==
           this->get_fe().n_dofs_per_cell(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());
  Assert(values.size() == this->get_dof_handler().n_dofs(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());

  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(this->get_fe().n_dofs_per_cell());
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  boost::container::small_vector<Number, 27> values_temp(local_values_end -
                                                         local_values_begin);
  auto view = make_array_view(values_temp.begin(), values_temp.end());
  values.extract_subvector_to(make_array_view(dof_indices.begin(),
                                              dof_indices.end()),
                              view);
  using view_type = std::remove_reference_t<decltype(*local_values_begin)>;
  ArrayView<view_type> values_view2(&*local_values_begin,
                                    local_values_end - local_values_begin);
  std::copy(values_temp.begin(), values_temp.end(), values_view2.begin());
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <class InputVector, typename ForwardIterator>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::get_dof_values(
  const AffineConstraints<typename InputVector::value_type> &constraints,
  const InputVector                                         &values,
  ForwardIterator                                            local_values_begin,
  ForwardIterator local_values_end) const
{
  Assert(this->is_artificial() == false,
         ExcMessage("Can't ask for DoF indices on artificial cells."));
  Assert(this->is_active(), ExcMessage("Cell must be active."));

  Assert(static_cast<unsigned int>(local_values_end - local_values_begin) ==
           this->get_fe().n_dofs_per_cell(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());
  Assert(values.size() == this->get_dof_handler().n_dofs(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());


  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(this->get_fe().n_dofs_per_cell());
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  constraints.get_dof_values(values,
                             dof_indices.data(),
                             local_values_begin,
                             local_values_end);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <class OutputVector, typename number>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::set_dof_values(
  const Vector<number> &local_values,
  OutputVector         &values) const
{
  Assert(this->is_artificial() == false,
         ExcMessage("Can't ask for DoF indices on artificial cells."));
  Assert(this->is_active(), ExcMessage("Cell must be active."));

  Assert(static_cast<unsigned int>(local_values.size()) ==
           this->get_fe().n_dofs_per_cell(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());
  Assert(values.size() == this->get_dof_handler().n_dofs(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());


  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());
  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(this->get_fe().n_dofs_per_cell());
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  for (unsigned int i = 0; i < this->get_fe().n_dofs_per_cell(); ++i)
    internal::ElementAccess<OutputVector>::set(local_values(i),
                                               dof_indices[i],
                                               values);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline const FiniteElement<dimension_, space_dimension_> &
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::get_fe() const
{
  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           this->is_active(),
         ExcMessage(
           "For DoFHandler objects in hp-mode, finite elements are only "
           "associated with active cells. Consequently, you can not ask "
           "for the active finite element on cells with children."));

  const auto &fe = this->dof_handler->get_fe(active_fe_index());

  Assert(this->reference_cell() == fe.reference_cell(),
         internal::ExcNonMatchingReferenceCellTypes(this->reference_cell(),
                                                    fe.reference_cell()));

  return fe;
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline types::fe_index
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  active_fe_index() const
{
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           this->is_active(),
         ExcMessage(
           "You can not ask for the active FE index on a cell that has "
           "children because no degrees of freedom are assigned "
           "to this cell and, consequently, no finite element "
           "is associated with it."));
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->is_locally_owned() || this->is_ghost()),
         ExcMessage("You can only query active FE index information on cells "
                    "that are either locally owned or (after distributing "
                    "degrees of freedom) are ghost cells."));

  return dealii::internal::DoFCellAccessorImplementation::Implementation::
    active_fe_index(*this);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  set_active_fe_index(const types::fe_index i) const
{
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           this->is_active(),
         ExcMessage("You can not set the active FE index on a cell that has "
                    "children because no degrees of freedom will be assigned "
                    "to this cell."));

  Assert((this->dof_handler->hp_capability_enabled == false) ||
           this->is_locally_owned(),
         ExcMessage("You can only set active FE index information on cells "
                    "that are locally owned. On ghost cells, this information "
                    "will automatically be propagated from the owning process "
                    "of that cell, and there is no information at all on "
                    "artificial cells."));

  dealii::internal::DoFCellAccessorImplementation::Implementation::
    set_active_fe_index(*this, i);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline const FiniteElement<dimension_, space_dimension_> &
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::get_future_fe()
  const
{
  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           this->is_active(),
         ExcMessage(
           "For DoFHandler objects in hp-mode, finite elements are only "
           "associated with active cells. Consequently, you can not ask "
           "for the future finite element on cells with children."));

  return this->dof_handler->get_fe(future_fe_index());
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline types::fe_index
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  future_fe_index() const
{
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->has_children() == false),
         ExcMessage(
           "You can not ask for the future FE index on a cell that has "
           "children because no degrees of freedom are assigned "
           "to this cell and, consequently, no finite element "
           "is associated with it."));
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->is_locally_owned()),
         ExcMessage("You can only query future FE index information on cells "
                    "that are locally owned."));

  return dealii::internal::DoFCellAccessorImplementation::Implementation::
    future_fe_index(*this);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  set_future_fe_index(const types::fe_index i) const
{
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->has_children() == false),
         ExcMessage("You can not set the future FE index on a cell that has "
                    "children because no degrees of freedom will be assigned "
                    "to this cell."));

  Assert((this->dof_handler->hp_capability_enabled == false) ||
           this->is_locally_owned(),
         ExcMessage("You can only set future FE index information on cells "
                    "that are locally owned."));

  dealii::internal::DoFCellAccessorImplementation::Implementation::
    set_future_fe_index(*this, i);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline bool
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  future_fe_index_set() const
{
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->has_children() == false),
         ExcMessage(
           "You can not ask for the future FE index on a cell that has "
           "children because no degrees of freedom are assigned "
           "to this cell and, consequently, no finite element "
           "is associated with it."));
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->is_locally_owned()),
         ExcMessage("You can only query future FE index information on cells "
                    "that are locally owned."));

  return dealii::internal::DoFCellAccessorImplementation::Implementation::
    future_fe_index_set(*this);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  clear_future_fe_index() const
{
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->has_children() == false),
         ExcMessage(
           "You can not ask for the future FE index on a cell that has "
           "children because no degrees of freedom are assigned "
           "to this cell and, consequently, no finite element "
           "is associated with it."));
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->is_locally_owned()),
         ExcMessage("You can only query future FE index information on cells "
                    "that are locally owned."));

  dealii::internal::DoFCellAccessorImplementation::Implementation::
    clear_future_fe_index(*this);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <typename number, typename OutputVector>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  distribute_local_to_global(const Vector<number> &local_source,
                             OutputVector         &global_destination) const
{
  this->distribute_local_to_global(local_source.begin(),
                                   local_source.end(),
                                   global_destination);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <typename ForwardIterator, typename OutputVector>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  distribute_local_to_global(ForwardIterator local_source_begin,
                             ForwardIterator local_source_end,
                             OutputVector   &global_destination) const
{
  Assert(this->dof_handler != nullptr,
         (typename std::decay_t<decltype(*this)>::ExcInvalidObject()));
  Assert(static_cast<unsigned int>(local_source_end - local_source_begin) ==
           this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_destination.size(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));

  Assert(!this->has_children(), ExcMessage("Cell must be active"));

  Assert(
    internal::ArrayViewHelper::is_contiguous(local_source_begin,
                                             local_source_end),
    ExcMessage(
      "This function can not be called with iterator types that do not point to contiguous memory."));

  const unsigned int n_dofs = local_source_end - local_source_begin;

  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(n_dofs);
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  // distribute cell vector
  global_destination.add(n_dofs, dof_indices.data(), &(*local_source_begin));
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <typename ForwardIterator, typename OutputVector>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  distribute_local_to_global(
    const AffineConstraints<typename OutputVector::value_type> &constraints,
    ForwardIterator local_source_begin,
    ForwardIterator local_source_end,
    OutputVector   &global_destination) const
{
  Assert(this->dof_handler != nullptr,
         (typename std::decay_t<decltype(*this)>::ExcInvalidObject()));
  Assert(local_source_end - local_source_begin ==
           this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_destination.size(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));

  Assert(!this->has_children(), ExcMessage("Cell must be active."));

  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(this->get_fe().n_dofs_per_cell());
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  // distribute cell vector
  constraints.distribute_local_to_global(local_source_begin,
                                         local_source_end,
                                         dof_indices.data(),
                                         global_destination);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <typename number, typename OutputMatrix>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  distribute_local_to_global(const FullMatrix<number> &local_source,
                             OutputMatrix             &global_destination) const
{
  Assert(this->dof_handler != nullptr,
         (typename std::decay_t<decltype(*this)>::ExcInvalidObject()));
  Assert(local_source.m() == this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));
  Assert(local_source.n() == this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_destination.m(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_destination.n(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));

  Assert(!this->has_children(), ExcMessage("Cell must be active."));

  const unsigned int n_dofs = local_source.m();

  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(n_dofs);
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  // distribute cell matrix
  for (unsigned int i = 0; i < n_dofs; ++i)
    global_destination.add(dof_indices[i],
                           n_dofs,
                           dof_indices.data(),
                           &local_source(i, 0));
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <typename number, typename OutputMatrix, typename OutputVector>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  distribute_local_to_global(const FullMatrix<number> &local_matrix,
                             const Vector<number>     &local_vector,
                             OutputMatrix             &global_matrix,
                             OutputVector             &global_vector) const
{
  Assert(this->dof_handler != nullptr,
         (typename std::decay_t<decltype(*this)>::ExcInvalidObject()));
  Assert(local_matrix.m() == this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));
  Assert(local_matrix.n() == this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_matrix.m(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_matrix.n(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));
  Assert(local_vector.size() == this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_vector.size(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));

  Assert(!this->has_children(), ExcMessage("Cell must be active."));

  const unsigned int n_dofs = this->get_fe().n_dofs_per_cell();
  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(n_dofs);
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  // distribute cell matrices
  for (unsigned int i = 0; i < n_dofs; ++i)
    {
      global_matrix.add(dof_indices[i],
                        n_dofs,
                        dof_indices.data(),
                        &local_matrix(i, 0));
      global_vector(dof_indices[i]) += local_vector(i);
    }
}



DEAL_II_NAMESPACE_CLOSE


#endif
