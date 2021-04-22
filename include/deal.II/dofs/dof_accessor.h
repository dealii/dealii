// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_dof_accessor_h
#define dealii_dof_accessor_h


#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_iterator_selector.h>

#include <deal.II/grid/tria_accessor.h>

#include <deal.II/hp/dof_handler.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/container/small_vector.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

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

// note: the file dof_accessor.templates.h is included at the end of
// this file.  this includes a lot of templates and thus makes
// compilation slower, but at the same time allows for more aggressive
// inlining and thus faster code.


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
 * detail in connection to @ref Iterators.
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
  static const unsigned int dimension = dim;

  /**
   * A static variable that allows users of this class to discover the value
   * of the third template argument.
   */
  static const unsigned int space_dimension = spacedim;

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
              const DoFHandler<dim, spacedim> *   dof_handler);

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
  get_dof_indices(std::vector<types::global_dof_index> &dof_indices,
                  const unsigned int                    fe_index =
                    DoFHandler<dim, spacedim>::invalid_fe_index) const;

  /**
   * Return the global multilevel indices of the degrees of freedom that live
   * on the current object with respect to the given level within the
   * multigrid hierarchy. The indices refer to the local numbering for the
   * level this line lives on.
   */
  void
  get_mg_dof_indices(const int                             level,
                     std::vector<types::global_dof_index> &dof_indices,
                     const unsigned int                    fe_index =
                       DoFHandler<dim, spacedim>::invalid_fe_index) const;

  /**
   * Set the level DoF indices that are returned by get_mg_dof_indices.
   */
  void
  set_mg_dof_indices(
    const int                                   level,
    const std::vector<types::global_dof_index> &dof_indices,
    const unsigned int fe_index = DoFHandler<dim, spacedim>::invalid_fe_index);

  /**
   * Global DoF index of the <i>i</i> degree associated with the @p vertexth
   * vertex of the present cell.
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
   */
  types::global_dof_index
  vertex_dof_index(const unsigned int vertex,
                   const unsigned int i,
                   const unsigned int fe_index =
                     DoFHandler<dim, spacedim>::invalid_fe_index) const;

  /**
   * Return the global DoF index of the <code>i</code>th degree of freedom
   * associated with the <code>vertex</code>th vertex on level @p level. Also
   * see vertex_dof_index().
   */
  types::global_dof_index
  mg_vertex_dof_index(const int          level,
                      const unsigned int vertex,
                      const unsigned int i,
                      const unsigned int fe_index =
                        DoFHandler<dim, spacedim>::invalid_fe_index) const;

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
  dof_index(const unsigned int i,
            const unsigned int fe_index =
              DoFHandler<dim, spacedim>::invalid_fe_index) const;

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
   * Return the @p n-th active FE index on this object. For cells and all non-
   * hp-objects, there is only a single active FE index, so the argument must
   * be equal to zero. For lower-dimensional hp-objects, there are
   * n_active_fe_indices() active finite elements, and this function can be
   * queried for their indices.
   */
  unsigned int
  nth_active_fe_index(const unsigned int n) const;

  /**
   * Returns all active FE indices on this object.
   *
   * The size of the returned set equals the number of finite elements that
   * are active on this object.
   */
  std::set<unsigned int>
  get_active_fe_indices() const;

  /**
   * Return true if the finite element with given index is active on the
   * present object. When the current DoFHandler does not have hp-
   * capabilities, this is of course the case only if @p fe_index equals
   * zero. For cells, it is the case if @p fe_index equals active_fe_index()
   * of this cell. For faces and other lower- dimensional objects, there may
   * be more than one @p fe_index that are active on any given object (see
   * n_active_fe_indices()).
   */
  bool
  fe_index_is_active(const unsigned int fe_index) const;

  /**
   * Return a reference to the finite element used on this object with the
   * given @p fe_index. @p fe_index must be used on this object, i.e.
   * <code>fe_index_is_active(fe_index)</code> must return true.
   */
  const FiniteElement<dim, spacedim> &
  get_fe(const unsigned int fe_index) const;

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
  set_dof_index(const unsigned int            i,
                const types::global_dof_index index,
                const unsigned int            fe_index =
                  DoFHandler<dim, spacedim>::invalid_fe_index) const;

  void
  set_mg_dof_index(const int                     level,
                   const unsigned int            i,
                   const types::global_dof_index index) const;

  /**
   * Set the global index of the <i>i</i> degree on the @p vertex-th vertex of
   * the present cell to @p index.
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
  set_vertex_dof_index(const unsigned int            vertex,
                       const unsigned int            i,
                       const types::global_dof_index index,
                       const unsigned int            fe_index =
                         DoFHandler<dim, spacedim>::invalid_fe_index) const;

  void
  set_mg_vertex_dof_index(const int                     level,
                          const unsigned int            vertex,
                          const unsigned int            i,
                          const types::global_dof_index index,
                          const unsigned int            fe_index =
                            DoFHandler<dim, spacedim>::invalid_fe_index) const;

  // Iterator classes need to be friends because they need to access
  // operator== and operator!=.
  template <typename>
  friend class TriaRawIterator;
  template <int, int, int, bool>
  friend class DoFAccessor;

private:
  // Make the DoFHandler class a friend so that it can call the set_xxx()
  // functions.
  template <int, int>
  friend class DoFHandler;

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
  static const unsigned int dimension = 1;

  /**
   * A static variable that allows users of this class to discover the value
   * of the third template argument.
   */
  static const unsigned int space_dimension = spacedim;

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
   * Constructor to be used if the object here refers to a vertex of a one-
   * dimensional triangulation, i.e. a face of the triangulation.
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
    const Triangulation<1, spacedim> *                      tria,
    const typename TriaAccessor<0, 1, spacedim>::VertexKind vertex_kind,
    const unsigned int                                      vertex_index,
    const DoFHandler<1, spacedim> *                         dof_handler);

  /**
   * Constructor. This constructor exists in order to maintain interface
   * compatibility with the other accessor classes. However, it doesn't do
   * anything useful here and so may not actually be called.
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
    const unsigned int fe_index = AccessorData::invalid_fe_index) const;

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
    const unsigned int fe_index = AccessorData::invalid_fe_index) const;

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
    const unsigned int vertex,
    const unsigned int i,
    const unsigned int fe_index = AccessorData::invalid_fe_index) const;

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
  dof_index(const unsigned int i,
            const unsigned int fe_index = AccessorData::invalid_fe_index) const;

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
  unsigned int
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
  fe_index_is_active(const unsigned int fe_index) const;

  /**
   * Return a reference to the finite element used on this object with the
   * given @p fe_index. @p fe_index must be used on this object, i.e.
   * <code>fe_index_is_active(fe_index)</code> must return true.
   */
  const FiniteElement<1, spacedim> &
  get_fe(const unsigned int fe_index) const;

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
  void set_dof_handler(DoFHandler<1, spacedim> *dh);

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
    const unsigned int fe_index = AccessorData::invalid_fe_index) const;

  /**
   * Set the global index of the <i>i</i> degree on the @p vertex-th vertex of
   * the present cell to @p index.
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
  set_vertex_dof_index(
    const unsigned int            vertex,
    const unsigned int            i,
    const types::global_dof_index index,
    const unsigned int fe_index = AccessorData::invalid_fe_index) const;

  // Iterator classes need to be friends because they need to access
  // operator== and operator!=.
  template <typename>
  friend class TriaRawIterator;


  // Make the DoFHandler class a friend so that it can call the set_xxx()
  // functions.
  template <int, int>
  friend class DoFHandler;

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
  DoFInvalidAccessor(const Triangulation<dim, spacedim> *parent     = 0,
                     const int                           level      = -1,
                     const int                           index      = -1,
                     const AccessorData *                local_data = 0);

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
  dof_index(const unsigned int i,
            const unsigned int fe_index =
              DoFHandler<dim, spacedim>::default_fe_index) const;

  /**
   * Set the index of the <i>i</i>th degree of freedom of this object to @p
   * index. Since the current object doesn't point to anything useful, like
   * all other functions in this class this function only throws an exception.
   */
  void
  set_dof_index(const unsigned int            i,
                const types::global_dof_index index,
                const unsigned int            fe_index =
                  DoFHandler<dim, spacedim>::invalid_fe_index) const;
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
   * 1D, a DoFAccessor::line in 2D, and a DoFAccessor::quad in 3d.
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
  template <class InputVector, typename ForwardIterator>
  void
  get_dof_values(const InputVector &values,
                 ForwardIterator    local_values_begin,
                 ForwardIterator    local_values_end) const;

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
    const InputVector &                                        values,
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
                 OutputVector &        values) const;

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
  template <class InputVector, typename number>
  void
  get_interpolated_dof_values(
    const InputVector &values,
    Vector<number> &   interpolated_values,
    const unsigned int fe_index =
      DoFHandler<dimension_, space_dimension_>::invalid_fe_index) const;

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
   * Consequently, it carries a third argument that can be used in the hp-
   * context that denotes the finite element space we are supposed to
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
   * @note Unlike the get_dof_values() function, this function is only
   * available on cells, rather than on lines, quads, and hexes, since
   * interpolation is presently only provided for cells by the finite element
   * classes.
   */
  template <class OutputVector, typename number>
  void
  set_dof_values_by_interpolation(
    const Vector<number> &local_values,
    OutputVector &        values,
    const unsigned int    fe_index =
      DoFHandler<dimension_, space_dimension_>::invalid_fe_index) const;

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
                             OutputVector &        global_destination) const;

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
                             OutputVector &  global_destination) const;

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
    OutputVector &  global_destination) const;

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
                             const Vector<number> &    local_vector,
                             OutputMatrix &            global_matrix,
                             OutputVector &            global_vector) const;

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
  unsigned int
  active_fe_index() const;

  /**
   * Set the index of the FiniteElement used for this cell. This determines
   * which element in an hp::FECollection to use. This function is only useful
   * if the DoF handler object associated with the current cell has hp-
   * capabilities enabled.
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
  set_active_fe_index(const unsigned int i) const;
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
   * Update the cache in which we store the dof indices of this cell.
   */
  void
  update_cell_dof_indices_cache() const;

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
   * the same element, independent of the cell we are presently on, but for hp-
   * DoFHandler objects this may change from cell to cell.
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
  unsigned int
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
  set_future_fe_index(const unsigned int i) const;

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
  // Make the DoFHandler class a friend so that it can call the
  // update_cell_dof_indices_cache() function
  template <int, int>
  friend class DoFHandler;
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



DEAL_II_NAMESPACE_CLOSE

// include more templates
#include "dof_accessor.templates.h"


#endif
