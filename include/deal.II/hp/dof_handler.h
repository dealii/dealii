// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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

#ifndef dealii_hp_dof_handler_h
#define dealii_hp_dof_handler_h



#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/dofs/deprecated_function_map.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_iterator_selector.h>
#include <deal.II/dofs/number_cache.h>

#include <deal.II/hp/dof_faces.h>
#include <deal.II/hp/dof_level.h>
#include <deal.II/hp/fe_collection.h>

#include <map>
#include <memory>
#include <set>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class Triangulation;

namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim, typename VectorType>
    class CellDataTransfer;
  }
} // namespace parallel

namespace internal
{
  namespace DoFHandlerImplementation
  {
    struct Implementation;

    namespace Policy
    {
      template <int dim, int spacedim>
      class PolicyBase;
      struct Implementation;
    } // namespace Policy
  }   // namespace DoFHandlerImplementation

  namespace hp
  {
    class DoFLevel;

    namespace DoFHandlerImplementation
    {
      struct Implementation;
    }
  } // namespace hp
} // namespace internal

namespace internal
{
  namespace DoFAccessorImplementation
  {
    struct Implementation;
  }

  namespace DoFCellAccessorImplementation
  {
    struct Implementation;
  }
} // namespace internal
#endif


namespace hp
{
  /**
   * Manage the distribution and numbering of the degrees of freedom for hp-
   * FEM algorithms. This class satisfies the
   * @ref ConceptMeshType "MeshType concept"
   * requirements.
   *
   * The purpose of this class is to allow for an enumeration of degrees of
   * freedom in the same way as the ::DoFHandler class, but it allows to use a
   * different finite element on every cell. To this end, one assigns an
   * <code>active_fe_index</code> to every cell that indicates which element
   * within a collection of finite elements (represented by an object of type
   * hp::FECollection) is the one that lives on this cell. The class then
   * enumerates the degree of freedom associated with these finite elements on
   * each cell of a triangulation and, if possible, identifies degrees of
   * freedom at the interfaces of cells if they match. If neighboring cells
   * have degrees of freedom along the common interface that do not immediate
   * match (for example, if you have $Q_2$ and $Q_3$ elements meeting at a
   * common face), then one needs to compute constraints to ensure that the
   * resulting finite element space on the mesh remains conforming.
   *
   * The whole process of working with objects of this type is explained in
   * step-27. Many of the algorithms this class implements are described in
   * the
   * @ref hp_paper "hp paper".
   *
   *
   * <h3>Active FE indices and their behavior under mesh refinement</h3>
   *
   * The typical workflow for using this class is to create a mesh, assign an
   * active FE index to every active cell, calls
   * hp::DoFHandler::distribute_dofs(), and then assemble a linear system and
   * solve a problem on this finite element space. However, one can skip
   * assigning active FE indices upon mesh refinement in certain
   * circumstances. In particular, the following rules apply:
   * - Upon mesh refinement, child cells inherit the active FE index of
   *   the parent.
   * - When coarsening cells, the (now active) parent cell will be assigned
   *   an active FE index that is determined from its (no longer active)
   *   children, following the FiniteElementDomination logic: Out of the set of
   *   elements previously assigned to the former children, we choose the one
   *   dominated by all children for the parent cell. If none was found, we pick
   *   the most dominant element in the whole collection that is dominated by
   *   all former children. See hp::FECollection::find_dominated_fe_extended()
   *   for further information on this topic.
   *
   * @note Finite elements need to be assigned to each cell by either calling
   * set_fe() or distribute_dofs() first to make this functionality available.
   *
   *
   * <h3>Active FE indices and parallel meshes</h3>
   *
   * When this class is used with either a parallel::shared::Triangulation
   * or a parallel::distributed::Triangulation, you can only set active
   * FE indices on cells that are locally owned,
   * using a call such as <code>cell-@>set_active_fe_index(...)</code>.
   * On the other hand, setting the active FE index on ghost
   * or artificial cells is not allowed.
   *
   * Ghost cells do acquire the information what element
   * is active on them, however: whenever
   * you call hp::DoFHandler::distribute_dofs(), all processors that
   * participate in the parallel mesh exchange information in such a way
   * that the active FE index on ghost cells equals the active FE index
   * that was set on that processor that owned that particular ghost cell.
   * Consequently, one can <i>query</i> the @p active_fe_index on ghost
   * cells, just not set it by hand.
   *
   * On artificial cells, no information is available about the
   * @p active_fe_index used there. That's because we don't even know
   * whether these cells exist at all, and even if they did, the
   * current processor does not know anything specific about them.
   * See
   * @ref GlossArtificialCell "the glossary entry on artificial cells"
   * for more information.
   *
   * During refinement and coarsening, information about the @p active_fe_index
   * of each cell will be automatically transferred.
   *
   * However, using a parallel::distributed::Triangulation with an
   * hp::DoFHandler requires additional attention during serialization, since no
   * information on active FE indices will be automatically transferred. This
   * has to be done manually using the
   * prepare_for_serialization_of_active_fe_indices() and
   * deserialize_active_fe_indices() functions. The former has to be called
   * before parallel::distributed::Triangulation::save() is invoked, and the
   * latter needs to be run after parallel::distributed::Triangulation::load().
   * If further data will be attached to the triangulation via the
   * parallel::distributed::CellDataTransfer,
   * parallel::distributed::SolutionTransfer, or Particles::ParticleHandler
   * classes, all corresponding preparation and deserialization function calls
   * need to happen in the same order. Consult the documentation of
   * parallel::distributed::SolutionTransfer for more information.
   *
   *
   * @ingroup dofs
   * @ingroup hp
   *
   * @author Wolfgang Bangerth, 2003, 2004, 2017, 2018
   * @author Oliver Kayser-Herold, 2003, 2004
   * @author Marc Fehling, 2018
   */
  template <int dim, int spacedim = dim>
  class DoFHandler : public Subscriptor
  {
    using ActiveSelector = dealii::internal::DoFHandlerImplementation::
      Iterators<DoFHandler<dim, spacedim>, false>;
    using LevelSelector = dealii::internal::DoFHandlerImplementation::
      Iterators<DoFHandler<dim, spacedim>, true>;

  public:
    /**
     * An alias that is used to identify cell iterators in DoFHandler objects.
     * The concept of iterators is discussed at length in the
     * @ref Iterators "iterators documentation module".
     *
     * The current alias works, in essence, like the corresponding
     * Triangulation::cell_accessor alias. However, it also makes available
     * the member functions of DoFCellAccessor, in addition to the ones
     * already available through the CellAccessor class.
     *
     * @ingroup Iterators
     */
    using cell_accessor = typename ActiveSelector::CellAccessor;

    /**
     * An alias that is used to identify iterators that point to faces.
     * The concept of iterators is discussed at length in the
     * @ref Iterators "iterators documentation module".
     *
     * The current alias works, in essence, like the corresponding
     * Triangulation::face_accessor alias. However, it also makes available
     * the member functions of DoFAccessor, in addition to the ones
     * already available through the TriaAccessor class.
     *
     * @ingroup Iterators
     */
    using face_accessor = typename ActiveSelector::FaceAccessor;

    /**
     * An alias that defines an iterator over the (one-dimensional) lines
     * of a mesh. In one-dimensional meshes, these are the cells of the mesh,
     * whereas in two-dimensional meshes the lines are the faces of cells.
     *
     * @ingroup Iterators
     */
    using line_iterator = typename ActiveSelector::line_iterator;

    /**
     * An alias that allows iterating over the <i>active</i> lines, i.e.,
     * that subset of lines that have no children. In one-dimensional meshes,
     * these are the cells of the mesh, whereas in two-dimensional
     * meshes the lines are the faces of cells.
     *
     * In two- or three-dimensional meshes, lines without children (i.e.,
     * the active lines) are part of at least one active cell. Each such line
     * may additionally be a child of a line of a coarser cell adjacent to a
     * cell that is active. (This coarser neighbor would then also be active.)
     *
     * @ingroup Iterators
     */
    using active_line_iterator = typename ActiveSelector::active_line_iterator;

    /**
     * An alias that defines an iterator over the (two-dimensional) quads
     * of a mesh. In two-dimensional meshes, these are the cells of the mesh,
     * whereas in three-dimensional meshes the quads are the faces of cells.
     *
     * @ingroup Iterators
     */
    using quad_iterator = typename ActiveSelector::quad_iterator;

    /**
     * An alias that allows iterating over the <i>active</i> quads, i.e.,
     * that subset of quads that have no children. In two-dimensional meshes,
     * these are the cells of the mesh, whereas in three-dimensional
     * meshes the quads are the faces of cells.
     *
     * In three-dimensional meshes, quads without children (i.e.,
     * the active quads) are faces of at least one active cell. Each such quad
     * may additionally be a child of a quad face of a coarser cell adjacent to
     * a cell that is active. (This coarser neighbor would then also be active.)
     *
     * @ingroup Iterators
     */
    using active_quad_iterator = typename ActiveSelector::active_quad_iterator;

    /**
     * An alias that defines an iterator over the (three-dimensional) hexes
     * of a mesh. This iterator only makes sense in three-dimensional meshes,
     * where hexes are the cells of the mesh.
     *
     * @ingroup Iterators
     */
    using hex_iterator = typename ActiveSelector::hex_iterator;

    /**
     * An alias that allows iterating over the <i>active</i> hexes of a mesh.
     * This iterator only makes sense in three-dimensional meshes,
     * where hexes are the cells of the mesh. Consequently, in these
     * three-dimensional meshes, this iterator is equivalent to the
     * @p active_cell_iterator alias.
     *
     * @ingroup Iterators
     */
    using active_hex_iterator = typename ActiveSelector::active_hex_iterator;

    /**
     * @copydoc ::DoFHandler::active_cell_iterator
     * @ingroup Iterators
     */
#ifndef _MSC_VER
    using active_cell_iterator = typename ActiveSelector::active_cell_iterator;
#else
    using active_cell_iterator = TriaActiveIterator<
      dealii::DoFCellAccessor<DoFHandler<dim, spacedim>, false>>;
#endif

    using level_cell_iterator = typename LevelSelector::cell_iterator;

    /**
     * @copydoc ::DoFHandler::cell_iterator
     * @ingroup Iterators
     */
#ifndef _MSC_VER
    using cell_iterator = typename ActiveSelector::cell_iterator;
#else
    using cell_iterator =
      TriaIterator<dealii::DoFCellAccessor<DoFHandler<dim, spacedim>, false>>;
#endif

    /**
     * @copydoc ::DoFHandler::face_iterator
     * @ingroup Iterators
     */
    using face_iterator = typename ActiveSelector::face_iterator;

    /**
     * @copydoc ::DoFHandler::active_face_iterator
     * @ingroup Iterators
     */
    using active_face_iterator = typename ActiveSelector::active_face_iterator;

    using level_cell_accessor = typename LevelSelector::CellAccessor;
    using level_face_accessor = typename LevelSelector::FaceAccessor;

    using level_face_iterator = typename LevelSelector::face_iterator;

    /**
     * Make the dimension available in function templates.
     */
    static const unsigned int dimension = dim;

    /**
     * Make the space dimension available in function templates.
     */
    static const unsigned int space_dimension = spacedim;

    /**
     * Make the type of this DoFHandler available in function templates.
     */
    static const bool is_hp_dof_handler = true;

    /**
     * The default index of the finite element to be used on a given cell. For
     * the usual, non-hp dealii::DoFHandler class that only supports the same
     * finite element to be used on all cells, the index of the finite element
     * needs to be the same on all cells anyway, and by convention we pick
     * zero for this value. The situation here is different, since the hp
     * classes support the case where different finite element indices may be
     * used on different cells. The default index consequently corresponds to
     * an invalid value.
     */
    static const unsigned int default_fe_index = numbers::invalid_unsigned_int;


    /**
     * Default Constructor.
     */
    DoFHandler();

    /**
     * Constructor. Take @p tria as the triangulation to work on.
     */
    DoFHandler(const Triangulation<dim, spacedim> &tria);

    /**
     * Copy constructor. DoFHandler objects are large and expensive.
     * They should not be copied, in particular not by accident, but
     * rather deliberately constructed. As a consequence, this constructor
     * is explicitly removed from the interface of this class.
     */
    DoFHandler(const DoFHandler &) = delete;

    /**
     * Destructor.
     */
    virtual ~DoFHandler() override;

    /**
     * Copy operator. DoFHandler objects are large and expensive.
     * They should not be copied, in particular not by accident, but
     * rather deliberately constructed. As a consequence, this operator
     * is explicitly removed from the interface of this class.
     */
    DoFHandler &
    operator=(const DoFHandler &) = delete;

    /**
     * Assign a Triangulation and a FECollection to the DoFHandler and compute
     * the distribution of degrees of freedom over the mesh.
     */
    void
    initialize(const Triangulation<dim, spacedim> &   tria,
               const hp::FECollection<dim, spacedim> &fe);

    /**
     * Assign a hp::FECollection @p fe to this object.
     *
     * In case a parallel::TriangulationBase is assigned to this object,
     * the active_fe_indices will be exchanged between processors so that
     * each one knows the indices on its own cells and all ghost cells.
     *
     * @note In accordance with dealii::DoFHandler::set_fe(),
     * this function also makes a copy of the object given as argument.
     *
     * @warning This function only sets a hp::FECollection. Degrees of freedom
     * have either not been distributed yet, or are distributed using a
     * previously set collection. In both cases, accessing degrees of freedom
     * will lead to invalid results. To restore consistency, call
     * distribute_dofs().
     */
    virtual void
    set_fe(const hp::FECollection<dim, spacedim> &fe);

    /**
     * Go through the triangulation and "distribute" the degrees of
     * freedom needed for the given finite element. "Distributing"
     * degrees of freedom involves allocating memory to store the
     * indices on all entities on which degrees of freedom can be
     * located (e.g., vertices, edges, faces, etc.) and to then enumerate
     * all degrees of freedom. In other words, while the mesh and the
     * finite element object by themselves simply define a finite
     * element space $V_h$, the process of distributing degrees of
     * freedom makes sure that there is a basis for this space and that
     * the shape functions of this basis are enumerated in an indexable,
     * predictable way.
     *
     * The exact order in which degrees of freedom on a mesh are
     * ordered, i.e., the order in which basis functions of the finite
     * element space are enumerated, is something that deal.II treats as
     * an implementation detail. By and large, degrees of freedom are
     * enumerated in the same order in which we traverse cells, but you
     * should not rely on any specific numbering. In contrast, if you
     * want a particular ordering, use the functions in namespace
     * DoFRenumbering.
     *
     * @note In accordance with dealii::DoFHandler::distribute_dofs(),
     * this function also makes a copy of the object given as argument.
     */
    virtual void
    distribute_dofs(const hp::FECollection<dim, spacedim> &fe);

    /**
     * Go through the triangulation and set the active FE indices of all
     * active cells to the values given in @p active_fe_indices.
     */
    void
    set_active_fe_indices(const std::vector<unsigned int> &active_fe_indices);

    /**
     * Go through the triangulation and store the active FE indices of all
     * active cells to the vector @p active_fe_indices. This vector is
     * resized, if necessary.
     */
    void
    get_active_fe_indices(std::vector<unsigned int> &active_fe_indices) const;

    /**
     * Clear all data of this object and especially delete the lock this
     * object has to the finite element used the last time when @p
     * distribute_dofs was called.
     */
    virtual void
    clear();

    /**
     * Renumber degrees of freedom based on a list of new DoF indices for each
     * of the degrees of freedom.
     *
     * This function is called by the functions in DoFRenumbering function after
     * computing a new ordering of the degree of freedom indices. However, it
     * can of course also be called from user code.
     *
     * @arg new_number This array must have a size equal to the number of
     * degrees of freedom owned by the current processor, i.e. the size must be
     * equal to what n_locally_owned_dofs() returns. If only one processor
     * participates in storing the current mesh, then this equals the total
     * number of degrees of freedom, i.e. the result of n_dofs(). The contents
     * of this array are the new global indices for each freedom listed in the
     * IndexSet returned by locally_owned_dofs(). In the case of a sequential
     * mesh this means that the array is a list of new indices for each of the
     * degrees of freedom on the current mesh. In the case that we have a
     * parallel::shared::Triangulation or
     * parallel::distributed::Triangulation underlying this DoFHandler object,
     * the array is a list of new indices for all the locally owned degrees of
     * freedom, enumerated in the same order as the currently locally owned
     * DoFs. In other words, assume that degree of freedom <code>i</code> is
     * currently locally owned, then
     * <code>new_numbers[locally_owned_dofs().index_within_set(i)]</code>
     * returns the new global DoF index of <code>i</code>. Since the IndexSet of
     * locally_owned_dofs() is complete in the sequential case, the latter
     * convention for the content of the array reduces to the former in the case
     * that only one processor participates in the mesh.
     *
     * @note While it follows from the above, it may be surprising to know that
     *   the <i>number</i> of locally owned degrees of freedom in a parallel
     *   computation is an invariant
     *   under renumbering, even if the <i>indices</i> associated with these
     *   locally owned degrees of freedom are not. At a fundamental level,
     *   this invariant exists because the <i>decision</i> whether a degree of
     *   freedom is locally owned or not has nothing to do with that
     *   degree of freedom's (old or new) index. Indeed, degrees of freedom
     *   are locally owned if they are on a locally owned cell and not on
     *   an interface between cells where the neighboring cell has a lower
     *   @ref GlossSubdomainId "subdomain id". Since both of these conditions
     *   are independent of the index associated with the DoF, a locally
     *   owned degree of freedom will also be locally owned after renumbering.
     *   On the other hand, properties such as whether the set of indices
     *   of locally owned DoFs forms a contiguous range or not
     *   (i.e., whether the locally_owned_dofs() returns an IndexSet object
     *   for which IndexSet::is_contiguous() returns @p true) are of
     *   course affected by the exact renumbering performed here. For example,
     *   while the initial numbering of DoF indices done in distribute_dofs()
     *   yields a contiguous numbering, the renumbering performed by
     *   DoFRenumbering::component_wise() will, in general, not yield
     *   contiguous locally owned DoF indices.
     */
    void
    renumber_dofs(const std::vector<types::global_dof_index> &new_numbers);

    /**
     * Return the maximum number of degrees of freedom a degree of freedom in
     * the given triangulation with the given finite element may couple with.
     * This is the maximum number of entries per line in the system matrix;
     * this information can therefore be used upon construction of the
     * SparsityPattern object.
     *
     * The returned number is not really the maximum number but an estimate
     * based on the finite element and the maximum number of cells meeting at
     * a vertex. The number holds for the constrained matrix also.
     *
     * As for ::DoFHandler::max_couplings_between_dofs(), the result of this
     * function is often not very accurate for 3d and/or high polynomial
     * degrees. The consequences are discussed in the documentation of the
     * module on
     * @ref Sparsity.
     */
    unsigned int
    max_couplings_between_dofs() const;

    /**
     * Return the number of degrees of freedom located on the boundary another
     * dof on the boundary can couple with.
     *
     * The number is the same as for @p max_coupling_between_dofs in one
     * dimension less.
     *
     * @note The same applies to this function as to max_couplings_per_dofs()
     * as regards the performance of this function. Think about one of the
     * dynamic sparsity pattern classes instead (see
     * @ref Sparsity).
     */
    unsigned int
    max_couplings_between_boundary_dofs() const;

    /**
     * @name Cell iterator functions
     */
    /*@{*/
    /**
     * Iterator to the first used cell on level @p level.
     */
    cell_iterator
    begin(const unsigned int level = 0) const;

    /**
     * Iterator to the first active cell on level @p level. If the given level
     * does not contain any active cells (i.e., all cells on this level are
     * further refined, then this function returns
     * <code>end_active(level)</code> so that loops of the kind
     * @code
     *   for (cell=dof_handler.begin_active(level);
     *        cell!=dof_handler.end_active(level);
     *        ++cell)
     *     {
     *       ...
     *     }
     * @endcode
     * have zero iterations, as may be expected if there are no active cells
     * on this level.
     */
    active_cell_iterator
    begin_active(const unsigned int level = 0) const;

    /**
     * Iterator past the end; this iterator serves for comparisons of
     * iterators with past-the-end or before-the-beginning states.
     */
    cell_iterator
    end() const;

    /**
     * Return an iterator which is the first iterator not on level. If @p
     * level is the last level, then this returns <tt>end()</tt>.
     */
    cell_iterator
    end(const unsigned int level) const;

    /**
     * Return an active iterator which is the first active iterator not on the
     * given level. If @p level is the last level, then this returns
     * <tt>end()</tt>.
     */
    active_cell_iterator
    end_active(const unsigned int level) const;

    /**
     * @name Cell iterator functions returning ranges of iterators
     */

    /**
     * Return an iterator range that contains all cells (active or not) that
     * make up this DoFHandler. Such a range is useful to initialize range-
     * based for loops as supported by C++11. See the example in the
     * documentation of active_cell_iterators().
     *
     * @return The half open range <code>[this->begin(), this->end())</code>
     *
     * @ingroup CPP11
     */
    IteratorRange<cell_iterator>
    cell_iterators() const;

    /**
     * Return an iterator range that contains all active cells that make up
     * this DoFHandler. Such a range is useful to initialize range-based for
     * loops as supported by C++11, see also
     * @ref CPP11 "C++11 standard".
     *
     * Range-based for loops are useful in that they require much less code
     * than traditional loops (see <a
     * href="http://en.wikipedia.org/wiki/C%2B%2B11#Range-based_for_loop">here</a>
     * for a discussion of how they work). An example is that without
     * range-based for loops, one often writes code such as the following:
     * @code
     *   DoFHandler<dim> dof_handler;
     *   ...
     *   typename DoFHandler<dim>::active_cell_iterator
     *     cell = dof_handler.begin_active(),
     *     endc = dof_handler.end();
     *   for (; cell!=endc; ++cell)
     *     {
     *       fe_values.reinit (cell);
     *       ...do the local integration on 'cell'...;
     *     }
     * @endcode
     * Using C++11's range-based for loops, this is now entirely equivalent to
     * the following:
     * @code
     *   DoFHandler<dim> dof_handler;
     *   ...
     *   for (const auto &cell : dof_handler.active_cell_iterators())
     *     {
     *       fe_values.reinit (cell);
     *       ...do the local integration on 'cell'...;
     *     }
     * @endcode
     *
     * @return The half open range <code>[this->begin_active(),
     * this->end())</code>
     *
     * @ingroup CPP11
     */
    IteratorRange<active_cell_iterator>
    active_cell_iterators() const;

    /**
     * Return an iterator range that contains all cells (active or not) that
     * make up the given level of this DoFHandler. Such a range is useful to
     * initialize range-based for loops as supported by C++11. See the example
     * in the documentation of active_cell_iterators().
     *
     * @param[in] level A given level in the refinement hierarchy of this
     * triangulation.
     * @return The half open range <code>[this->begin(level),
     * this->end(level))</code>
     *
     * @pre level must be less than this->n_levels().
     *
     * @ingroup CPP11
     */
    IteratorRange<cell_iterator>
    cell_iterators_on_level(const unsigned int level) const;

    /**
     * Return an iterator range that contains all active cells that make up
     * the given level of this DoFHandler. Such a range is useful to
     * initialize range-based for loops as supported by C++11. See the example
     * in the documentation of active_cell_iterators().
     *
     * @param[in] level A given level in the refinement hierarchy of this
     * triangulation.
     * @return The half open range <code>[this->begin_active(level),
     * this->end(level))</code>
     *
     * @pre level must be less than this->n_levels().
     *
     * @ingroup CPP11
     */
    IteratorRange<active_cell_iterator>
    active_cell_iterators_on_level(const unsigned int level) const;

    /*
     * @}
     */

    /*---------------------------------------*/


    /**
     * Return the global number of degrees of freedom. If the current object
     * handles all degrees of freedom itself (even if you may intend to solve
     * your linear system in parallel, such as in step-17 or step-18), then
     * this number equals the number of locally owned degrees of freedom since
     * this object doesn't know anything about what you want to do with it and
     * believes that it owns every degree of freedom it knows about.
     *
     * On the other hand, if this object operates on a
     * parallel::distributed::Triangulation object, then this function returns
     * the global number of degrees of freedom, accumulated over all
     * processors.
     *
     * In either case, included in the returned number are those DoFs which
     * are constrained by hanging nodes, see
     * @ref constraints.
     *
     * Mathematically speaking, the number returned by this function equals the
     * dimension of the finite element space (without taking into account
     * constraints) that corresponds to (i) the mesh on which it is defined,
     * and (ii) the finite element that is used by the current object. It
     * also, of course, equals the number of shape functions that span this
     * space.
     */
    types::global_dof_index
    n_dofs() const;

    /**
     * The number of multilevel dofs on given level. Since hp::DoFHandler does
     * not support multilevel methods yet, this function throws an exception
     * ExcNotImplemented() independent of its argument.
     */
    types::global_dof_index
    n_dofs(const unsigned int level) const;

    /**
     * Return the number of locally owned degrees of freedom located on the
     * boundary.
     */
    types::global_dof_index
    n_boundary_dofs() const;

    /**
     * Return the number of degrees of freedom located on those parts of the
     * boundary which have a boundary indicator listed in the given set. The
     * reason that a @p map rather than a @p set is used is the same as
     * described in the documentation of that variant of
     * DoFTools::make_boundary_sparsity_pattern() that takes a map.
     *
     * There is, however, another overload of this function that takes
     * a @p set argument (see below).
     */
    template <typename number>
    types::global_dof_index
    n_boundary_dofs(
      const std::map<types::boundary_id, const Function<spacedim, number> *>
        &boundary_ids) const;

    /**
     * Return the number of locally owned degrees of freedom located on those
     * parts of the boundary which have a boundary indicator listed in the given
     * set.
     */
    types::global_dof_index
    n_boundary_dofs(const std::set<types::boundary_id> &boundary_ids) const;

    /**
     * Return the number of degrees of freedom that belong to this process.
     *
     * If this is a sequential DoFHandler, then the result equals that produced
     * by n_dofs(). (Here, "sequential" means that either the whole program does
     * not use MPI, or that it uses MPI but only uses a single MPI process, or
     * that there are multiple MPI processes but the Triangulation on which this
     * DoFHandler builds works only on one MPI process.) On the other hand, if
     * we are operating on a parallel::distributed::Triangulation or
     * parallel::shared::Triangulation, then it includes only the degrees of
     * freedom that the current processor owns. Note that in this case this does
     * not include all degrees of freedom that have been distributed on the
     * current processor's image of the mesh: in particular, some of the degrees
     * of freedom on the interface between the cells owned by this processor and
     * cells owned by other processors may be theirs, and degrees of freedom on
     * ghost cells are also not necessarily included.
     */
    types::global_dof_index
    n_locally_owned_dofs() const;

    /**
     * Return an IndexSet describing the set of locally owned DoFs as a subset
     * of 0..n_dofs(). The number of elements of this set equals
     * n_locally_owned_dofs().
     */
    const IndexSet &
    locally_owned_dofs() const;

    /**
     * Compute a vector with the locally owned DoFs of each processor.
     *
     * This function involves global communication via the @p MPI_Allgather
     * function, so it must be called on all processors participating in the MPI
     * communicator underlying the triangulation.
     *
     * If you are only interested in the number of elements each processor owns
     * then compute_n_locally_owned_dofs_per_processor() is a better choice.
     *
     * If this is a sequential DoFHandler, then the vector has a single element
     * that equals the IndexSet representing the entire range [0,n_dofs()].
     * (Here, "sequential" means that either the whole program does not use MPI,
     * or that it uses MPI but only uses a single MPI process, or that there are
     * multiple MPI processes but the Triangulation on which this DoFHandler
     * builds works only on one MPI process.)
     */
    std::vector<IndexSet>
    compute_locally_owned_dofs_per_processor() const;

    /**
     * Compute a vector with the number of degrees of freedom each
     * processor that participates in this triangulation owns locally. The sum
     * of all these numbers equals the number of degrees of freedom that exist
     * globally, i.e. what n_dofs() returns.
     *
     * This function involves global communication via the @p MPI_Allgather
     * function, so it must be called on all processors participating in the MPI
     * communicator underlying the triangulation.
     *
     * Each element of the vector returned by this function equals the number of
     * elements of the corresponding sets returned by
     * compute_locally_owned_dofs_per_processor().
     *
     * If this is a sequential DoFHandler, then the vector has a single element
     * equal to n_dofs(). (Here, "sequential" means that either the whole
     * program does not use MPI, or that it uses MPI but only uses a single MPI
     * process, or that there are multiple MPI processes but the Triangulation
     * on which this DoFHandler builds works only on one MPI process.)
     */
    std::vector<types::global_dof_index>
    compute_n_locally_owned_dofs_per_processor() const;

    /**
     * Return a vector that stores the locally owned DoFs of each processor.
     *
     * @deprecated As of deal.II version 9.2, we do not populate a vector with
     * the index sets of all processors by default any more due to a possibly
     * large memory footprint on many processors. As a consequence, this
     * function needs to call compute_locally_owned_dofs_per_processor() upon
     * the first invocation, including global communication. Use
     * compute_locally_owned_dofs_per_processor() instead if using up to a few
     * thousands of MPI ranks or some variant involving local communication with
     * more processors.
     */
    DEAL_II_DEPRECATED const std::vector<IndexSet> &
                             locally_owned_dofs_per_processor() const;

    /**
     * Return a vector that stores the number of degrees of freedom each
     * processor that participates in this triangulation owns locally. The sum
     * of all these numbers equals the number of degrees of freedom that exist
     * globally, i.e. what n_dofs() returns.
     *
     * @deprecated As of deal.II version 9.2, we do not populate a vector with
     * the numbers of dofs of all processors by default any more due to a
     * possibly large memory footprint on many processors. As a consequence,
     * this function needs to call compute_n_locally_owned_dofs_per_processor()
     * upon the first invocation, including global communication. Use
     * compute_n_locally_owned_dofs_per_processor() instead if using up to a few
     * thousands of MPI ranks or some variant involving local communication with
     * more processors.
     */
    DEAL_II_DEPRECATED const std::vector<types::global_dof_index> &
                             n_locally_owned_dofs_per_processor() const;

    /**
     * Return an IndexSet describing the set of locally owned DoFs used for
     * the given multigrid level. Since hp::DoFHandler does not support
     * multilevel methods yet, this function throws an exception
     * ExcNotImplemented() independent of its argument.
     */
    const IndexSet &
    locally_owned_mg_dofs(const unsigned int level) const;

    /**
     * Compute a vector with the locally owned DoFs of each processor on
     * the given level @p level for geometric multigrid.
     *
     * This function involves global communication via the @p MPI_Allgather
     * function, so it must be called on all processors participating in the MPI
     * communicator underlying the triangulation.
     *
     * If this is a sequential DoFHandler, then the vector has a single element
     * that equals the IndexSet representing the entire range [0,n_dofs()].
     * (Here, "sequential" means that either the whole program does not use MPI,
     * or that it uses MPI but only uses a single MPI process, or that there are
     * multiple MPI processes but the Triangulation on which this DoFHandler
     * builds works only on one MPI process.)
     */
    std::vector<IndexSet>
    compute_locally_owned_mg_dofs_per_processor(const unsigned int level) const;

    /**
     * Return a vector that stores the locally owned DoFs of each processor on
     * the given level @p level.
     *
     * @deprecated As of deal.II version 9.2, we do not populate a vector with
     * the index sets of all processors by default any more due to a possibly
     * large memory footprint on many processors. As a consequence, this
     * function needs to call compute_locally_owned_dofs_mg_per_processor() upon
     * the first invocation, including global communication. Use
     * compute_locally_owned_mg_dofs_per_processor() instead if using up to a
     * few thousands of MPI ranks or some variant involving local communication
     * with more processors.
     */
    DEAL_II_DEPRECATED const std::vector<IndexSet> &
                             locally_owned_mg_dofs_per_processor(const unsigned int level) const;

    /**
     * Return a constant reference to the indexth finite element object that is
     * used by this @p DoFHandler.
     */
    const FiniteElement<dim, spacedim> &
    get_fe(const unsigned int index) const;

    /**
     * Return a constant reference to the set of finite element objects that
     * are used by this @p DoFHandler.
     */
    const hp::FECollection<dim, spacedim> &
    get_fe_collection() const;

    /**
     * Return a constant reference to the triangulation underlying this
     * object.
     */
    const Triangulation<dim, spacedim> &
    get_triangulation() const;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     *
     * This function is made virtual, since a dof handler object might be
     * accessed through a pointers to this base class, although the actual
     * object might be a derived class.
     */
    virtual std::size_t
    memory_consumption() const;

    /**
     * Whenever serialization with a parallel::distributed::Triangulation as the
     * underlying triangulation is considered, we also need to consider storing
     * the active_fe_indices on all active cells as well.
     *
     * This function registers that these indices are to be stored whenever the
     * parallel::distributed::Triangulation::save() function is called on the
     * underlying triangulation.
     *
     * @note Currently only implemented for triangulations of type
     *   parallel::distributed::Triangulation. An assertion will be triggered if
     *   a different type is registered.
     *
     * @see The documentation of parallel::distributed::SolutionTransfer has further
     *   information on serialization.
     */
    void
    prepare_for_serialization_of_active_fe_indices();

    /**
     * Whenever serialization with a parallel::distributed::Triangulation as the
     * underlying triangulation is considered, we also need to consider storing
     * the active_fe_indices on all active cells as well.
     *
     * This function deserializes and distributes the previously stored
     * active_fe_indices on all active cells.
     *
     * @note Currently only implemented for triangulations of type
     *   parallel::distributed::Triangulation. An assertion will be triggered if
     *   a different type is registered.
     *
     * @see The documentation of parallel::distributed::SolutionTransfer has further
     *   information on serialization.
     */
    void
    deserialize_active_fe_indices();

    /**
     * Write the data of this object to a stream for the purpose of
     * serialization.
     */
    template <class Archive>
    void
    save(Archive &ar, const unsigned int version) const;

    /**
     * Read the data of this object from a stream for the purpose of
     * serialization.
     */
    template <class Archive>
    void
    load(Archive &ar, const unsigned int version);

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    /**
     * Exception
     */
    DeclException0(ExcNoFESelected);
    /**
     * Exception
     */
    DeclException0(ExcGridsDoNotMatch);
    /**
     * Exception
     */
    DeclException0(ExcInvalidBoundaryIndicator);
    /**
     * Exception
     */
    DeclException1(ExcMatrixHasWrongSize,
                   int,
                   << "The matrix has the wrong dimension " << arg1);
    /**
     * Exception
     */
    DeclException0(ExcFunctionNotUseful);
    /**
     * Exception
     */
    DeclException1(ExcNewNumbersNotConsecutive,
                   types::global_dof_index,
                   << "The given list of new dof indices is not consecutive: "
                   << "the index " << arg1 << " does not exist.");
    /**
     * Exception
     */
    DeclException2(ExcInvalidFEIndex,
                   int,
                   int,
                   << "The mesh contains a cell with an active_fe_index of "
                   << arg1 << ", but the finite element collection only has "
                   << arg2 << " elements");
    /**
     * Exception
     */
    DeclException1(ExcInvalidLevel,
                   int,
                   << "The given level " << arg1
                   << " is not in the valid range!");
    /**
     * Exception
     */
    DeclException0(ExcFacesHaveNoLevel);
    /**
     * The triangulation level you accessed is empty.
     */
    DeclException1(ExcEmptyLevel,
                   int,
                   << "You tried to do something on level " << arg1
                   << ", but this level is empty.");

  private:
    /**
     * Address of the triangulation to work on.
     */
    SmartPointer<const Triangulation<dim, spacedim>, DoFHandler<dim, spacedim>>
      tria;

    /**
     * Store a copy of the finite element set given latest to distribute_dofs().
     */
    hp::FECollection<dim, spacedim> fe_collection;

    /**
     * An object that describes how degrees of freedom should be distributed and
     * renumbered.
     */
    std::unique_ptr<dealii::internal::DoFHandlerImplementation::Policy::
                      PolicyBase<dim, spacedim>>
      policy;


    /**
     * Setup policy and listeners based on the underlying Triangulation.
     */
    void
    setup_policy_and_listeners();

    /**
     * Free all used memory.
     */
    void
    clear_space();

    template <int structdim>
    types::global_dof_index
    get_dof_index(const unsigned int obj_level,
                  const unsigned int obj_index,
                  const unsigned int fe_index,
                  const unsigned int local_index) const;

    template <int structdim>
    void
    set_dof_index(const unsigned int            obj_level,
                  const unsigned int            obj_index,
                  const unsigned int            fe_index,
                  const unsigned int            local_index,
                  const types::global_dof_index global_index) const;

    /**
     * Create default tables for the active_fe_indices in the
     * dealii::internal::hp::DoFLevel. They are initialized with a zero
     * indicator, meaning that fe[0] is going to be used by default. This
     * method is called before refinement and while setting the finite elements
     * via set_fe(). It ensures each cell has a valid active_fe_index.
     */
    void
    create_active_fe_table();

    /**
     * A function that will be triggered through a triangulation
     * signal just before the triangulation is modified.
     *
     * The function that stores the active_fe_flags of all cells that will
     * be refined or coarsened before the refinement happens, so that
     * they can be set again after refinement.
     */
    void
    pre_refinement_action();

    /**
     * A function that will be triggered through a triangulation
     * signal just after the triangulation is modified.
     *
     * The function that restores the active_fe_flags of all cells that
     * were refined.
     */
    void
    post_refinement_action();

    /**
     * A function that will be triggered through a triangulation
     * signal just before the associated Triangulation or
     * parallel::shared::Triangulation is modified.
     *
     * The function that stores the active_fe_indices of all cells that will
     * be refined or coarsened before the refinement happens, so that
     * they can be set again after refinement.
     */
    void
    pre_active_fe_index_transfer();

    /**
     * A function that will be triggered through a triangulation
     * signal just before the associated parallel::distributed::Triangulation is
     * modified.
     *
     * The function that stores all active_fe_indices on locally owned cells for
     * distribution over all participating processors.
     */
    void
    pre_distributed_active_fe_index_transfer();

    /**
     * A function that will be triggered through a triangulation
     * signal just after the associated Triangulation or
     * parallel::shared::Triangulation is modified.
     *
     * The function that restores the active_fe_indices of all cells that
     * were refined or coarsened.
     */
    void
    post_active_fe_index_transfer();

    /**
     * A function that will be triggered through a triangulation
     * signal just after the associated parallel::distributed::Triangulation is
     * modified.
     *
     * The function that restores all active_fe_indices on locally owned cells
     * that have been communicated.
     */
    void
    post_distributed_active_fe_index_transfer();

    /**
     * A function that will be triggered through a triangulation
     * signal just after the associated parallel::distributed::Triangulation has
     * been saved.
     *
     * The function frees all memory related to the transfer of
     * active_fe_indices.
     */
    void
    post_distributed_serialization_of_active_fe_indices();

    /**
     * Space to store the DoF numbers for the different levels. Analogous to
     * the <tt>levels[]</tt> tree of the Triangulation objects.
     */
    std::vector<std::unique_ptr<dealii::internal::hp::DoFLevel>> levels;

    /**
     * Space to store the DoF numbers for the faces. Analogous to the
     * <tt>faces</tt> pointer of the Triangulation objects.
     */
    std::unique_ptr<dealii::internal::hp::DoFIndicesOnFaces<dim>> faces;

    /**
     * A structure that contains all sorts of numbers that characterize the
     * degrees of freedom this object works on.
     *
     * For most members of this structure, there is an accessor function in
     * this class that returns its value.
     */
    dealii::internal::DoFHandlerImplementation::NumberCache number_cache;

    /**
     * A structure that contains all sorts of numbers that characterize the
     * degrees of freedom on multigrid levels. Since multigrid is not currently
     * supported, this table is not filled with valid entries.
     */
    std::vector<dealii::internal::DoFHandlerImplementation::NumberCache>
      mg_number_cache;

    /**
     * Array to store the indices for degrees of freedom located at vertices.
     *
     * The format used here, in the form of a linked list, is the same as used
     * for the arrays used in the internal::hp::DoFLevel hierarchy. Starting
     * indices into this array are provided by the vertex_dof_offsets field.
     *
     * Access to this field is generally through the
     * DoFAccessor::get_vertex_dof_index() and
     * DoFAccessor::set_vertex_dof_index() functions, encapsulating the actual
     * data format used to the present class.
     */
    std::vector<types::global_dof_index> vertex_dofs;

    /**
     * For each vertex in the triangulation, store the offset within the
     * vertex_dofs array where the dofs for this vertex start.
     *
     * As for that array, the format is the same as described in the
     * documentation of hp::DoFLevel.
     *
     * Access to this field is generally through the
     * Accessor::get_vertex_dof_index() and Accessor::set_vertex_dof_index()
     * functions, encapsulating the actual data format used to the present
     * class.
     */
    std::vector<unsigned int> vertex_dof_offsets;

    /**
     * Whenever the underlying triangulation changes by either
     * h/p refinement/coarsening and serialization, the active_fe_index of cells
     * needs to be transferred. This structure stores all temporary information
     * required during that process.
     */
    struct ActiveFEIndexTransfer
    {
      /**
       * Container to temporarily store the iterator and future active FE index
       * of cells that persist.
       */
      std::map<const cell_iterator, const unsigned int>
        persisting_cells_fe_index;

      /**
       * Container to temporarily store the iterator and future active FE index
       * of cells that will be refined.
       */
      std::map<const cell_iterator, const unsigned int> refined_cells_fe_index;

      /**
       * Container to temporarily store the iterator and future active FE index
       * of parent cells that will remain after coarsening.
       */
      std::map<const cell_iterator, const unsigned int>
        coarsened_cells_fe_index;

      /**
       * Container to temporarily store the active_fe_index of every locally
       * owned cell for transfer across parallel::distributed::Triangulation
       * objects.
       */
      std::vector<unsigned int> active_fe_indices;

      /**
       * Helper object to transfer all active_fe_indices on
       * parallel::distributed::Triangulation objects during
       * refinement/coarsening and serialization.
       */
      std::unique_ptr<
        parallel::distributed::
          CellDataTransfer<dim, spacedim, std::vector<unsigned int>>>
        cell_data_transfer;
    };

    /**
     * We embed our data structure into a pointer to control that
     * all transfer related data only exists during the actual transfer process.
     */
    std::unique_ptr<ActiveFEIndexTransfer> active_fe_index_transfer;

    /**
     * A list of connections with which this object connects to the
     * triangulation to get information about when the triangulation changes.
     */
    std::vector<boost::signals2::connection> tria_listeners;

    // Make accessor objects friends.
    template <int, class, bool>
    friend class dealii::DoFAccessor;
    template <class, bool>
    friend class dealii::DoFCellAccessor;
    friend struct dealii::internal::DoFAccessorImplementation::Implementation;
    friend struct dealii::internal::DoFCellAccessorImplementation::
      Implementation;

    // Likewise for DoFLevel objects since they need to access the vertex dofs
    // in the functions that set and retrieve vertex dof indices.
    template <int>
    friend class dealii::internal::hp::DoFIndicesOnFacesOrEdges;
    friend struct dealii::internal::hp::DoFHandlerImplementation::
      Implementation;
    friend struct dealii::internal::DoFHandlerImplementation::Policy::
      Implementation;
  };



#ifndef DOXYGEN


  /* ----------------------- Inline functions ----------------------------------
   */


  template <int dim, int spacedim>
  template <typename number>
  types::global_dof_index
  DoFHandler<dim, spacedim>::n_boundary_dofs(
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &boundary_ids) const
  {
    // extract the set of boundary ids and forget about the function object
    // pointers
    std::set<types::boundary_id> boundary_ids_only;
    for (typename std::map<types::boundary_id,
                           const Function<spacedim, number> *>::const_iterator
           p = boundary_ids.begin();
         p != boundary_ids.end();
         ++p)
      boundary_ids_only.insert(p->first);

    // then just hand everything over to the other function that does the work
    return n_boundary_dofs(boundary_ids_only);
  }



  template <>
  inline types::global_dof_index
  DoFHandler<2, 3>::n_boundary_dofs() const
  {
    Assert(false, ExcNotImplemented());
    return 0;
  }



  template <>
  template <typename number>
  inline types::global_dof_index
  DoFHandler<2, 3>::n_boundary_dofs(
    const std::map<types::boundary_id, const Function<3, number> *> &) const
  {
    Assert(false, ExcNotImplemented());
    return 0;
  }



  template <>
  inline types::global_dof_index
  DoFHandler<2, 3>::n_boundary_dofs(const std::set<types::boundary_id> &) const
  {
    Assert(false, ExcNotImplemented());
    return 0;
  }
}


namespace internal
{
  /**
   * Return a string representing the dynamic type of the given argument.
   * This is basically the same what typeid(...).name() does, but it turns out
   * this is broken on Intel 13+.
   *
   * Defined in source/dofs/dof_handler.cc.
   */
  template <int dim, int spacedim>
  std::string
  policy_to_string(const dealii::internal::DoFHandlerImplementation::Policy::
                     PolicyBase<dim, spacedim> &policy);
} // namespace internal


namespace hp
{
  template <int dim, int spacedim>
  template <class Archive>
  void
  DoFHandler<dim, spacedim>::save(Archive &ar, unsigned int) const
  {
    ar &vertex_dofs;
    ar &vertex_dof_offsets;
    ar &number_cache;
    ar &mg_number_cache;

    // some versions of gcc have trouble with loading vectors of
    // std::unique_ptr objects because std::unique_ptr does not
    // have a copy constructor. do it one level at a time
    const unsigned int n_levels = levels.size();
    ar &               n_levels;
    for (unsigned int i = 0; i < n_levels; ++i)
      ar &levels[i];

    // boost dereferences a nullptr when serializing a nullptr
    // at least up to 1.65.1. This causes problems with clang-5.
    // Therefore, work around it.
    bool faces_is_nullptr = (faces.get() == nullptr);
    ar & faces_is_nullptr;
    if (!faces_is_nullptr)
      ar &faces;

    // write out the number of triangulation cells and later check during
    // loading that this number is indeed correct; same with something that
    // identifies the policy
    const unsigned int n_cells = tria->n_cells();
    std::string policy_name    = dealii::internal::policy_to_string(*policy);

    ar &n_cells &policy_name;
  }



  template <int dim, int spacedim>
  template <class Archive>
  void
  DoFHandler<dim, spacedim>::load(Archive &ar, unsigned int)
  {
    ar &vertex_dofs;
    ar &vertex_dof_offsets;
    ar &number_cache;
    ar &mg_number_cache;

    // boost::serialization can restore pointers just fine, but if the
    // pointer object still points to something useful, that object is not
    // destroyed and we end up with a memory leak. consequently, first delete
    // previous content before re-loading stuff
    levels.clear();
    faces.reset();

    // some versions of gcc have trouble with loading vectors of
    // std::unique_ptr objects because std::unique_ptr does not
    // have a copy constructor. do it one level at a time
    unsigned int size;
    ar &         size;
    levels.resize(size);
    for (unsigned int i = 0; i < size; ++i)
      {
        std::unique_ptr<dealii::internal::hp::DoFLevel> level;
        ar &                                            level;
        levels[i] = std::move(level);
      }

    // Workaround for nullptr, see in save().
    bool faces_is_nullptr = true;
    ar & faces_is_nullptr;
    if (!faces_is_nullptr)
      ar &faces;

    // these are the checks that correspond to the last block in the save()
    // function
    unsigned int n_cells;
    std::string  policy_name;

    ar &n_cells &policy_name;

    AssertThrow(
      n_cells == tria->n_cells(),
      ExcMessage(
        "The object being loaded into does not match the triangulation "
        "that has been stored previously."));
    AssertThrow(policy_name == dealii::internal::policy_to_string(*policy),
                ExcMessage(
                  "The policy currently associated with this DoFHandler (" +
                  dealii::internal::policy_to_string(*policy) +
                  ") does not match the one that was associated with the "
                  "DoFHandler previously stored (" +
                  policy_name + ")."));
  }

  template <int dim, int spacedim>
  inline types::global_dof_index
  DoFHandler<dim, spacedim>::n_dofs() const
  {
    return number_cache.n_global_dofs;
  }



  template <int dim, int spacedim>
  inline types::global_dof_index
  DoFHandler<dim, spacedim>::n_dofs(const unsigned int) const
  {
    Assert(false, ExcNotImplemented());
    return numbers::invalid_dof_index;
  }



  template <int dim, int spacedim>
  types::global_dof_index
  DoFHandler<dim, spacedim>::n_locally_owned_dofs() const
  {
    return number_cache.n_locally_owned_dofs;
  }



  template <int dim, int spacedim>
  const IndexSet &
  DoFHandler<dim, spacedim>::locally_owned_dofs() const
  {
    return number_cache.locally_owned_dofs;
  }



  template <int dim, int spacedim>
  const std::vector<types::global_dof_index> &
  DoFHandler<dim, spacedim>::n_locally_owned_dofs_per_processor() const
  {
    if (number_cache.n_locally_owned_dofs_per_processor.empty() &&
        number_cache.n_global_dofs > 0)
      {
        const_cast<dealii::internal::DoFHandlerImplementation::NumberCache &>(
          number_cache)
          .n_locally_owned_dofs_per_processor =
          compute_n_locally_owned_dofs_per_processor();
      }
    return number_cache.n_locally_owned_dofs_per_processor;
  }



  template <int dim, int spacedim>
  const std::vector<IndexSet> &
  DoFHandler<dim, spacedim>::locally_owned_dofs_per_processor() const
  {
    if (number_cache.locally_owned_dofs_per_processor.empty() &&
        number_cache.n_global_dofs > 0)
      {
        const_cast<dealii::internal::DoFHandlerImplementation::NumberCache &>(
          number_cache)
          .locally_owned_dofs_per_processor =
          compute_locally_owned_dofs_per_processor();
      }
    return number_cache.locally_owned_dofs_per_processor;
  }



  template <int dim, int spacedim>
  std::vector<types::global_dof_index>
  DoFHandler<dim, spacedim>::compute_n_locally_owned_dofs_per_processor() const
  {
    const parallel::TriangulationBase<dim, spacedim> *tr =
      (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &this->get_triangulation()));
    if (tr != nullptr)
      return number_cache.get_n_locally_owned_dofs_per_processor(
        tr->get_communicator());
    else
      return number_cache.get_n_locally_owned_dofs_per_processor(MPI_COMM_SELF);
  }



  template <int dim, int spacedim>
  std::vector<IndexSet>
  DoFHandler<dim, spacedim>::compute_locally_owned_dofs_per_processor() const
  {
    const parallel::TriangulationBase<dim, spacedim> *tr =
      (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &this->get_triangulation()));
    if (tr != nullptr)
      return number_cache.get_locally_owned_dofs_per_processor(
        tr->get_communicator());
    else
      return number_cache.get_locally_owned_dofs_per_processor(MPI_COMM_SELF);
  }



  template <int dim, int spacedim>
  const IndexSet &
  DoFHandler<dim, spacedim>::locally_owned_mg_dofs(
    const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
    (void)level;
    Assert(level < this->get_triangulation().n_global_levels(),
           ExcMessage("The given level index exceeds the number of levels "
                      "present in the triangulation"));
    return mg_number_cache[0].locally_owned_dofs;
  }



  template <int dim, int spacedim>
  const std::vector<IndexSet> &
  DoFHandler<dim, spacedim>::locally_owned_mg_dofs_per_processor(
    const unsigned int level) const
  {
    Assert(level < this->get_triangulation().n_global_levels(),
           ExcMessage("The given level index exceeds the number of levels "
                      "present in the triangulation"));
    Assert(
      mg_number_cache.size() == this->get_triangulation().n_global_levels(),
      ExcMessage(
        "The level dofs are not set up properly! Did you call distribute_mg_dofs()?"));
    if (mg_number_cache[level].locally_owned_dofs_per_processor.empty() &&
        mg_number_cache[level].n_global_dofs > 0)
      {
        const_cast<dealii::internal::DoFHandlerImplementation::NumberCache &>(
          mg_number_cache[level])
          .locally_owned_dofs_per_processor =
          compute_locally_owned_mg_dofs_per_processor(level);
      }
    return mg_number_cache[level].locally_owned_dofs_per_processor;
  }



  template <int dim, int spacedim>
  std::vector<IndexSet>
  DoFHandler<dim, spacedim>::compute_locally_owned_mg_dofs_per_processor(
    const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
    (void)level;
    Assert(level < this->get_triangulation().n_global_levels(),
           ExcMessage("The given level index exceeds the number of levels "
                      "present in the triangulation"));
    const parallel::TriangulationBase<dim, spacedim> *tr =
      (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &this->get_triangulation()));
    if (tr != nullptr)
      return mg_number_cache[level].get_locally_owned_dofs_per_processor(
        tr->get_communicator());
    else
      return mg_number_cache[level].get_locally_owned_dofs_per_processor(
        MPI_COMM_SELF);
  }



  template <int dim, int spacedim>
  inline const FiniteElement<dim, spacedim> &
  DoFHandler<dim, spacedim>::get_fe(const unsigned int number) const
  {
    Assert(fe_collection.size() > 0,
           ExcMessage("No finite element collection is associated with "
                      "this DoFHandler"));
    return fe_collection[number];
  }



  template <int dim, int spacedim>
  inline const hp::FECollection<dim, spacedim> &
  DoFHandler<dim, spacedim>::get_fe_collection() const
  {
    Assert(fe_collection.size() > 0,
           ExcMessage("No finite element collection is associated with "
                      "this DoFHandler"));
    return fe_collection;
  }



  template <int dim, int spacedim>
  inline const Triangulation<dim, spacedim> &
  DoFHandler<dim, spacedim>::get_triangulation() const
  {
    Assert(tria != nullptr,
           ExcMessage("This DoFHandler object has not been associated "
                      "with a triangulation."));
    return *tria;
  }

#endif

} // namespace hp

DEAL_II_NAMESPACE_CLOSE

#endif
