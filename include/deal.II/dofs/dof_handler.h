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

#ifndef dealii_dof_handler_h
#define dealii_dof_handler_h



#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/base/observer_pointer.h>
#include <deal.II/base/types.h>

#include <deal.II/dofs/block_info.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_faces.h>
#include <deal.II/dofs/dof_iterator_selector.h>
#include <deal.II/dofs/dof_levels.h>
#include <deal.II/dofs/number_cache.h>

#include <deal.II/hp/fe_collection.h>

#include <boost/serialization/split_member.hpp>
#include <boost/signals2/connection.hpp>

#include <map>
#include <memory>
#include <set>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class FiniteElement;
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class Triangulation;

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

  namespace DoFAccessorImplementation
  {
    struct Implementation;
  }

  namespace DoFCellAccessorImplementation
  {
    struct Implementation;
  }

  namespace hp
  {
    namespace DoFHandlerImplementation
    {
      struct Implementation;
    }
  } // namespace hp
} // namespace internal

namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim, typename VectorType>
    class CellDataTransfer;
  }
} // namespace parallel

template <int structdim, int dim, int spacedim, bool level_dof_access>
class DoFAccessor;

template <int dimension_, int space_dimension_, bool level_dof_access>
class DoFCellAccessor;

#endif

/**
 * Given a triangulation and a description of a finite element, this
 * class enumerates degrees of freedom on all vertices, edges, faces,
 * and cells of the triangulation. As a result, it also provides a
 * <i>basis</i> for a discrete space $V_h$ whose elements are finite
 * element functions defined on each cell by a FiniteElement object.
 * This class satisfies the
 * @ref ConceptMeshType "MeshType concept"
 * requirements.
 *
 * It is first used in the step-2 tutorial program.
 *
 * For each 0d, 1d, 2d, and 3d subobject, this class stores a list of the
 * indices of degrees of freedom defined on this DoFHandler. These indices
 * refer to the unconstrained degrees of freedom, i.e. constrained degrees of
 * freedom are numbered in the same way as unconstrained ones, and are only
 * later eliminated.  This leads to the fact that indices in global vectors
 * and matrices also refer to all degrees of freedom and some kind of
 * condensation is needed to restrict the systems of equations to the
 * unconstrained degrees of freedom only. The actual layout of storage of the
 * indices is described in the
 * dealii::internal::DoFHandlerImplementation::DoFLevel class documentation.
 *
 * The class offers iterators to traverse all cells, in much the same way as
 * the Triangulation class does. Using the begin() and end() functions (and
 * companions, like begin_active()), one can obtain iterators to walk over
 * cells, and query the degree of freedom structures as well as the
 * triangulation data. These iterators are built on top of those of the
 * Triangulation class, but offer the additional information on degrees of
 * freedom functionality compared to pure triangulation iterators. The order
 * in which dof iterators are presented by the <tt>++</tt> and <tt>\--</tt>
 * operators is the same as that for the corresponding iterators traversing
 * the triangulation on which this DoFHandler is constructed.
 *
 * The <tt>spacedim</tt> parameter has to be used if one wants to solve
 * problems on surfaces. If not specified, this parameter takes the default
 * value <tt>=dim</tt> implying that we want to solve problems in a domain
 * whose dimension equals the dimension of the space in which it is embedded.
 *
 *
 * <h3>Distribution of indices for degrees of freedom</h3>
 *
 * The degrees of freedom (`dofs') are distributed on the given triangulation
 * by the function distribute_dofs(). It gets passed a finite element object
 * describing how many degrees of freedom are located on vertices, lines, etc.
 * It traverses the triangulation cell by cell and numbers the dofs of that
 * cell if not yet numbered. For non-multigrid algorithms, only active cells
 * are considered. Active cells are defined to be those cells which have no
 * children, i.e. they are the most refined ones.
 *
 * Since the triangulation is traversed starting with the cells of the
 * coarsest active level and going to more refined levels, the lowest numbers
 * for dofs are given to the largest cells as well as their bounding lines and
 * vertices, with the dofs of more refined cells getting higher numbers.
 *
 * This numbering implies very large bandwidths of the resulting matrices and
 * is thus vastly suboptimal for some solution algorithms. For this reason,
 * the DoFRenumbering class offers several algorithms to reorder the dof
 * numbering according. See there for a discussion of the implemented
 * algorithms.
 *
 *
 * <h3>Interaction with distributed meshes</h3>
 *
 * Upon construction, this class takes a reference to a triangulation object.
 * In most cases, this will be a reference to an object of type Triangulation,
 * i.e. the class that represents triangulations that entirely reside on a
 * single processor. However, it can also be of type
 * parallel::distributed::Triangulation (see, for example, step-32, step-40
 * and in particular the
 * @ref distributed
 * topic) in which case the DoFHandler object will proceed to only manage
 * degrees of freedom on locally owned and ghost cells. This process is
 * entirely transparent to the used.
 *
 *
 * <h3>User defined renumbering schemes</h3>
 *
 * The DoFRenumbering class offers a number of renumbering schemes like the
 * Cuthill-McKee scheme. Basically, the function sets up an array in which for
 * each degree of freedom we store the new index this DoF should have after
 * renumbering. Using this array, the renumber_dofs() function of the present
 * class is called, which actually performs the change from old DoF indices to
 * the ones given in the array. In some cases, however, a user may want to
 * compute their own renumbering order; in this case, one can allocate an array
 * with one element per degree of freedom and fill it with the number that the
 * respective degree of freedom shall be assigned. This number may, for
 * example, be obtained by sorting the support points of the degrees of
 * freedom in downwind direction.  Then call the
 * <tt>renumber_dofs(vector<types::global_dof_index>)</tt> function with the
 * array, which converts old into new degree of freedom indices.
 *
 *
 * <h3>Serializing (loading or storing) DoFHandler objects</h3>
 *
 * Like many other classes in deal.II, the DoFHandler class can stream its
 * contents to an archive using BOOST's serialization facilities. The data so
 * stored can later be retrieved again from the archive to restore the
 * contents of this object. This facility is frequently used to save the state
 * of a program to disk for possible later resurrection, often in the context
 * of checkpoint/restart strategies for long running computations or on
 * computers that aren't very reliable (e.g. on very large clusters where
 * individual nodes occasionally fail and then bring down an entire MPI job).
 *
 * The model for doing so is similar for the DoFHandler class as it is for the
 * Triangulation class (see the section in the general documentation of that
 * class). In particular, the load() function does not exactly restore the
 * same state as was stored previously using the save() function. Rather, the
 * function assumes that you load data into a DoFHandler object that is
 * already associated with a triangulation that has a content that matches the
 * one that was used when the data was saved. Likewise, the load() function
 * assumes that the current object is already associated with a finite element
 * object that matches the one that was associated with it when data was
 * saved; the latter can be achieved by calling DoFHandler::distribute_dofs()
 * using the same kind of finite element before re-loading data from the
 * serialization archive.
 *
 *
 * <h3>hp-adaptive finite element methods</h3>
 *
 * Instead of only using one particular FiniteElement on all cells, this class
 * also allows for an enumeration of degrees of freedom on different finite
 * elements on every cells. To this end, one assigns an
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
 * @ref hp_paper "hp-paper".
 *
 *
 * <h3>Active FE indices and their behavior under mesh refinement</h3>
 *
 * The typical workflow for using this class is to create a mesh, assign an
 * active FE index to every active cell, call DoFHandler::distribute_dofs(),
 * and then assemble a linear system and solve a problem on this finite element
 * space.
 *
 * Active FE indices will be automatically transferred during mesh adaptation
 * from the old to the new mesh. Future FE indices are meant to determine the
 * active FE index after mesh adaptation, and are used to prepare data on the
 * old mesh for the new one. If no future FE index is specified, the finite
 * element prevails.
 *
 * In particular, the following rules apply during adaptation:
 * - Upon mesh refinement, child cells inherit the future FE index of
 *   the parent.
 * - When coarsening cells, the (now active) parent cell will be assigned
 *   a future FE index that is determined from its (no longer active)
 *   children, following the FiniteElementDomination logic: Out of the set of
 *   elements previously assigned to the former children, we choose the one
 *   dominated by all children for the parent cell. If none was found, we pick
 *   the most dominant element in the whole collection that is dominated by
 *   all former children. See hp::FECollection::find_dominated_fe_extended()
 *   for further information on this topic.
 *
 * Strategies for automatic hp-adaptation which will set future FE indices based
 * on criteria are available in the hp::Refinement namespace.
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
 * is active on them, however: whenever you call DoFHandler::distribute_dofs(),
 * all processors that participate in the parallel mesh exchange information in
 * such a way that the active FE index on ghost cells equals the active FE index
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
 * However, using a parallel::distributed::Triangulation with a DoFHandler
 * in hp-mode requires additional attention during serialization, since no
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
 * @ingroup dofs
 *
 * @dealiiConceptRequires{(concepts::is_valid_dim_spacedim<dim, spacedim>)}
 */
template <int dim, int spacedim = dim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler : public EnableObserverPointer
{
  using ActiveSelector =
    dealii::internal::DoFHandlerImplementation::Iterators<dim, spacedim, false>;
  using LevelSelector =
    dealii::internal::DoFHandlerImplementation::Iterators<dim, spacedim, true>;

public:
  /**
   * An alias that is used to identify cell iterators in DoFHandler objects.
   * The concept of iterators is discussed at length in the
   * @ref Iterators "iterators documentation topic".
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
   * @ref Iterators "iterators documentation topic".
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
   * the active lines) are part of at least one active cell. Each such line may
   * additionally be a child of a line of a coarser cell adjacent to a cell
   * that is active. (This coarser neighbor would then also be active.)
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
   * the active quads) are faces of at least one active cell. Each such quad may
   * additionally be a child of a quad face of a coarser cell adjacent to a cell
   * that is active. (This coarser neighbor would then also be active.)
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
   * An alias that is used to identify
   * @ref GlossActive "active cell iterators".
   * The concept of iterators is discussed at length in the
   * @ref Iterators "iterators documentation topic".
   *
   * The current alias identifies active cells in a DoFHandler object. While
   * the actual data type of the alias is hidden behind a few layers of
   * (unfortunately necessary) indirections, it is in essence
   * TriaActiveIterator<DoFCellAccessor>. The TriaActiveIterator class works
   * like a pointer to active objects that when you dereference it yields an
   * object of type DoFCellAccessor. DoFCellAccessor is a class that
   * identifies properties that are specific to cells in a DoFHandler, but it
   * is derived (and consequently inherits) from both DoFAccessor,
   * TriaCellAccessor and TriaAccessor that describe what you can ask of more
   * general objects (lines, faces, as well as cells) in a triangulation and
   * DoFHandler objects.
   *
   * @ingroup Iterators
   */
  using active_cell_iterator = typename ActiveSelector::active_cell_iterator;

  /**
   * An alias that is used to identify cell iterators. The concept of
   * iterators is discussed at length in the
   * @ref Iterators "iterators documentation topic".
   *
   * The current alias identifies cells in a DoFHandler object. Some of
   * these cells may in fact be active (see
   * @ref GlossActive "active cell iterators")
   * in which case they can in fact be asked for the degrees of freedom that
   * live on them. On the other hand, if the cell is not active, any such
   * query will result in an error. Note that this is what distinguishes this
   * alias from the level_cell_iterator alias.
   *
   * While the actual data type of the alias is hidden behind a few layers
   * of (unfortunately necessary) indirections, it is in essence
   * TriaIterator<DoFCellAccessor>. The TriaIterator class works like a
   * pointer to objects that when you dereference it yields an object of type
   * DoFCellAccessor. DoFCellAccessor is a class that identifies properties
   * that are specific to cells in a DoFHandler, but it is derived (and
   * consequently inherits) from both DoFAccessor, TriaCellAccessor and
   * TriaAccessor that describe what you can ask of more general objects
   * (lines, faces, as well as cells) in a triangulation and DoFHandler
   * objects.
   *
   * @ingroup Iterators
   */
  using cell_iterator = typename ActiveSelector::cell_iterator;

  /**
   * An alias that is used to identify iterators that point to faces.
   * The concept of iterators is discussed at length in the
   * @ref Iterators "iterators documentation topic".
   *
   * While the actual data type of the alias is hidden behind a few layers
   * of (unfortunately necessary) indirections, it is in essence
   * TriaIterator<DoFAccessor>. The
   * TriaIterator class works like a pointer to objects that when
   * you dereference it yields an object of type DoFAccessor. DoFAccessor,
   * in turn, is a class that can be used to query DoF indices on faces,
   * but it is also derived from TriaAccessor and consequently can be used
   * to query geometric properties such as vertices of faces, their area, etc.
   *
   * @ingroup Iterators
   */
  using face_iterator = typename ActiveSelector::face_iterator;

  /**
   * An alias that is used to identify iterators that point to active faces,
   * i.e., to faces that have no children. Active faces must be faces of at
   * least one active cell.
   *
   * Other than the "active" qualification, this alias is identical to the
   * @p face_iterator alias. In particular, dereferencing either yields
   * the same kind of object.
   *
   * @ingroup Iterators
   */
  using active_face_iterator = typename ActiveSelector::active_face_iterator;

  using level_cell_accessor = typename LevelSelector::CellAccessor;
  using level_face_accessor = typename LevelSelector::FaceAccessor;

  using level_cell_iterator = typename LevelSelector::cell_iterator;
  using level_face_iterator = typename LevelSelector::face_iterator;


  /**
   * Make the dimension available in function templates.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Make the space dimension available in function templates.
   */
  static constexpr unsigned int space_dimension = spacedim;

  /**
   * The default index of the finite element to be used on a given cell.
   */
  static const types::fe_index default_fe_index = 0;

  /**
   * The type in which we store the offsets in the CRS data structures.
   */
  using offset_type = unsigned int;

  /**
   * Standard constructor, not initializing any data. After constructing an
   * object with this constructor, use reinit() to get a valid DoFHandler.
   */
  DoFHandler();

  /**
   * Constructor. Take @p tria as the triangulation to work on.
   */
  explicit DoFHandler(const Triangulation<dim, spacedim> &tria);

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
   * For each locally owned cell, set the active finite element index to the
   * corresponding value given in @p active_fe_indices.
   *
   * The vector @p active_fe_indices needs to have as many entries as there
   * are active cells. The FE indices must be in the order in which we iterate
   * over active cells. Vector entries corresponding to active cells that are
   * not locally owned are ignored.
   *
   * Active FE indices will only be set for locally owned cells. Ghost and
   * artificial cells will be ignored; no active FE index will be assigned to
   * them. To exchange active FE indices on ghost cells, call distribute_dofs()
   * afterwards.
   */
  void
  set_active_fe_indices(const std::vector<types::fe_index> &active_fe_indices);

  /**
   * @copydoc set_active_fe_indices()
   *
   * @deprecated Use set_active_fe_indices() with the types::fe_index datatype.
   */
  DEAL_II_DEPRECATED
  void
  set_active_fe_indices(const std::vector<unsigned int> &active_fe_indices);

  /**
   * For each locally relevant cell, extract the active finite element index and
   * return them in the order in which we iterate over active cells.
   *
   * As we do not know the active FE index on artificial cells, they are set to
   * the invalid value numbers::invalid_fe_index.
   *
   * For DoFHandler objects without hp-capabilities, the vector will consist of
   * zeros, indicating that all cells use the same finite element. In hp-mode,
   * the values may be different, though.
   *
   * The returned vector has as many entries as there are active cells.
   */
  std::vector<types::fe_index>
  get_active_fe_indices() const;

  /**
   * For each locally relevant cell, extract the active finite element index and
   * fill the vector @p active_fe_indices in the order in which we iterate over
   * active cells. This vector is resized, if necessary.
   *
   * As we do not know the active FE index on artificial cells, they are set to
   * the invalid value numbers::invalid_fe_index.
   *
   * For DoFHandler objects without hp-capabilities, the vector will consist of
   * zeros, indicating that all cells use the same finite element. In hp-mode,
   * the values may be different, though.
   *
   * The returned vector has as many entries as there are active cells.
   *
   * @deprecated Use get_active_fe_indices() that returns the result vector.
   */
  DEAL_II_DEPRECATED
  void
  get_active_fe_indices(std::vector<unsigned int> &active_fe_indices) const;

  /**
   * For each locally owned cell, set the future finite element index to the
   * corresponding value given in @p future_fe_indices.
   *
   * The vector @p future_fe_indices needs to have as many entries as there
   * are active cells. The FE indices must be in the order in which we iterate
   * over active cells. Vector entries corresponding to active cells that are
   * not locally owned are ignored.
   *
   * Future FE indices will only be set for locally owned cells. Ghost and
   * artificial cells will be ignored; no future FE index will be assigned to
   * them.
   */
  void
  set_future_fe_indices(const std::vector<types::fe_index> &future_fe_indices);

  /**
   * For each locally owned cell, extract the future finite element index and
   * return them in the order in which we iterate over active cells.
   *
   * As we do not know the future FE index on ghost and artificial cells, they
   * are set to the invalid value numbers::invalid_fe_index. The same applies to
   * locally owned cells that have no future FE index assigned.
   *
   * The returned vector has as many entries as there are active cells.
   */
  std::vector<types::fe_index>
  get_future_fe_indices() const;

  /**
   * Assign a Triangulation to the DoFHandler.
   *
   * Remove all associations with the previous Triangulation object and
   * establish connections with the new one. All information about previous
   * degrees of freedom will be removed. Activates hp-mode.
   */
  void
  reinit(const Triangulation<dim, spacedim> &tria);

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
   * This function is first discussed in the introduction to the
   * step-2 tutorial program.
   *
   * @note This function makes a copy of the finite element given as
   * argument, and stores it as a member variable, similarly to the above
   * function set_fe().
   */
  void
  distribute_dofs(const FiniteElement<dim, spacedim> &fe);

  /**
   * Same as above but taking an hp::FECollection object.
   */
  void
  distribute_dofs(const hp::FECollection<dim, spacedim> &fe);

  /**
   * Distribute level degrees of freedom on each level for geometric
   * multigrid. The active DoFs need to be distributed using distribute_dofs()
   * before calling this function.
   */
  void
  distribute_mg_dofs();

  /**
   * Returns whether this DoFHandler has hp-capabilities.
   */
  bool
  has_hp_capabilities() const;

  /**
   * This function returns whether this DoFHandler has DoFs distributed on
   * each multigrid level or in other words if distribute_mg_dofs() has been
   * called.
   */
  bool
  has_level_dofs() const;

  /**
   * This function returns whether this DoFHandler has active DoFs. This is
   * equivalent to asking whether (i) distribute_dofs() has been called and
   * (ii) the finite element for which degrees of freedom have been
   * distributed actually has degrees of freedom (which is not the case for
   * FE_Nothing, for example).
   *
   * If this object is based on a parallel::distributed::Triangulation, then
   * the current function returns true if <i>any</i> partition of the parallel
   * DoFHandler object has any degrees of freedom. In other words, the
   * function returns true even if the Triangulation does not own any active
   * cells on the current MPI process, but at least one process owns cells and
   * at least this one process has any degrees of freedom associated with it.
   */
  bool
  has_active_dofs() const;

  /**
   * After distribute_dofs() with an FESystem element, the block structure of
   * global and level vectors is stored in a BlockInfo object accessible with
   * block_info(). This function initializes the local block structure on each
   * cell in the same object.
   */
  void
  initialize_local_block_info();

  /**
   * Clear all data of this object.
   */
  void
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
   * The same function as above, but renumber the degrees of freedom of a
   * single level of a multigrid hierarchy.
   */
  void
  renumber_dofs(const unsigned int                          level,
                const std::vector<types::global_dof_index> &new_numbers);

  /**
   * Return the maximum number of degrees of freedom a degree of freedom in
   * the given triangulation with the given finite element may couple with.
   * This is the maximum number of entries per line in the system matrix; this
   * information can therefore be used upon construction of the
   * SparsityPattern object.
   *
   * The returned number is not really the maximum number but an estimate
   * based on the finite element and the maximum number of cells meeting at a
   * vertex. The number holds for the constrained matrix as well.
   *
   * The determination of the number of couplings can be done by simple
   * picture drawing. An example can be found in the implementation of this
   * function.
   *
   * @note This function is most often used to determine the maximal row
   * length for sparsity patterns. Unfortunately, while the estimates returned
   * by this function are rather accurate in 1d and 2d, they are often
   * significantly too high in 3d, leading the SparsityPattern class to
   * allocate much too much memory in some cases. Unless someone comes around
   * to improving the present function for 3d, there is not very much one can
   * do about these cases. The typical way to work around this problem is to
   * use an intermediate compressed sparsity pattern that only allocates
   * memory on demand. Refer to the step-2 and step-11 example programs on how
   * to do this. The problem is also discussed in the documentation of the
   * topic on
   * @ref Sparsity.
   */
  unsigned int
  max_couplings_between_dofs() const;

  /**
   * Return the number of degrees of freedom located on the boundary another
   * dof on the boundary can couple with.
   *
   * The number is the same as for max_couplings_between_dofs() in one
   * dimension less.
   *
   * @note The same applies to this function as to max_couplings_between_dofs() as
   * regards the performance of this function. Think about one of the dynamic
   * sparsity pattern classes instead (see
   * @ref Sparsity).
   */
  unsigned int
  max_couplings_between_boundary_dofs() const;

  /*--------------------------------------*/

  /**
   * @name Cell iterator functions
   */

  /** @{ */

  /**
   * Iterator to the first used cell on level @p level.
   */
  cell_iterator
  begin(const unsigned int level = 0) const;

  /**
   * Iterator to the first active cell on level @p level. If the given level
   * does not contain any active cells (i.e., all cells on this level are
   * further refined), then this function returns
   * <code>end_active(level)</code> so that loops of the kind
   * @code
   *   for (cell=dof_handler.begin_active(level);
   *        cell!=dof_handler.end_active(level);
   *        ++cell)
   *     {
   *       ...
   *     }
   * @endcode
   * have zero iterations, as may be expected if there are no active cells on
   * this level.
   */
  active_cell_iterator
  begin_active(const unsigned int level = 0) const;

  /**
   * Iterator past the end; this iterator serves for comparisons of iterators
   * with past-the-end or before-the-beginning states.
   */
  cell_iterator
  end() const;

  /**
   * Return an iterator which is the first iterator not on the given level. If
   * @p level is the last level, then this returns <tt>end()</tt>.
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
   * Iterator to the first used cell on level @p level. This returns a
   * level_cell_iterator that returns level dofs when dof_indices() is called.
   */
  level_cell_iterator
  begin_mg(const unsigned int level = 0) const;

  /**
   * Iterator past the last cell on level @p level. This returns a
   * level_cell_iterator that returns level dofs when dof_indices() is called.
   */
  level_cell_iterator
  end_mg(const unsigned int level) const;

  /**
   * Iterator past the end; this iterator serves for comparisons of iterators
   * with past-the-end or before-the-beginning states.
   */
  level_cell_iterator
  end_mg() const;
  /** @} */

  /**
   * @name Cell iterator functions returning ranges of iterators
   */

  /** @{ */
  /**
   * Return an iterator range that contains all cells (active or not) that
   * make up this DoFHandler. Such a range is useful to initialize range-based
   * for loops as supported by C++11. See the example in the documentation of
   * active_cell_iterators().
   *
   * @return The half open range <code>[this->begin(), this->end())</code>
   *
   * @ingroup CPP11
   */
  IteratorRange<cell_iterator>
  cell_iterators() const;

  /**
   * Return an iterator range that contains all active cells that make up this
   * DoFHandler. Such a range is useful to initialize range-based for loops as
   * supported by C++11, see also
   * @ref CPP11 "C++11 standard".
   *
   * Range-based for loops are useful in that they require much less code than
   * traditional loops (see <a href="http://en.wikipedia.org/wiki/C%2B%2B11
   * #Range-based_for_loop">here</a> for a discussion of how they work). An
   * example is that without range-based for loops, one often writes code such
   * as the following:
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
   * make up this DoFHandler in their level-cell form. Such a range is useful
   * to initialize range-based for loops as supported by C++11. See the
   * example in the documentation of active_cell_iterators().
   *
   * @return The half open range <code>[this->begin_mg(),
   * this->end_mg())</code>
   *
   * @ingroup CPP11
   */
  IteratorRange<level_cell_iterator>
  mg_cell_iterators() const;

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
   * Return an iterator range that contains all active cells that make up the
   * given level of this DoFHandler. Such a range is useful to initialize
   * range-based for loops as supported by C++11. See the example in the
   * documentation of active_cell_iterators().
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

  /**
   * Return an iterator range that contains all cells (active or not) that
   * make up the given level of this DoFHandler in their level-cell form. Such
   * a range is useful to initialize range-based for loops as supported by
   * C++11. See the example in the documentation of active_cell_iterators().
   *
   * @param[in] level A given level in the refinement hierarchy of this
   * triangulation.
   * @return The half open range <code>[this->begin_mg(level),
   * this->end_mg(level))</code>
   *
   * @pre level must be less than this->n_levels().
   *
   * @ingroup CPP11
   */
  IteratorRange<level_cell_iterator>
  mg_cell_iterators_on_level(const unsigned int level) const;

  /** @} */


  /*---------------------------------------*/


  /**
   * Return the global number of degrees of freedom. If the current object
   * handles all degrees of freedom itself (even if you may intend to solve
   * your linear system in parallel, such as in step-17 or step-18), then this
   * number equals the number of locally owned degrees of freedom since this
   * object doesn't know anything about what you want to do with it and
   * believes that it owns every degree of freedom it knows about.
   *
   * On the other hand, if this object operates on a
   * parallel::distributed::Triangulation object, then this function returns
   * the global number of degrees of freedom, accumulated over all processors.
   *
   * In either case, included in the returned number are those DoFs which are
   * constrained by hanging nodes, see
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
   * The (global) number of multilevel degrees of freedom on a given level.
   *
   * If no level degrees of freedom have been assigned to this level, returns
   * numbers::invalid_dof_index. Else returns the number of degrees of freedom
   * on this level.
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
   * Return the number of locally owned degrees of freedom located on those
   * parts of the boundary which have a boundary indicator listed in the given
   * set.
   * The reason that a @p map rather than a @p set is used is the same as
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
   * Return the number of degrees of freedom located on those parts of the
   * boundary which have a boundary indicator listed in the given set. The
   */
  types::global_dof_index
  n_boundary_dofs(const std::set<types::boundary_id> &boundary_ids) const;

  /**
   * Access to an object informing of the block structure of the dof handler.
   *
   * If an FESystem is used in distribute_dofs(), degrees of freedom naturally
   * split into several
   * @ref GlossBlock "blocks".
   * For each base element as many blocks appear as its multiplicity.
   *
   * At the end of distribute_dofs(), the number of degrees of freedom in each
   * block is counted, and stored in a BlockInfo object, which can be accessed
   * here. If you have previously called distribute_mg_dofs(), the same is
   * done on each level of the multigrid hierarchy. Additionally, the block
   * structure on each cell can be generated in this object by calling
   * initialize_local_block_info().
   */
  const BlockInfo &
  block_info() const;

  /**
   * Return the number of degrees of freedom that belong to this process.
   *
   * If this is a sequential DoFHandler, then the result equals that produced by
   * n_dofs(). (Here, "sequential" means that either
   * the whole program does not use MPI, or that it uses MPI
   * but only uses a single MPI process, or that there are multiple MPI
   * processes but the Triangulation on which this DoFHandler builds
   * works only on one MPI process.)
   * On the other hand, if we are operating on a
   * parallel::distributed::Triangulation or parallel::shared::Triangulation,
   * then it includes only the degrees
   * of freedom that the current processor owns. Note that in this case this
   * does not include all degrees of freedom that have been distributed on the
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
   * Return an IndexSet describing the set of locally owned DoFs used for the
   * given multigrid level as a subset of 0..n_dofs(level).
   */
  const IndexSet &
  locally_owned_mg_dofs(const unsigned int level) const;

  /**
   * Return a constant reference to the indexth finite element object that is
   * used by this object.
   */
  const FiniteElement<dim, spacedim> &
  get_fe(const types::fe_index index = 0) const;

  /**
   * Return a constant reference to the set of finite element objects that
   * are used by this object.
   */
  const hp::FECollection<dim, spacedim> &
  get_fe_collection() const;

  /**
   * Return a constant reference to the triangulation underlying this object.
   */
  const Triangulation<dim, spacedim> &
  get_triangulation() const;

  /**
   * Return MPI communicator used by the underlying triangulation.
   */
  MPI_Comm
  get_mpi_communicator() const;

  /**
   * Return MPI communicator used by the underlying triangulation.
   *
   * @deprecated Use get_mpi_communicator() instead.
   */
  DEAL_II_DEPRECATED_EARLY_WITH_COMMENT(
    "Access the MPI communicator with get_mpi_communicator() instead.")
  MPI_Comm
  get_communicator() const;

  /**
   * Whenever serialization with a parallel::distributed::Triangulation as the
   * underlying triangulation is considered, we also need to consider storing
   * the active FE indices on all active cells as well.
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
   * the active FE indices on all active cells as well.
   *
   * This function deserializes and distributes the previously stored
   * active FE indices on all active cells.
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
   * Write the data of this object to a stream for the purpose of
   * serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int version) const;

  /**
   * Read the data of this object from a stream for the purpose of
   * serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
  /**
   * Write and read the data of this object from a stream for the purpose
   * of serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int version);
#else
  // This macro defines the serialize() method that is compatible with
  // the templated save() and load() method that have been implemented.
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  /**
   * Exception
   */
  DeclException0(ExcNoFESelected);
  /**
   * Exception
   * @ingroup Exceptions
   */
  DeclException0(ExcInvalidBoundaryIndicator);
  /**
   * Exception
   * @ingroup Exceptions
   */
  DeclException1(ExcInvalidLevel,
                 int,
                 << "The given level " << arg1
                 << " is not in the valid range!");
  /**
   * Exception
   * @ingroup Exceptions
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
                 << "The mesh contains a cell with an active FE index of "
                 << arg1 << ", but the finite element collection only has "
                 << arg2 << " elements");

  /**
   * Exception used when a certain feature doesn't make sense when
   * DoFHandler does not hp-capabilities.
   */
  DeclExceptionMsg(ExcOnlyAvailableWithHP,
                   "The current function doesn't make sense when used with a "
                   "DoFHandler without hp-capabilities.");

  /**
   * Exception used when a certain feature is not implemented when the
   * DoFHandler has hp-capabilities.
   */
  DeclExceptionMsg(ExcNotImplementedWithHP,
                   "The current function has not yet been implemented for a "
                   "DoFHandler with hp-capabilities.");

private:
  /**
   * A data structure that is used to store the DoF indices associated with a
   * particular vertex. Unlike cells, vertices live on several levels of a
   * multigrid hierarchy; consequently, we need to store DoF indices for each
   * vertex for each of the levels it lives on. This class does this.
   */
  class MGVertexDoFs
  {
  public:
    /**
     * Constructor.
     */
    MGVertexDoFs();

    /**
     * A function that is called to allocate the necessary amount of memory to
     * store the indices of the DoFs that live on this vertex for the given
     * (inclusive) range of levels.
     */
    void
    init(const unsigned int coarsest_level,
         const unsigned int finest_level,
         const unsigned int dofs_per_vertex);

    /**
     * Return the coarsest level for which this structure stores data.
     */
    unsigned int
    get_coarsest_level() const;

    /**
     * Return the finest level for which this structure stores data.
     */
    unsigned int
    get_finest_level() const;

    /**
     * Return the index of the <code>dof_number</code>th degree of freedom for
     * the given level stored for the current vertex.
     */
    types::global_dof_index &
    access_index(const unsigned int level,
                 const unsigned int dof_number,
                 const unsigned int dofs_per_vertex);

  private:
    /**
     * Coarsest level for which this object stores DoF indices.
     */
    unsigned int coarsest_level;

    /**
     * Finest level for which this object stores DoF indices.
     */
    unsigned int finest_level;

    /**
     * A pointer to an array where we store the indices of the DoFs that live
     * on the various levels this vertex exists on.
     *
     * The starting offset of the DoFs that belong to a @p level are given by
     * <code>n_dofs_per_vertex() * (level-coarsest_level)</code>. @p n_dofs_per_vertex()
     * must therefore be passed as an argument to the functions that set or
     * read an index.
     */
    std::unique_ptr<types::global_dof_index[]> indices;
  };

  /**
   * Whenever the underlying triangulation changes by either
   * h/p-refinement/coarsening and serialization, the active FE index of cells
   * needs to be transferred. This structure stores all temporary information
   * required during that process.
   */
  struct ActiveFEIndexTransfer
  {
    /**
     * Container to temporarily store the iterator and future active FE index
     * of cells that persist.
     */
    std::map<const cell_iterator, const types::fe_index>
      persisting_cells_fe_index;

    /**
     * Container to temporarily store the iterator and future active FE index
     * of cells that will be refined.
     */
    std::map<const cell_iterator, const types::fe_index> refined_cells_fe_index;

    /**
     * Container to temporarily store the iterator and future active FE index
     * of parent cells that will remain after coarsening.
     */
    std::map<const cell_iterator, const types::fe_index>
      coarsened_cells_fe_index;

    /**
     * Container to temporarily store the active FE index of every locally
     * owned cell for transfer across parallel::distributed::Triangulation
     * objects.
     */
    std::vector<types::fe_index> active_fe_indices;

    /**
     * Helper object to transfer all active FE indices on
     * parallel::distributed::Triangulation objects during
     * refinement/coarsening and serialization.
     */
    std::unique_ptr<
      parallel::distributed::
        CellDataTransfer<dim, spacedim, std::vector<types::fe_index>>>
      cell_data_transfer;
  };

  /**
   * An object containing information on the block structure.
   */
  BlockInfo block_info_object;

  /**
   * Boolean indicating whether or not the current DoFHandler has
   * hp-capabilities.
   */
  bool hp_capability_enabled;

  /**
   * Address of the triangulation to work on.
   */
  ObserverPointer<const Triangulation<dim, spacedim>, DoFHandler<dim, spacedim>>
    tria;

  /**
   * Store a hp::FECollection object. If only a single FiniteElement is
   * used during initialization of this object, it contains the (one)
   * FiniteElement.
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
   * A structure that contains all sorts of numbers that characterize the
   * degrees of freedom this object works on.
   *
   * For most members of this structure, there is an accessor function in this
   * class that returns its value.
   */
  dealii::internal::DoFHandlerImplementation::NumberCache number_cache;

  /**
   * Data structure like number_cache, but for each multigrid level.
   */
  std::vector<dealii::internal::DoFHandlerImplementation::NumberCache>
    mg_number_cache;

  /**
   * Indices of degree of freedom of each d+1 geometric object (3d: vertex,
   * line, quad, hex) for all relevant active finite elements. Identification
   * of the appropriate position is done via object_dof_ptr (CRS scheme).
   */
  mutable std::vector<std::array<std::vector<types::global_dof_index>, dim + 1>>
    object_dof_indices;

  /**
   * Pointer to the first cached degree of freedom of a geometric object for all
   * relevant active finite elements.
   *
   * @note In normal mode it is possible to access this data structure directly.
   *   In hp-mode, an indirection via hp_object_fe_indices/hp_object_fe_ptr is
   * necessary.
   */
  mutable std::vector<std::array<std::vector<offset_type>, dim + 1>>
    object_dof_ptr;

  /**
   * Active FE indices of each geometric object. Identification
   * of the appropriate position of a cell in the vectors is done via
   * hp_object_fe_ptr (CRS scheme).
   */
  mutable std::array<std::vector<types::fe_index>, dim + 1>
    hp_object_fe_indices;

  /**
   * Pointer to the first FE index of a geometric object.
   */
  mutable std::array<std::vector<offset_type>, dim + 1> hp_object_fe_ptr;

  /**
   * Active FE index of an active cell (identified by level and level index).
   * This vector is only used in hp-mode.
   */
  mutable std::vector<std::vector<types::fe_index>> hp_cell_active_fe_indices;

  /**
   * Future FE index of an active cell (identified by level and level index).
   * This vector is only used in hp-mode.
   */
  mutable std::vector<std::vector<types::fe_index>> hp_cell_future_fe_indices;

  /**
   * An array to store the indices for level degrees of freedom located at
   * vertices.
   */
  std::vector<MGVertexDoFs> mg_vertex_dofs;

  /**
   * Space to store the DoF numbers for the different multigrid levels.
   */
  std::vector<
    std::unique_ptr<dealii::internal::DoFHandlerImplementation::DoFLevel<dim>>>
    mg_levels;

  /**
   * Space to store DoF numbers of faces in the multigrid context.
   */
  std::unique_ptr<dealii::internal::DoFHandlerImplementation::DoFFaces<dim>>
    mg_faces;

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

  /**
   * A list of connections with which this object connects to the
   * triangulation. They get triggered specifically when data needs to be
   * transferred due to refinement or repartitioning. Only active in hp-mode.
   */
  std::vector<boost::signals2::connection> tria_listeners_for_transfer;

  /**
   * Free all memory used for non-multigrid data structures.
   */
  void
  clear_space();

  /**
   * Free all memory used for multigrid data structures.
   */
  void
  clear_mg_space();

  /**
   * Set up DoFHandler policy.
   */
  void
  setup_policy();

  /**
   * Set up connections to signals of the underlying triangulation.
   */
  void
  connect_to_triangulation_signals();

  /**
   * Create default tables for the active and future fe_indices.
   *
   * Active indices are initialized with a zero indicator, meaning that fe[0] is
   * going to be used by default. Future indices are initialized with an invalid
   * indicator, meaning that no p-adaptation is scheduled by default.
   *
   * This method is called upon construction and whenever the underlying
   * triangulation gets created. This ensures that each cell has a valid active
   * and future fe_index.
   */
  void
  create_active_fe_table();

  /**
   * Update tables for active and future fe_indices.
   *
   * Whenever the underlying triangulation changes (either by adaptation or
   * deserialization), active and future FE index tables will be adjusted to the
   * current structure of the triangulation. Missing values of active and future
   * indices will be initialized with their defaults (see
   * create_active_fe_table()).
   *
   * This method is called post refinement and post deserialization. This
   * ensures that each cell has a valid active and future fe_index.
   */
  void
  update_active_fe_table();

  /**
   * A function that will be triggered through a triangulation
   * signal just before the associated Triangulation or
   * parallel::shared::Triangulation is modified.
   *
   * The function that stores the active FE indices of all cells that will
   * be refined or coarsened before the refinement happens, so that
   * they can be set again after refinement.
   */
  void
  pre_transfer_action();

  /**
   * A function that will be triggered through a triangulation
   * signal just after the associated Triangulation or
   * parallel::shared::Triangulation is modified.
   *
   * The function that restores the active FE indices of all cells that
   * were refined or coarsened.
   */
  void
  post_transfer_action();

  /**
   * A function that will be triggered through a triangulation
   * signal just before the associated parallel::distributed::Triangulation is
   * modified.
   *
   * The function that stores all active FE indices on locally owned cells for
   * distribution over all participating processors.
   */
  void
  pre_distributed_transfer_action();

  /**
   * A function that will be triggered through a triangulation
   * signal just after the associated parallel::distributed::Triangulation is
   * modified.
   *
   * The function that restores all active FE indices on locally owned cells
   * that have been communicated.
   */
  void
  post_distributed_transfer_action();


  // Make accessor objects friends.
  template <int, int, int, bool>
  friend class dealii::DoFAccessor;
  template <int, int, bool>
  friend class dealii::DoFCellAccessor;
  friend struct dealii::internal::DoFAccessorImplementation::Implementation;
  friend struct dealii::internal::DoFCellAccessorImplementation::Implementation;

  // Likewise for DoFLevel objects since they need to access the vertex dofs
  // in the functions that set and retrieve vertex dof indices.
  friend struct dealii::internal::DoFHandlerImplementation::Implementation;
  friend struct dealii::internal::hp::DoFHandlerImplementation::Implementation;
  friend struct dealii::internal::DoFHandlerImplementation::Policy::
    Implementation;

  // explicitly check for sensible template arguments, but not on windows
  // because MSVC creates bogus warnings during normal compilation
#ifndef DEAL_II_MSVC
  static_assert(dim <= spacedim,
                "The dimension <dim> of a DoFHandler must be less than or "
                "equal to the space dimension <spacedim> in which it lives.");
#endif
};

namespace internal
{
  namespace hp
  {
    namespace DoFHandlerImplementation
    {
      /**
       * Given a DoFHandler object in hp-mode, make sure that the
       * future FE indices that a user has set for locally owned cells are
       * communicated to all other relevant cells as well.
       *
       * For parallel::shared::Triangulation objects,
       * this information is distributed on both ghost and artificial cells.
       *
       * In case a parallel::distributed::Triangulation is used,
       * indices are communicated only to ghost cells.
       */
      template <int dim, int spacedim>
      void
      communicate_future_fe_indices(DoFHandler<dim, spacedim> &dof_handler);

      /**
       * Return the index of the finite element from the entire hp::FECollection
       * that is dominated by those assigned as future finite elements to the
       * children of @p parent.
       *
       * We find the corresponding finite element among the future finite
       * elements on the children of this cell. If none of them qualify, we
       * extend our search on the whole hp::FECollection, which is the element
       * that describes the smallest finite element space that includes all
       * future finite elements assigned to the children. If the function is not
       * able to find a finite element at all, an assertion will be triggered.
       *
       * In this way, we determine the finite element of the parent cell in case
       * of h-coarsening in the hp-context.
       *
       * @note This function can only be called on direct parent cells, i.e.,
       * non-active cells whose children are all active.
       *
       * @note On parallel Triangulation objects where sibling cells
       * can be ghost cells, make sure that future FE indices have been properly
       * communicated with communicate_future_fe_indices() first. Otherwise,
       * results might differ on different processors. There is no check for
       * consistency of future FE indices.
       */
      template <int dim, int spacedim = dim>
      unsigned int
      dominated_future_fe_on_children(
        const typename DoFHandler<dim, spacedim>::cell_iterator &parent);

      /**
       * Exception
       */
      DeclExceptionMsg(
        ExcNoDominatedFiniteElementOnChildren,
        "No FiniteElement has been found in your FECollection that is "
        "dominated by all children of a cell you are trying to coarsen!");
    } // namespace DoFHandlerImplementation
  }   // namespace hp
} // namespace internal

#ifndef DOXYGEN

/* ----------------------- Inline functions ----------------------------------
 */


template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
inline bool DoFHandler<dim, spacedim>::has_hp_capabilities() const
{
  return hp_capability_enabled;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
inline bool DoFHandler<dim, spacedim>::has_level_dofs() const
{
  return mg_number_cache.size() > 0;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
inline bool DoFHandler<dim, spacedim>::has_active_dofs() const
{
  return number_cache.n_global_dofs > 0;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
inline types::global_dof_index DoFHandler<dim, spacedim>::n_dofs() const
{
  return number_cache.n_global_dofs;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
inline types::global_dof_index
  DoFHandler<dim, spacedim>::n_dofs(const unsigned int level) const
{
  Assert(has_level_dofs(),
         ExcMessage(
           "n_dofs(level) can only be called after distribute_mg_dofs()"));
  Assert(level < mg_number_cache.size(), ExcInvalidLevel(level));
  return mg_number_cache[level].n_global_dofs;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
types::global_dof_index DoFHandler<dim, spacedim>::n_locally_owned_dofs() const
{
  return number_cache.n_locally_owned_dofs;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
const IndexSet &DoFHandler<dim, spacedim>::locally_owned_dofs() const
{
  return number_cache.locally_owned_dofs;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
const IndexSet &DoFHandler<dim, spacedim>::locally_owned_mg_dofs(
  const unsigned int level) const
{
  Assert(level < this->get_triangulation().n_global_levels(),
         ExcMessage("The given level index exceeds the number of levels "
                    "present in the triangulation"));
  Assert(
    mg_number_cache.size() == this->get_triangulation().n_global_levels(),
    ExcMessage(
      "The level dofs are not set up properly! Did you call distribute_mg_dofs()?"));
  return mg_number_cache[level].locally_owned_dofs;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
inline const FiniteElement<dim, spacedim> &DoFHandler<dim, spacedim>::get_fe(
  const types::fe_index number) const
{
  Assert(fe_collection.size() > 0,
         ExcMessage("No finite element collection is associated with "
                    "this DoFHandler"));
  return fe_collection[number];
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
inline const hp::FECollection<dim, spacedim>
  &DoFHandler<dim, spacedim>::get_fe_collection() const
{
  return fe_collection;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
inline const Triangulation<dim, spacedim>
  &DoFHandler<dim, spacedim>::get_triangulation() const
{
  Assert(tria != nullptr,
         ExcMessage("This DoFHandler object has not been associated "
                    "with a triangulation."));
  return *tria;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
inline MPI_Comm DoFHandler<dim, spacedim>::get_mpi_communicator() const
{
  Assert(tria != nullptr,
         ExcMessage("This DoFHandler object has not been associated "
                    "with a triangulation."));
  return tria->get_mpi_communicator();
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
inline MPI_Comm DoFHandler<dim, spacedim>::get_communicator() const
{
  return get_mpi_communicator();
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
inline const BlockInfo &DoFHandler<dim, spacedim>::block_info() const
{
  Assert(this->hp_capability_enabled == false, ExcNotImplementedWithHP());

  return block_info_object;
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
template <typename number>
types::global_dof_index DoFHandler<dim, spacedim>::n_boundary_dofs(
  const std::map<types::boundary_id, const Function<spacedim, number> *>
    &boundary_ids) const
{
  Assert(!(dim == 2 && spacedim == 3) || this->hp_capability_enabled == false,
         ExcNotImplementedWithHP());

  // extract the set of boundary ids and forget about the function object
  // pointers
  std::set<types::boundary_id> boundary_ids_only;
  for (typename std::map<types::boundary_id,
                         const Function<spacedim, number> *>::const_iterator p =
         boundary_ids.begin();
       p != boundary_ids.end();
       ++p)
    boundary_ids_only.insert(p->first);

  // then just hand everything over to the other function that does the work
  return n_boundary_dofs(boundary_ids_only);
}



namespace internal
{
  /**
   * Return a string representing the dynamic type of the given argument.
   * This is basically the same what typeid(...).name() does, but it turns out
   * this is broken on Intel 13+.
   *
   * Defined in dof_handler.cc.
   */
  template <int dim, int spacedim>
  std::string
  policy_to_string(const dealii::internal::DoFHandlerImplementation::Policy::
                     PolicyBase<dim, spacedim> &policy);
} // namespace internal



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
template <class Archive>
void DoFHandler<dim, spacedim>::save(Archive &ar, const unsigned int) const
{
  if (this->hp_capability_enabled)
    {
      ar &this->object_dof_indices;
      ar &this->object_dof_ptr;

      ar &this->hp_cell_active_fe_indices;
      ar &this->hp_cell_future_fe_indices;

      ar &hp_object_fe_ptr;
      ar &hp_object_fe_indices;

      ar &number_cache;

      ar &mg_number_cache;

      // write out the number of triangulation cells and later check during
      // loading that this number is indeed correct; same with something that
      // identifies the policy
      const unsigned int n_cells = this->tria->n_cells();
      std::string        policy_name =
        dealii::internal::policy_to_string(*this->policy);

      ar &n_cells &policy_name;
    }
  else
    {
      ar &this->block_info_object;
      ar &number_cache;

      ar &this->object_dof_indices;
      ar &this->object_dof_ptr;

      // write out the number of triangulation cells and later check during
      // loading that this number is indeed correct; same with something that
      // identifies the FE and the policy
      unsigned int n_cells     = this->tria->n_cells();
      std::string  fe_name     = this->get_fe(0).get_name();
      std::string  policy_name = internal::policy_to_string(*this->policy);

      ar &n_cells &fe_name &policy_name;
    }
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
template <class Archive>
void DoFHandler<dim, spacedim>::load(Archive &ar, const unsigned int)
{
  if (this->hp_capability_enabled)
    {
      ar &this->object_dof_indices;
      ar &this->object_dof_ptr;

      ar &this->hp_cell_active_fe_indices;
      ar &this->hp_cell_future_fe_indices;

      ar &hp_object_fe_ptr;
      ar &hp_object_fe_indices;

      ar &number_cache;

      ar &mg_number_cache;

      // these are the checks that correspond to the last block in the save()
      // function
      unsigned int n_cells;
      std::string  policy_name;

      ar &n_cells &policy_name;

      AssertThrow(
        n_cells == this->tria->n_cells(),
        ExcMessage(
          "The object being loaded into does not match the triangulation "
          "that has been stored previously."));
      AssertThrow(
        policy_name == dealii::internal::policy_to_string(*this->policy),
        ExcMessage("The policy currently associated with this DoFHandler (" +
                   dealii::internal::policy_to_string(*this->policy) +
                   ") does not match the one that was associated with the "
                   "DoFHandler previously stored (" +
                   policy_name + ")."));
    }
  else
    {
      ar &this->block_info_object;
      ar &number_cache;

      object_dof_indices.clear();

      object_dof_ptr.clear();

      ar &this->object_dof_indices;
      ar &this->object_dof_ptr;

      // these are the checks that correspond to the last block in the save()
      // function
      unsigned int n_cells;
      std::string  fe_name;
      std::string  policy_name;

      ar &n_cells &fe_name &policy_name;

      AssertThrow(
        n_cells == this->tria->n_cells(),
        ExcMessage(
          "The object being loaded into does not match the triangulation "
          "that has been stored previously."));
      AssertThrow(
        fe_name == this->get_fe(0).get_name(),
        ExcMessage(
          "The finite element associated with this DoFHandler does not match "
          "the one that was associated with the DoFHandler previously stored."));
      AssertThrow(policy_name == internal::policy_to_string(*this->policy),
                  ExcMessage(
                    "The policy currently associated with this DoFHandler (" +
                    internal::policy_to_string(*this->policy) +
                    ") does not match the one that was associated with the "
                    "DoFHandler previously stored (" +
                    policy_name + ")."));
    }
}



template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
inline types::global_dof_index
  &DoFHandler<dim, spacedim>::MGVertexDoFs::access_index(
    const unsigned int level,
    const unsigned int dof_number,
    const unsigned int dofs_per_vertex)
{
  Assert((level >= coarsest_level) && (level <= finest_level),
         ExcInvalidLevel(level));
  return indices[dofs_per_vertex * (level - coarsest_level) + dof_number];
}



// Declare the existence of explicit instantiations of the class
// above. This is not strictly necessary, but tells the compiler to
// avoid instantiating templates that we know are instantiated in
// .cc files and so can be referenced without implicit
// instantiations.
//
// Unfortunately, this does not seem to work when building modules
// because the compiler (well, Clang at least) then just doesn't
// instantiate these classes at all, even though their members are
// defined and explicitly instantiated in a .cc file.
#  ifndef DEAL_II_BUILDING_CXX20_MODULE
extern template class DoFHandler<1, 1>;
extern template class DoFHandler<1, 2>;
extern template class DoFHandler<1, 3>;
extern template class DoFHandler<2, 2>;
extern template class DoFHandler<2, 3>;
extern template class DoFHandler<3, 3>;
#  endif

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
