// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_distributed_tria_base_h
#define dealii_distributed_tria_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/mpi_stub.h>
#include <deal.II/base/observer_pointer.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/grid/tria.h>

#include <functional>
#include <list>
#include <set>
#include <utility>
#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  /**
   * This class describes the interface for all triangulation classes that
   * work in parallel, namely parallel::distributed::Triangulation,
   * parallel::fullydistributed::Triangulation, and
   * parallel::shared::Triangulation.
   *
   * It is, consequently, a class that can be used to test whether a
   * pointer of reference to a triangulation object refers to a
   * sequential triangulation, or whether the triangulation is in fact
   * parallel. In other words, one could write a function like this:
   * @code
   *   template <int dim, int spacedim>
   *   bool is_parallel (const dealii::Triangulation<dim,spacedim> &tria)
   *   {
   *     if (dynamic_cast<const parallel::TriangulationBase<dim,spacedim>*>
   *                     (&tria)
   *         != nullptr)
   *       return true;
   *     else
   *       return false;
   *   }
   * @endcode
   *
   * All parallel triangulations share certain traits, such as the fact that
   * they communicate via
   * @ref GlossMPICommunicator "MPI communicators"
   * or that they have
   * @ref GlossLocallyOwnedCell "locally owned",
   * @ref GlossGhostCell "ghost",
   * and possibly
   * @ref GlossArtificialCell "artificial cells".
   * This class provides
   * a number of member functions that allows querying some information
   * about the triangulation that is independent of how exactly a
   * parallel triangulation is implemented (i.e., which of the various
   * classes derived from the current one it actually is).
   */
  template <int dim, int spacedim = dim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  class TriangulationBase : public dealii::Triangulation<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    TriangulationBase(
      const MPI_Comm mpi_communicator,
      const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
                 smooth_grid = (dealii::Triangulation<dim, spacedim>::none),
      const bool check_for_distorted_cells = false);

    /**
     * Destructor.
     */
    virtual ~TriangulationBase() override;

    /**
     * Return MPI communicator used by this triangulation.
     */
    virtual MPI_Comm
    get_mpi_communicator() const override;

    /**
     * Return if multilevel hierarchy is supported and has been constructed.
     */
    virtual bool
    is_multilevel_hierarchy_constructed() const = 0;

    /**
     * Implementation of the same function as in the base class.
     *
     * @note This function copies the cells, but not the communicator,
     * of the source triangulation. In other words, the resulting
     * triangulation will operate on the communicator it was constructed
     * with.
     */
    virtual void
    copy_triangulation(
      const dealii::Triangulation<dim, spacedim> &old_tria) override;

    /**
     * Return the number of active cells in the triangulation that are locally
     * owned, i.e. that have a subdomain_id equal to
     * locally_owned_subdomain(). Note that there may be more active cells in
     * the triangulation stored on the present processor, such as for example
     * ghost cells, or cells further away from the locally owned block of
     * cells but that are needed to ensure that the triangulation that stores
     * this processor's set of active cells still remains balanced with
     * respect to the 2:1 size ratio of adjacent cells.
     *
     * As a consequence of the remark above, the result of this function is
     * always smaller or equal to the result of the function with the same
     * name in the ::Triangulation base class, which includes the active ghost
     * and artificial cells (see also
     * @ref GlossArtificialCell
     * and
     * @ref GlossGhostCell).
     */
    unsigned int
    n_locally_owned_active_cells() const;

    /**
     * Return the sum over all processors of the number of active cells owned
     * by each processor. This equals the overall number of active cells in
     * the triangulation.
     */
    virtual types::global_cell_index
    n_global_active_cells() const override;

    /**
     * Return the local memory consumption in bytes.
     */
    virtual std::size_t
    memory_consumption() const override;


    /**
     * Return the global maximum level. This may be bigger than the number
     * dealii::Triangulation::n_levels() (a function in this class's base
     * class) returns if the current processor only stores cells in parts of
     * the domain that are not very refined, but if other processors store
     * cells in more deeply refined parts of the domain.
     */
    virtual unsigned int
    n_global_levels() const override;

    /**
     * Return the subdomain id of those cells that are owned by the current
     * processor. All cells in the triangulation that do not have this
     * subdomain id are either owned by another processor or have children
     * that only exist on other processors.
     */
    types::subdomain_id
    locally_owned_subdomain() const override;

    /**
     * Return a set of MPI ranks of the processors that have at least one
     * ghost cell adjacent to the cells of the local processor. In other
     * words, this is the set of subdomain_id() for all ghost cells.
     *
     * The returned sets are symmetric, that is if @p i is contained in the
     * list of processor @p j, then @p j will also be contained in the list of
     * processor @p i.
     */
    const std::set<types::subdomain_id> &
    ghost_owners() const;

    /**
     * Return a set of MPI ranks of the processors that have at least one
     * level ghost cell adjacent to our cells used in geometric multigrid. In
     * other words, this is the set of level_subdomain_id() for all level
     * ghost cells.
     *
     * The returned sets are symmetric, that is if @p i is contained in the
     * list of processor @p j, then @p j will also be contained in the list of
     * processor @p i.
     *
     * @note The level ghost owners can only be determined if the multigrid
     * ownership has been assigned (by setting the
     * construct_multigrid_hierarchy flag at construction time), otherwise the
     * returned set will be empty.
     */
    const std::set<types::subdomain_id> &
    level_ghost_owners() const;

    std::weak_ptr<const Utilities::MPI::Partitioner>
    global_active_cell_index_partitioner() const override;

    std::weak_ptr<const Utilities::MPI::Partitioner>
    global_level_cell_index_partitioner(
      const unsigned int level) const override;

    /**
     * Return a vector containing all boundary indicators assigned to boundary
     * faces of active cells of this Triangulation object. Note, that each
     * boundary indicator is reported only once. The size of the return vector
     * will represent the number of different indicators (which is greater or
     * equal one).
     *
     * @see
     * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
     *
     * @note This function involves a global communication gathering all current
     *   IDs from all processes.
     */
    virtual std::vector<types::boundary_id>
    get_boundary_ids() const override;

    /**
     * Return a vector containing all manifold indicators assigned to the
     * objects of the active cells of this Triangulation. Note, that each
     * manifold indicator is reported only once. The size of the return vector
     * will represent the number of different indicators (which is greater or
     * equal one).
     *
     * @ingroup manifold
     *
     * @see
     * @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
     *
     * @note This function involves a global communication gathering all current
     *   IDs from all processes.
     */
    virtual std::vector<types::manifold_id>
    get_manifold_ids() const override;

    /**
     * When vertices have been moved locally, for example using code like
     * @code
     *   cell->vertex(0) = new_location;
     * @endcode
     * then this function can be used to update the location of vertices
     * between MPI processes.
     *
     * All the vertices that have been moved and might be in the ghost layer
     * of a process have to be reported in the @p vertex_locally_moved
     * argument. This ensures that that part of the information that has to
     * be send between processes is actually sent. Additionally, it is quite
     * important that vertices on the boundary between processes are
     * reported on exactly one process (e.g. the one with the highest id).
     * Otherwise we could expect undesirable results if multiple processes
     * move a vertex differently. A typical strategy is to let processor $i$
     * move those vertices that are adjacent to cells whose owners include
     * processor $i$ but no other processor $j$ with $j<i$; in other words,
     * for vertices at the boundary of a subdomain, the processor with the
     * lowest subdomain id "owns" a vertex.
     *
     * @note It only makes sense to move vertices that are either located on
     * locally owned cells or on cells in the ghost layer. This is because
     * you can be sure that these vertices indeed exist on the finest mesh
     * aggregated over all processors, whereas vertices on artificial cells
     * but not at least in the ghost layer may or may not exist on the
     * globally finest mesh. Consequently, the @p vertex_locally_moved
     * argument may not contain vertices that aren't at least on ghost
     * cells.
     *
     * @note This function moves vertices in such a way that on every
     * processor, the vertices of every locally owned and ghost cell is
     * consistent with the corresponding location of these cells on other
     * processors. On the other hand, the locations of artificial cells will
     * in general be wrong since artificial cells may or may not exist on
     * other processors and consequently it is not possible to determine
     * their location in any way. This is not usually a problem since one
     * never does anything on artificial cells. However, it may lead to
     * problems if the mesh with moved vertices is refined in a later step.
     * If that's what you want to do, the right way to do it is to save the
     * offset applied to every vertex, call this function, and before
     * refining or coarsening the mesh apply the opposite offset and call
     * this function again.
     *
     * @param vertex_locally_moved A bitmap indicating which vertices have
     * been moved. The size of this array must be equal to
     * Triangulation::n_vertices() and must be a subset of those vertices
     * flagged by GridTools::get_locally_owned_vertices().
     *
     * @see This function is used, for example, in
     * GridTools::distort_random().
     */
    void
    communicate_locally_moved_vertices(
      const std::vector<bool> &vertex_locally_moved);

    virtual types::coarse_cell_id
    n_global_coarse_cells() const override;

    /**
     * Reset this triangulation to an empty state by deleting all data.
     *
     * Note that this operation is only allowed if no subscriptions to this
     * object exist any more, such as DoFHandler objects using it.
     */
    virtual void
    clear() override;

  protected:
    /**
     * MPI communicator to be used for the triangulation. We create a unique
     * communicator for this class, which is a duplicate of the one passed to
     * the constructor.
     */
    const MPI_Comm mpi_communicator;

    /**
     * The subdomain id to be used for the current processor. This is the MPI
     * rank.
     */
    types::subdomain_id my_subdomain;

    /**
     * The total number of subdomains (or the size of the MPI communicator).
     */
    types::subdomain_id n_subdomains;

    /**
     * A structure that contains information about the distributed
     * triangulation.
     */
    struct NumberCache
    {
      /**
       * Number of locally owned active cells of this MPI rank.
       */
      unsigned int n_locally_owned_active_cells;

      /**
       * The total number of active cells (sum of @p
       * n_locally_owned_active_cells).
       */
      types::global_cell_index n_global_active_cells;

      /**
       * Number of global coarse cells.
       */
      types::coarse_cell_id number_of_global_coarse_cells;

      /**
       * The global number of levels computed as the maximum number of levels
       * taken over all MPI ranks, so <tt>n_levels()<=n_global_levels =
       * max(n_levels() on proc i)</tt>.
       */
      unsigned int n_global_levels;

      /**
       * A set containing the subdomain_id (MPI rank) of the owners of the
       * ghost cells on this processor.
       */
      std::set<types::subdomain_id> ghost_owners;

      /**
       * A set containing the MPI ranks of the owners of the level ghost cells
       * on this processor (for all levels).
       */
      std::set<types::subdomain_id> level_ghost_owners;

      /**
       * Partitioner for the global active cell indices.
       */
      std::shared_ptr<const Utilities::MPI::Partitioner>
        active_cell_index_partitioner;

      /**
       * Partitioner for the global level cell indices for each level.
       */
      std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
        level_cell_index_partitioners;

      NumberCache();
    };

    NumberCache number_cache;

    /**
     * Update the number_cache variable after mesh creation or refinement.
     */
    virtual void
    update_number_cache();

    /**
     * Reset global active cell indices and global level cell indices.
     */
    void
    reset_global_cell_indices();

    void
    update_reference_cells() override;
  };


  /**
   * A base class for distributed triangulations, i.e., triangulations that
   * do not store all cells on all processors. This implies that not
   * every detail of a triangulation may be known on each processor.
   * In particular, you have to expect that triangulations of classes
   * derived from this one only store some of the active cells (namely,
   * the
   * @ref GlossLocallyOwnedCell "locally owned cells"),
   * along with
   * @ref GlossGhostCell "ghost cells"
   * and possibly
   * @ref GlossArtificialCell "artificial cells".
   * In contrast to the classes
   * derived from parallel::TriangulationBase, it is certain that the
   * classes derived from the current class will not store the entire
   * triangulation as long as it has a large enough number of cells. (The
   * difference to parallel::TriangulationBase is that the
   * parallel::shared::Triangulation is derived from
   * parallel::TriangulationBase, but not from the current class.) The
   * distinction is not large in practice: Everything that is difficult for
   * parallel distributed triangulation is generally also difficult for any
   * other kind of parallel triangulation classes; however, this intermediate
   * base class allows to further differentiate between the different kinds of
   * classes providing parallel mesh functionality.
   *
   * This class can, then, be used to test whether a
   * pointer or reference to a triangulation object refers to any kind of
   * parallel triangulation, or whether the triangulation is in fact
   * parallel distributed. In other words, one could write a function like
   * this:
   * @code
   *   template <int dim, int spacedim>
   *   bool
   *   is_parallel_distributed(const dealii::Triangulation<dim,spacedim> &tria)
   *   {
   *     if(dynamic_cast<const
   *                     parallel::DistributedTriangulationBase<dim,spacedim>*>
   *                    (&tria)
   *        != nullptr)
   *       return true;
   *     else
   *       return false;
   *   }
   * @endcode
   *
   * @dealiiConceptRequires{(concepts::is_valid_dim_spacedim<dim, spacedim>)}
   */
  template <int dim, int spacedim = dim>
  DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
  class DistributedTriangulationBase
    : public dealii::parallel::TriangulationBase<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    DistributedTriangulationBase(
      const MPI_Comm mpi_communicator,
      const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
                 smooth_grid = (dealii::Triangulation<dim, spacedim>::none),
      const bool check_for_distorted_cells = false);

    using cell_iterator =
      typename dealii::Triangulation<dim, spacedim>::cell_iterator;

    /**
     * Return true if the triangulation has hanging nodes.
     *
     * In the context of parallel distributed triangulations, every
     * processor stores only that part of the triangulation it owns locally.
     * However, it also stores coarser levels, and to guarantee the
     * 2:1 relationship between cells, this may mean that there are hanging
     * nodes between cells that are not locally owned or ghost cells (i.e.,
     * between ghost cells and artificial cells, or between artificial and
     * artificial cells; see
     * @ref GlossArtificialCell "the glossary").
     * One is not typically interested in this case, so the function returns
     * whether there are hanging nodes between any two cells of the "global"
     * mesh, i.e., the union of locally owned cells on all processors.
     */
    virtual bool
    has_hanging_nodes() const override;

    /**
     * Save the triangulation into the given file. This file needs to be
     * reachable from all nodes in the computation on a shared network file
     * system. See the SolutionTransfer class on how to store solution vectors
     * into this file. Additional cell-based data can be saved using
     * register_data_attach().
     */
    using Triangulation<dim, spacedim>::save;


    /**
     * Load the triangulation saved with save() back in. Cell-based data that
     * was saved with register_data_attach() can be read in with
     * notify_ready_to_unpack() after calling Triangulation<dim,
     * spacedim>::load.
     */
    using Triangulation<dim, spacedim>::load;
  };

} // namespace parallel

DEAL_II_NAMESPACE_CLOSE

#endif
