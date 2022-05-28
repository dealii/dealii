// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2022 by the deal.II authors
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

#ifndef dealii_distributed_tria_base_h
#define dealii_distributed_tria_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>
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
  class TriangulationBase : public dealii::Triangulation<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    TriangulationBase(
      const MPI_Comm &mpi_communicator,
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
    get_communicator() const override;

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

    /**
     * Return partitioner for the global indices of the cells on the active
     * level of the triangulation.
     */
    const std::weak_ptr<const Utilities::MPI::Partitioner>
    global_active_cell_index_partitioner() const;

    /**
     * Return partitioner for the global indices of the cells on the given @p
     * level of the triangulation.
     */
    const std::weak_ptr<const Utilities::MPI::Partitioner>
    global_level_cell_index_partitioner(const unsigned int level) const;

    /**
     * @copydoc dealii::Triangulation::get_boundary_ids()
     *
     * @note This function involves a global communication gathering all current
     *   IDs from all processes.
     */
    virtual std::vector<types::boundary_id>
    get_boundary_ids() const override;

    /**
     * @copydoc dealii::Triangulation::get_manifold_ids()
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
     * @copydoc dealii::Triangulation::update_reference_cells()
     */
    void
    update_reference_cells() override;

    /**
     * Reset global active cell indices and global level cell indices.
     */
    void
    reset_global_cell_indices();
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
   */
  template <int dim, int spacedim = dim>
  class DistributedTriangulationBase
    : public dealii::parallel::TriangulationBase<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    DistributedTriangulationBase(
      const MPI_Comm &mpi_communicator,
      const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
                 smooth_grid = (dealii::Triangulation<dim, spacedim>::none),
      const bool check_for_distorted_cells = false);

    /**
     * Reset this triangulation into a virgin state by deleting all data.
     *
     * Note that this operation is only allowed if no subscriptions to this
     * object exist any more, such as DoFHandler objects using it.
     */
    virtual void
    clear() override;

    using cell_iterator =
      typename dealii::Triangulation<dim, spacedim>::cell_iterator;

    using CellStatus =
      typename dealii::Triangulation<dim, spacedim>::CellStatus;

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
    virtual void
    save(const std::string &filename) const = 0;

    /**
     * Load the triangulation saved with save() back in. Cell-based data that
     * was saved with register_data_attach() can be read in with
     * notify_ready_to_unpack() after calling load().
     */
    virtual void
    load(const std::string &filename) = 0;

    /**
     * @copydoc load()
     *
     * @deprecated The autopartition parameter has been removed.
     */
    DEAL_II_DEPRECATED_EARLY
    virtual void
    load(const std::string &filename, const bool autopartition) = 0;

    /**
     * Register a function that can be used to attach data of fixed size
     * to cells. This is useful for two purposes: (i) Upon refinement and
     * coarsening of a triangulation (@a e.g. in
     * parallel::distributed::Triangulation::execute_coarsening_and_refinement()),
     * one needs to be able to store one or more data vectors per cell that
     * characterizes the solution values on the cell so that this data can
     * then be transferred to the new owning processor of the cell (or
     * its parent/children) when the mesh is re-partitioned; (ii) when
     * serializing a computation to a file, it is necessary to attach
     * data to cells so that it can be saved (@a e.g. in
     * parallel::distributed::Triangulation::save()) along with the cell's
     * other information and, if necessary, later be reloaded from disk
     * with a different subdivision of cells among the processors.
     *
     * The way this function works is that it allows any number of interest
     * parties to register their intent to attach data to cells. One example
     * of classes that do this is parallel::distributed::SolutionTransfer
     * where each parallel::distributed::SolutionTransfer object that works
     * on the current Triangulation object then needs to register its intent.
     * Each of these parties registers a callback function (the first
     * argument here, @p pack_callback) that will be called whenever the
     * triangulation's execute_coarsening_and_refinement() or save()
     * functions are called.
     *
     * The current function then returns an integer handle that corresponds
     * to the number of data set that the callback provided here will attach.
     * While this number could be given a precise meaning, this is
     * not important: You will never actually have to do anything with
     * this number except return it to the notify_ready_to_unpack() function.
     * In other words, each interested party (i.e., the caller of the current
     * function) needs to store their respective returned handle for later use
     * when unpacking data in the callback provided to
     * notify_ready_to_unpack().
     *
     * Whenever @p pack_callback is then called by
     * execute_coarsening_and_refinement() or load() on a given cell, it
     * receives a number of arguments. In particular, the first
     * argument passed to the callback indicates the cell for which
     * it is supposed to attach data. This is always an active cell.
     *
     * The second, CellStatus, argument provided to the callback function
     * will tell you if the given cell will be coarsened, refined, or will
     * persist as is. (This status may be different than the refinement
     * or coarsening flags set on that cell, to accommodate things such as
     * the "one hanging node per edge" rule.). These flags need to be
     * read in context with the p4est quadrant they belong to, as their
     * relations are gathered in local_cell_relations.
     *
     * Specifically, the values for this argument mean the following:
     *
     * - `CELL_PERSIST`: The cell won't be refined/coarsened, but might be
     * moved to a different processor. If this is the case, the callback
     * will want to pack up the data on this cell into an array and store
     * it at the provided address for later unpacking wherever this cell
     * may land.
     * - `CELL_REFINE`: This cell will be refined into 4 or 8 cells (in 2d
     * and 3d, respectively). However, because these children don't exist
     * yet, you cannot access them at the time when the callback is
     * called. Thus, in local_cell_relations, the corresponding
     * p4est quadrants of the children cells are linked to the deal.II
     * cell which is going to be refined. To be specific, only the very
     * first child is marked with `CELL_REFINE`, whereas the others will be
     * marked with `CELL_INVALID`, which indicates that these cells will be
     * ignored by default during the packing or unpacking process. This
     * ensures that data is only transferred once onto or from the parent
     * cell. If the callback is called with `CELL_REFINE`, the callback
     * will want to pack up the data on this cell into an array and store
     * it at the provided address for later unpacking in a way so that
     * it can then be transferred to the children of the cell that will
     * then be available. In other words, if the data the callback
     * will want to pack up corresponds to a finite element field, then
     * the prolongation from parent to (new) children will have to happen
     * during unpacking.
     * - `CELL_COARSEN`: The children of this cell will be coarsened into the
     * given cell. These children still exist, so if this is the value
     * given to the callback as second argument, the callback will want
     * to transfer data from the children to the current parent cell and
     * pack it up so that it can later be unpacked again on a cell that
     * then no longer has any children (and may also be located on a
     * different processor). In other words, if the data the callback
     * will want to pack up corresponds to a finite element field, then
     * it will need to do the restriction from children to parent at
     * this point.
     * - `CELL_INVALID`: See `CELL_REFINE`.
     *
     * @note If this function is used for serialization of data
     *   using save() and load(), then the cell status argument with which
     *   the callback is called will always be `CELL_PERSIST`.
     *
     * The callback function is expected to return a memory chunk of the
     * format `std::vector<char>`, representing the packed data on a
     * certain cell.
     *
     * The second parameter @p returns_variable_size_data indicates whether
     * the returned size of the memory region from the callback function
     * varies by cell (<tt>=true</tt>) or stays constant on each one
     * throughout the whole domain (<tt>=false</tt>).
     *
     * @note The purpose of this function is to register intent to
     *   attach data for a single, subsequent call to
     *   execute_coarsening_and_refinement() and notify_ready_to_unpack(),
     *   save(), load(). Consequently, notify_ready_to_unpack(), save(),
     *   and load() all forget the registered callbacks once these
     *   callbacks have been called, and you will have to re-register
     *   them with a triangulation if you want them to be active for
     *   another call to these functions.
     */
    unsigned int
    register_data_attach(
      const std::function<std::vector<char>(const cell_iterator &,
                                            const CellStatus)> &pack_callback,
      const bool returns_variable_size_data);

    /**
     * This function is the opposite of register_data_attach(). It is called
     * <i>after</i> the execute_coarsening_and_refinement() or save()/load()
     * functions are done when classes and functions that have previously
     * attached data to a triangulation for either transfer to other
     * processors, across mesh refinement, or serialization of data to
     * a file are ready to receive that data back. The important part about
     * this process is that the triangulation cannot do this right away from
     * the end of execute_coarsening_and_refinement() or load() via a
     * previously attached callback function (as the register_data_attach()
     * function does) because the classes that eventually want the data
     * back may need to do some setup between the point in time where the
     * mesh has been recreated and when the data can actually be received.
     * An example is the parallel::distributed::SolutionTransfer class
     * that can really only receive the data once not only the mesh is
     * completely available again on the current processor, but only
     * after a DoFHandler has been reinitialized and distributed
     * degrees of freedom. In other words, there is typically a significant
     * amount of set up that needs to happen in user space before the classes
     * that can receive data attached to cell are ready to actually do so.
     * When they are, they use the current function to tell the triangulation
     * object that now is the time when they are ready by calling the
     * current function.
     *
     * The supplied callback function is then called for each newly locally
     * owned cell. The first argument to the callback is an iterator that
     * designates the cell; the second argument indicates the status of the
     * cell in question; and the third argument localizes a memory area by
     * two iterators that contains the data that was previously saved from
     * the callback provided to register_data_attach().
     *
     * The CellStatus will indicate if the cell was refined, coarsened, or
     * persisted unchanged. The @p cell_iterator argument to the callback
     * will then either be an active,
     * locally owned cell (if the cell was not refined), or the immediate
     * parent if it was refined during execute_coarsening_and_refinement().
     * Therefore, contrary to during register_data_attach(), you can now
     * access the children if the status is `CELL_REFINE` but no longer for
     * callbacks with status `CELL_COARSEN`.
     *
     * The first argument to this function, `handle`, corresponds to
     * the return value of register_data_attach(). (The precise
     * meaning of what the numeric value of this handle is supposed
     * to represent is neither important, nor should you try to use
     * it for anything other than transmit information between a
     * call to register_data_attach() to the corresponding call to
     * notify_ready_to_unpack().)
     */
    void
    notify_ready_to_unpack(
      const unsigned int handle,
      const std::function<
        void(const cell_iterator &,
             const CellStatus,
             const boost::iterator_range<std::vector<char>::const_iterator> &)>
        &unpack_callback);

  protected:
    /**
     * Save additional cell-attached data into the given file. The first
     * arguments are used to determine the offsets where to write buffers to.
     *
     * Called by
     * @ref save.
     */
    void
    save_attached_data(const unsigned int global_first_cell,
                       const unsigned int global_num_cells,
                       const std::string &filename) const;

    /**
     * Load additional cell-attached data from the given file, if any was saved.
     * The first arguments are used to determine the offsets where to read
     * buffers from.
     *
     * Called by
     * @ref load.
     */
    void
    load_attached_data(const unsigned int global_first_cell,
                       const unsigned int global_num_cells,
                       const unsigned int local_num_cells,
                       const std::string &filename,
                       const unsigned int n_attached_deserialize_fixed,
                       const unsigned int n_attached_deserialize_variable);

    /**
     * A function to record the CellStatus of currently active cells that
     * are locally owned. This information is mandatory to transfer data
     * between meshes during adaptation or serialization, e.g., using
     * parallel::distributed::SolutionTransfer.
     *
     * Relations will be stored in the private member local_cell_relations. For
     * an extensive description of CellStatus, see the documentation for the
     * member function register_data_attach().
     */
    virtual void
    update_cell_relations() = 0;

    /**
     * Auxiliary data structure for assigning a CellStatus to a deal.II cell
     * iterator. For an extensive description of the former, see the
     * documentation for the member function register_data_attach().
     */
    using cell_relation_t = typename std::pair<cell_iterator, CellStatus>;

    /**
     * Vector of pairs, each containing a deal.II cell iterator and its
     * respective CellStatus. To update its contents, use the
     * update_cell_relations() member function.
     */
    std::vector<cell_relation_t> local_cell_relations;

    /**
     * A structure that stores information about the data that has been, or
     * will be, attached to cells via the register_data_attach() function
     * and later retrieved via notify_ready_to_unpack().
     */
    struct CellAttachedData
    {
      /**
       * number of functions that get attached to the Triangulation through
       * register_data_attach() for example SolutionTransfer.
       */
      unsigned int n_attached_data_sets;

      /**
       * number of functions that need to unpack their data after a call from
       * load()
       */
      unsigned int n_attached_deserialize;

      using pack_callback_t = std::function<std::vector<char>(
        typename dealii::Triangulation<dim, spacedim>::cell_iterator,
        typename dealii::Triangulation<dim, spacedim>::CellStatus)>;

      /**
       * These callback functions will be stored in the order in which they
       * have been registered with the register_data_attach() function.
       */
      std::vector<pack_callback_t> pack_callbacks_fixed;
      std::vector<pack_callback_t> pack_callbacks_variable;
    };

    CellAttachedData cell_attached_data;

    /**
     * This class in the private scope of parallel::DistributedTriangulationBase
     * is dedicated to the data transfer across repartitioned meshes
     * and to/from the file system.
     *
     * It is designed to store all data buffers intended for transfer.
     */
    class DataTransfer
    {
    public:
      DataTransfer(const MPI_Comm &mpi_communicator);

      /**
       * Prepare data transfer by calling the pack callback functions on each
       * cell
       * in @p cell_relations.
       *
       * All registered callback functions in @p pack_callbacks_fixed will write
       * into the fixed size buffer, whereas each entry of @p pack_callbacks_variable
       * will write its data into the variable size buffer.
       */
      void
      pack_data(const std::vector<cell_relation_t> &cell_relations,
                const std::vector<typename CellAttachedData::pack_callback_t>
                  &pack_callbacks_fixed,
                const std::vector<typename CellAttachedData::pack_callback_t>
                  &pack_callbacks_variable);



      /**
       * Unpack the CellStatus information on each entry of
       * @p cell_relations.
       *
       * Data has to be previously transferred with execute_transfer()
       * or read from the file system via load().
       */
      void
      unpack_cell_status(std::vector<cell_relation_t> &cell_relations) const;

      /**
       * Unpack previously transferred data on each cell registered in
       * @p cell_relations with the provided @p unpack_callback function.
       *
       * The parameter @p handle corresponds to the position where the
       * @p unpack_callback function is allowed to read from the memory. Its
       * value needs to be in accordance with the corresponding pack_callback
       * function that has been registered previously.
       *
       * Data has to be previously transferred with execute_transfer()
       * or read from the file system via load().
       */
      void
      unpack_data(
        const std::vector<cell_relation_t> &cell_relations,
        const unsigned int                  handle,
        const std::function<void(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
          const typename dealii::Triangulation<dim, spacedim>::CellStatus &,
          const boost::iterator_range<std::vector<char>::const_iterator> &)>
          &unpack_callback) const;

      /**
       * Transfer data to file system.
       *
       * The data will be written in a separate file, whose name
       * consists of the stem @p filename and an attached identifier
       * <tt>_fixed.data</tt> for fixed size data and <tt>_variable.data</tt>
       * for variable size data.
       *
       * All processors write into these files simultaneously via MPIIO.
       * Each processor's position to write to will be determined
       * from the provided input parameters.
       *
       * Data has to be previously packed with pack_data().
       */
      void
      save(const unsigned int global_first_cell,
           const unsigned int global_num_cells,
           const std::string &filename) const;

      /**
       * Transfer data from file system.
       *
       * The data will be read from separate file, whose name
       * consists of the stem @p filename and an attached identifier
       * <tt>_fixed.data</tt> for fixed size data and <tt>_variable.data</tt>
       * for variable size data.
       * The @p n_attached_deserialize_fixed and @p n_attached_deserialize_variable
       * parameters are required to gather the memory offsets for each
       * callback.
       *
       * All processors read from these files simultaneously via MPIIO.
       * Each processor's position to read from will be determined
       * from the provided input arguments.
       *
       * After loading, unpack_data() needs to be called to finally
       * distribute data across the associated triangulation.
       */
      void
      load(const unsigned int global_first_cell,
           const unsigned int global_num_cells,
           const unsigned int local_num_cells,
           const std::string &filename,
           const unsigned int n_attached_deserialize_fixed,
           const unsigned int n_attached_deserialize_variable);

      /**
       * Clears all containers and associated data, and resets member
       * values to their default state.
       *
       * Frees memory completely.
       */
      void
      clear();

      /**
       * Flag that denotes if variable size data has been packed.
       */
      bool variable_size_data_stored;

      /**
       * Cumulative size in bytes that those functions that have called
       * register_data_attach() want to attach to each cell. This number
       * only pertains to fixed-sized buffers where the data attached to
       * each cell has exactly the same size.
       *
       * The last entry of this container corresponds to the data size
       * packed per cell in the fixed size buffer (which can be accessed
       * calling <tt>sizes_fixed_cumulative.back()</tt>).
       */
      std::vector<unsigned int> sizes_fixed_cumulative;

      /**
       * Consecutive buffers designed for the fixed size transfer
       * functions of p4est.
       */
      std::vector<char> src_data_fixed;
      std::vector<char> dest_data_fixed;

      /**
       * Consecutive buffers designed for the variable size transfer
       * functions of p4est.
       */
      std::vector<int>  src_sizes_variable;
      std::vector<int>  dest_sizes_variable;
      std::vector<char> src_data_variable;
      std::vector<char> dest_data_variable;

    private:
      MPI_Comm mpi_communicator;
    };

    DataTransfer data_transfer;
  };

} // namespace parallel

DEAL_II_NAMESPACE_CLOSE

#endif
