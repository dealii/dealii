// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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

#ifndef __deal2__distributed_tria_h
#define __deal2__distributed_tria_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/grid/tria.h>

#include <deal.II/base/std_cxx11/function.h>
#include <deal.II/base/std_cxx11/tuple.h>

#include <set>
#include <vector>
#include <list>
#include <utility>

#ifdef DEAL_II_WITH_MPI
#  include <mpi.h>
#endif

#ifdef DEAL_II_WITH_P4EST
#include <p4est_connectivity.h>
#include <p4est.h>
#include <p4est_ghost.h>

#include <p8est_connectivity.h>
#include <p8est.h>
#include <p8est_ghost.h>
#endif


DEAL_II_NAMESPACE_OPEN

template <int, int> class Triangulation;

#ifdef DEAL_II_WITH_P4EST

namespace internal
{
  namespace DoFHandler
  {
    namespace Policy
    {
      template <int, int> class ParallelDistributed;
    }
  }
}


namespace internal
{
  namespace p4est
  {
    /**
     * A structure whose explicit specializations contain typedefs to
     * the relevant p4est_* and p8est_* types. Using this structure,
     * for example by saying <code>types@<dim@>::connectivity</code>
     * we can write code in a dimension independent way, either
     * referring to p4est_connectivity_t or p8est_connectivity_t,
     * depending on template argument.
     */
    template <int> struct types;

    template <>
    struct types<2>
    {
      typedef p4est_connectivity_t connectivity;
      typedef p4est_t              forest;
      typedef p4est_tree_t         tree;
      typedef p4est_quadrant_t     quadrant;
      typedef p4est_topidx_t       topidx;
      typedef p4est_locidx_t       locidx;
#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
      typedef p4est_connect_type_t balance_type;
#else
      typedef p4est_balance_type_t balance_type;
#endif
      typedef p4est_ghost_t        ghost;
    };

    template <>
    struct types<3>
    {
      typedef p8est_connectivity_t connectivity;
      typedef p8est_t              forest;
      typedef p8est_tree_t         tree;
      typedef p8est_quadrant_t     quadrant;
      typedef p4est_topidx_t       topidx;
      typedef p4est_locidx_t       locidx;
#if DEAL_II_P4EST_VERSION_GTE(0,3,4,3)
      typedef p8est_connect_type_t balance_type;
#else
      typedef p8est_balance_type_t balance_type;
#endif
      typedef p8est_ghost_t        ghost;
    };


    /**
     * Initialize the GeometryInfo<dim>::max_children_per_cell
     * children of the cell p4est_cell.
     */
    template <int dim>
    void
    init_quadrant_children
    (const typename types<dim>::quadrant &p4est_cell,
     typename types<dim>::quadrant (&p4est_children)[GeometryInfo<dim>::max_children_per_cell]);


    /**
     * Initialize quadrant to represent a coarse cell.
     */
    template <int dim>
    void
    init_coarse_quadrant(typename types<dim>::quadrant &quad);



    /**
     * Returns whether q1 and q2 are equal
     */
    template <int dim>
    bool
    quadrant_is_equal (const typename types<dim>::quadrant &q1,
                       const typename types<dim>::quadrant &q2);

    //TODO: remove these functions from
    //public interface somehow? [TH]

    /**
     * returns whether q1 is an ancestor of q2
     */
    template <int dim>
    bool
    quadrant_is_ancestor (const typename types<dim>::quadrant &q1,
                          const typename types<dim>::quadrant &q2);
  }
}

//forward declaration of the data type for periodic face pairs
namespace GridTools
{
  template <typename CellIterator> struct PeriodicFacePair;
}

namespace parallel
{
  namespace distributed
  {


    /**
     * This class acts like the dealii::Triangulation class, but it
     * distributes the mesh across a number of different processors when
     * using MPI. The class's interface does not add a lot to the
     * dealii::Triangulation class but there are a number of difficult
     * algorithms under the hood that ensure we always have a
     * load-balanced, fully distributed mesh. Use of this class is
     * explained in step-40, step-32, the @ref distributed documentation
     * module, as well as the @ref distributed_paper . See there for more
     * information. This class satisfies the requirements outlined in
     * @ref GlossMeshAsAContainer "Meshes as containers".
     *
     * @note This class does not support anisotropic refinement, because
     * it relies on the p4est library that does not support this. Attempts
     * to refine cells anisotropically will result in errors.
     * @note There is currently no support for distributing 1d triangulations.
     *
     *
     * <h3> Interaction with boundary description </h3>
     *
     * Refining and coarsening a distributed triangulation is a complicated
     * process because cells may have to be migrated from one processor to
     * another. On a single processor, materializing that part of the global
     * mesh that we want to store here from what we have stored before therefore
     * may involve several cycles of refining and coarsening the locally stored
     * set of cells until we have finally gotten from the previous to the next
     * triangulation. (This process is described in more detail in the
     * @ref distributed_paper.) Unfortunately, in this process, some information
     * can get lost relating to flags that are set by user code and that are
     * inherited from mother to child cell but that are not moved along with
     * a cell if that cell is migrated from one processor to another.
     *
     * An example are boundary indicators. Assume, for example, that you start
     * with a single cell that is refined once globally, yielding four children.
     * If you have four processors, each one owns one cell. Assume now that processor
     * 1 sets the boundary indicators of the external boundaries of the cell it owns
     * to 42. Since processor 0 does not own this cell, it doesn't set the boundary
     * indicators of its ghost cell copy of this cell. Now, assume we do several mesh
     * refinement cycles and end up with a configuration where this processor suddenly finds itself
     * as the owner of this cell. If boundary indicator 42 means that we need to
     * integrate Neumann boundary conditions along this boundary, then processor 0
     * will forget to do so because it has never set the boundary indicator along
     * this cell's boundary to 42.
     *
     * The way to avoid this dilemma is to make sure that things like setting
     * boundary indicators or material ids is done immediately
     * every time a parallel triangulation is refined. This is not necessary
     * for sequential triangulations because, there, these flags are inherited
     * from mother to child cell and remain with a cell even if it is refined
     * and the children are later coarsened again, but this does not hold for
     * distributed triangulations. It is made even more difficult by the fact
     * that in the process of refining a parallel distributed triangulation,
     * the triangulation may call dealii::Triangulation::execute_coarsening_and_refinement
     * multiple times and this function needs to know about boundaries. In
     * other words, it is <i>not</i> enough to just set boundary indicators on
     * newly created faces only <i>after</i> calling
     * distributed::parallel::Triangulation::execute_coarsening_and_refinement:
     * it actually has to happen while that function is still running.
     *
     * The way to do this is by writing a function that sets boundary
     * indicators and that will be called by the dealii::Triangulation class. The
     * triangulation does not provide a pointer to itself to the function being
     * called, nor any other information, so the trick is to get this information
     * into the function. C++ provides a nice mechanism for this that is best
     * explained using an example:
     * @code
     *     #include <deal.II/base/std_cxx11/bind.h>
     *
     *     template <int dim>
     *     void set_boundary_indicators (parallel::distributed::Triangulation<dim> &triangulation)
     *     {
     *       ... set boundary indicators on the triangulation object ...
     *     }
     *
     *     template <int dim>
     *     void
     *     MyClass<dim>::
     *     create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
     *     {
     *       ... create the coarse mesh ...
     *
     *       coarse_grid.signals.post_refinement.connect
     *         (std_cxx11::bind (&set_boundary_indicators<dim>,
     *                           std_cxx11::ref(coarse_grid)));
     *
     *     }
     * @endcode
     *
     * What the call to <code>std_cxx11::bind</code> does is to produce an object that
     * can be called like a function with no arguments. It does so by taking the
     * address of a function that does, in fact, take an argument but permanently fix
     * this one argument to a reference to the coarse grid triangulation. After each
     * refinement step, the triangulation will then call the object so created which
     * will in turn call <code>set_boundary_indicators<dim></code> with the reference
     * to the coarse grid as argument.
     *
     * This approach can be generalized. In the example above, we have used a global
     * function that will be called. However, sometimes it is necessary that this
     * function is in fact a member function of the class that generates the mesh,
     * for example because it needs to access run-time parameters. This can be
     * achieved as follows: assuming the <code>set_boundary_indicators()</code>
     * function has been declared as a (non-static, but possibly private) member
     * function of the <code>MyClass</code> class, then the following will work:
     * @code
     *     #include <deal.II/base/std_cxx11/bind.h>
     *
     *     template <int dim>
     *     void
     *     MyClass<dim>::
     *     set_boundary_indicators (parallel::distributed::Triangulation<dim> &triangulation) const
     *     {
     *       ... set boundary indicators on the triangulation object ...
     *     }
     *
     *     template <int dim>
     *     void
     *     MyClass<dim>::
     *     create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
     *     {
     *       ... create the coarse mesh ...
     *
     *       coarse_grid.signals.post_refinement.connect
     *         (std_cxx11::bind (&MyGeometry<dim>::set_boundary_indicators,
     *                           std_cxx11::cref(*this),
     *                           std_cxx11::ref(coarse_grid)));
     *     }
     * @endcode
     * Here, like any other member function, <code>set_boundary_indicators</code>
     * implicitly takes a pointer or reference to the object it belongs to as first
     * argument. <code>std::bind</code> again creates an object that can be called like a
     * global function with no arguments, and this object in turn calls
     * <code>set_boundary_indicators</code> with a pointer to the current object and a
     * reference to the triangulation to work on. Note that because the
     * <code>create_coarse_mesh</code> function is declared as <code>const</code>, it is
     * necessary that the <code>set_boundary_indicators</code> function is also
     * declared <code>const</code>.
     *
     * <b>Note:</b>For reasons that have to do with the way the
     *   parallel::distributed::Triangulation is implemented, functions that
     *   have been attached to the post-refinement signal of the triangulation are
     *   called more than once, sometimes several times, every time the triangulation
     *  is actually refined.
     *
     *
     * @author Wolfgang Bangerth, Timo Heister 2008, 2009, 2010, 2011
     * @ingroup distributed
     */
    template <int dim, int spacedim = dim>
    class Triangulation : public dealii::Triangulation<dim,spacedim>
    {
    public:
      /**
       * A typedef that is used to to identify cell iterators. The
       * concept of iterators is discussed at length in the
       * @ref Iterators "iterators documentation module".
       *
       * The current typedef identifies cells in a triangulation. You
       * can find the exact type it refers to in the base class's own
       * typedef, but it should be TriaIterator<CellAccessor<dim,spacedim> >. The
       * TriaIterator class works like a pointer that when you
       * dereference it yields an object of type CellAccessor.
       * CellAccessor is a class that identifies properties that
       * are specific to cells in a triangulation, but it is derived
       * (and consequently inherits) from TriaAccessor that describes
       * what you can ask of more general objects (lines, faces, as
       * well as cells) in a triangulation.
       *
       * @ingroup Iterators
       */
      typedef typename dealii::Triangulation<dim,spacedim>::cell_iterator        cell_iterator;

      /**
       * A typedef that is used to to identify
       * @ref GlossActive "active cell iterators". The
       * concept of iterators is discussed at length in the
       * @ref Iterators "iterators documentation module".
       *
       * The current typedef identifies active cells in a triangulation. You
       * can find the exact type it refers to in the base class's own
       * typedef, but it should be TriaActiveIterator<CellAccessor<dim,spacedim> >. The
       * TriaActiveIterator class works like a pointer to active objects that when you
       * dereference it yields an object of type CellAccessor.
       * CellAccessor is a class that identifies properties that
       * are specific to cells in a triangulation, but it is derived
       * (and consequently inherits) from TriaAccessor that describes
       * what you can ask of more general objects (lines, faces, as
       * well as cells) in a triangulation.
       *
       * @ingroup Iterators
       */
      typedef typename dealii::Triangulation<dim,spacedim>::active_cell_iterator active_cell_iterator;

      /**
       * Generic settings for distributed Triangulations. If
       * mesh_reconstruction_after_repartitioning is set, the deal.II
       * mesh will be reconstructed from the coarse mesh every time a
       * repartioning in p4est happens. This can be a bit more
       * expensive, but guarantees the same memory layout and
       * therefore cell ordering in the deal.II mesh. As assembly is
       * done in the deal.II cell ordering, this flag is required to
       * get reproducible behaviour after snapshot/resume.
       *
       * The flag construct_multigrid_hierarchy needs to be set to use
       * the geometric multigrid functionality. This option requires
       * additional computation and communication. Note: geometric
       * multigrid is still a work in progress.
       */
      enum Settings
      {
        default_setting = 0x0,
        mesh_reconstruction_after_repartitioning = 0x1,
        construct_multigrid_hierarchy = 0x2
      };



      /**
       * Constructor.
       *
       * @param mpi_communicator denotes the MPI communicator to be
       * used for the triangulation.
       *
       * @param smooth_grid Degree and kind of mesh smoothing to be
       * applied to the mesh. See the dealii::Triangulation class for
       * a description of the kinds of smoothing operations that can
       * be applied.
       *
       * @param settings See the description of the Settings
       * enumerator.
       *
       * @note This class does not currently support the
       * <code>check_for_distorted_cells</code> argument provided by
       * the base class.
       *
       * @note While it is possible to pass all of the mesh smoothing
       * flags listed in the base class to objects of this type, it is
       * not always possible to honor all of these smoothing options
       * if they would require knowledge of refinement/coarsening
       * flags on cells not locally owned by this processor. As a
       * consequence, for some of these flags, the ultimate number of
       * cells of the parallel triangulation may depend on the number
       * of processors into which it is partitioned. On the other
       * hand, if no smoothing flags are passed, if you always mark
       * the same cells of the mesh, you will always get the exact
       * same refined mesh independent of the number of processors
       * into which the triangulation is partitioned.
       */
      Triangulation (MPI_Comm mpi_communicator,
                     const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing
                     smooth_grid = (dealii::Triangulation<dim,spacedim>::none),
                     const Settings settings = default_setting);

      /**
       * Destructor.
       */
      virtual ~Triangulation ();

      /**
       * Reset this triangulation into a virgin state by deleting all
       * data.
       *
       * Note that this operation is only allowed if no subscriptions
       * to this object exist any more, such as DoFHandler objects
       * using it.
       */
      virtual void clear ();

      /**
       * Implementation of the same function as in the base class.
       */
      virtual void copy_triangulation (const dealii::Triangulation<dim, spacedim> &old_tria);

      /**
       * Create a triangulation as documented in the base class.
       *
       * This function also sets up the various data structures
       * necessary to distribute a mesh across a number of
       * processors. This will be necessary once the mesh is being
       * refined, though we will always keep the entire coarse mesh
       * that is generated by this function on all processors.
       */
      virtual void create_triangulation (const std::vector<Point<spacedim> >    &vertices,
                                         const std::vector<CellData<dim> > &cells,
                                         const SubCellData                 &subcelldata);

      /**
       * Coarsen and refine the mesh according to refinement and
       * coarsening flags set.
       *
       * Since the current processor only has control over those cells
       * it owns (i.e. the ones for which <code>cell-@>subdomain_id() ==
       * this-@>locally_owned_subdomain()</code>), refinement and
       * coarsening flags are only respected for those locally owned
       * cells. Flags may be set on other cells as well (and may
       * often, in fact, if you call
       * dealii::Triangulation::prepare_coarsening_and_refinement) but
       * will be largely ignored: the decision to refine the global
       * mesh will only be affected by flags set on locally owned
       * cells.
       */
      virtual void execute_coarsening_and_refinement ();

      /**
       * Return the subdomain id of those cells that are owned by the
       * current processor. All cells in the triangulation that do not
       * have this subdomain id are either owned by another processor
       * or have children that only exist on other processors.
       */
      types::subdomain_id locally_owned_subdomain () const;

      /**
       * Return the number of active cells in the triangulation that
       * are locally owned, i.e. that have a subdomain_id equal to
       * locally_owned_subdomain(). Note that there may be more active
       * cells in the triangulation stored on the present processor,
       * such as for example ghost cells, or cells further away from
       * the locally owned block of cells but that are needed to
       * ensure that the triangulation that stores this processor's
       * set of active cells still remains balanced with respect to
       * the 2:1 size ratio of adjacent cells.
       *
       * As a consequence of the remark above, the result of this
       * function is always smaller or equal to the result of the
       * function with the same name in the ::Triangulation base
       * class, which includes the active ghost and artificial cells
       * (see also @ref GlossArtificialCell and @ref GlossGhostCell).
       */
      unsigned int n_locally_owned_active_cells () const;

      /**
       * Return the sum over all processors of the number of active
       * cells owned by each processor. This equals the overall number
       * of active cells in the distributed triangulation.
       */
      virtual types::global_dof_index n_global_active_cells () const;

      /**
       * Returns the global maximum level. This may be bigger than
       * the number dealii::Triangulation::n_levels() (a function in this
       * class's base class) returns if the current processor only stores
       * cells in parts of the domain that are not very refined, but
       * if other processors store cells in more deeply refined parts of
       * the domain.
       */
      virtual unsigned int n_global_levels () const;

      /**
       * Returns true if the triangulation has hanging nodes.
       *
       * In the context of parallel distributed triangulations, every
       * processor stores only that part of the triangulation it
       * locally owns. However, it also stores the entire coarse
       * mesh, and to guarantee the 2:1 relationship between cells,
       * this may mean that there are hanging nodes between cells that
       * are not locally owned or ghost cells (i.e., between ghost cells
       * and artificial cells, or between artificial and artificial cells;
       * see @ref GlossArtificialCell "the glossary").
       * One is not typically interested in this case, so the function
       * returns whether there are hanging nodes between any two cells
       * of the "global" mesh, i.e., the union of locally owned cells
       * on all processors.
       */
      virtual
      bool has_hanging_nodes() const;

      /**
       * Return the number of active cells owned by each of the MPI
       * processes that contribute to this triangulation. The element
       * of this vector indexed by locally_owned_subdomain() equals
       * the result of n_locally_owned_active_cells().
       */
      const std::vector<unsigned int> &
      n_locally_owned_active_cells_per_processor () const;

      /**
       * Return the MPI communicator used by this triangulation.
       */
      MPI_Comm get_communicator () const;

      /**
       * Return the local memory consumption in bytes.
       */
      virtual std::size_t memory_consumption () const;

      /**
       * Return the local memory consumption contained in the p4est
       * data structures alone. This is already contained in
       * memory_consumption() but made available separately for
       * debugging purposes.
       */
      virtual std::size_t memory_consumption_p4est () const;

      /**
       * A collective operation that produces a sequence of output
       * files with the given file base name that contain the mesh in
       * VTK format.
       *
       * More than anything else, this function is useful for
       * debugging the interface between deal.II and p4est.
       */
      void write_mesh_vtk (const char *file_basename) const;

      /**
       * Produce a check sum of the triangulation.  This is a
       * collective operation and is mostly useful for debugging
       * purposes.
       */
      unsigned int get_checksum () const;

      /**
       * Save the refinement information from the coarse mesh into the given
       * file. This file needs to be reachable from all nodes in the computation
       * on a shared network file system. See the SolutionTransfer class
       * on how to store solution vectors into this file. Additional cell-based data can be saved
       * using register_data_attach().
       */
      void save(const char *filename) const;

      /**
       * Load the refinement information saved with save() back in. The
       * mesh must contain the same coarse mesh that was used in save()
       * before calling this function.
       *
       * You do not need to load with the same number of MPI processes that
       * you saved with. Rather, if a mesh is loaded with a different
       * number of MPI processes than used at the time of saving, the mesh
       * is repartitioned appropriately. Cell-based data that was saved
       * with register_data_attach() can be read in with
       * notify_ready_to_unpack() after calling load().
       *
       * If you use p4est version > 0.3.4.2 the @p autopartition flag tells
       * p4est to ignore the partitioning that the triangulation had when
       * it was saved and make it uniform upon loading. If @p autopartition
       * is set to false, the triangulation is only repartitioned if needed
       * (i.e. if a different number of MPI processes is encountered).
       */
      void load(const char *filename,
                const bool autopartition = true);

      /**
       * Used to inform in the callbacks of register_data_attach() and
       * notify_ready_to_unpack() how the cell with the given
       * cell_iterator is going to change.  Note that this may me
       * different than the refine_flag() and coarsen_flag() in the
       * cell_iterator because of refinement constraints that this
       * machine does not see.
       */
      enum CellStatus
      {
        /**
         * The cell will not be refined or coarsened and might or might
         * not move to a different processor.
         */
        CELL_PERSIST,
        /**
         * The cell will be or was refined.
         */
        CELL_REFINE,
        /**
         * The children of this cell will be or were coarsened into this cell.
         */
        CELL_COARSEN,
        /**
         * Invalid status. Will not occur for the user.
         */
        CELL_INVALID
      };

      /**
       * Register a function with the current Triangulation object
       * that will be used to attach data to active cells before
       * execute_coarsening_and_refinement(). In
       * execute_coarsening_and_refinement() the Triangulation will
       * call the given function pointer and provide @p size bytes to
       * store data. If necessary, this data will be transferred to
       * the new owner of that cell during repartitioning the
       * tree. See notify_ready_to_unpack() on how to retrieve the
       * data.
       *
       * Callers need to store the return value.  It specifies an
       * offset of the position at which data can later be retrieved
       * during a call to notify_ready_to_unpack().
       *
       * The CellStatus argument in the callback function will tell you if the
       * given cell will be coarsened, refined, or will persist as is (this
       * can be different than the coarsen and refine flags set by you). If it
       * is
       *
       * - CELL_PERIST: the cell won't be refined/coarsened, but might be
       *   moved to a different processor
       * - CELL_REFINE: this cell will be refined into 4/8 cells, you can not
       *   access the children (because they don't exist yet)
       * - CELL_COARSEN: the children of this cell will be coarsened into the
       *   given cell (you can access the active children!)
       *
       * When unpacking the data with notify_ready_to_unpack() you can access
       * the children of the cell if the status is CELL_REFINE but not for
       * CELL_COARSEN. As a consequence you need to handle coarsening while
       * packing and refinement during unpacking.
       *
       * @note The two functions can also be used for serialization of data
       * using save() and load() in the same way. Then the status will always
       * be CELL_PERSIST.
       */
      unsigned int
      register_data_attach (const std::size_t size,
                            const std_cxx11::function<void (const cell_iterator &,
                                                            const CellStatus,
                                                            void *)> &pack_callback);

      /**
       * The supplied callback function is called for each newly locally owned
       * cell and corresponding data saved with register_data_attach().  This
       * function needs to be called after execute_coarsening_and_refinement()
       * with the offset returned by register_data_attach().
       *
       * The CellStatus will indicate if the cell was refined, coarsened, or
       * persisted unchanged. The cell_iterator will either by an active,
       * locally owned cell (if the cell was not refined), or the immediate
       * parent if it was refined during
       * execute_coarsening_and_refinement(). Therefore, contrary to during
       * register_data_attach(), you can now access the children if the status
       * is CELL_REFINE but no longer for callbacks with status CELL_COARSEN.
       */
      void
      notify_ready_to_unpack (const unsigned int offset,
                              const std_cxx11::function<void (const cell_iterator &,
                                                              const CellStatus,
                                                              const void *)> &unpack_callback);

      /**
       * Returns a permutation vector for the order the coarse cells
       * are handed of to p4est. For example the first element i in
       * this vector denotes that the first cell in hierarchical
       * ordering is the ith deal cell starting from begin(0).
       */
      const std::vector<types::global_dof_index> &
      get_p4est_tree_to_coarse_cell_permutation() const;



      /**
       * Join faces in the p4est forest for periodic boundary
       * conditions. As a result, each pair of faces will differ by at
       * most one refinement level and ghost neighbors will be
       * available across these faces.
       *
       * The vector can be filled by the function
       * GridTools::collect_periodic_faces.
       *
       * @todo At the moment just default orientation is implemented.
       *
       * @note Before this function can be used the Triangulation has to be
       * initialized and must not be refined.
       * Calling this function more than once is possible, but not recommended:
       * The function destroys and rebuilds the p4est forest each time it is
       * called.
       */
      void
      add_periodicity
      (const std::vector<GridTools::PeriodicFacePair<cell_iterator> > &);



    private:
      /**
       * MPI communicator to be used for the triangulation. We create
       * a unique communicator for this class, which is a duplicate of
       * the one passed to the constructor.
       */
      MPI_Comm mpi_communicator;

      /**
       * store the Settings.
       */
      Settings settings;

      /**
       * The subdomain id to be used for the current processor.
       */
      types::subdomain_id my_subdomain;

      /**
       * A flag that indicates whether the triangulation has actual
       * content.
       */
      bool triangulation_has_content;

      /**
       * A structure that contains some numbers about the distributed
       * triangulation.
       */
      struct NumberCache
      {
        std::vector<unsigned int> n_locally_owned_active_cells;
        types::global_dof_index   n_global_active_cells;
        unsigned int              n_global_levels;

        NumberCache();
      };

      NumberCache number_cache;

      /**
       * A data structure that holds the connectivity between
       * trees. Since each tree is rooted in a coarse grid cell, this
       * data structure holds the connectivity between the cells of
       * the coarse grid.
       */
      typename dealii::internal::p4est::types<dim>::connectivity *connectivity;

      /**
       * A data structure that holds the local part of the global
       * triangulation.
       */
      typename dealii::internal::p4est::types<dim>::forest *parallel_forest;
      /**
       * A data structure that holds some information about the ghost
       * cells of the triangulation.
       */
      typename dealii::internal::p4est::types<dim>::ghost  *parallel_ghost;

      /**
       * A flag that indicates whether refinement of a triangulation
       * is currently in progress. This flag is used to disambiguate
       * whether a call to execute_coarsening_and_triangulation came
       * from the outside or through a recursive call. While the first
       * time we want to take over work to copy things from a refined
       * p4est, the other times we don't want to get in the way as
       * these latter calls to
       * Triangulation::execute_coarsening_and_refinement() are simply
       * there in order to re-create a triangulation that matches the
       * p4est.
       */
      bool refinement_in_progress;


      /**
       * number of bytes that get attached to the Triangulation
       * through register_data_attach() for example SolutionTransfer.
       */
      unsigned int attached_data_size;

      /**
       * number of functions that get attached to the Triangulation
       * through register_data_attach() for example SolutionTransfer.
       */
      unsigned int n_attached_datas;

      /**
       * number of functions that need to unpack their data after a
       * call from load()
       */
      unsigned int n_attached_deserialize;

      typedef  std_cxx11::function<
      void(typename Triangulation<dim,spacedim>::cell_iterator, CellStatus, void *)
      > pack_callback_t;

      typedef std::pair<unsigned int, pack_callback_t> callback_pair_t;

      typedef std::list<callback_pair_t> callback_list_t;

      /**
       * List of callback functions registered by
       * register_data_attach() that are going to be called for
       * packing data.
       */
      callback_list_t attached_data_pack_callbacks;


      /**
       * Two arrays that store which p4est tree corresponds to which
       * coarse grid cell and vice versa. We need these arrays because
       * p4est goes with the original order of coarse cells when it
       * sets up its forest, and then applies the Morton ordering
       * within each tree. But if coarse grid cells are badly ordered
       * this may mean that individual parts of the forest stored on a
       * local machine may be split across coarse grid cells that are
       * not geometrically close. Consequently, we apply a
       * Cuthill-McKee preordering to ensure that the part of the
       * forest stored by p4est is located on geometrically close
       * coarse grid cells.
       */
      std::vector<types::global_dof_index> coarse_cell_to_p4est_tree_permutation;
      std::vector<types::global_dof_index> p4est_tree_to_coarse_cell_permutation;

      /**
       * Return a pointer to the p4est tree that belongs to the given
       * dealii_coarse_cell_index()
       */
      typename dealii::internal::p4est::types<dim>::tree *
      init_tree(const int dealii_coarse_cell_index) const;

      /**
       * The function that computes the permutation between the two
       * data storage schemes.
       */
      void setup_coarse_cell_to_p4est_tree_permutation ();

      /**
       * Take the contents of a newly created triangulation we are
       * attached to and copy it to p4est data structures.
       *
       * This function exists in 2d and 3d variants.
       */
      void copy_new_triangulation_to_p4est (dealii::internal::int2type<2>);
      void copy_new_triangulation_to_p4est (dealii::internal::int2type<3>);

      /**
       * Copy the local part of the refined forest from p4est into the
       * attached triangulation.
       */
      void copy_local_forest_to_triangulation ();


      /**
       * Update the number_cache variable after mesh creation or
       * refinement.
       */
      void update_number_cache ();

      /**
       * Internal function notifying all registered classes to attach
       * their data before repartitioning occurs. Called from
       * execute_coarsening_and_refinement().
       */
      void attach_mesh_data();

      /**
       * fills a map that, for each vertex, lists all the processors whose
       * subdomains are adjacent to that vertex.  Used by
       * DoFHandler::Policy::ParallelDistributed.
       */
      void
      fill_vertices_with_ghost_neighbors
      (std::map<unsigned int, std::set<dealii::types::subdomain_id> >
       &vertices_with_ghost_neighbors);

      template <int, int> friend class dealii::internal::DoFHandler::Policy::ParallelDistributed;
    };


    /**
     * Specialization of the general template for the 1d case. There
     * is currently no support for distributing 1d
     * triangulations. Consequently, all this class does is throw an
     * exception.
     */
    template <int spacedim>
    class Triangulation<1,spacedim> : public dealii::Triangulation<1,spacedim>
    {
    public:
      /**
       * Constructor. The argument denotes the MPI communicator to be
       * used for the triangulation.
       */
      Triangulation (MPI_Comm mpi_communicator);

      /**
       * Destructor.
       */
      virtual ~Triangulation ();

      /**
       * Return the MPI communicator used by this triangulation.
       */
      MPI_Comm get_communicator () const;

      /**
       * Return the sum over all processors of the number of active
       * cells owned by each processor. This equals the overall number
       * of active cells in the distributed triangulation.
       */
      types::global_dof_index n_global_active_cells () const;
      virtual unsigned int n_global_levels () const;

      /**
       * Returns a permutation vector for the order the coarse cells
       * are handed of to p4est. For example the first element i in
       * this vector denotes that the first cell in hierarchical
       * ordering is the ith deal cell starting from begin(0).
       */
      const std::vector<types::global_dof_index> &
      get_p4est_tree_to_coarse_cell_permutation() const;

      /**
       * Return the subdomain id of those cells that are owned by the
       * current processor. All cells in the triangulation that do not
       * have this subdomain id are either owned by another processor
       * or have children that only exist on other processors.
       */
      types::subdomain_id locally_owned_subdomain () const;

      /**
       * Dummy arrays. This class isn't usable but the compiler wants
       * to see these variables at a couple places anyway.
       */
      std::vector<types::global_dof_index> coarse_cell_to_p4est_tree_permutation;
      std::vector<types::global_dof_index> p4est_tree_to_coarse_cell_permutation;

      /**
       * dummy settings
       */
      enum Settings
      {
        default_setting = 0x0,
        mesh_reconstruction_after_repartitioning = 0x1,
        construct_multigrid_hierarchy = 0x2
      };


//TODO: The following variable should really be private, but it is used in dof_handler_policy.cc ...
      /**
       * dummy settings object
       */
      Settings settings;

      /**
       * Like above, this method, which is only implemented for dim = 2 or 3,
       * needs a stub because it is used in dof_handler_policy.cc
       */
      void
      fill_vertices_with_ghost_neighbors
      (std::map<unsigned int, std::set<dealii::types::subdomain_id> >
       &vertices_with_ghost_neighbors);

    };
  }
}


#else // DEAL_II_WITH_P4EST

namespace parallel
{
  namespace distributed
  {
    /**
     * Dummy class the compiler chooses for parallel distributed
     * triangulations if we didn't actually configure deal.II with the
     * p4est library. The existence of this class allows us to refer
     * to parallel::distributed::Triangulation objects throughout the
     * library even if it is disabled.
     *
     * Since the constructor of this class is private, no such objects
     * can actually be created if we don't have p4est available.
     */
    template <int dim, int spacedim = dim>
    class Triangulation : public dealii::Triangulation<dim,spacedim>
    {
    private:
      /**
       * Constructor.
       */
      Triangulation ();

    public:

      /**
       * Destructor.
       */
      virtual ~Triangulation ();

      /**
       * Return the subdomain id of those cells that are owned by the
       * current processor. All cells in the triangulation that do not
       * have this subdomain id are either owned by another processor
       * or have children that only exist on other processors.
       */
      types::subdomain_id locally_owned_subdomain () const;

      /**
       * Return the MPI communicator used by this triangulation.
       */
#ifdef DEAL_II_WITH_MPI
      MPI_Comm get_communicator () const;
#endif
    };
  }
}


#endif


DEAL_II_NAMESPACE_CLOSE

#endif
