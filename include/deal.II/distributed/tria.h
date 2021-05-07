// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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

#ifndef dealii_distributed_tria_h
#define dealii_distributed_tria_h


#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/distributed/p4est_wrappers.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/tria.h>

#include <boost/range/iterator_range.hpp>

#include <functional>
#include <list>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef DEAL_II_WITH_MPI
#  include <mpi.h>
#endif

#ifdef DEAL_II_WITH_P4EST
#  include <p4est.h>
#  include <p4est_connectivity.h>
#  include <p4est_ghost.h>
#  include <p8est.h>
#  include <p8est_connectivity.h>
#  include <p8est_ghost.h>
#endif


DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_P4EST

// Forward declarations
#  ifndef DOXYGEN

namespace FETools
{
  namespace internal
  {
    template <int, int, class>
    class ExtrapolateImplementation;
  }
} // namespace FETools

// forward declaration of the data type for periodic face pairs
namespace GridTools
{
  template <typename CellIterator>
  struct PeriodicFacePair;
}

namespace parallel
{
  namespace distributed
  {
    template <int, int>
    class TemporarilyMatchRefineFlags;
  }
} // namespace parallel
#  endif

namespace parallel
{
  namespace distributed
  {
    /**
     * This class acts like the dealii::Triangulation class, but it
     * distributes the mesh across a number of different processors when using
     * MPI. The class's interface does not add a lot to the
     * dealii::Triangulation class but there are a number of difficult
     * algorithms under the hood that ensure we always have a load-balanced,
     * fully distributed mesh. Use of this class is explained in step-40,
     * step-32, the
     * @ref distributed
     * documentation module, as well as the
     * @ref distributed_paper.
     * See there for more information. This class satisfies the
     * @ref ConceptMeshType "MeshType concept".
     *
     * @note This class does not support anisotropic refinement, because it
     * relies on the p4est library that does not support this. Attempts to
     * refine cells anisotropically will result in errors.
     *
     * @note There is currently no support for distributing 1d triangulations.
     *
     *
     * <h3> Interaction with boundary description </h3>
     *
     * Refining and coarsening a distributed triangulation is a complicated
     * process because cells may have to be migrated from one processor to
     * another. On a single processor, materializing that part of the global
     * mesh that we want to store here from what we have stored before
     * therefore may involve several cycles of refining and coarsening the
     * locally stored set of cells until we have finally gotten from the
     * previous to the next triangulation. This process is described in more
     * detail in the
     * @ref distributed_paper.
     * Unfortunately, in this process, some information can get lost relating
     * to flags that are set by user code and that are inherited from mother
     * to child cell but that are not moved along with a cell if that cell is
     * migrated from one processor to another.
     *
     * An example are boundary indicators. Assume, for example, that you start
     * with a single cell that is refined once globally, yielding four
     * children. If you have four processors, each one owns one cell. Assume
     * now that processor 1 sets the boundary indicators of the external
     * boundaries of the cell it owns to 42. Since processor 0 does not own
     * this cell, it doesn't set the boundary indicators of its ghost cell
     * copy of this cell. Now, assume we do several mesh refinement cycles and
     * end up with a configuration where this processor suddenly finds itself
     * as the owner of this cell. If boundary indicator 42 means that we need
     * to integrate Neumann boundary conditions along this boundary, then
     * processor 0 will forget to do so because it has never set the boundary
     * indicator along this cell's boundary to 42.
     *
     * The way to avoid this dilemma is to make sure that things like setting
     * boundary indicators or material ids is done immediately every time a
     * parallel triangulation is refined. This is not necessary for sequential
     * triangulations because, there, these flags are inherited from mother to
     * child cell and remain with a cell even if it is refined and the
     * children are later coarsened again, but this does not hold for
     * distributed triangulations. It is made even more difficult by the fact
     * that in the process of refining a parallel distributed triangulation,
     * the triangulation may call
     * dealii::Triangulation::execute_coarsening_and_refinement multiple times
     * and this function needs to know about boundaries. In other words, it is
     * <i>not</i> enough to just set boundary indicators on newly created
     * faces only <i>after</i> calling
     * <tt>distributed::parallel::TriangulationBase::execute_coarsening_and_refinement</tt>:
     * it actually has to happen while that function is still running.
     *
     * The way to do this is by writing a function that sets boundary
     * indicators and that will be called by the dealii::Triangulation class.
     * The triangulation does not provide a pointer to itself to the function
     * being called, nor any other information, so the trick is to get this
     * information into the function. C++ provides a nice mechanism for this
     * that is best explained using an example:
     * @code
     * #include <functional>
     *
     * template <int dim>
     * void set_boundary_ids (
     *   parallel::distributed::Triangulation<dim> &triangulation)
     * {
     *   ... set boundary indicators on the triangulation object ...
     * }
     *
     * template <int dim>
     * void
     * MyClass<dim>::create_coarse_mesh (
     *   parallel::distributed::Triangulation<dim> &coarse_grid) const
     * {
     *   ... create the coarse mesh ...
     *
     *   coarse_grid.signals.post_refinement.connect(
     *     [&coarse_grid](){
     *       set_boundary_ids<dim>(coarse_grid);
     *     });
     * }
     * @endcode
     *
     * The object passed as argument to <code>connect</code> is an object
     * that can be called like a function with no arguments. It does so by
     * wrapping a function that does, in fact, take an argument but this one
     * argument is stored as a reference to the coarse grid triangulation when
     * the lambda function is created. After each refinement step, the
     * triangulation will then call the object so created which will in turn
     * call <code>set_boundary_ids<dim></code> with the reference to the coarse
     * grid as argument.
     *
     * This approach can be generalized. In the example above, we have used a
     * global function that will be called. However, sometimes it is necessary
     * that this function is in fact a member function of the class that
     * generates the mesh, for example because it needs to access run-time
     * parameters. This can be achieved as follows: assuming the
     * <code>set_boundary_ids()</code> function has been declared as a (non-
     * static, but possibly private) member function of the
     * <code>MyClass</code> class, then the following will work:
     * @code
     * #include <functional>
     *
     * template <int dim>
     * void
     * MyClass<dim>::set_boundary_ids (
     *   parallel::distributed::Triangulation<dim> &triangulation) const
     * {
     *   ... set boundary indicators on the triangulation object ...
     * }
     *
     * template <int dim>
     * void
     * MyClass<dim>::create_coarse_mesh (
     *   parallel::distributed::Triangulation<dim> &coarse_grid) const
     * {
     *   ... create the coarse mesh ...
     *
     *   coarse_grid.signals.post_refinement.connect(
     *     [this, &coarse_grid]()
     *     {
     *       this->set_boundary_ids(coarse_grid);
     *     });
     * }
     * @endcode
     * The lambda function above again is an object that can
     * be called like a global function with no arguments, and this object in
     * turn calls the current object's member function
     * <code>set_boundary_ids</code> with a reference to the triangulation to
     * work on. Note that
     * because the <code>create_coarse_mesh</code> function is declared as
     * <code>const</code>, it is necessary that the
     * <code>set_boundary_ids</code> function is also declared
     * <code>const</code>.
     *
     * <b>Note:</b>For reasons that have to do with the way the
     * parallel::distributed::Triangulation is implemented, functions that
     * have been attached to the post-refinement signal of the triangulation
     * are called more than once, sometimes several times, every time the
     * triangulation is actually refined.
     *
     *
     * @ingroup distributed
     */
    template <int dim, int spacedim = dim>
    class Triangulation
      : public dealii::parallel::DistributedTriangulationBase<dim, spacedim>
    {
    public:
      /**
       * An alias that is used to identify cell iterators. The concept of
       * iterators is discussed at length in the
       * @ref Iterators "iterators documentation module".
       *
       * The current alias identifies cells in a triangulation. You can find
       * the exact type it refers to in the base class's own alias, but it
       * should be TriaIterator<CellAccessor<dim,spacedim> >. The TriaIterator
       * class works like a pointer that when you dereference it yields an
       * object of type CellAccessor. CellAccessor is a class that identifies
       * properties that are specific to cells in a triangulation, but it is
       * derived (and consequently inherits) from TriaAccessor that describes
       * what you can ask of more general objects (lines, faces, as well as
       * cells) in a triangulation.
       *
       * @ingroup Iterators
       */
      using cell_iterator =
        typename dealii::Triangulation<dim, spacedim>::cell_iterator;

      /**
       * An alias that is used to identify
       * @ref GlossActive "active cell iterators".
       * The concept of iterators is discussed at length in the
       * @ref Iterators "iterators documentation module".
       *
       * The current alias identifies active cells in a triangulation. You
       * can find the exact type it refers to in the base class's own alias,
       * but it should be TriaActiveIterator<CellAccessor<dim,spacedim> >. The
       * TriaActiveIterator class works like a pointer to active objects that
       * when you dereference it yields an object of type CellAccessor.
       * CellAccessor is a class that identifies properties that are specific
       * to cells in a triangulation, but it is derived (and consequently
       * inherits) from TriaAccessor that describes what you can ask of more
       * general objects (lines, faces, as well as cells) in a triangulation.
       *
       * @ingroup Iterators
       */
      using active_cell_iterator =
        typename dealii::Triangulation<dim, spacedim>::active_cell_iterator;

      using CellStatus =
        typename dealii::Triangulation<dim, spacedim>::CellStatus;

      /**
       * Configuration flags for distributed Triangulations to be set in the
       * constructor. Settings can be combined using bitwise OR.
       */
      enum Settings
      {
        /**
         * Default settings, other options are disabled.
         */
        default_setting = 0x0,
        /**
         * If set, the deal.II mesh will be reconstructed from the coarse mesh
         * every time a repartitioning in p4est happens. This can be a bit more
         * expensive, but guarantees the same memory layout and therefore cell
         * ordering in the deal.II mesh. As assembly is done in the deal.II
         * cell ordering, this flag is required to get reproducible behavior
         * after snapshot/resume.
         */
        mesh_reconstruction_after_repartitioning = 0x1,
        /**
         * This flags needs to be set to use the geometric multigrid
         * functionality. This option requires additional computation and
         * communication.
         */
        construct_multigrid_hierarchy = 0x2,
        /**
         * Setting this flag will disable automatic repartitioning of the cells
         * after a refinement cycle. It can be executed manually by calling
         * repartition().
         */
        no_automatic_repartitioning = 0x4
      };



      /**
       * Constructor.
       *
       * @param mpi_communicator The MPI communicator to be used for
       * the triangulation.
       *
       * @param smooth_grid Degree and kind of mesh smoothing to be applied to
       * the mesh. See the dealii::Triangulation class for a description of
       * the kinds of smoothing operations that can be applied.
       *
       * @param settings See the description of the Settings enumerator.
       * Providing <code>construct_multigrid_hierarchy</code> enforces
       * <code>Triangulation::limit_level_difference_at_vertices</code>
       * for smooth_grid.
       *
       * @note This class does not currently support the
       * <code>check_for_distorted_cells</code> argument provided by the base
       * class.
       *
       * @note While it is possible to pass all of the mesh smoothing flags
       * listed in the base class to objects of this type, it is not always
       * possible to honor all of these smoothing options if they would
       * require knowledge of refinement/coarsening flags on cells not locally
       * owned by this processor. As a consequence, for some of these flags,
       * the ultimate number of cells of the parallel triangulation may depend
       * on the number of processors into which it is partitioned. On the
       * other hand, if no smoothing flags are passed, if you always mark the
       * same cells of the mesh, you will always get the exact same refined
       * mesh independent of the number of processors into which the
       * triangulation is partitioned.
       */
      explicit Triangulation(
        const MPI_Comm &mpi_communicator,
        const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
                       smooth_grid = (dealii::Triangulation<dim, spacedim>::none),
        const Settings settings    = default_setting);

      /**
       * Destructor.
       */
      virtual ~Triangulation() override;

      /**
       * Reset this triangulation into a virgin state by deleting all data.
       *
       * Note that this operation is only allowed if no subscriptions to this
       * object exist any more, such as DoFHandler objects using it.
       */
      virtual void
      clear() override;

      /**
       * Return if multilevel hierarchy is supported and has been constructed.
       */
      bool
      is_multilevel_hierarchy_constructed() const override;

      /**
       * Transfer data across forests.
       *
       * Besides the actual @p parallel_forest, which has been already refined
       * and repartitioned, this function also needs information about its
       * previous state, i.e. the locally owned intervals in p4est's
       * sc_array of each processor. This information needs to be memcopyied
       * out of the old p4est object and has to be provided via the parameter
       * @p previous_global_first_quadrant.
       *
       * Data has to be previously packed with
       * DistributedTriangulationBase::DataTransfer::pack_data().
       */
      void
      execute_transfer(
        const typename dealii::internal::p4est::types<dim>::forest
          *parallel_forest,
        const typename dealii::internal::p4est::types<dim>::gloidx
          *previous_global_first_quadrant);

      /**
       * Implementation of the same function as in the base class.
       *
       * @note This function can be used to copy a serial Triangulation to a
       * parallel::distributed::Triangulation but only if the serial
       * Triangulation has never been refined.
       */
      virtual void
      copy_triangulation(
        const dealii::Triangulation<dim, spacedim> &other_tria) override;

      /**
       * Create a triangulation as documented in the base class.
       *
       * This function also sets up the various data structures necessary to
       * distribute a mesh across a number of processors. This will be
       * necessary once the mesh is being refined, though we will always keep
       * the entire coarse mesh that is generated by this function on all
       * processors.
       */
      virtual void
      create_triangulation(const std::vector<Point<spacedim>> &vertices,
                           const std::vector<CellData<dim>> &  cells,
                           const SubCellData &subcelldata) override;

      /**
       * @copydoc Triangulation::create_triangulation()
       *
       * @note Not implemented yet.
       */
      virtual void
      create_triangulation(
        const TriangulationDescription::Description<dim, spacedim>
          &construction_data) override;

      /**
       * Coarsen and refine the mesh according to refinement and coarsening
       * flags set.
       *
       * Since the current processor only has control over those cells it owns
       * (i.e. the ones for which <code>cell-@>subdomain_id() ==
       * this-@>locally_owned_subdomain()</code>), refinement and coarsening
       * flags are only respected for those locally owned cells. Flags may be
       * set on other cells as well (and may often, in fact, if you call
       * dealii::Triangulation::prepare_coarsening_and_refinement()) but will
       * be largely ignored: the decision to refine the global mesh will only
       * be affected by flags set on locally owned cells.
       *
       * @note This function by default partitions the mesh in such a way that
       * the number of cells on all processors is roughly equal. If you want
       * to set weights for partitioning, e.g. because some cells are more
       * expensive to compute than others, you can use the signal cell_weight
       * as documented in the dealii::Triangulation class. This function will
       * check whether a function is connected to the signal and if so use it.
       * If you prefer to repartition the mesh yourself at user-defined
       * intervals only, you can create your triangulation object by passing
       * the parallel::distributed::Triangulation::no_automatic_repartitioning
       * flag to the constructor, which ensures that calling the current
       * function only refines and coarsens the triangulation, but doesn't
       * partition it. You can then call the repartition() function manually.
       * The usage of the cell_weights signal is identical in both cases, if a
       * function is connected to the signal it will be used to balance the
       * calculated weights, otherwise the number of cells is balanced.
       */
      virtual void
      execute_coarsening_and_refinement() override;

      /**
       * Override the implementation of prepare_coarsening_and_refinement from
       * the base class. This is necessary if periodic boundaries are enabled
       * and the level difference over vertices over the periodic boundary
       * must not be more than 2:1.
       */
      virtual bool
      prepare_coarsening_and_refinement() override;

      /**
       * Manually repartition the active cells between processors. Normally
       * this repartitioning will happen automatically when calling
       * execute_coarsening_and_refinement() (or refine_global()) unless the
       * @p no_automatic_repartitioning is set in the constructor. Setting the
       * flag and then calling repartition() gives the same result.
       *
       * If you want to transfer data (using SolutionTransfer or manually with
       * register_data_attach() and notify_ready_to_unpack()), you need to set
       * it up twice: once when calling execute_coarsening_and_refinement(),
       * which will handle coarsening and refinement but obviously won't ship
       * any data between processors, and a second time when calling
       * repartition().  Here, no coarsening and refinement will be done but
       * information will be packed and shipped to different processors. In
       * other words, you probably want to treat a call to repartition() in
       * the same way as execute_coarsening_and_refinement() with respect to
       * dealing with data movement (SolutionTransfer, etc.).
       *
       * @note If no function is connected to the cell_weight signal described
       * in the dealii::Triangulation class, this function will balance the
       * number of cells on each processor. If one or more functions are
       * connected, it will calculate the sum of the weights and balance the
       * weights across processors. The only requirement on the weights is
       * that every cell's weight is positive and that the sum over all
       * weights on all processors can be formed using a 64-bit integer.
       * Beyond that, it is your choice how you want to interpret the weights.
       * A common approach is to consider the weights proportional to the cost
       * of doing computations on a cell, e.g., by summing the time for
       * assembly and solving. In practice, determining this cost is of course
       * not trivial since we don't solve on isolated cells, but on the entire
       * mesh. In such cases, one could, for example, choose the weight equal
       * to the number of unknowns per cell (in the context of hp-finite
       * element methods), or using a heuristic that estimates the cost on
       * each cell depending on whether, for example, one has to run some
       * expensive algorithm on some cells but not others (such as forming
       * boundary integrals during the assembly only on cells that are
       * actually at the boundary, or computing expensive nonlinear terms only
       * on some cells but not others, e.g., in the elasto-plastic problem in
       * step-42).
       */
      void
      repartition();


      /**
       * Return true if the triangulation has hanging nodes.
       *
       * In the context of parallel distributed triangulations, every
       * processor stores only that part of the triangulation it locally owns.
       * However, it also stores the entire coarse mesh, and to guarantee the
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
       * Return the local memory consumption in bytes.
       */
      virtual std::size_t
      memory_consumption() const override;

      /**
       * Return the local memory consumption contained in the p4est data
       * structures alone. This is already contained in memory_consumption()
       * but made available separately for debugging purposes.
       */
      virtual std::size_t
      memory_consumption_p4est() const;

      /**
       * A collective operation that produces a sequence of output files with
       * the given file base name that contain the mesh in VTK format.
       *
       * More than anything else, this function is useful for debugging the
       * interface between deal.II and p4est.
       */
      void
      write_mesh_vtk(const std::string &file_basename) const;

      /**
       * Produce a check sum of the triangulation.  This is a collective
       * operation and is mostly useful for debugging purposes.
       */
      unsigned int
      get_checksum() const;

      /**
       * Save the refinement information from the coarse mesh into the given
       * file. This file needs to be reachable from all nodes in the
       * computation on a shared network file system. See the SolutionTransfer
       * class on how to store solution vectors into this file. Additional
       * cell-based data can be saved using
       * DistributedTriangulationBase::DataTransfer::register_data_attach().
       */
      virtual void
      save(const std::string &filename) const override;

      /**
       * Load the refinement information saved with save() back in. The mesh
       * must contain the same coarse mesh that was used in save() before
       * calling this function.
       *
       * You do not need to load with the same number of MPI processes that
       * you saved with. Rather, if a mesh is loaded with a different number
       * of MPI processes than used at the time of saving, the mesh is
       * repartitioned appropriately. Cell-based data that was saved with
       * DistributedTriangulationBase::DataTransfer::register_data_attach() can
       * be read in with
       * DistributedTriangulationBase::DataTransfer::notify_ready_to_unpack()
       * after calling load().
       *
       * If you use p4est version > 0.3.4.2 the @p autopartition flag tells
       * p4est to ignore the partitioning that the triangulation had when it
       * was saved and make it uniform upon loading. If @p autopartition is
       * set to false, the triangulation is only repartitioned if needed (i.e.
       * if a different number of MPI processes is encountered).
       */
      virtual void
      load(const std::string &filename,
           const bool         autopartition = true) override;

      /**
       * Load the refinement information from a given parallel forest. This
       * forest might be obtained from the function call to
       * parallel::distributed::Triangulation::get_p4est().
       */
      void
      load(const typename dealii::internal::p4est::types<dim>::forest *forest);

      /**
       * Return a permutation vector for the order the coarse cells are handed
       * off to p4est. For example the value of the $i$th element in this
       * vector is the index of the deal.II coarse cell (counting from
       * begin(0)) that corresponds to the $i$th tree managed by p4est.
       */
      const std::vector<types::global_dof_index> &
      get_p4est_tree_to_coarse_cell_permutation() const;

      /**
       * Return a permutation vector for the mapping from the coarse deal
       * cells to the p4est trees. This is the inverse of
       * get_p4est_tree_to_coarse_cell_permutation.
       */
      const std::vector<types::global_dof_index> &
      get_coarse_cell_to_p4est_tree_permutation() const;

      /**
       * This returns a pointer to the internally stored p4est object (of type
       * p4est_t or p8est_t depending on @p dim).
       *
       * @warning If you modify the p4est object, internal data structures
       * can become inconsistent.
       */
      const typename dealii::internal::p4est::types<dim>::forest *
      get_p4est() const;

      /**
       * In addition to the action in the base class Triangulation, this
       * function joins faces in the p4est forest for periodic boundary
       * conditions. As a result, each pair of faces will differ by at most one
       * refinement level and ghost neighbors will be available across these
       * faces.
       *
       * The vector can be filled by the function
       * GridTools::collect_periodic_faces.
       *
       * For more information on periodic boundary conditions see
       * GridTools::collect_periodic_faces,
       * DoFTools::make_periodicity_constraints and step-45.
       *
       * @note Before this function can be used the Triangulation has to be
       * initialized and must not be refined. Calling this function more than
       * once is possible, but not recommended: The function destroys and
       * rebuilds the p4est forest each time it is called.
       */
      virtual void
      add_periodicity(
        const std::vector<dealii::GridTools::PeriodicFacePair<cell_iterator>> &)
        override;


    private:
      /**
       * store the Settings.
       */
      Settings settings;

      /**
       * A flag that indicates whether the triangulation has actual content.
       */
      bool triangulation_has_content;

      /**
       * A data structure that holds the connectivity between trees. Since
       * each tree is rooted in a coarse grid cell, this data structure holds
       * the connectivity between the cells of the coarse grid.
       */
      typename dealii::internal::p4est::types<dim>::connectivity *connectivity;

      /**
       * A data structure that holds the local part of the global
       * triangulation.
       */
      typename dealii::internal::p4est::types<dim>::forest *parallel_forest;

      /**
       * A data structure that holds some information about the ghost cells of
       * the triangulation.
       */
      typename dealii::internal::p4est::types<dim>::ghost *parallel_ghost;

      /**
       * Go through all p4est trees and record the relations between locally
       * owned p4est quadrants and active deal.II cells in the private member
       * vector local_cell_relations.
       *
       * The vector contains an active cell iterator for every locally owned
       * p4est quadrant, as well as a CellStatus flag to describe their
       * relation.
       *
       * The stored vector will be ordered by the occurrence of quadrants in
       * the corresponding local sc_array of the parallel_forest. p4est requires
       * this specific ordering for its transfer functions. Therefore, the size
       * of this vector will be equal to the number of locally owned quadrants
       * in the parallel_forest object.
       *
       * These relations will be established for example in the mesh refinement
       * process: after adapting the parallel_forest, but before applying these
       * changes to this triangulation, we will record how cells will change in
       * the refinement process. With this information, we can prepare all
       * buffers for data transfer accordingly.
       */
      virtual void
      update_cell_relations() override;

      /**
       * Two arrays that store which p4est tree corresponds to which coarse
       * grid cell and vice versa. We need these arrays because p4est goes
       * with the original order of coarse cells when it sets up its forest,
       * and then applies the Morton ordering within each tree. But if coarse
       * grid cells are badly ordered this may mean that individual parts of
       * the forest stored on a local machine may be split across coarse grid
       * cells that are not geometrically close. Consequently, we apply a
       * hierarchical preordering according to
       * SparsityTools::reorder_hierarchical() to ensure that the part of the
       * forest stored by p4est is located on geometrically close coarse grid
       * cells.
       */
      std::vector<types::global_dof_index>
        coarse_cell_to_p4est_tree_permutation;
      std::vector<types::global_dof_index>
        p4est_tree_to_coarse_cell_permutation;

      /**
       * Return a pointer to the p4est tree that belongs to the given
       * dealii_coarse_cell_index()
       */
      typename dealii::internal::p4est::types<dim>::tree *
      init_tree(const int dealii_coarse_cell_index) const;

      /**
       * The function that computes the permutation between the two data
       * storage schemes.
       */
      void
      setup_coarse_cell_to_p4est_tree_permutation();

      /**
       * Take the contents of a newly created triangulation we are attached to
       * and copy it to p4est data structures.
       *
       * This function exists in 2d and 3d variants.
       */
      void
      copy_new_triangulation_to_p4est(std::integral_constant<int, 2>);
      void
      copy_new_triangulation_to_p4est(std::integral_constant<int, 3>);

      /**
       * Copy the local part of the refined forest from p4est into the
       * attached triangulation.
       */
      void
      copy_local_forest_to_triangulation();

      /**
       * Internal function notifying all registered slots to provide their
       * weights before repartitioning occurs. Called from
       * execute_coarsening_and_refinement() and repartition().
       *
       * @return A vector of unsigned integers representing the weight or
       * computational load of every cell after the refinement/coarsening/
       * repartition cycle. Note that the number of entries does not need to
       * be equal to either n_active_cells() or n_locally_owned_active_cells(),
       * because the triangulation is not updated yet. The weights are sorted
       * in the order that p4est will encounter them while iterating over
       * them.
       */
      std::vector<unsigned int>
      get_cell_weights() const;

      /**
       * This method returns a bit vector of length tria.n_vertices()
       * indicating the locally active vertices on a level, i.e., the vertices
       * touched by the locally owned level cells for use in geometric
       * multigrid (possibly including the vertices due to periodic boundary
       * conditions) are marked by true.
       *
       * Used by DoFHandler::Policy::ParallelDistributed.
       */
      std::vector<bool>
      mark_locally_active_vertices_on_level(const int level) const;

      virtual unsigned int
      coarse_cell_id_to_coarse_cell_index(
        const types::coarse_cell_id coarse_cell_id) const override;

      virtual types::coarse_cell_id
      coarse_cell_index_to_coarse_cell_id(
        const unsigned int coarse_cell_index) const override;

      template <int, int, class>
      friend class dealii::FETools::internal::ExtrapolateImplementation;

      template <int, int>
      friend class TemporarilyMatchRefineFlags;
    };


    /**
     * Specialization of the general template for the 1d case. There is
     * currently no support for distributing 1d triangulations. Consequently,
     * all this class does is throw an exception.
     */
    template <int spacedim>
    class Triangulation<1, spacedim>
      : public dealii::parallel::DistributedTriangulationBase<1, spacedim>
    {
    public:
      /**
       * dummy settings
       */
      enum Settings
      {
        default_setting                          = 0x0,
        mesh_reconstruction_after_repartitioning = 0x1,
        construct_multigrid_hierarchy            = 0x2
      };

      /**
       * Constructor. The argument denotes the MPI communicator to be used for
       * the triangulation.
       */
      Triangulation(
        const MPI_Comm &mpi_communicator,
        const typename dealii::Triangulation<1, spacedim>::MeshSmoothing
                       smooth_grid = (dealii::Triangulation<1, spacedim>::none),
        const Settings settings    = default_setting);

      /**
       * Destructor.
       */
      virtual ~Triangulation() override;

      /**
       * Return a permutation vector for the order the coarse cells are
       * handed of to p4est. For example the first element i in this vector
       * denotes that the first cell in hierarchical ordering is the ith deal
       * cell starting from begin(0).
       */
      const std::vector<types::global_dof_index> &
      get_p4est_tree_to_coarse_cell_permutation() const;

      /**
       * This function is not implemented, but needs to be present for the
       * compiler.
       */
      virtual void
      load(const std::string &filename,
           const bool         autopartition = true) override;

      /**
       * This function is not implemented, but needs to be present for the
       * compiler.
       */
      virtual void
      save(const std::string &filename) const override;

      /**
       * This function is not implemented, but needs to be present for the
       * compiler.
       */
      virtual bool
      is_multilevel_hierarchy_constructed() const override;

      /**
       * This function is not implemented, but needs to be present for the
       * compiler.
       */
      virtual void
      update_cell_relations() override;

      /**
       * Dummy arrays. This class isn't usable but the compiler wants to see
       * these variables at a couple places anyway.
       */
      std::vector<types::global_dof_index>
        coarse_cell_to_p4est_tree_permutation;
      std::vector<types::global_dof_index>
        p4est_tree_to_coarse_cell_permutation;

      /**
       * This method, which is only implemented for dim = 2 or 3,
       * needs a stub because it is used in dof_handler_policy.cc
       */
      virtual std::map<unsigned int, std::set<dealii::types::subdomain_id>>
      compute_level_vertices_with_ghost_neighbors(
        const unsigned int level) const;

      /**
       * Like above, this method, which is only implemented for dim = 2 or 3,
       * needs a stub because it is used in dof_handler_policy.cc
       */
      virtual std::vector<bool>
      mark_locally_active_vertices_on_level(const unsigned int level) const;

      virtual unsigned int
      coarse_cell_id_to_coarse_cell_index(
        const types::coarse_cell_id coarse_cell_id) const override;

      virtual types::coarse_cell_id
      coarse_cell_index_to_coarse_cell_id(
        const unsigned int coarse_cell_index) const override;

      template <int, int>
      friend class TemporarilyMatchRefineFlags;
    };
  } // namespace distributed
} // namespace parallel


#else // DEAL_II_WITH_P4EST

namespace parallel
{
  namespace distributed
  {
    /**
     * Dummy class the compiler chooses for parallel distributed
     * triangulations if we didn't actually configure deal.II with the p4est
     * library. The existence of this class allows us to refer to
     * parallel::distributed::Triangulation objects throughout the library
     * even if it is disabled.
     *
     * Since the constructor of this class is deleted, no such objects
     * can actually be created as this would be pointless given that
     * p4est is not available.
     */
    template <int dim, int spacedim = dim>
    class Triangulation
      : public dealii::parallel::DistributedTriangulationBase<dim, spacedim>
    {
    public:
      /**
       * Dummy settings to allow defining the deleted constructor.
       */
      enum Settings
      {
        default_setting                          = 0x0,
        mesh_reconstruction_after_repartitioning = 0x1,
        construct_multigrid_hierarchy            = 0x2,
        no_automatic_repartitioning              = 0x4
      };

      /**
       * Constructor. Deleted to make sure that objects of this type cannot be
       * constructed (see also the class documentation).
       */
      explicit Triangulation(
        const MPI_Comm & /*mpi_communicator*/,
        const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
        /*smooth_grid*/
        = (dealii::Triangulation<dim, spacedim>::none),
        const Settings /*settings*/ = default_setting) = delete;

      /**
       * Dummy replacement to allow for better error messages when compiling
       * this class.
       */
      virtual bool
      is_multilevel_hierarchy_constructed() const override
      {
        return false;
      }

      /**
       * Dummy replacement to allow for better error messages when compiling
       * this class.
       */
      virtual void
      save(const std::string & /*filename*/) const override
      {}

      /**
       * Dummy replacement to allow for better error messages when compiling
       * this class.
       */
      virtual void
      load(const std::string & /*filename*/,
           const bool /*autopartition*/ = true) override
      {}

      /**
       * Dummy replacement to allow for better error messages when compiling
       * this class.
       */
      virtual void
      update_cell_relations() override
      {}
    };
  } // namespace distributed
} // namespace parallel


#endif


namespace parallel
{
  namespace distributed
  {
    /**
     * This class temporarily modifies the refine and coarsen flags of all
     * active cells to match the p4est oracle.
     *
     * The modification only happens on parallel::distributed::Triangulation
     * objects, and persists for the lifetime of an instantiation of this
     * class.
     *
     * The TemporarilyMatchRefineFlags class should only be used in
     * combination with the Triangulation::Signals::post_p4est_refinement
     * signal. At this stage, the p4est orcale already has been refined, but
     * the triangulation is still unchanged. After the modification, all
     * refine and coarsen flags describe how the traingulation will acutally
     * be refined.
     */
    template <int dim, int spacedim = dim>
    class TemporarilyMatchRefineFlags : public Subscriptor
    {
    public:
      /**
       * Constructor.
       *
       * Stores the refine and coarsen flags of all active cells if the
       * provided Triangulation is of type
       * parallel::distributed::Triangulation.
       *
       * Adjusts them to be consistent with the p4est oracle.
       */
      TemporarilyMatchRefineFlags(dealii::Triangulation<dim, spacedim> &tria);

      /**
       * Destructor.
       *
       * Returns the refine and coarsen flags of all active cells on the
       * parallel::distributed::Triangulation into their previous state.
       */
      ~TemporarilyMatchRefineFlags();

    private:
      /**
       * The modified parallel::distributed::Triangulation.
       */
      const SmartPointer<
        dealii::parallel::distributed::Triangulation<dim, spacedim>>
        distributed_tria;

      /**
       * A vector that temporarily stores the refine flags before they have
       * been modified on the parallel::distributed::Triangulation.
       */
      std::vector<bool> saved_refine_flags;

      /**
       * A vector that temporarily stores the coarsen flags before they have
       * been modified on the parallel::distributed::Triangulation.
       */
      std::vector<bool> saved_coarsen_flags;
    };
  } // namespace distributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif
