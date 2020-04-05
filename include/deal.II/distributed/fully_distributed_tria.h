// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_fully_distributed_tria_h
#define dealii_fully_distributed_tria_h


#include <deal.II/base/config.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/grid_tools.h>

#include <vector>

#ifdef DEAL_II_WITH_MPI
#  include <mpi.h>
#endif

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
// forward declaration of the data type for periodic face pairs
namespace GridTools
{
  template <typename CellIterator>
  struct PeriodicFacePair;
}
#endif

namespace parallel
{
  /**
   * A namespace for the fully distributed triangulation.
   *
   * @ingroup parallel
   */
  namespace fullydistributed
  {
    /**
     * A distributed triangulation with a distributed coarse grid.
     *
     * The motivation for parallel::fullydistributed::Triangulation has its
     * origins in the following observations about complex geometries and/or
     * about given meshes created by an external mesh generator. We regard
     * complex geometries as geometries that can be meshed only with a
     * non-negligible number of coarse cells (>10,000):
     * - storing the coarse-grid information on every process is too expensive
     *   from a memory point of view (as done by
     *   parallel::distributed::Triangulation). Normally, a process only needs a
     *   small section of the global triangulation, i.e., a small section of the
     *   coarse grid such that a partitioning of the coarse grid is indeed
     *   essential. The cells stored on each process consist of the
     *   @ref GlossLocallyOwnedCell "locally owned cells" and the
     *   @ref GlossGhostCell "ghost cells".
     * - the distribution of the active cells - on the finest level - among all
     *   processes by simply partitioning a space-filling curve might not lead
     *   to an optimal result for triangulations that originate from large
     *   coarse grids: e.g. partitions that belong to the same process might
     *   be discontinuous, leading to increased communication (within a
     *   node and beyond). Graph-based partitioning algorithms might be a sound
     *   alternative to the space filling curve used by
     *   parallel::distributed::Triangulation.
     *
     * To be able to construct a fully partitioned triangulation that
     * distributes the coarse grid and gives flexibility regarding partitioning,
     * the following ingredients are required:
     * - a locally relevant coarse-grid triangulation
     *   (vertices, cell definition; including a layer of ghost cells)
     * - a mapping of the locally relevant coarse-grid triangulation into the
     *   global coarse-grid triangulation
     * - information about which cell should be refined as well as information
     *   regarding the subdomain_id, the level_subdomain_id, manifold_id,
     *   and boundary_id of each cell.
     *
     * The ingredients listed above are bundled in the struct
     * TriangulationDescription::Description. The user has to fill this
     * data structure - in a pre-processing step - before actually creating the
     * triangulation. Predefined functions to create
     * TriangulationDescription::Description can be found in the namespace
     * TriangulationDescription::Utilities.
     *
     * Once the TriangulationDescription::Description `construction_data` has
     * been constructed, the triangulation `tria` can be created by calling
     * `tria.create_triangulation(construction_data);`.
     *
     * @note This triangulation supports: 1D/2D/3D, hanging nodes,
     *       geometric multigrid, and periodicity.
     *
     * @note You can create a triangulation with hanging nodes and multigrid
     *       levels with create_triangulation(). However, once it has been
     *       created, it cannot be altered anymore, i.e. you cannot coarsen or
     *       refine afterwards.
     *
     * @note Currently only simple periodicity conditions (i.e. without offsets
     *       and rotation matrices - see also the documentation of
     *       GridTools::collect_periodic_faces()) are supported.
     *
     * @author Peter Munch, 2019
     */
    template <int dim, int spacedim = dim>
    class Triangulation
      : public parallel::DistributedTriangulationBase<dim, spacedim>
    {
    public:
      using cell_iterator =
        typename dealii::Triangulation<dim, spacedim>::cell_iterator;

      using active_cell_iterator =
        typename dealii::Triangulation<dim, spacedim>::active_cell_iterator;

      using CellStatus =
        typename dealii::Triangulation<dim, spacedim>::CellStatus;

      /**
       * Constructor.
       *
       * @param mpi_communicator The MPI communicator to be used for the
       *                         triangulation.
       */
      explicit Triangulation(MPI_Comm mpi_communicator);

      /**
       * Destructor.
       */
      virtual ~Triangulation() = default;

      /**
       * @copydoc dealii::Triangulation::create_triangulation()
       *
       * @note This is the function to be used instead of
       * Triangulation::create_triangulation() for some of the other
       * triangulations of deal.II.
       */
      void
      create_triangulation(
        const TriangulationDescription::Description<dim, spacedim>
          &construction_data) override;

      /**
       * @note This function is not implemented for this class  and throws
       *       an assertion. Instead, use
       *       the other create_triangulation() function to create the
       *       triangulation.
       */
      virtual void
      create_triangulation(const std::vector<Point<spacedim>> &      vertices,
                           const std::vector<dealii::CellData<dim>> &cells,
                           const SubCellData &subcelldata) override;

      /**
       * Implementation of the same function as in the base class.
       *
       * @param other_tria The triangulation to be copied. It can be a serial
       *        Triangulation or a parallel::distributed::Triangulation. Both
       *        can have been refined already.
       *
       * @note This function uses the partitioner registered with
       *       set_partitioner().
       */
      void
      copy_triangulation(
        const dealii::Triangulation<dim, spacedim> &other_tria) override;

      /**
       * Register a partitioner, which is used within the method
       * copy_triangulation.
       *
       * @param partitioner A partitioning function, which takes as input argument
       *                    a reference to the triangulation to be partitioned
       *                    and the number of partitions to be created.
       *                    The function needs to set subdomain
       *                    ids for each active cell of the given triangulation,
       *                    with values between zero (inclusive)
       *                    and the second argument to the function (exclusive).
       * @param settings See the description of the Settings enumerator.
       *
       * @note As a default, GridTools::partition_triangulation_zorder() is used
       *       as partitioner and data structures on multigrid levels are not
       *       set up.
       */
      void
      set_partitioner(
        const std::function<void(dealii::Triangulation<dim, spacedim> &,
                                 const unsigned int)> &partitioner,
        const TriangulationDescription::Settings &     settings);

      /**
       * Coarsen and refine the mesh according to refinement and coarsening
       * flags set.
       *
       * @note Not implemented yet.
       */
      virtual void
      execute_coarsening_and_refinement() override;

      /**
       * Override the implementation of prepare_coarsening_and_refinement from
       * the base class.
       *
       * @note Not implemented yet.
       */
      virtual bool
      prepare_coarsening_and_refinement() override;

      /**
       * Return true if the triangulation has hanging nodes.
       *
       * @note Not implemented yet.
       */
      virtual bool
      has_hanging_nodes() const override;

      /**
       * Return the local memory consumption in bytes.
       */
      virtual std::size_t
      memory_consumption() const override;

      virtual bool
      is_multilevel_hierarchy_constructed() const override;

      virtual unsigned int
      coarse_cell_id_to_coarse_cell_index(
        const types::coarse_cell_id coarse_cell_id) const override;

      virtual types::coarse_cell_id
      coarse_cell_index_to_coarse_cell_id(
        const unsigned int coarse_cell_index) const override;

    private:
      /**
       * Override the function to update the number cache so we can fill data
       * like @p level_ghost_owners.
       */
      virtual void
      update_number_cache() override;

      /**
       * store the Settings.
       */
      TriangulationDescription::Settings settings;

      /**
       * Partitioner used in copy_triangulation().
       */
      std::function<void(dealii::Triangulation<dim, spacedim> &,
                         const unsigned int)>
        partitioner;

      /**
       * Sorted list of pairs of coarse-cell ids and their indices.
       */
      std::vector<std::pair<types::coarse_cell_id, unsigned int>>
        coarse_cell_id_to_coarse_cell_index_vector;

      /**
       * List of the coarse-cell id for each coarse cell (stored at
       * cell->index()).
       */
      std::vector<types::coarse_cell_id>
        coarse_cell_index_to_coarse_cell_id_vector;

      /**
       * Boolean indicating that the function create_triangulation() was called
       * for internal usage.
       */
      bool currently_processing_create_triangulation_for_internal_usage;

      /**
       * Boolean indicating that the function
       * prepare_coarsening_and_refinement() was called for internal usage.
       */
      bool
        currently_processing_prepare_coarsening_and_refinement_for_internal_usage;
    };

  } // namespace fullydistributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif
