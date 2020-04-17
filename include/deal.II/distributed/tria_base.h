// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2019 by the deal.II authors
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
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/grid/tria.h>

#include <functional>
#include <list>
#include <set>
#include <tuple>
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
   * they communicate via @ref GlossMPICommunicator "MPI communicators" or
   * that they have
   * @ref GlossLocallyOwnedCell "locally owned",
   * @ref GlossGhostCell "ghost", and possibly
   * @ref GlossArtificialCell "artificial cells". This class provides
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
      MPI_Comm mpi_communicator,
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
    virtual const MPI_Comm &
    get_communicator() const;

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
     * Return the number of active cells owned by each of the MPI processes
     * that contribute to this triangulation. The element of this vector
     * indexed by locally_owned_subdomain() equals the result of
     * n_locally_owned_active_cells().
     *
     * @note This function involves global communication!
     */
    std::vector<unsigned int>
    compute_n_locally_owned_active_cells_per_processor() const;

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
     * Return a map that, for each vertex, lists all the processors whose
     * subdomains are adjacent to that vertex.
     */
    virtual std::map<unsigned int, std::set<dealii::types::subdomain_id>>
    compute_vertices_with_ghost_neighbors() const;

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

  protected:
    /**
     * MPI communicator to be used for the triangulation. We create a unique
     * communicator for this class, which is a duplicate of the one passed to
     * the constructor.
     */
    MPI_Comm mpi_communicator;

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

      NumberCache();
    };

    NumberCache number_cache;

    /**
     * Update the number_cache variable after mesh creation or refinement.
     */
    virtual void
    update_number_cache();

    /**
     * Store MPI ranks of level ghost owners of this processor on all levels.
     */
    void
    fill_level_ghost_owners();
  };

  /**
   *  Using directive for backwards-compatibility.
   *  @deprecated Use TriangulationBase instead of Triangulation.
   */
  template <int dim, int spacedim = dim>
  using Triangulation DEAL_II_DEPRECATED = TriangulationBase<dim, spacedim>;



  /**
   * A base class for distributed triangulations, i.e., triangulations that
   * do not store all cells on all processors. This implies that not
   * every detail of a triangulation may be known on each processor.
   * In particular, you have to expect that triangulations of classes
   * derived from this one only store some of the active cells (namely,
   * the @ref GlossLocallyOwnedCell "locally owned cells"), along
   * with @ref GlossGhostCell "ghost cells" and possibly
   * @ref GlossArtificialCell "artificial cells". In contrast to the classes
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
      MPI_Comm mpi_communicator,
      const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
                 smooth_grid = (dealii::Triangulation<dim, spacedim>::none),
      const bool check_for_distorted_cells = false);
  };

} // namespace parallel

DEAL_II_NAMESPACE_CLOSE

#endif
