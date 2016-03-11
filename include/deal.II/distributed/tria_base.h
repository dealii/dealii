// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2016 by the deal.II authors
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

#ifndef dealii__distributed__tria_base_h
#define dealii__distributed__tria_base_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/mpi.h>
#include <deal.II/grid/tria.h>

#include <deal.II/base/std_cxx1x/function.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <set>
#include <vector>
#include <list>
#include <utility>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  /**
   * This class describes the interface for all triangulation classes that
   * work in parallel, namely parallel::distributed::Triangulation and
   * parallel::shared::Triangulation.
   */
  template <int dim, int spacedim = dim>
  class Triangulation : public dealii::Triangulation<dim,spacedim>
  {
  public:

    /**
     * Constructor.
     */
    Triangulation (MPI_Comm mpi_communicator,
                   const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing smooth_grid = (dealii::Triangulation<dim,spacedim>::none),
                   const bool check_for_distorted_cells = false);

    /**
     * Destructor.
     */
    virtual ~Triangulation ();

    /**
     * Return MPI communicator used by this triangulation.
     */
    virtual MPI_Comm get_communicator () const;

    /**
     * Implementation of the same function as in the base class.
     */
    virtual void copy_triangulation (const dealii::Triangulation<dim, spacedim> &old_tria);

    /**
     * Return the number of active cells owned by each of the MPI processes
     * that contribute to this triangulation. The element of this vector
     * indexed by locally_owned_subdomain() equals the result of
     * n_locally_owned_active_cells().
     */
    const std::vector<unsigned int> &
    n_locally_owned_active_cells_per_processor () const;


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
    unsigned int n_locally_owned_active_cells () const;

    /**
     * Return the sum over all processors of the number of active cells owned
     * by each processor. This equals the overall number of active cells in
     * the triangulation.
     */
    virtual types::global_dof_index n_global_active_cells () const;

    /**
     * Return the local memory consumption in bytes.
     */
    virtual std::size_t memory_consumption () const;


    /**
     * Returns the global maximum level. This may be bigger than the number
     * dealii::Triangulation::n_levels() (a function in this class's base
     * class) returns if the current processor only stores cells in parts of
     * the domain that are not very refined, but if other processors store
     * cells in more deeply refined parts of the domain.
     */
    virtual unsigned int n_global_levels () const;

    /**
     * Return the subdomain id of those cells that are owned by the current
     * processor. All cells in the triangulation that do not have this
     * subdomain id are either owned by another processor or have children
     * that only exist on other processors.
     */
    types::subdomain_id locally_owned_subdomain () const;

    /**
     * Returns a set of MPI ranks of the processors that have at least one
     * ghost cell adjacent to the cells of the local processor. In other
     * words, this is the set of subdomain_id() for all ghost cells.
     *
     * @note: If @p i is contained in the list of processor @p j, then @p j
     * will also be contained in the list of processor @p i.
     */
    const std::set<unsigned int> &ghost_owners () const;

    /**
     * Returns a set of MPI ranks of the processors that have at least one
     * level ghost cell adjacent to our cells used in geometric multigrid. In
     * other words, this is the set of level_subdomain_id() for all level
     * ghost cells.
     *
     * @note: If @p i is contained in the list of processor @p j, then @p j
     * will also be contained in the list of processor @p i.
     */
    const std::set<unsigned int> &level_ghost_owners () const;

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
       * This vector stores the number of locally owned active cells per MPI
       * rank.
       */
      std::vector<unsigned int> n_locally_owned_active_cells;
      /**
       * The total number of active cells (sum of @p
       * n_locally_owned_active_cells).
       */
      types::global_dof_index   n_global_active_cells;
      /**
       * The global number of levels computed as the maximum number of levels
       * taken over all MPI ranks, so <tt>n_levels()<=n_global_levels =
       * max(n_levels() on proc i)</tt>.
       */
      unsigned int              n_global_levels;
      /**
       * A set containing the subdomain_id (MPI rank) of the owners of the
       * ghost cells on this processor.
       */
      std::set<unsigned int> ghost_owners;
      /**
       * A set containing the MPI ranks of the owners of the level ghost cells
       * on this processor (for all levels).
       */
      std::set<unsigned int> level_ghost_owners;

      NumberCache();
    };

    NumberCache number_cache;

    /**
     * Update the number_cache variable after mesh creation or refinement.
     */
    virtual void update_number_cache ();


  };

} // namespace parallel

DEAL_II_NAMESPACE_CLOSE

#endif
