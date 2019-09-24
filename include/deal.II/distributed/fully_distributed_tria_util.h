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

#ifndef dealii_fully_distributed_tria_util_h
#define dealii_fully_distributed_tria_util_h

#include <deal.II/distributed/fully_distributed_tria.h>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace fullydistributed
  {
    /**
     * A namespace for fully distributed triangulation utility functions.
     *
     * @ingroup fullydistributed
     */
    namespace Utilities
    {
      /**
       * Construct parallel::fullydistributed::ConstructionData from a given
       * partitioned triangulation `tria` and a specified process.
       * The input triangulation can be either
       * a serial triangulation of type dealii::Triangulation which has been
       * colored (subdomain_id and/or level_subdomain_id has been set) or a
       * distributed triangulation of type dealii::distributed::Triangulation,
       * where the partitioning is adopted unaltered.
       *
       * @param tria Partitioned input triangulation.
       * @param comm MPI_Communicator to be used. In the case
       *             of dealii::distributed::Triangulation, the communicators
       *             have to match.
       * @param construct_multilevel_hierarchy Signal if the multigrid levels
       *        should be constructed.
       * @param my_rank_in Construct ConstructionData for this rank (only
       *                   working for serial triangulations).
       * @return ConstructionData to be used to setup
       *         parallel::fullydistributed::Triangulation.
       *
       * @note Multilevel hierarchies are supported if it is enabled in
       *       parallel::fullydistributed::Triangulation.
       *
       * @note Hanging nodes in the input triangulation are supported. However,
       *       to be able to use this
       *       feature, the user has to enable multilevel hierarchy support in
       *       parallel::fullydistributed::Triangulation.
       *
       * @author Peter Munch, 2019
       */
      template <int dim, int spacedim = dim>
      ConstructionData<dim, spacedim>
      create_construction_data_from_triangulation(
        const dealii::Triangulation<dim, spacedim> &tria,
        const MPI_Comm                              comm,
        const bool         construct_multilevel_hierarchy = false,
        const unsigned int my_rank_in = numbers::invalid_unsigned_int);


      /**
       * Construct a parallel::fullydistributed::ConstructionData. In contrast
       * to the function above, this function is also responsible for creating
       * a serial triangulation and for its partitioning (by calling the
       * provided std::functions). Internally only selected processes (every
       * n-th/each root of a group of size group_size) create a serial
       * triangulation and the ConstructionData for all processes in its group,
       * which is communicated.
       *
       * @note A reasonable group size is the size of a NUMA domain or the
       * size of a compute node.
       *
       * @param serial_grid_generator A function, which creates a serial triangulation.
       * @param serial_grid_partitioner A function, which can partition a serial
       *        triangulation, i.e., sets the sudomain_ids of the active cells.
       *        The function takes as the first argument a serial triangulation,
       *        as the second argument the MPI communicator, and as the third
       *        argument the group size.
       * @param comm MPI communicator
       * @param group_size The size of each group.
       * @param construct_multilevel_hierarchy Construct multigrid levels.
       * @return ConstructionData to be used to setup
       *         parallel::fullydistributed::Triangulation.
       *
       * @author Peter Munch, 2019
       */
      template <int dim, int spacedim = dim>
      ConstructionData<dim, spacedim>
      create_construction_data_from_triangulation_in_groups(
        std::function<void(dealii::Triangulation<dim, spacedim> &)>
                                                serial_grid_generator,
        std::function<void(dealii::Triangulation<dim, spacedim> &,
                           const MPI_Comm,
                           const unsigned int)> serial_grid_partitioner,
        const MPI_Comm                          comm,
        const int                               group_size = 1,
        const bool construct_multilevel_hierarchy          = false);

    } // namespace Utilities
  }   // namespace fullydistributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif
