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

    } // namespace Utilities
  }   // namespace fullydistributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif
