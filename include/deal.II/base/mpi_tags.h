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

#ifndef dealii_mpi_tags_h
#define dealii_mpi_tags_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{
  namespace MPI
  {
    namespace internal
    {
      /**
       * This enum holds all MPI tags used in collective MPI communications with
       * using MPI_ANY_SOURCE inside the deal.II library.
       *
       * We keep these tags in a central location so that they are unique within
       * the library. Otherwise, communication that receives packages with
       * MPI_ANY_SOURCE might pick up packets from a different algorithm.
       */
      struct Tags
      {
        /**
         * The enum with the tags.
         */
        enum enumeration
        {
          /// Utilities::MPI::some_to_some()
          mpi_some_to_some = 300,

          /// Utilities::MPI::compute_point_to_point_communication_pattern()
          compute_point_to_point_communication_pattern,

          /// GridTools::exchange_cell_data_to_ghosts():
          exchange_cell_data_to_ghosts,

          /// Triangulation<dim, spacedim>::communicate_locally_moved_vertices()
          triangulation_communicate_locally_moved_vertices,

          /// dof_handler_policy.cc: communicate_mg_ghost_cells()
          dofhandler_communicate_mg_ghost_cells,

          /// dof_handler_policy.cc: communicate_mg_ghost_cells()
          dofhandler_communicate_mg_ghost_cells_reply,

          /// mg_transfer_internal.cc: fill_copy_indices()
          mg_transfer_fill_copy_indices,

          /// SparsityTools::sparsity_tools_distribute_sparsity_pattern()
          sparsity_tools_distribute_sparsity_pattern,

          /// Dictionary::reinit()
          dictionary_reinit,

          /// ConsensusAlgorithmPayload::get_requesters()
          consensus_algorithm_payload_get_requesters,

          /// FETools::extrapolate()
          fe_tools_extrapolate,
          /// FETools::extrapolate(), allocate space for 10 rounds:
          fe_tools_extrapolate_end = fe_tools_extrapolate + 10,

          /// ConsensusAlgorithm_NBX::process
          consensus_algorithm_nbx_process_request,
          /// ConsensusAlgorithm_NBX::process
          consensus_algorithm_nbx_process_deliver,

          /// ConsensusAlgorithm_PEX::process
          consensus_algorithm_pex_process_request,
          /// ConsensusAlgorithm_PEX::process
          consensus_algorithm_pex_process_deliver,

        };
      };
    } // namespace internal
  }   // namespace MPI
} // namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
