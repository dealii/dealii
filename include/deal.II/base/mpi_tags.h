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
       * This enum holds all MPI tags used in point to point MPI communications
       * inside the deal.II library.
       *
       * We keep these tags in a central location so that they are unique within
       * the library. Otherwise, communication that receives packages might pick
       * up packets from a different algorithm. This is especially true if
       * MPI_ANY_SOURCE is used.
       *
       * The list of MPI functions that use an MPI tag is:
       * - MPI_Send, MPI_Recv, MPI_Sendrecv
       * - MPI_Isend, MPI_Irecv
       * - MPI_Probe, MPI_Iprobe
       * - MPI_Comm_create_group, MPI_Intercomm_create,
       * Utilities::MPI::create_group
       */
      namespace Tags
      {
        /**
         * The enum with the tags.
         */
        enum enumeration : std::uint16_t
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

          /// ConsensusAlgorithms::Payload::get_requesters()
          consensus_algorithm_payload_get_requesters,

          /// FETools::extrapolate()
          fe_tools_extrapolate,
          /// FETools::extrapolate(), allocate space for 10 rounds:
          fe_tools_extrapolate_end = fe_tools_extrapolate + 10,

          /// ConsensusAlgorithms::NBX::process
          consensus_algorithm_nbx_answer_request,
          /// ConsensusAlgorithms::NBX::process
          consensus_algorithm_nbx_process_deliver,

          /// ConsensusAlgorithms::PEX::process
          consensus_algorithm_pex_answer_request,
          /// ConsensusAlgorithms::PEX::process
          consensus_algorithm_pex_process_deliver,

          /// TriangulationDescription::Utilities::create_description_from_triangulation()
          fully_distributed_create,

          /// TriangulationBase<dim, spacedim>::fill_level_ghost_owners()
          triangulation_base_fill_level_ghost_owners,

          /// GridTools::compute_local_to_global_vertex_index_map
          grid_tools_compute_local_to_global_vertex_index_map,
          /// GridTools::compute_local_to_global_vertex_index_map second tag
          grid_tools_compute_local_to_global_vertex_index_map2,

          /// ParticleHandler<dim, spacedim>::send_recv_particles
          particle_handler_send_recv_particles_setup,
          /// ParticleHandler<dim, spacedim>::send_recv_particles
          particle_handler_send_recv_particles_send,

          /// ScaLAPACKMatrix<NumberType>::copy_to
          scalapack_copy_to,
          /// ScaLAPACKMatrix<NumberType>::copy_to
          scalapack_copy_to2,
          /// ScaLAPACKMatrix<NumberType>::copy_from
          scalapack_copy_from,

          /// ProcessGrid::ProcessGrid
          process_grid_constructor,

          /// 200 tags for Partitioner::import_from_ghosted_array_start
          partitioner_import_start,
          partitioner_import_end = partitioner_import_start + 200,

          /// 200 tags for Partitioner::export_to_ghosted_array_start
          partitioner_export_start,
          partitioner_export_end = partitioner_export_start + 200,

          /// NoncontiguousPartitioner::update_values
          noncontiguous_partitioner_update_values,

          // Utilities::MPI::compute_union
          compute_union,

        };
      } // namespace Tags
    }   // namespace internal
  }     // namespace MPI
} // namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
