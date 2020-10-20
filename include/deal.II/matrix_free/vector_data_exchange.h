// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


#ifndef dealii_matrix_free_vector_data_exchange_h
#define dealii_matrix_free_vector_data_exchange_h


#include <deal.II/base/config.h>

#include <deal.II/base/partitioner.h>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    namespace VectorDataExchange
    {
      class Base
      {
      public:
        virtual unsigned int
        local_size() const = 0;

        virtual unsigned int
        n_ghost_indices() const = 0;

        virtual unsigned int
        n_import_indices() const = 0;

        virtual unsigned int
        n_import_sm_procs() const
        {
          return 0; // TODO
        }

        virtual types::global_dof_index
        size() const = 0;

        virtual void
        export_to_ghosted_array_start(
          const unsigned int                          communication_channel,
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const = 0;

        virtual void
        export_to_ghosted_array_finish(
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          std::vector<MPI_Request> &                  requests) const = 0;

        virtual void
        import_from_ghosted_array_start(
          const VectorOperation::values               vector_operation,
          const unsigned int                          communication_channel,
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const = 0;

        virtual void
        import_from_ghosted_array_finish(
          const VectorOperation::values               vector_operation,
          const ArrayView<double> &                   locally_owned_storage,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<const double> &             temporary_storage,
          std::vector<MPI_Request> &                  requests) const = 0;

        virtual void
        reset_ghost_values(const ArrayView<double> &ghost_array) const = 0;

        virtual void
        export_to_ghosted_array_start(
          const unsigned int                         communication_channel,
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &                 requests) const = 0;

        virtual void
        export_to_ghosted_array_finish(
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          std::vector<MPI_Request> &                 requests) const = 0;

        virtual void
        import_from_ghosted_array_start(
          const VectorOperation::values              vector_operation,
          const unsigned int                         communication_channel,
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &                 requests) const = 0;

        virtual void
        import_from_ghosted_array_finish(
          const VectorOperation::values              vector_operation,
          const ArrayView<float> &                   locally_owned_storage,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<const float> &             temporary_storage,
          std::vector<MPI_Request> &                 requests) const = 0;

        virtual void
        reset_ghost_values(const ArrayView<float> &ghost_array) const = 0;
      };



      class PartitionerWrapper : public Base
      {
      public:
        PartitionerWrapper(
          const std::shared_ptr<const Utilities::MPI::Partitioner> partitioner)
          : partitioner(partitioner)
        {}


        unsigned int
        local_size() const override
        {
          return partitioner->local_size();
        }

        unsigned int
        n_ghost_indices() const override
        {
          return partitioner->n_ghost_indices();
        }

        unsigned int
        n_import_indices() const override
        {
          return partitioner->n_import_indices();
        }

        virtual types::global_dof_index
        size() const override
        {
          return partitioner->size();
        }

        void
        export_to_ghosted_array_start(
          const unsigned int                          communication_channel,
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const override
        {
          (void)shared_arrays;
          partitioner->export_to_ghosted_array_start(communication_channel,
                                                     locally_owned_array,
                                                     temporary_storage,
                                                     ghost_array,
                                                     requests);
        }

        void
        export_to_ghosted_array_finish(
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          std::vector<MPI_Request> &                  requests) const override
        {
          (void)locally_owned_array;
          (void)shared_arrays;
          partitioner->export_to_ghosted_array_finish(ghost_array, requests);
        }

        void
        import_from_ghosted_array_start(
          const VectorOperation::values               vector_operation,
          const unsigned int                          communication_channel,
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const override
        {
          (void)locally_owned_array;
          (void)shared_arrays;
          partitioner->import_from_ghosted_array_start(vector_operation,
                                                       communication_channel,
                                                       ghost_array,
                                                       temporary_storage,
                                                       requests);
        }

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values               vector_operation,
          const ArrayView<double> &                   locally_owned_storage,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<const double> &             temporary_storage,
          std::vector<MPI_Request> &                  requests) const override
        {
          (void)shared_arrays;
          partitioner->import_from_ghosted_array_finish(vector_operation,
                                                        temporary_storage,
                                                        locally_owned_storage,
                                                        ghost_array,
                                                        requests);
        }

        void
        reset_ghost_values(const ArrayView<double> &ghost_array) const override
        {
          reset_ghost_values_impl(ghost_array);
        }

        void
        export_to_ghosted_array_start(
          const unsigned int                         communication_channel,
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &                 requests) const override
        {
          (void)shared_arrays;
          partitioner->export_to_ghosted_array_start(communication_channel,
                                                     locally_owned_array,
                                                     temporary_storage,
                                                     ghost_array,
                                                     requests);
        }

        void
        export_to_ghosted_array_finish(
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          std::vector<MPI_Request> &                 requests) const override
        {
          (void)locally_owned_array;
          (void)shared_arrays;
          partitioner->export_to_ghosted_array_finish(ghost_array, requests);
        }

        void
        import_from_ghosted_array_start(
          const VectorOperation::values              vector_operation,
          const unsigned int                         communication_channel,
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &                 requests) const override
        {
          (void)locally_owned_array;
          (void)shared_arrays;
          partitioner->import_from_ghosted_array_start(vector_operation,
                                                       communication_channel,
                                                       ghost_array,
                                                       temporary_storage,
                                                       requests);
        }

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values              vector_operation,
          const ArrayView<float> &                   locally_owned_storage,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<const float> &             temporary_storage,
          std::vector<MPI_Request> &                 requests) const override
        {
          (void)shared_arrays;
          partitioner->import_from_ghosted_array_finish(vector_operation,
                                                        temporary_storage,
                                                        locally_owned_storage,
                                                        ghost_array,
                                                        requests);
        }

        void
        reset_ghost_values(const ArrayView<float> &ghost_array) const override
        {
          reset_ghost_values_impl(ghost_array);
        }

      private:
        template <typename Number>
        void
        reset_ghost_values_impl(const ArrayView<Number> &ghost_array) const
        {
          for (const auto &my_ghosts :
               partitioner->ghost_indices_within_larger_ghost_set())
            for (unsigned int j = my_ghosts.first; j < my_ghosts.second; ++j)
              ghost_array[j] = 0.;
        }

        const std::shared_ptr<const Utilities::MPI::Partitioner> partitioner;
      };



      class Full : public Base
      {
      public:
        Full(const IndexSet &is_locally_owned,
             const IndexSet &is_locally_ghost,
             const MPI_Comm &comm,
             const MPI_Comm &comm_sm,
             const std::vector<std::pair<unsigned int, unsigned int>>
               &ghost_indices_within_larger_ghost_set);

        unsigned int
        local_size() const override;

        unsigned int
        n_ghost_indices() const override;

        unsigned int
        n_import_indices() const override;

        virtual unsigned int
        n_import_sm_procs() const override
        {
          return send_sm_ranks.size() + recv_sm_ranks.size(); // TODO
        }

        virtual types::global_dof_index
        size() const override
        {
          return n_global_elements;
        }

        void
        export_to_ghosted_array_start(
          const unsigned int                          communication_channel,
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const override;

        void
        export_to_ghosted_array_finish(
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          std::vector<MPI_Request> &                  requests) const override;

        void
        import_from_ghosted_array_start(
          const VectorOperation::values               vector_operation,
          const unsigned int                          communication_channel,
          const ArrayView<const double> &             locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const override;

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values               vector_operation,
          const ArrayView<double> &                   locally_owned_storage,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double> &                   ghost_array,
          const ArrayView<const double> &             temporary_storage,
          std::vector<MPI_Request> &                  requests) const override;

        void
        reset_ghost_values(const ArrayView<double> &ghost_array) const override;

        void
        export_to_ghosted_array_start(
          const unsigned int                         communication_channel,
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &                 requests) const override;

        void
        export_to_ghosted_array_finish(
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          std::vector<MPI_Request> &                 requests) const override;

        void
        import_from_ghosted_array_start(
          const VectorOperation::values              vector_operation,
          const unsigned int                         communication_channel,
          const ArrayView<const float> &             locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &                 requests) const override;

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values              vector_operation,
          const ArrayView<float> &                   locally_owned_storage,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float> &                   ghost_array,
          const ArrayView<const float> &             temporary_storage,
          std::vector<MPI_Request> &                 requests) const override;

        void
        reset_ghost_values(const ArrayView<float> &ghost_array) const override;

      private:
        template <typename Number>
        void
        export_to_ghosted_array_start_impl(
          const unsigned int                          communication_channel,
          const ArrayView<const Number> &             locally_owned_array,
          const std::vector<ArrayView<const Number>> &shared_arrays,
          const ArrayView<Number> &                   ghost_array,
          const ArrayView<Number> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const;

        template <typename Number>
        void
        export_to_ghosted_array_finish_impl(
          const ArrayView<const Number> &             locally_owned_array,
          const std::vector<ArrayView<const Number>> &shared_arrays,
          const ArrayView<Number> &                   ghost_array,
          std::vector<MPI_Request> &                  requests) const;

        template <typename Number>
        void
        import_from_ghosted_array_start_impl(
          const VectorOperation::values               vector_operation,
          const unsigned int                          communication_channel,
          const ArrayView<const Number> &             locally_owned_array,
          const std::vector<ArrayView<const Number>> &shared_arrays,
          const ArrayView<Number> &                   ghost_array,
          const ArrayView<Number> &                   temporary_storage,
          std::vector<MPI_Request> &                  requests) const;

        template <typename Number>
        void
        import_from_ghosted_array_finish_impl(
          const VectorOperation::values               vector_operation,
          const ArrayView<Number> &                   locally_owned_storage,
          const std::vector<ArrayView<const Number>> &shared_arrays,
          const ArrayView<Number> &                   ghost_array,
          const ArrayView<const Number> &             temporary_storage,
          std::vector<MPI_Request> &                  requests) const;

        template <typename Number>
        void
        reset_ghost_values_impl(const ArrayView<Number> &ghost_array) const;

      private:
        /**
         * Global communicator.
         */
        MPI_Comm comm;

        /**
         * Shared-memory sub-communicator.
         */
        MPI_Comm comm_sm;

        /**
         * Number of processes in comm.
         */
        unsigned int n_mpi_processes_;

        /**
         * Number of locally-owned vector entries.
         */
        unsigned int n_local_elements;
        
        /**
         * Number of global vector entries.
         */
        types::global_dof_index n_global_elements;

        /**
         * Number of ghost vector entries.
         */
        unsigned int n_ghost_elements;

        /**
         * TODO
         */
        std::vector<unsigned int> recv_remote_ranks;

        /**
         * TODO
         */
        std::vector<std::pair<types::global_dof_index, types::global_dof_index>>
          recv_remote_ptr;

        std::vector<unsigned int> shifts;
        std::vector<unsigned int> shifts_ptr;

        /**
         * TODO
         */
        std::vector<unsigned int> recv_sm_ranks;

        /**
         * TODO
         */
        std::vector<unsigned int> recv_sm_ptr = {0};

        /**
         * TODO
         */
        std::vector<unsigned int> recv_sm_indices;

        /**
         * TODO
         */
        std::vector<unsigned int> recv_sm_len;

        /**
         * TODO
         */
        std::vector<unsigned int> recv_sm_offset;

        /**
         * TODO
         */
        std::vector<unsigned int> send_remote_ranks;

        /**
         * TODO
         */
        std::vector<unsigned int> send_remote_ptr = {0};

        /**
         * TODO
         */
        std::vector<unsigned int> send_remote_indices;

        /**
         * TODO
         */
        std::vector<unsigned int> send_remote_len;

        /**
         * TODO
         */
        std::vector<unsigned int> send_remote_offset = {0};

        /**
         * TODO
         */
        std::vector<unsigned int> send_sm_ranks;

        /**
         * TODO
         */
        std::vector<unsigned int> send_sm_ptr = {0};

        /**
         * TODO
         */
        std::vector<unsigned int> send_sm_indices;

        /**
         * TODO
         */
        std::vector<unsigned int> send_sm_len;

        /**
         * TODO
         */
        std::vector<unsigned int> send_sm_offset;
      };



      /**
       * Partitioner for discontinuous Galerkin discretizations, exploiting
       * shared memory.
       */
      class Contiguous : public Base
      {
        const dealii::types::global_dof_index         dofs_per_cell;
        const dealii::types::global_dof_index         dofs_per_face;
        const std::vector<std::vector<unsigned int>> &face_to_cell_index_nodal;

      public:
        using RankType     = unsigned int;
        using LocalDoFType = unsigned int;
        using CellIdType   = dealii::types::global_dof_index;
        using FaceIdType   = std::pair<CellIdType, unsigned int>;

        /**
         * Constructor.
         */
        template <typename ShapeInfo>
        Contiguous(const ShapeInfo &shape_info)
          : dofs_per_cell(shape_info.dofs_per_cell)
          , dofs_per_face(shape_info.dofs_per_face)
          , face_to_cell_index_nodal(shape_info.face_to_cell_index_nodal)
        {}

        /**
         * Initialize partitioner with a list of locally owned cells and
         * a list of ghost faces (cell and face no).
         */
        void
        reinit(const std::vector<dealii::types::global_dof_index> local_cells,
               const std::vector<std::pair<dealii::types::global_dof_index,
                                           std::vector<unsigned int>>>
                              local_ghost_faces,
               const MPI_Comm comm,
               const MPI_Comm comm_sm,
               const bool     do_buffering);

        unsigned int
        local_size() const override;

        unsigned int
        n_ghost_indices() const override;

        unsigned int
        n_import_indices() const override;

        void
        reset_ghost_values(
          const dealii::ArrayView<double> &ghost_array) const override;

        void
        reset_ghost_values(
          const dealii::ArrayView<float> &ghost_array) const override;

        /**
         * Start to export to ghost array.
         */
        void
        export_to_ghosted_array_start(
          const unsigned int                     communication_channel,
          const dealii::ArrayView<const double> &locally_owned_array,
          const std::vector<dealii::ArrayView<const double>> &shared_arrays,
          const dealii::ArrayView<double> &                   ghost_array,
          const dealii::ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &requests) const override;

        /**
         * Finish to export to ghost array.
         */
        void
        export_to_ghosted_array_finish(
          const dealii::ArrayView<const double> &locally_owned_array,
          const std::vector<dealii::ArrayView<const double>> &shared_arrays,
          const dealii::ArrayView<double> &                   ghost_array,
          std::vector<MPI_Request> &requests) const override;

        /**
         * Start to import from ghost array.
         */
        void
        import_from_ghosted_array_start(
          const dealii::VectorOperation::values  vector_operation,
          const unsigned int                     communication_channel,
          const dealii::ArrayView<const double> &locally_owned_array,
          const std::vector<dealii::ArrayView<const double>> &shared_arrays,
          const dealii::ArrayView<double> &                   ghost_array,
          const dealii::ArrayView<double> &                   temporary_storage,
          std::vector<MPI_Request> &requests) const override;

        /**
         * Finish to import from ghost array.
         */
        void
        import_from_ghosted_array_finish(
          const dealii::VectorOperation::values vector_operation,
          const dealii::ArrayView<double> &     locally_owned_storage,
          const std::vector<dealii::ArrayView<const double>> &shared_arrays,
          const dealii::ArrayView<double> &                   ghost_array,
          const dealii::ArrayView<const double> &             temporary_storage,
          std::vector<MPI_Request> &requests) const override;

        /**
         * Start to export to ghost array.
         */
        void
        export_to_ghosted_array_start(
          const unsigned int                    communication_channel,
          const dealii::ArrayView<const float> &locally_owned_array,
          const std::vector<dealii::ArrayView<const float>> &shared_arrays,
          const dealii::ArrayView<float> &                   ghost_array,
          const dealii::ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &requests) const override;

        /**
         * Finish to export to ghost array.
         */
        void
        export_to_ghosted_array_finish(
          const dealii::ArrayView<const float> &locally_owned_array,
          const std::vector<dealii::ArrayView<const float>> &shared_arrays,
          const dealii::ArrayView<float> &                   ghost_array,
          std::vector<MPI_Request> &requests) const override;

        /**
         * Start to import from ghost array.
         */
        void
        import_from_ghosted_array_start(
          const dealii::VectorOperation::values vector_operation,
          const unsigned int                    communication_channel,
          const dealii::ArrayView<const float> &locally_owned_array,
          const std::vector<dealii::ArrayView<const float>> &shared_arrays,
          const dealii::ArrayView<float> &                   ghost_array,
          const dealii::ArrayView<float> &                   temporary_storage,
          std::vector<MPI_Request> &requests) const override;

        /**
         * Finish to import from ghost array.
         */
        void
        import_from_ghosted_array_finish(
          const dealii::VectorOperation::values vector_operation,
          const dealii::ArrayView<float> &      locally_owned_storage,
          const std::vector<dealii::ArrayView<const float>> &shared_arrays,
          const dealii::ArrayView<float> &                   ghost_array,
          const dealii::ArrayView<const float> &             temporary_storage,
          std::vector<MPI_Request> &requests) const override;

        /**
         * TODO.
         */
        template <typename Number>
        void
        export_to_ghosted_array_finish_0(
          const dealii::ArrayView<const Number> &locally_owned_array,
          const std::vector<dealii::ArrayView<const Number>> &shared_arrays,
          const dealii::ArrayView<Number> &                   ghost_array,
          std::vector<MPI_Request> &                          requests) const;

        /**
         * TODO.
         */
        template <typename Number>
        void
        export_to_ghosted_array_finish_1(
          const dealii::ArrayView<const Number> &locally_owned_array,
          const std::vector<dealii::ArrayView<const Number>> &shared_arrays,
          const dealii::ArrayView<Number> &                   ghost_array,
          std::vector<MPI_Request> &                          requests) const;

      private:
        /**
         * Actual type-independent implementation of
         * export_to_ghosted_array_start().
         */
        template <typename Number>
        void
        export_to_ghosted_array_start_impl(
          const unsigned int                     communication_channel,
          const dealii::ArrayView<const Number> &locally_owned_array,
          const std::vector<dealii::ArrayView<const Number>> &shared_arrays,
          const dealii::ArrayView<Number> &                   ghost_array,
          const dealii::ArrayView<Number> &                   temporary_storage,
          std::vector<MPI_Request> &                          requests) const;

        /**
         * Actual type-independent implementation of
         * export_to_ghosted_array_finish().
         */
        template <typename Number>
        void
        export_to_ghosted_array_finish_impl(
          const dealii::ArrayView<const Number> &locally_owned_array,
          const std::vector<dealii::ArrayView<const Number>> &shared_arrays,
          const dealii::ArrayView<Number> &                   ghost_array,
          std::vector<MPI_Request> &                          requests) const;

        /**
         * Actual type-independent implementation of
         * import_from_ghosted_array_start().
         */
        template <typename Number>
        void
        import_from_ghosted_array_start_impl(
          const dealii::VectorOperation::values  vector_operation,
          const unsigned int                     communication_channel,
          const dealii::ArrayView<const Number> &locally_owned_array,
          const std::vector<dealii::ArrayView<const Number>> &shared_arrays,
          const dealii::ArrayView<Number> &                   ghost_array,
          const dealii::ArrayView<Number> &                   temporary_storage,
          std::vector<MPI_Request> &                          requests) const;

        /**
         * Actual type-independent implementation of
         * import_from_ghosted_array_finish().
         */
        template <typename Number>
        void
        import_from_ghosted_array_finish_impl(
          const dealii::VectorOperation::values vector_operation,
          const dealii::ArrayView<Number> &     locally_owned_storage,
          const std::vector<dealii::ArrayView<const Number>> &shared_arrays,
          const dealii::ArrayView<Number> &                   ghost_array,
          const dealii::ArrayView<const Number> &             temporary_storage,
          std::vector<MPI_Request> &                          requests) const;

      public:
        /**
         * Return position of shared cell: cell -> (owner, offset)
         */
        const std::map<dealii::types::global_dof_index,
                       std::pair<unsigned int, unsigned int>> &
        get_maps() const;


        /**
         * Return position of ghost face: (cell, no) -> (owner, offset)
         */
        const std::map<std::pair<dealii::types::global_dof_index, unsigned int>,
                       std::pair<unsigned int, unsigned int>> &
        get_maps_ghost() const;

        /**
         * Return memory consumption.
         *
         * @note: Only counts the buffers [TODO].
         */
        std::size_t
        memory_consumption() const;

        /**
         * Synchronize.
         */
        void
        sync(const unsigned int tag = 0) const;

      private:
        /**
         * Global communicator.
         */
        MPI_Comm comm;

        /**
         * Shared-memory sub-communicator.
         */
        MPI_Comm comm_sm;

        /**
         * Number of processes in comm.
         */
        unsigned int n_mpi_processes_;

        /**
         * Number of locally-owned vector entries.
         */
        unsigned int n_local_elements;

        /**
         * Number of ghost vector entries.
         */
        unsigned int n_ghost_elements;

        // I) configuration parameters
        bool         do_buffering;   // buffering vs. non-buffering modus
        unsigned int dofs_per_ghost; // ghost face or ghost cell

        // II) MPI-communicator related stuff
        unsigned int sm_size;
        unsigned int sm_rank;

        // III) access cells and ghost faces
        std::map<CellIdType, std::pair<RankType, LocalDoFType>> maps;
        std::map<FaceIdType, std::pair<RankType, LocalDoFType>> maps_ghost;

        // III) information to pack/unpack buffers
        std::vector<unsigned int>                    send_ranks;
        std::vector<dealii::types::global_dof_index> send_ptr;
        std::vector<dealii::types::global_dof_index> send_data_id;
        std::vector<unsigned int>                    send_data_face_no;

        std::vector<unsigned int>                    recv_ranks;
        std::vector<dealii::types::global_dof_index> recv_ptr;
        std::vector<dealii::types::global_dof_index> recv_size;

        std::vector<unsigned int> sm_targets;
        std::vector<unsigned int> sm_sources;


        std::vector<dealii::types::global_dof_index> sm_send_ptr;
        std::vector<unsigned int>                    sm_send_rank;
        std::vector<unsigned int>                    sm_send_offset_1;
        std::vector<unsigned int>                    sm_send_offset_2;
        std::vector<unsigned int>                    sm_send_no;

        std::vector<dealii::types::global_dof_index> sm_recv_ptr;
        std::vector<unsigned int>                    sm_recv_rank;
        std::vector<unsigned int>                    sm_recv_offset_1;
        std::vector<unsigned int>                    sm_recv_offset_2;
        std::vector<unsigned int>                    sm_recv_no;
      };

    } // namespace VectorDataExchange
  }   // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
