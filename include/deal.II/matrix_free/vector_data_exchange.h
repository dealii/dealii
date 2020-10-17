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

        virtual const std::vector<std::pair<unsigned int, unsigned int>> &
        ghost_indices_within_larger_ghost_set() const = 0;

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

        const std::vector<std::pair<unsigned int, unsigned int>> &
        ghost_indices_within_larger_ghost_set() const override
        {
          return partitioner->ghost_indices_within_larger_ghost_set();
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

      private:
        const std::shared_ptr<const Utilities::MPI::Partitioner> partitioner;
      };

      class Full : public Base
      {
      public:
        Full(const IndexSet &is_locally_owned,
             const IndexSet &is_locally_ghost,
             const MPI_Comm &comm,
             const MPI_Comm &comm_sm);

        unsigned int
        local_size() const override;

        unsigned int
        n_ghost_indices() const override;

        unsigned int
        n_import_indices() const override
        {
          return 0; // TODO
        }

        const std::vector<std::pair<unsigned int, unsigned int>> &
        ghost_indices_within_larger_ghost_set() const override
        {
          return dummy; // TODO
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

        /**
         * TODO
         */
        std::vector<unsigned int> recv_remote_ranks;

        /**
         * TODO
         */
        std::vector<types::global_dof_index> recv_remote_ptr = {0};

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
        std::vector<unsigned int> send_remote_offset;

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

        std::vector<std::pair<unsigned int, unsigned int>> dummy;
      };

    } // namespace VectorDataExchange
  }   // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
