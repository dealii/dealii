// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_vector_data_exchange_h
#define dealii_matrix_free_vector_data_exchange_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/mpi_stub.h>
#include <deal.II/base/partitioner.h>

#include <deal.II/lac/vector_operation.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * Namespace containing classes for inter-process data exchange (i.e.,
     * for update_ghost_values and compress) in MatrixFree.
     */
    namespace VectorDataExchange
    {
      /**
       * Interface needed by MatrixFree.
       */
      class Base
      {
      public:
        virtual ~Base() = default;

        virtual unsigned int
        locally_owned_size() const = 0;

        virtual unsigned int
        n_ghost_indices() const = 0;

        virtual unsigned int
        n_import_indices() const = 0;

        virtual unsigned int
        n_import_sm_procs() const = 0;

        virtual types::global_dof_index
        size() const = 0;

        virtual void
        export_to_ghosted_array_start(
          const unsigned int                          communication_channel,
          const ArrayView<const double>              &locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double>                    &ghost_array,
          const ArrayView<double>                    &temporary_storage,
          std::vector<MPI_Request>                   &requests) const = 0;

        virtual void
        export_to_ghosted_array_finish(
          const ArrayView<const double>              &locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double>                    &ghost_array,
          std::vector<MPI_Request>                   &requests) const = 0;

        virtual void
        import_from_ghosted_array_start(
          const VectorOperation::values               vector_operation,
          const unsigned int                          communication_channel,
          const ArrayView<const double>              &locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double>                    &ghost_array,
          const ArrayView<double>                    &temporary_storage,
          std::vector<MPI_Request>                   &requests) const = 0;

        virtual void
        import_from_ghosted_array_finish(
          const VectorOperation::values               vector_operation,
          const ArrayView<double>                    &locally_owned_storage,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double>                    &ghost_array,
          const ArrayView<const double>              &temporary_storage,
          std::vector<MPI_Request>                   &requests) const = 0;

        virtual void
        reset_ghost_values(const ArrayView<double> &ghost_array) const = 0;

        virtual void
        export_to_ghosted_array_start(
          const unsigned int                         communication_channel,
          const ArrayView<const float>              &locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float>                    &ghost_array,
          const ArrayView<float>                    &temporary_storage,
          std::vector<MPI_Request>                  &requests) const = 0;

        virtual void
        export_to_ghosted_array_finish(
          const ArrayView<const float>              &locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float>                    &ghost_array,
          std::vector<MPI_Request>                  &requests) const = 0;

        virtual void
        import_from_ghosted_array_start(
          const VectorOperation::values              vector_operation,
          const unsigned int                         communication_channel,
          const ArrayView<const float>              &locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float>                    &ghost_array,
          const ArrayView<float>                    &temporary_storage,
          std::vector<MPI_Request>                  &requests) const = 0;

        virtual void
        import_from_ghosted_array_finish(
          const VectorOperation::values              vector_operation,
          const ArrayView<float>                    &locally_owned_storage,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float>                    &ghost_array,
          const ArrayView<const float>              &temporary_storage,
          std::vector<MPI_Request>                  &requests) const = 0;

        virtual void
        reset_ghost_values(const ArrayView<float> &ghost_array) const = 0;
      };


      /**
       * Class that simply delegates the task to a Utilities::MPI::Partitioner.
       */
      class PartitionerWrapper : public Base
      {
      public:
        PartitionerWrapper(
          const std::shared_ptr<const Utilities::MPI::Partitioner>
            &partitioner);

        virtual ~PartitionerWrapper() = default;

        unsigned int
        locally_owned_size() const override;

        unsigned int
        n_ghost_indices() const override;

        unsigned int
        n_import_indices() const override;

        unsigned int
        n_import_sm_procs() const override;

        types::global_dof_index
        size() const override;

        void
        export_to_ghosted_array_start(
          const unsigned int                          communication_channel,
          const ArrayView<const double>              &locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double>                    &ghost_array,
          const ArrayView<double>                    &temporary_storage,
          std::vector<MPI_Request>                   &requests) const override;

        void
        export_to_ghosted_array_finish(
          const ArrayView<const double>              &locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double>                    &ghost_array,
          std::vector<MPI_Request>                   &requests) const override;

        void
        import_from_ghosted_array_start(
          const VectorOperation::values               vector_operation,
          const unsigned int                          communication_channel,
          const ArrayView<const double>              &locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double>                    &ghost_array,
          const ArrayView<double>                    &temporary_storage,
          std::vector<MPI_Request>                   &requests) const override;

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values               vector_operation,
          const ArrayView<double>                    &locally_owned_storage,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double>                    &ghost_array,
          const ArrayView<const double>              &temporary_storage,
          std::vector<MPI_Request>                   &requests) const override;

        void
        reset_ghost_values(const ArrayView<double> &ghost_array) const override;

        void
        export_to_ghosted_array_start(
          const unsigned int                         communication_channel,
          const ArrayView<const float>              &locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float>                    &ghost_array,
          const ArrayView<float>                    &temporary_storage,
          std::vector<MPI_Request>                  &requests) const override;

        void
        export_to_ghosted_array_finish(
          const ArrayView<const float>              &locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float>                    &ghost_array,
          std::vector<MPI_Request>                  &requests) const override;

        void
        import_from_ghosted_array_start(
          const VectorOperation::values              vector_operation,
          const unsigned int                         communication_channel,
          const ArrayView<const float>              &locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float>                    &ghost_array,
          const ArrayView<float>                    &temporary_storage,
          std::vector<MPI_Request>                  &requests) const override;

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values              vector_operation,
          const ArrayView<float>                    &locally_owned_storage,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float>                    &ghost_array,
          const ArrayView<const float>              &temporary_storage,
          std::vector<MPI_Request>                  &requests) const override;

        void
        reset_ghost_values(const ArrayView<float> &ghost_array) const override;

      private:
        template <typename Number>
        void
        reset_ghost_values_impl(const ArrayView<Number> &ghost_array) const;

        const std::shared_ptr<const Utilities::MPI::Partitioner> partitioner;
      };



      /**
       * Similar to the above but using the internal data structures in the
       * partitioner in order to identify indices of degrees of freedom that are
       * in the same shared memory region.
       */
      class Full : public Base
      {
      public:
        Full(
          const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
          const MPI_Comm communicator_sm);

        unsigned int
        locally_owned_size() const override;

        unsigned int
        n_ghost_indices() const override;

        unsigned int
        n_import_indices() const override;

        virtual unsigned int
        n_import_sm_procs() const override;

        virtual types::global_dof_index
        size() const override;

        MPI_Comm
        get_sm_mpi_communicator() const;

        void
        export_to_ghosted_array_start(
          const unsigned int                          communication_channel,
          const ArrayView<const double>              &locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double>                    &ghost_array,
          const ArrayView<double>                    &temporary_storage,
          std::vector<MPI_Request>                   &requests) const override;

        void
        export_to_ghosted_array_finish(
          const ArrayView<const double>              &locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double>                    &ghost_array,
          std::vector<MPI_Request>                   &requests) const override;

        void
        import_from_ghosted_array_start(
          const VectorOperation::values               vector_operation,
          const unsigned int                          communication_channel,
          const ArrayView<const double>              &locally_owned_array,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double>                    &ghost_array,
          const ArrayView<double>                    &temporary_storage,
          std::vector<MPI_Request>                   &requests) const override;

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values               vector_operation,
          const ArrayView<double>                    &locally_owned_storage,
          const std::vector<ArrayView<const double>> &shared_arrays,
          const ArrayView<double>                    &ghost_array,
          const ArrayView<const double>              &temporary_storage,
          std::vector<MPI_Request>                   &requests) const override;

        void
        reset_ghost_values(const ArrayView<double> &ghost_array) const override;

        void
        export_to_ghosted_array_start(
          const unsigned int                         communication_channel,
          const ArrayView<const float>              &locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float>                    &ghost_array,
          const ArrayView<float>                    &temporary_storage,
          std::vector<MPI_Request>                  &requests) const override;

        void
        export_to_ghosted_array_finish(
          const ArrayView<const float>              &locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float>                    &ghost_array,
          std::vector<MPI_Request>                  &requests) const override;

        void
        import_from_ghosted_array_start(
          const VectorOperation::values              vector_operation,
          const unsigned int                         communication_channel,
          const ArrayView<const float>              &locally_owned_array,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float>                    &ghost_array,
          const ArrayView<float>                    &temporary_storage,
          std::vector<MPI_Request>                  &requests) const override;

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values              vector_operation,
          const ArrayView<float>                    &locally_owned_storage,
          const std::vector<ArrayView<const float>> &shared_arrays,
          const ArrayView<float>                    &ghost_array,
          const ArrayView<const float>              &temporary_storage,
          std::vector<MPI_Request>                  &requests) const override;

        void
        reset_ghost_values(const ArrayView<float> &ghost_array) const override;

      private:
        template <typename Number>
        void
        export_to_ghosted_array_start_impl(
          const unsigned int                          communication_channel,
          const ArrayView<const Number>              &locally_owned_array,
          const std::vector<ArrayView<const Number>> &shared_arrays,
          const ArrayView<Number>                    &ghost_array,
          const ArrayView<Number>                    &temporary_storage,
          std::vector<MPI_Request>                   &requests) const;

        template <typename Number>
        void
        export_to_ghosted_array_finish_impl(
          const ArrayView<const Number>              &locally_owned_array,
          const std::vector<ArrayView<const Number>> &shared_arrays,
          const ArrayView<Number>                    &ghost_array,
          std::vector<MPI_Request>                   &requests) const;

        template <typename Number>
        void
        import_from_ghosted_array_start_impl(
          const VectorOperation::values               vector_operation,
          const unsigned int                          communication_channel,
          const ArrayView<const Number>              &locally_owned_array,
          const std::vector<ArrayView<const Number>> &shared_arrays,
          const ArrayView<Number>                    &ghost_array,
          const ArrayView<Number>                    &temporary_storage,
          std::vector<MPI_Request>                   &requests) const;

        template <typename Number>
        void
        import_from_ghosted_array_finish_impl(
          const VectorOperation::values               vector_operation,
          const ArrayView<Number>                    &locally_owned_storage,
          const std::vector<ArrayView<const Number>> &shared_arrays,
          const ArrayView<Number>                    &ghost_array,
          const ArrayView<const Number>              &temporary_storage,
          std::vector<MPI_Request>                   &requests) const;

        template <typename Number>
        void
        reset_ghost_values_impl(const ArrayView<Number> &ghost_array) const;

      private:
        /**
         * Global communicator.
         */
        const MPI_Comm comm;

        /**
         * Shared-memory sub-communicator.
         */
        const MPI_Comm comm_sm;

        /**
         * Number of locally-owned vector entries.
         */
        const unsigned int n_local_elements;

        /**
         * Number of ghost vector entries.
         */
        const unsigned int n_ghost_elements;

        /**
         * Number of global vector entries.
         */
        const types::global_dof_index n_global_elements;

        /**
         * A variable caching the number of ghost indices in a larger set of
         * indices by rank.
         */
        std::vector<unsigned int> n_ghost_indices_in_larger_set_by_remote_rank;

        /**
         * The set of indices that appear for an IndexSet that is a subset of a
         * larger set for each rank in a compressed manner.
         */
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
          ghost_indices_subset_data;

        /**
         * An array that contains information which processors my ghost indices
         * belong to, at which offset and how many those indices are
         */
        std::vector<std::array<unsigned int, 3>> ghost_targets_data;

        /**
         * The set of processors and length of data field which send us their
         * ghost data.
         *
         * @note Structured as ghost_targets_data.
         */
        std::vector<std::array<unsigned int, 3>> import_targets_data;

        /**
         * An array that caches the number of chunks in the import indices per
         * MPI rank. The length is import_indices_data.size()+1.
         *
         * The set of (local) indices that we are importing during compress()
         * from remote processes, i.e., others' ghosts that belong to the local
         * range.
         */
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
          import_indices_data;

        /**
         * Shared-memory ranks from which data is copied from during
         * export_to_ghosted_array_finish().
         */
        std::vector<unsigned int> sm_ghost_ranks;

        /**
         * Indices from where to copy data from during
         * export_to_ghosted_array_finish().
         */
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
          sm_export_data;

        /**
         * Indices where to copy data to during
         * export_to_ghosted_array_finish().
         */
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
          sm_export_data_this;

        /**
         * Shared-memory ranks from where to copy data from during
         * import_from_ghosted_array_finish().
         */
        std::vector<unsigned int> sm_import_ranks;

        /**
         * Indices from where to copy data from during
         * import_from_ghosted_array_finish().
         */
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
          sm_import_data;

        /**
         * Indices where to copy data to during
         * import_from_ghosted_array_finish().
         */
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
          sm_import_data_this;
      };

    } // namespace VectorDataExchange
  }   // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
