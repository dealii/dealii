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
        local_size() const = 0;

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
        local_size() const override;

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
        reset_ghost_values_impl(const ArrayView<Number> &ghost_array) const;

        const std::shared_ptr<const Utilities::MPI::Partitioner> partitioner;
      };

    } // namespace VectorDataExchange
  }   // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
