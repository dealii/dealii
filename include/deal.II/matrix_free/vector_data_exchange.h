// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
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
          const unsigned int             communication_channel,
          const ArrayView<const double> &locally_owned_array,
          const ArrayView<double> &      temporary_storage,
          const ArrayView<double> &      ghost_array,
          std::vector<MPI_Request> &     requests) const = 0;

        virtual void
        export_to_ghosted_array_finish(
          const ArrayView<double> & ghost_array,
          std::vector<MPI_Request> &requests) const = 0;

        virtual void
        import_from_ghosted_array_start(
          const VectorOperation::values vector_operation,
          const unsigned int            communication_channel,
          const ArrayView<double> &     ghost_array,
          const ArrayView<double> &     temporary_storage,
          std::vector<MPI_Request> &    requests) const = 0;

        virtual void
        import_from_ghosted_array_finish(
          const VectorOperation::values  vector_operation,
          const ArrayView<const double> &temporary_storage,
          const ArrayView<double> &      locally_owned_storage,
          const ArrayView<double> &      ghost_array,
          std::vector<MPI_Request> &     requests) const = 0;

        virtual void
        export_to_ghosted_array_start(
          const unsigned int            communication_channel,
          const ArrayView<const float> &locally_owned_array,
          const ArrayView<float> &      temporary_storage,
          const ArrayView<float> &      ghost_array,
          std::vector<MPI_Request> &    requests) const = 0;

        virtual void
        export_to_ghosted_array_finish(
          const ArrayView<float> &  ghost_array,
          std::vector<MPI_Request> &requests) const = 0;

        virtual void
        import_from_ghosted_array_start(
          const VectorOperation::values vector_operation,
          const unsigned int            communication_channel,
          const ArrayView<float> &      ghost_array,
          const ArrayView<float> &      temporary_storage,
          std::vector<MPI_Request> &    requests) const = 0;

        virtual void
        import_from_ghosted_array_finish(
          const VectorOperation::values vector_operation,
          const ArrayView<const float> &temporary_storage,
          const ArrayView<float> &      locally_owned_storage,
          const ArrayView<float> &      ghost_array,
          std::vector<MPI_Request> &    requests) const = 0;
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
          const unsigned int             communication_channel,
          const ArrayView<const double> &locally_owned_array,
          const ArrayView<double> &      temporary_storage,
          const ArrayView<double> &      ghost_array,
          std::vector<MPI_Request> &     requests) const override
        {
          partitioner->export_to_ghosted_array_start(communication_channel,
                                                     locally_owned_array,
                                                     temporary_storage,
                                                     ghost_array,
                                                     requests);
        }

        void
        export_to_ghosted_array_finish(
          const ArrayView<double> & ghost_array,
          std::vector<MPI_Request> &requests) const override
        {
          partitioner->export_to_ghosted_array_finish(ghost_array, requests);
        }

        void
        import_from_ghosted_array_start(
          const VectorOperation::values vector_operation,
          const unsigned int            communication_channel,
          const ArrayView<double> &     ghost_array,
          const ArrayView<double> &     temporary_storage,
          std::vector<MPI_Request> &    requests) const override
        {
          partitioner->import_from_ghosted_array_start(vector_operation,
                                                       communication_channel,
                                                       ghost_array,
                                                       temporary_storage,
                                                       requests);
        }

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values  vector_operation,
          const ArrayView<const double> &temporary_storage,
          const ArrayView<double> &      locally_owned_storage,
          const ArrayView<double> &      ghost_array,
          std::vector<MPI_Request> &     requests) const override
        {
          partitioner->import_from_ghosted_array_finish(vector_operation,
                                                        temporary_storage,
                                                        locally_owned_storage,
                                                        ghost_array,
                                                        requests);
        }
        void
        export_to_ghosted_array_start(
          const unsigned int            communication_channel,
          const ArrayView<const float> &locally_owned_array,
          const ArrayView<float> &      temporary_storage,
          const ArrayView<float> &      ghost_array,
          std::vector<MPI_Request> &    requests) const override
        {
          partitioner->export_to_ghosted_array_start(communication_channel,
                                                     locally_owned_array,
                                                     temporary_storage,
                                                     ghost_array,
                                                     requests);
        }

        void
        export_to_ghosted_array_finish(
          const ArrayView<float> &  ghost_array,
          std::vector<MPI_Request> &requests) const override
        {
          partitioner->export_to_ghosted_array_finish(ghost_array, requests);
        }

        void
        import_from_ghosted_array_start(
          const VectorOperation::values vector_operation,
          const unsigned int            communication_channel,
          const ArrayView<float> &      ghost_array,
          const ArrayView<float> &      temporary_storage,
          std::vector<MPI_Request> &    requests) const override
        {
          partitioner->import_from_ghosted_array_start(vector_operation,
                                                       communication_channel,
                                                       ghost_array,
                                                       temporary_storage,
                                                       requests);
        }

        void
        import_from_ghosted_array_finish(
          const VectorOperation::values vector_operation,
          const ArrayView<const float> &temporary_storage,
          const ArrayView<float> &      locally_owned_storage,
          const ArrayView<float> &      ghost_array,
          std::vector<MPI_Request> &    requests) const override
        {
          partitioner->import_from_ghosted_array_finish(vector_operation,
                                                        temporary_storage,
                                                        locally_owned_storage,
                                                        ghost_array,
                                                        requests);
        }

      private:
        const std::shared_ptr<const Utilities::MPI::Partitioner> partitioner;
      };

    } // namespace VectorDataExchange
  }   // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
