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

#ifndef dealii_la_sm_partitioner_h
#define dealii_la_sm_partitioner_h


#include <deal.II/base/config.h>

#include <deal.II/base/mpi_compute_index_owner_internal.h>

#include <deal.II/lac/communication_pattern_base.h>

DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  namespace SharedMPI
  {
    /**
     * Partitioner base class to be used in context of
     * LinearAlgebra::SharedMPI::Vector.
     */
    class PartitionerBase : public LinearAlgebra::CommunicationPatternBase
    {
    public:
      /**
       * Constructor.
       */
      PartitionerBase(const bool contiguous_allocation);

      /**
       * Return global communicator.
       */
      virtual const MPI_Comm &
      get_mpi_communicator() const override final;

      /**
       * Return shared-memory communicator.
       */
      virtual const MPI_Comm &
      get_sm_mpi_communicator() const final;

      /**
       * Return if the partitioner implementation has been enabled with
       * contiguous or non-noncontiguous shared-memory allocation.
       */
      bool
      contiguous_allocation_enabled() const;

      /**
       * Return number of locally-owned vector entries.
       */
      virtual std::size_t
      local_size() const final;

      /**
       * Return number of ghost vector entries.
       */
      virtual std::size_t
      n_ghost_indices() const final;

      /**
       * Return number of processes in the global communicator.
       */
      virtual std::size_t
      n_mpi_processes() const final;

      /**
       * Return an estimate for the memory consumption, in bytes, of this
       * object.
       */
      virtual std::size_t
      memory_consumption() const = 0;

      /**
       * Start to export to ghost array.
       */
      virtual void
      export_to_ghosted_array_start(
        const unsigned int             communication_channel,
        double *const                  data_this,
        const std::vector<double *> &  data_others,
        dealii::AlignedVector<double> &buffer,
        std::vector<MPI_Request> &     requests) const = 0;

      /**
       * Finish to export to ghost array.
       */
      virtual void
      export_to_ghosted_array_finish(
        double *const                data_this,
        const std::vector<double *> &data_others,
        std::vector<MPI_Request> &   requests) const = 0;

      /**
       * Start to import from ghost array.
       */
      virtual void
      import_from_ghosted_array_start(
        const VectorOperation::values  operation,
        const unsigned int             communication_channel,
        double *const                  data_this,
        const std::vector<double *> &  data_others,
        dealii::AlignedVector<double> &buffer,
        std::vector<MPI_Request> &     requests) const = 0;

      /**
       * Finish to import from ghost array.
       */
      virtual void
      import_from_ghosted_array_finish(
        const VectorOperation::values        operation,
        double *const                        data_this,
        const std::vector<double *> &        data_others,
        const dealii::AlignedVector<double> &buffer,
        std::vector<MPI_Request> &           requests) const = 0;

      /**
       * Start to export to ghost array.
       */
      virtual void
      export_to_ghosted_array_start(
        const unsigned int            communication_channel,
        float *const                  data_this,
        const std::vector<float *> &  data_others,
        dealii::AlignedVector<float> &buffer,
        std::vector<MPI_Request> &    requests) const = 0;

      /**
       * Finish to export to ghost array.
       */
      virtual void
      export_to_ghosted_array_finish(
        float *const                data_this,
        const std::vector<float *> &data_others,
        std::vector<MPI_Request> &  requests) const = 0;

      /**
       * Start to import from ghost array.
       */
      virtual void
      import_from_ghosted_array_start(
        const VectorOperation::values operation,
        const unsigned int            communication_channel,
        float *const                  data_this,
        const std::vector<float *> &  data_others,
        dealii::AlignedVector<float> &buffer,
        std::vector<MPI_Request> &    requests) const = 0;

      /**
       * Finish to import from ghost array.
       */
      virtual void
      import_from_ghosted_array_finish(
        const VectorOperation::values       operation,
        float *const                        data_this,
        const std::vector<float *> &        data_others,
        const dealii::AlignedVector<float> &buffer,
        std::vector<MPI_Request> &          requests) const = 0;

    private:
      /**
       * Flag indicating if the partitioner implementation has been enabled
       * with contiguous or non-noncontiguous shared-memory allocation.
       */
      bool contiguous_allocation;

    protected:
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
    };

    /**
     * Partitioner implementation that is built around index-sets similar to
     * Utilities::MPI::Partitioner, however, distinguishes between remote and
     * shared ghost vector entries.
     */
    class Partitioner : public PartitionerBase
    {
    public:
      /**
       * Constructor.
       */
      Partitioner(const IndexSet &is_locally_owned,
                  const IndexSet &is_locally_ghost,
                  const MPI_Comm &comm,
                  const MPI_Comm &comm_sm);

      /**
       * Set up internal data structures.
       *
       * @note Not implemented. Be explicit and use the other one.
       */
      void
      reinit(const IndexSet &is_locally_owned,
             const IndexSet &is_locally_ghost,
             const MPI_Comm &communicator) override;

      /**
       * Set up internal data structures with a global and a shared-memory
       * communicator.
       */
      void
      reinit(const IndexSet &is_locally_owned,
             const IndexSet &is_locally_ghost,
             const MPI_Comm &communicator,
             const MPI_Comm &communicator_sm);

      /**
       * @copydoc PartitionerBase::export_to_ghosted_array_start()
       */
      void
      export_to_ghosted_array_start(
        const unsigned int             communication_channel,
        double *const                  data_this,
        const std::vector<double *> &  data_others,
        dealii::AlignedVector<double> &buffer,
        std::vector<MPI_Request> &     requests) const override;

      /**
       * @copydoc PartitionerBase::export_to_ghosted_array_finish()
       */
      void
      export_to_ghosted_array_finish(
        double *const                data_this,
        const std::vector<double *> &data_others,
        std::vector<MPI_Request> &   requests) const override;

      /**
       * @copydoc PartitionerBase::import_from_ghosted_array_start()
       */
      void
      import_from_ghosted_array_start(
        const VectorOperation::values  operation,
        const unsigned int             communication_channel,
        double *const                  data_this,
        const std::vector<double *> &  data_others,
        dealii::AlignedVector<double> &buffer,
        std::vector<MPI_Request> &     requests) const override;

      /**
       * @copydoc PartitionerBase::import_from_ghosted_array_finish()
       */
      void
      import_from_ghosted_array_finish(
        const VectorOperation::values        operation,
        double *const                        data_this,
        const std::vector<double *> &        data_others,
        const dealii::AlignedVector<double> &buffer,
        std::vector<MPI_Request> &           requests) const override;

      /**
       * @copydoc PartitionerBase::export_to_ghosted_array_start()
       */
      void
      export_to_ghosted_array_start(
        const unsigned int            communication_channel,
        float *const                  data_this,
        const std::vector<float *> &  data_others,
        dealii::AlignedVector<float> &buffer,
        std::vector<MPI_Request> &    requests) const override;

      /**
       * @copydoc PartitionerBase::export_to_ghosted_array_finish()
       */
      void
      export_to_ghosted_array_finish(
        float *const                data_this,
        const std::vector<float *> &data_others,
        std::vector<MPI_Request> &  requests) const override;

      /**
       * @copydoc PartitionerBase::import_from_ghosted_array_start()
       */
      void
      import_from_ghosted_array_start(
        const VectorOperation::values operation,
        const unsigned int            communication_channel,
        float *const                  data_this,
        const std::vector<float *> &  data_others,
        dealii::AlignedVector<float> &buffer,
        std::vector<MPI_Request> &    requests) const override;

      /**
       * @copydoc PartitionerBase::import_from_ghosted_array_finish()
       */
      void
      import_from_ghosted_array_finish(
        const VectorOperation::values       operation,
        float *const                        data_this,
        const std::vector<float *> &        data_others,
        const dealii::AlignedVector<float> &buffer,
        std::vector<MPI_Request> &          requests) const override;

      /**
       * @copydoc PartitionerBase::memory_consumption()
       */
      std::size_t
      memory_consumption() const override;

    private:
      /**
       * Actual type-independent implementation of
       * export_to_ghosted_array_start().
       */
      template <typename Number>
      void
      export_to_ghosted_array_start_impl(
        const unsigned int             communication_channel,
        Number *const                  data_this,
        const std::vector<Number *> &  data_others,
        dealii::AlignedVector<Number> &buffer,
        std::vector<MPI_Request> &     requests) const;

      /**
       * Actual type-independent implementation of
       * export_to_ghosted_array_finish().
       */
      template <typename Number>
      void
      export_to_ghosted_array_finish_impl(
        Number *const                data_this,
        const std::vector<Number *> &data_others,
        std::vector<MPI_Request> &   requests) const;

      /**
       * Actual type-independent implementation of
       * import_from_ghosted_array_start().
       */
      template <typename Number>
      void
      import_from_ghosted_array_start_impl(
        const VectorOperation::values  operation,
        const unsigned int             communication_channel,
        Number *const                  data_this,
        const std::vector<Number *> &  data_others,
        dealii::AlignedVector<Number> &buffer,
        std::vector<MPI_Request> &     requests) const;

      /**
       * Actual type-independent implementation of
       * import_from_ghosted_array_finish().
       */
      template <typename Number>
      void
      import_from_ghosted_array_finish_impl(
        const VectorOperation::values        operation,
        Number *const                        data_this,
        const std::vector<Number *> &        data_others,
        const dealii::AlignedVector<Number> &buffer,
        std::vector<MPI_Request> &           requests) const;

    private:
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
    };

  } // end of namespace SharedMPI
} // end of namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#endif
