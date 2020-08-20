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
    class PartitionerBase : public LinearAlgebra::CommunicationPatternBase
    {
    public:
      virtual const MPI_Comm &
      get_mpi_communicator() const override = 0;

      virtual const MPI_Comm &
      get_sm_mpi_communicator() const = 0;

      // double versions

      virtual void
      export_to_ghosted_array_start(
        double *                       data_this,
        std::vector<double *> &        data_others,
        dealii::AlignedVector<double> &buffer,
        const unsigned int             communication_channel = 0) const = 0;

      virtual void
      export_to_ghosted_array_finish(
        double *                       data_this,
        std::vector<double *> &        data_others,
        dealii::AlignedVector<double> &buffer) const = 0;

      virtual void
      import_from_ghosted_array_start(
        double *                       data_this,
        std::vector<double *> &        data_others,
        dealii::AlignedVector<double> &buffer,
        const unsigned int             communication_channel = 0) const = 0;

      virtual void
      import_from_ghosted_array_finish(
        double *                       data_this,
        std::vector<double *> &        data_others,
        dealii::AlignedVector<double> &buffer) const = 0;

      // float versions

      virtual void
      export_to_ghosted_array_start(
        float *                       data_this,
        std::vector<float *> &        data_others,
        dealii::AlignedVector<float> &buffer,
        const unsigned int            communication_channel = 0) const = 0;

      virtual void
      export_to_ghosted_array_finish(
        float *                       data_this,
        std::vector<float *> &        data_others,
        dealii::AlignedVector<float> &buffer) const = 0;

      virtual void
      import_from_ghosted_array_start(
        float *                       data_this,
        std::vector<float *> &        data_others,
        dealii::AlignedVector<float> &buffer,
        const unsigned int            communication_channel = 0) const = 0;

      virtual void
      import_from_ghosted_array_finish(
        float *                       data_this,
        std::vector<float *> &        data_others,
        dealii::AlignedVector<float> &buffer) const = 0;

      virtual std::size_t
      local_size() const = 0;

      virtual std::size_t
      n_ghost_indices() const = 0;

      virtual std::size_t
      n_mpi_processes() const = 0;

      virtual std::size_t
      memory_consumption() const = 0;

      virtual std::vector<unsigned int>
      get_sm_view() const = 0;
    };

    class Partitioner : public PartitionerBase
    {
    public:
      Partitioner(const MPI_Comm &comm,
                  const MPI_Comm &comm_sm,
                  const IndexSet &is_locally_owned,
                  const IndexSet &is_locally_ghost);

      const MPI_Comm &
      get_mpi_communicator() const override;

      const MPI_Comm &
      get_sm_mpi_communicator() const override;

      void
      reinit(const IndexSet &is_locally_owned,
             const IndexSet &is_locally_ghost,
             const MPI_Comm &communicator) override;

      void
      reinit(const IndexSet &is_locally_owned,
             const IndexSet &is_locally_ghost);

      void
      export_to_ghosted_array_start(
        double *                       data_this,
        std::vector<double *> &        data_others,
        dealii::AlignedVector<double> &buffer,
        const unsigned int communication_channel = 0) const override;

      void
      export_to_ghosted_array_finish(
        double *                       data_this,
        std::vector<double *> &        data_others,
        dealii::AlignedVector<double> &buffer) const override;

      void
      import_from_ghosted_array_start(
        double *                       data_this,
        std::vector<double *> &        data_others,
        dealii::AlignedVector<double> &buffer,
        const unsigned int communication_channel = 0) const override;

      void
      import_from_ghosted_array_finish(
        double *                       data_this,
        std::vector<double *> &        data_others,
        dealii::AlignedVector<double> &buffer) const override;

      void
      export_to_ghosted_array_start(
        float *                       data_this,
        std::vector<float *> &        data_others,
        dealii::AlignedVector<float> &buffer,
        const unsigned int            communication_channel = 0) const override;

      void
      export_to_ghosted_array_finish(
        float *                       data_this,
        std::vector<float *> &        data_others,
        dealii::AlignedVector<float> &buffer) const override;

      void
      import_from_ghosted_array_start(
        float *                       data_this,
        std::vector<float *> &        data_others,
        dealii::AlignedVector<float> &buffer,
        const unsigned int            communication_channel = 0) const override;

      void
      import_from_ghosted_array_finish(
        float *                       data_this,
        std::vector<float *> &        data_others,
        dealii::AlignedVector<float> &buffer) const override;

      std::size_t
      local_size() const override;

      std::size_t
      n_ghost_indices() const override;

      std::size_t
      n_mpi_processes() const override;

      std::size_t
      memory_consumption() const override;

      std::vector<unsigned int>
      get_sm_view() const override;

    private:
      template <typename Number>
      void
      export_to_ghosted_array_start_impl(
        Number *                       data_this,
        std::vector<Number *> &        data_others,
        dealii::AlignedVector<Number> &buffer,
        const unsigned int             communication_channel = 0) const;

      template <typename Number>
      void
      export_to_ghosted_array_finish_impl(
        Number *                       data_this,
        std::vector<Number *> &        data_others,
        dealii::AlignedVector<Number> &buffer) const;

      template <typename Number>
      void
      import_from_ghosted_array_start_impl(
        Number *                       data_this,
        std::vector<Number *> &        data_others,
        dealii::AlignedVector<Number> &buffer,
        const unsigned int             communication_channel = 0) const;

      template <typename Number>
      void
      import_from_ghosted_array_finish_impl(
        Number *                       data_this,
        std::vector<Number *> &        data_others,
        dealii::AlignedVector<Number> &buffer) const;

    private:
      const MPI_Comm &comm;
      const MPI_Comm &comm_sm;

      unsigned int n_mpi_processes_;

      unsigned int n_local_elements;
      unsigned int n_ghost_elements;

      std::vector<unsigned int>            recv_remote_ranks;
      std::vector<types::global_dof_index> recv_remote_ptr = {0};
      mutable std::vector<MPI_Request>     recv_remote_req; // TODO: move

      std::vector<unsigned int>        recv_sm_ranks;
      std::vector<unsigned int>        recv_sm_ptr = {0};
      mutable std::vector<MPI_Request> recv_sm_req; // TODO: move
      std::vector<unsigned int>        recv_sm_indices;
      std::vector<unsigned int>        recv_sm_len;
      std::vector<unsigned int>        recv_sm_offset;

      std::vector<unsigned int>        send_remote_ranks;
      std::vector<unsigned int>        send_remote_ptr = {0};
      std::vector<unsigned int>        send_remote_indices;
      std::vector<unsigned int>        send_remote_len;
      std::vector<unsigned int>        send_remote_offset;
      mutable std::vector<MPI_Request> send_remote_req; // TODO: move

      std::vector<unsigned int>        send_sm_ranks;
      std::vector<unsigned int>        send_sm_ptr = {0};
      std::vector<unsigned int>        send_sm_indices;
      std::vector<unsigned int>        send_sm_len;
      mutable std::vector<MPI_Request> send_sm_req; // TODO: move
      std::vector<unsigned int>        send_sm_offset;

      std::vector<unsigned int> sm_view;
    };

  } // end of namespace SharedMPI
} // end of namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#endif