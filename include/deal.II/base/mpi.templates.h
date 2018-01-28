// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_mpi_templates_h
#define dealii_mpi_templates_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    namespace internal
    {
#ifdef DEAL_II_WITH_MPI
      /**
       * Return the corresponding MPI data type id for the argument given.
       */
      inline MPI_Datatype mpi_type_id (const int *)
      {
        return MPI_INT;
      }



      inline MPI_Datatype mpi_type_id (const long int *)
      {
        return MPI_LONG;
      }



      inline MPI_Datatype mpi_type_id (const unsigned int *)
      {
        return MPI_UNSIGNED;
      }



      inline MPI_Datatype mpi_type_id (const unsigned long int *)
      {
        return MPI_UNSIGNED_LONG;
      }



      inline MPI_Datatype mpi_type_id (const unsigned long long int *)
      {
        return MPI_UNSIGNED_LONG_LONG;
      }



      inline MPI_Datatype mpi_type_id (const float *)
      {
        return MPI_FLOAT;
      }



      inline MPI_Datatype mpi_type_id (const double *)
      {
        return MPI_DOUBLE;
      }



      inline MPI_Datatype mpi_type_id (const long double *)
      {
        return MPI_LONG_DOUBLE;
      }



      inline MPI_Datatype mpi_type_id (const std::complex<float> *)
      {
        return MPI_COMPLEX;
      }



      inline MPI_Datatype mpi_type_id (const std::complex<double> *)
      {
        return MPI_DOUBLE_COMPLEX;
      }
#endif


      template <typename T>
      void all_reduce (const MPI_Op             &mpi_op,
                       const ArrayView<const T> &values,
                       const MPI_Comm           &mpi_communicator,
                       const ArrayView<T>       &output)
      {
        AssertDimension(values.size(), output.size());
#ifdef DEAL_II_WITH_MPI
        if (job_supports_mpi())
          {
            const int ierr = MPI_Allreduce
                             (values != output
                              ?
                              // TODO This const_cast is only needed for older
                              // (e.g., openMPI 1.6, released in 2012)
                              // implementations of MPI-2. It is not needed as
                              // of MPI-3 and we should remove it at some
                              // point in the future.
                              const_cast<void *>(static_cast<const void *>(values.data()))
                              :
                              MPI_IN_PLACE,
                              static_cast<void *>(output.data()),
                              static_cast<int>(values.size()),
                              internal::mpi_type_id(values.data()),
                              mpi_op,
                              mpi_communicator);
            AssertThrowMPI(ierr);
          }
        else
#endif
          {
            (void)mpi_op;
            (void)mpi_communicator;
            if (values != output)
              std::copy(values.begin(), values.end(), output.begin());
          }
      }



      template <typename T>
      void all_reduce (const MPI_Op                           &mpi_op,
                       const ArrayView<const std::complex<T>> &values,
                       const MPI_Comm                         &mpi_communicator,
                       const ArrayView<std::complex<T>>       &output)
      {
        AssertDimension(values.size(), output.size());
#ifdef DEAL_II_WITH_MPI
        if (job_supports_mpi())
          {
            const int ierr = MPI_Allreduce
                             (values != output
                              ?
                              // TODO This const_cast is only needed for older
                              // (e.g., openMPI 1.6, released in 2012)
                              // implementations of MPI-2. It is not needed as
                              // of MPI-3 and we should remove it at some
                              // point in the future.
                              const_cast<void *>(static_cast<const void *>(values.data()))
                              :
                              MPI_IN_PLACE,
                              static_cast<void *>(output.data()),
                              static_cast<int>(values.size()*2),
                              internal::mpi_type_id(static_cast<T *>(nullptr)),
                              mpi_op,
                              mpi_communicator);
            AssertThrowMPI(ierr);
          }
        else
#endif
          {
            (void)mpi_op;
            (void)mpi_communicator;
            if (values != output)
              std::copy(values.begin(), values.end(), output.begin());
          }
      }
    }



    template<typename T>
    std::map<unsigned int, T>
    some_to_some(const MPI_Comm &comm,
                 const std::map<unsigned int, T> &objects_to_send)
    {
#ifndef DEAL_II_WITH_MPI
      (void)comm;
      Assert(objects_to_send.size() == 0, ExcMessage("Cannot send to more than one processor."));
      Assert(objects_to_send.find(0) != objects_to_send.end() || objects_to_send.size() == 0,
             ExcMessage("Can only send to myself or to nobody."));
      return objects_to_send;
#else

      std::vector<unsigned int> send_to(objects_to_send.size());
      {
        unsigned int i=0;
        for (const auto &m: objects_to_send)
          send_to[i++] = m.first;
      }
      AssertDimension(send_to.size(), objects_to_send.size());

      const auto receive_from =
        Utilities::MPI::compute_point_to_point_communication_pattern(comm, send_to);

      // Sending buffers
      std::vector<std::vector<char> > buffers_to_send(send_to.size());
      std::vector<MPI_Request> buffer_send_requests(send_to.size());
      {
        unsigned int i = 0;
        for (const auto &rank_obj : objects_to_send)
          {
            const auto &rank = rank_obj.first;
            buffers_to_send[i] = Utilities::pack(rank_obj.second);
            const int ierr = MPI_Isend(buffers_to_send[i].data(),
                                       buffers_to_send[i].size(), MPI_CHAR,
                                       rank, 21, comm, &buffer_send_requests[i]);
            AssertThrowMPI(ierr);
            ++i;
          }
      }

      // Receiving buffers
      std::map<unsigned int, T> received_objects;
      {
        std::vector<char> buffer;
        // We do this on a first come/first served basis
        for (unsigned int i = 0; i<receive_from.size(); ++i)
          {
            // Probe what's going on. Take data from the first available sender
            MPI_Status status;
            int ierr = MPI_Probe(MPI_ANY_SOURCE, 21, comm, &status);
            AssertThrowMPI(ierr);

            // Length of the message
            int len;
            ierr = MPI_Get_count(&status, MPI_CHAR, &len);
            AssertThrowMPI(ierr);
            buffer.resize(len);

            // Source rank
            const unsigned int rank = status.MPI_SOURCE;

            // Actually receive the message
            ierr = MPI_Recv(buffer.data(), len, MPI_CHAR,
                            rank, 21, comm, MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);
            Assert(received_objects.find(rank) == received_objects.end(),
                   ExcInternalError("I should not receive again from this rank"));
            received_objects[rank] = Utilities::unpack<T>(buffer);
          }
      }

      // Wait to have sent all objects.
      MPI_Waitall(send_to.size(), buffer_send_requests.data(),MPI_STATUSES_IGNORE);

      return received_objects;
#endif // deal.II with MPI
    }



    template<typename T>
    std::vector<T> all_gather(const MPI_Comm       &comm,
                              const T &object)
    {
#ifndef DEAL_II_WITH_MPI
      (void)comm;
      std::vector<T> v(1, object);
      return v;
#else
      const auto n_procs = dealii::Utilities::MPI::n_mpi_processes(comm);

      std::vector<char> buffer = Utilities::pack(object);

      int n_local_data = buffer.size();

      // Vector to store the size of loc_data_array for every process
      std::vector<int> size_all_data(n_procs,0);

      // Exchanging the size of each buffer
      MPI_Allgather(&n_local_data, 1, MPI_INT,
                    &(size_all_data[0]), 1, MPI_INT,
                    comm);

      // Now computing the the displacement, relative to recvbuf,
      // at which to store the incoming buffer
      std::vector<int> rdispls(n_procs);
      rdispls[0] = 0;
      for (unsigned int i=1; i < n_procs; ++i)
        rdispls[i] = rdispls[i-1] + size_all_data[i-1];

      // Step 3: exchange the buffer:
      std::vector<char> received_unrolled_buffer(rdispls.back() + size_all_data.back());

      MPI_Allgatherv(buffer.data(), n_local_data, MPI_CHAR,
                     received_unrolled_buffer.data(), size_all_data.data(),
                     rdispls.data(), MPI_CHAR, comm);

      std::vector<T> received_objects(n_procs);
      for (unsigned int i= 0; i < n_procs; ++i)
        {
          std::vector<char> local_buffer(received_unrolled_buffer.begin()+rdispls[i],
                                         received_unrolled_buffer.begin()+rdispls[i]+size_all_data[i]);
          received_objects[i] = Utilities::unpack<T>(local_buffer);
        }

      return received_objects;
#endif
    }



    template <typename T>
    T sum (const T &t,
           const MPI_Comm &mpi_communicator)
    {
      T return_value;
      internal::all_reduce(MPI_SUM, ArrayView<const T>(&t,1),
                           mpi_communicator, ArrayView<T>(&return_value,1));
      return return_value;
    }



    template <typename T, typename U>
    void sum (const T        &values,
              const MPI_Comm &mpi_communicator,
              U              &sums)
    {
      static_assert(std::is_same<typename std::decay<T>::type,
                    typename std::decay<U>::type>::value,
                    "Input and output arguments must have the same type!");
      const auto array_view_values = make_array_view(values);
      using const_type = ArrayView<const typename decltype(array_view_values)::value_type>;
      sum(static_cast<const_type> (array_view_values),
          mpi_communicator,
          make_array_view(sums));
    }



    template <typename T>
    void sum (const ArrayView<const T> &values,
              const MPI_Comm           &mpi_communicator,
              const ArrayView<T>       &sums)
    {
      internal::all_reduce(MPI_SUM, values, mpi_communicator, sums);
    }



    template <int rank, int dim, typename Number>
    Tensor<rank,dim,Number>
    sum (const Tensor<rank,dim,Number> &local,
         const MPI_Comm &mpi_communicator)
    {
      Tensor<rank, dim, Number> sums;
      dealii::Utilities::MPI::sum(local, mpi_communicator, sums);
      return sums;
    }



    template <int rank, int dim, typename Number>
    SymmetricTensor<rank,dim,Number>
    sum (const SymmetricTensor<rank,dim,Number> &local,
         const MPI_Comm &mpi_communicator)
    {
      const unsigned int n_entries = SymmetricTensor<rank,dim,Number>::n_independent_components;
      Number entries[ SymmetricTensor<rank,dim,Number>::n_independent_components ];

      for (unsigned int i=0; i< n_entries; ++i)
        entries[i] = local[ local.unrolled_to_component_indices(i) ];

      Number global_entries[ SymmetricTensor<rank,dim,Number>::n_independent_components ];
      dealii::Utilities::MPI::sum( entries, mpi_communicator, global_entries );

      SymmetricTensor<rank,dim,Number> global;
      for (unsigned int i=0; i< n_entries; ++i)
        global[ global.unrolled_to_component_indices(i) ] = global_entries[i];

      return global;
    }



    template <typename T>
    T max (const T &t,
           const MPI_Comm &mpi_communicator)
    {
      T return_value;
      internal::all_reduce(MPI_MAX, ArrayView<const T>(&t,1),
                           mpi_communicator, ArrayView<T> (&return_value,1));
      return return_value;
    }



    template <typename T, typename U>
    void max (const T         &values,
              const MPI_Comm  &mpi_communicator,
              U               &maxima)
    {
      static_assert(std::is_same<typename std::decay<T>::type,
                    typename std::decay<U>::type>::value,
                    "Input and output arguments must have the same type!");
      const auto array_view_values = make_array_view(values);
      using const_type = ArrayView<const typename decltype(array_view_values)::value_type>;
      max(static_cast<const_type> (array_view_values),
          mpi_communicator,
          make_array_view(maxima));
    }



    template <typename T>
    void max (const ArrayView<const T> &values,
              const MPI_Comm           &mpi_communicator,
              const ArrayView<T>       &maxima)
    {
      internal::all_reduce(MPI_MAX, values, mpi_communicator, maxima);
    }



    template <typename T>
    T min (const T &t,
           const MPI_Comm &mpi_communicator)
    {
      T return_value;
      internal::all_reduce(MPI_MIN, ArrayView<const T>(&t,1),
                           mpi_communicator, ArrayView<T>(&return_value,1));
      return return_value;
    }



    template <typename T, typename U>
    void min (const T        &values,
              const MPI_Comm &mpi_communicator,
              U              &minima)
    {
      static_assert(std::is_same<typename std::decay<T>::type,
                    typename std::decay<U>::type>::value,
                    "Input and output arguments must have the same type!");
      const auto array_view_values = make_array_view(values);
      using const_type = ArrayView<const typename decltype(array_view_values)::value_type>;
      min(static_cast<const_type> (array_view_values),
          mpi_communicator,
          make_array_view(minima));
    }



    template <typename T>
    void min (const ArrayView<const T> &values,
              const MPI_Comm           &mpi_communicator,
              const ArrayView<T>       &minima)
    {
      internal::all_reduce(MPI_MIN, values, mpi_communicator, minima);
    }
  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
