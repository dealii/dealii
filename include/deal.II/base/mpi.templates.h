// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2017 by the deal.II authors
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

      template<typename T>
      std::map<unsigned int, std::vector<T> >
      send_and_receive(const MPI_Comm                                         &comm,
                       const std::map<unsigned int, typename std::vector<T> > &objects_to_send)
      {
#ifndef DEAL_II_WITH_MPI
        (void)comm;
        (void)objects_to_send;
        Assert (false, ExcMessage ("The function send_and_receive doesn't make"
                                   "any sense if you do not have MPI enabled. "));
#else
        auto n_procs = dealii::Utilities::MPI::n_mpi_processes(comm);
        auto my_proc = dealii::Utilities::MPI::this_mpi_process(comm);

        // Collective communication telling each process
        // how many objects I'm sending: this value shall be contained
        // inside send_amount
        std::vector<unsigned int> send_amount(n_procs, 0);
        // Vector to store, in each entry, the number of objects the corresponding
        // process shall send to this process
        std::vector<unsigned int> rec_amount(n_procs, 0);

        for (const auto &m: objects_to_send)
          send_amount[m.first] = m.second.size();

        MPI_Alltoall(
          &send_amount[0], 1, MPI_UNSIGNED,
          &rec_amount[0], 1, MPI_UNSIGNED, comm);
        unsigned int number_p2p_comm = 0;
        for (unsigned int rank=0; rank< n_procs; ++rank)
          if (rec_amount[rank] != 0)
            ++number_p2p_comm;
        number_p2p_comm += objects_to_send.size();
        // Creating the MPI data type to send the object
        unsigned int size_of_T = sizeof(T);
        MPI_Datatype ptype;
        MPI_Type_contiguous(size_of_T,MPI_CHAR,&ptype);
        MPI_Type_commit(&ptype);
        std::map<unsigned int, std::vector<T> > received_objects;
        std::vector<MPI_Request> requests(number_p2p_comm);
        unsigned int req_n = 0; // index running on the request vector
        // Sending part
        for (const auto &rank_vec : objects_to_send)
          {
            const auto &rank = rank_vec.first;
            const auto &vec  = rank_vec.second;
            MPI_Isend(&(vec[0]), vec.size(), ptype,
                      rank, n_procs*rank+my_proc, comm, &requests[req_n++]);
          }

        // Receiving part
        for (unsigned int rank=0; rank < n_procs; ++rank)
          if (rec_amount[rank] != 0)
            {
              received_objects[rank].resize(rec_amount[rank]);
              auto &vec = received_objects[rank];
              MPI_Irecv(&(vec[0]), rec_amount[rank], ptype,
                        rank, n_procs*my_proc+rank, comm, &requests[req_n++]);
            }

        MPI_Waitall(number_p2p_comm,&requests[0],MPI_STATUSES_IGNORE);
        return received_objects;
#endif // deal.II with MPI
      }



      template<typename T>
      std::vector<std::vector<T> >
      send_and_receive(const MPI_Comm       &comm,
                       const std::vector<T> &objects_to_send)
      {
#ifndef DEAL_II_WITH_MPI
        (void)comm;
        (void)objects_to_send;
        Assert (false, ExcMessage ("The function send_and_receive doesn't make"
                                   "any sense if you do not have MPI enabled. "));
#else
        auto n_procs = dealii::Utilities::MPI::n_mpi_processes(comm);

        std::vector<unsigned int> rec_amount(n_procs, 0);
        const int n_local_data = objects_to_send.size();

        unsigned int size_of_T = sizeof(T);
        MPI_Datatype ptype;
        MPI_Type_contiguous(size_of_T,MPI_CHAR,&ptype);
        MPI_Type_commit(&ptype);

        // Vector to store the size of loc_data_array for every process
        std::vector<int> size_all_data(n_procs,0);

        // Exchanging the number of bboxes
        MPI_Allgather(&n_local_data, 1, MPI_INT,
                      &(size_all_data[0]), 1, MPI_INT,
                      comm);

        // Now computing the the displacement, relative to recvbuf,
        // at which to store the incoming data
        std::vector<int> rdispls(n_procs);
        rdispls[0] = 0;
        for (unsigned int i=1; i < n_procs; ++i)
          rdispls[i] = rdispls[i-1] + size_all_data[i-1];

        // Step 3: exchange the data and bounding boxes:
        // Allocating a vector to contain all the received data
        std::vector<T> data_array(rdispls.back() + size_all_data.back());

        MPI_Allgatherv(&(objects_to_send[0]), n_local_data, ptype,
                       &(data_array[0]), &(size_all_data[0]),
                       &(rdispls[0]), ptype, comm);

        std::vector<std::vector<T> > received_objects(n_procs);
        for (unsigned int i= 0; i < n_procs; ++i)
          {
            received_objects[i].insert(received_objects[i].begin(),
                                       data_array.begin()+rdispls[i],
                                       data_array.begin()+rdispls[i]+size_all_data[i]);
          }

        return received_objects;
#endif
      }
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



    template <typename T>
    void sum (const std::vector<T> &values,
              const MPI_Comm       &mpi_communicator,
              std::vector<T>       &sums)
    {
      internal::all_reduce(MPI_SUM, ArrayView<const T>(values),
                           mpi_communicator, ArrayView<T>(sums));
    }



    template <typename T>
    void sum (const ArrayView<const T> &values,
              const MPI_Comm           &mpi_communicator,
              const ArrayView<T>       &sums)
    {
      internal::all_reduce(MPI_SUM, values, mpi_communicator, sums);
    }



    template <typename T>
    void sum (const Vector<T> &values,
              const MPI_Comm &mpi_communicator,
              Vector<T> &sums)
    {
      const auto &size = values.size();
      internal::all_reduce(MPI_SUM, ArrayView<const T>(values.begin(), size),
                           mpi_communicator, ArrayView<T>(sums.begin(), size));
    }



    template <typename T>
    void sum (const FullMatrix<T> &values,
              const MPI_Comm &mpi_communicator,
              FullMatrix<T> &sums)
    {
      const auto size_values = values.n()*values.m();
      const auto size_sums = sums.n()*sums.m();
      internal::all_reduce(MPI_SUM, ArrayView<const T>(&values[0][0], size_values),
                           mpi_communicator, ArrayView<T>(&sums[0][0], size_sums));
    }



    template <typename T>
    void sum (const LAPACKFullMatrix<T> &values,
              const MPI_Comm &mpi_communicator,
              LAPACKFullMatrix<T> &sums)
    {
      const auto size_values = values.n()*values.m();
      const auto size_sums = sums.n()*sums.m();
      internal::all_reduce(MPI_SUM, ArrayView<const T>(&values(0,0), size_values),
                           mpi_communicator, ArrayView<T>(&sums(0,0), size_sums));
    }



    template <int rank, int dim, typename Number>
    Tensor<rank,dim,Number>
    sum (const Tensor<rank,dim,Number> &local,
         const MPI_Comm &mpi_communicator)
    {
      const unsigned int n_entries = Tensor<rank,dim,Number>::n_independent_components;
      Number entries[ Tensor<rank,dim,Number>::n_independent_components ];

      for (unsigned int i=0; i< n_entries; ++i)
        entries[i] = local[ local.unrolled_to_component_indices(i) ];

      Number global_entries[ Tensor<rank,dim,Number>::n_independent_components ];
      dealii::Utilities::MPI::sum( entries, mpi_communicator, global_entries );

      Tensor<rank,dim,Number> global;
      for (unsigned int i=0; i< n_entries; ++i)
        global[ global.unrolled_to_component_indices(i) ] = global_entries[i];

      return global;
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



    template <typename T>
    void max (const std::vector<T> &values,
              const MPI_Comm       &mpi_communicator,
              std::vector<T>       &maxima)
    {
      internal::all_reduce(MPI_MAX, ArrayView<const T>(values),
                           mpi_communicator, ArrayView<T> (maxima));
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



    template <typename T>
    void min (const std::vector<T> &values,
              const MPI_Comm       &mpi_communicator,
              std::vector<T>       &minima)
    {
      internal::all_reduce(MPI_MIN, ArrayView<const T>(values),
                           mpi_communicator, ArrayView<T> (minima));
    }
  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
