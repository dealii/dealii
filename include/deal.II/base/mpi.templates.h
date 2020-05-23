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

#ifndef dealii_mpi_templates_h
#define dealii_mpi_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <set>
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
      inline MPI_Datatype
      mpi_type_id(const int *)
      {
        return MPI_INT;
      }



      inline MPI_Datatype
      mpi_type_id(const long int *)
      {
        return MPI_LONG;
      }



      inline MPI_Datatype
      mpi_type_id(const unsigned int *)
      {
        return MPI_UNSIGNED;
      }



      inline MPI_Datatype
      mpi_type_id(const unsigned long int *)
      {
        return MPI_UNSIGNED_LONG;
      }



      inline MPI_Datatype
      mpi_type_id(const unsigned long long int *)
      {
        return MPI_UNSIGNED_LONG_LONG;
      }



      inline MPI_Datatype
      mpi_type_id(const float *)
      {
        return MPI_FLOAT;
      }



      inline MPI_Datatype
      mpi_type_id(const double *)
      {
        return MPI_DOUBLE;
      }



      inline MPI_Datatype
      mpi_type_id(const long double *)
      {
        return MPI_LONG_DOUBLE;
      }



      inline MPI_Datatype
      mpi_type_id(const std::complex<float> *)
      {
        return MPI_COMPLEX;
      }



      inline MPI_Datatype
      mpi_type_id(const std::complex<double> *)
      {
        return MPI_DOUBLE_COMPLEX;
      }
#endif


      template <typename T>
      void
      all_reduce(const MPI_Op &            mpi_op,
                 const ArrayView<const T> &values,
                 const MPI_Comm &          mpi_communicator,
                 const ArrayView<T> &      output)
      {
        AssertDimension(values.size(), output.size());
#ifdef DEAL_II_WITH_MPI
        if (job_supports_mpi())
          {
#  ifdef DEBUG
            {
              const unsigned int rank     = this_mpi_process(mpi_communicator);
              unsigned int       size     = values.size();
              unsigned int       size_min = 0;
              unsigned int       size_max = 0;
              int                ierr2    = 0;
              ierr2                       = MPI_Reduce(&size,
                                 &size_min,
                                 1,
                                 MPI_UNSIGNED,
                                 MPI_MIN,
                                 0,
                                 mpi_communicator);
              AssertThrowMPI(ierr2);
              ierr2 = MPI_Reduce(&size,
                                 &size_max,
                                 1,
                                 MPI_UNSIGNED,
                                 MPI_MAX,
                                 0,
                                 mpi_communicator);
              AssertThrowMPI(ierr2);
              if (rank == 0)
                Assert(size_min == size_max,
                       ExcMessage(
                         "values has different size across MPI processes."));
            }
#  endif
            const int ierr =
              MPI_Allreduce(values != output ?
                              DEAL_II_MPI_CONST_CAST(values.data()) :
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
      void
      all_reduce(const MPI_Op &                          mpi_op,
                 const ArrayView<const std::complex<T>> &values,
                 const MPI_Comm &                        mpi_communicator,
                 const ArrayView<std::complex<T>> &      output)
      {
        AssertDimension(values.size(), output.size());
#ifdef DEAL_II_WITH_MPI
        if (job_supports_mpi())
          {
            const int ierr =
              MPI_Allreduce(values != output ?
                              // TODO This const_cast is only needed for older
                              // (e.g., openMPI 1.6, released in 2012)
                              // implementations of MPI-2. It is not needed as
                              // of MPI-3 and we should remove it at some
                              // point in the future.
                              const_cast<void *>(
                                static_cast<const void *>(values.data())) :
                              MPI_IN_PLACE,
                            static_cast<void *>(output.data()),
                            static_cast<int>(values.size() * 2),
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
    } // namespace internal



    template <typename T>
    T
    sum(const T &t, const MPI_Comm &mpi_communicator)
    {
      T return_value{};
      internal::all_reduce(MPI_SUM,
                           ArrayView<const T>(&t, 1),
                           mpi_communicator,
                           ArrayView<T>(&return_value, 1));
      return return_value;
    }



    template <typename T, typename U>
    void
    sum(const T &values, const MPI_Comm &mpi_communicator, U &sums)
    {
      static_assert(std::is_same<typename std::decay<T>::type,
                                 typename std::decay<U>::type>::value,
                    "Input and output arguments must have the same type!");
      const auto array_view_values = make_array_view(values);
      using const_type =
        ArrayView<const typename decltype(array_view_values)::value_type>;
      sum(static_cast<const_type>(array_view_values),
          mpi_communicator,
          make_array_view(sums));
    }



    template <typename T>
    void
    sum(const ArrayView<const T> &values,
        const MPI_Comm &          mpi_communicator,
        const ArrayView<T> &      sums)
    {
      internal::all_reduce(MPI_SUM, values, mpi_communicator, sums);
    }



    template <int rank, int dim, typename Number>
    Tensor<rank, dim, Number>
    sum(const Tensor<rank, dim, Number> &local,
        const MPI_Comm &                 mpi_communicator)
    {
      Tensor<rank, dim, Number> sums;
      sum(local, mpi_communicator, sums);
      return sums;
    }



    template <int rank, int dim, typename Number>
    SymmetricTensor<rank, dim, Number>
    sum(const SymmetricTensor<rank, dim, Number> &local,
        const MPI_Comm &                          mpi_communicator)
    {
      const unsigned int n_entries =
        SymmetricTensor<rank, dim, Number>::n_independent_components;
      Number
        entries[SymmetricTensor<rank, dim, Number>::n_independent_components];

      for (unsigned int i = 0; i < n_entries; ++i)
        entries[i] = local[local.unrolled_to_component_indices(i)];

      Number global_entries
        [SymmetricTensor<rank, dim, Number>::n_independent_components];
      sum(entries, mpi_communicator, global_entries);

      SymmetricTensor<rank, dim, Number> global;
      for (unsigned int i = 0; i < n_entries; ++i)
        global[global.unrolled_to_component_indices(i)] = global_entries[i];

      return global;
    }



    template <typename Number>
    void
    sum(const SparseMatrix<Number> &local,
        const MPI_Comm &            mpi_communicator,
        SparseMatrix<Number> &      global)
    {
      Assert(
        local.get_sparsity_pattern() == global.get_sparsity_pattern(),
        ExcMessage(
          "The sparsity pattern of the local and the global matrices should match."));
#ifdef DEAL_II_WITH_MPI
      // makes use of the fact that the matrix stores its data in a
      // contiguous array.
      sum(ArrayView<const Number>(local.val.get(), local.n_nonzero_elements()),
          mpi_communicator,
          ArrayView<Number>(global.val.get(), global.n_nonzero_elements()));
#else
      (void)mpi_communicator;
      if (!PointerComparison::equal(&local, &global))
        global = local;
#endif
    }



    template <typename T>
    T
    max(const T &t, const MPI_Comm &mpi_communicator)
    {
      T return_value{};
      internal::all_reduce(MPI_MAX,
                           ArrayView<const T>(&t, 1),
                           mpi_communicator,
                           ArrayView<T>(&return_value, 1));
      return return_value;
    }



    template <typename T, typename U>
    void
    max(const T &values, const MPI_Comm &mpi_communicator, U &maxima)
    {
      static_assert(std::is_same<typename std::decay<T>::type,
                                 typename std::decay<U>::type>::value,
                    "Input and output arguments must have the same type!");
      const auto array_view_values = make_array_view(values);
      using const_type =
        ArrayView<const typename decltype(array_view_values)::value_type>;
      max(static_cast<const_type>(array_view_values),
          mpi_communicator,
          make_array_view(maxima));
    }



    template <typename T>
    void
    max(const ArrayView<const T> &values,
        const MPI_Comm &          mpi_communicator,
        const ArrayView<T> &      maxima)
    {
      internal::all_reduce(MPI_MAX, values, mpi_communicator, maxima);
    }



    template <typename T>
    T
    min(const T &t, const MPI_Comm &mpi_communicator)
    {
      T return_value{};
      internal::all_reduce(MPI_MIN,
                           ArrayView<const T>(&t, 1),
                           mpi_communicator,
                           ArrayView<T>(&return_value, 1));
      return return_value;
    }



    template <typename T, typename U>
    void
    min(const T &values, const MPI_Comm &mpi_communicator, U &minima)
    {
      static_assert(std::is_same<typename std::decay<T>::type,
                                 typename std::decay<U>::type>::value,
                    "Input and output arguments must have the same type!");
      const auto array_view_values = make_array_view(values);
      using const_type =
        ArrayView<const typename decltype(array_view_values)::value_type>;
      min(static_cast<const_type>(array_view_values),
          mpi_communicator,
          make_array_view(minima));
    }



    template <typename T>
    void
    min(const ArrayView<const T> &values,
        const MPI_Comm &          mpi_communicator,
        const ArrayView<T> &      minima)
    {
      internal::all_reduce(MPI_MIN, values, mpi_communicator, minima);
    }



    template <typename T>
    std::vector<T>
    compute_set_union(const std::vector<T> &vec, const MPI_Comm &comm)
    {
#ifdef DEAL_II_WITH_MPI
      // 1) collect vector entries and create union
      std::vector<T> result = vec;

      if (this_mpi_process(comm) == 0)
        {
          for (unsigned int i = 1; i < n_mpi_processes(comm); i++)
            {
              MPI_Status status;
              const auto ierr_1 = MPI_Probe(MPI_ANY_SOURCE,
                                            internal::Tags::compute_union,
                                            comm,
                                            &status);
              AssertThrowMPI(ierr_1);

              int        amount;
              const auto ierr_2 =
                MPI_Get_count(&status,
                              internal::mpi_type_id(vec.data()),
                              &amount);
              AssertThrowMPI(ierr_2);

              const unsigned int old_size = result.size();
              result.resize(old_size + amount);

              const auto ierr_3 = MPI_Recv(result.data() + old_size,
                                           amount,
                                           internal::mpi_type_id(vec.data()),
                                           status.MPI_SOURCE,
                                           status.MPI_TAG,
                                           comm,
                                           MPI_STATUS_IGNORE);
              AssertThrowMPI(ierr_3);

              std::sort(result.begin(), result.end());
              result.erase(std::unique(result.begin(), result.end()),
                           result.end());
            }
        }
      else
        {
          const auto ierr = MPI_Send(vec.data(),
                                     vec.size(),
                                     internal::mpi_type_id(vec.data()),
                                     0,
                                     internal::Tags::compute_union,
                                     comm);
          AssertThrowMPI(ierr);
        }

      // 2) broadcast result
      int size = result.size();
      MPI_Bcast(&size, 1, MPI_INT, 0, comm);
      result.resize(size);
      MPI_Bcast(
        result.data(), size, internal::mpi_type_id(vec.data()), 0, comm);

      return result;
#else
      (void)comm;
      return vec;
#endif
    }



    template <typename T>
    std::set<T>
    compute_set_union(const std::set<T> &set_in, const MPI_Comm &comm)
    {
      // convert vector to set
      std::vector<T> vector_in(set_in.begin(), set_in.end());

      // perform operation to vector
      const std::vector<T> vector_out = compute_set_union(vector_in, comm);

      // convert vector to set
      return std::set<T>(vector_out.begin(), vector_out.end());
    }

  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
