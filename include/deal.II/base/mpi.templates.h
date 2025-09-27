// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mpi_templates_h
#define dealii_mpi_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

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
      template <typename T>
      void
      all_reduce(const MPI_Op             &mpi_op,
                 const ArrayView<const T> &values,
                 const MPI_Comm            mpi_communicator,
                 const ArrayView<T>       &output)
      {
        AssertDimension(values.size(), output.size());
#ifdef DEAL_II_WITH_MPI
        if (job_supports_mpi())
          {
            if constexpr (running_in_debug_mode())
              {
                {
                  const unsigned int rank = this_mpi_process(mpi_communicator);
                  unsigned int       size = values.size();
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
                    Assert(
                      size_min == size_max,
                      ExcMessage(
                        "values has different size across MPI processes."));
                }
              }
            const int ierr =
              MPI_Allreduce(values != output ? values.data() : MPI_IN_PLACE,
                            static_cast<void *>(output.data()),
                            static_cast<int>(values.size()),
                            mpi_type_id_for_type<decltype(*values.data())>,
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
      all_reduce(const MPI_Op                           &mpi_op,
                 const ArrayView<const std::complex<T>> &values,
                 const MPI_Comm                          mpi_communicator,
                 const ArrayView<std::complex<T>>       &output)
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
                            mpi_type_id_for_type<T>,
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
    sum(const T &t, const MPI_Comm mpi_communicator)
    {
      if (mpi_communicator == MPI_COMM_SELF)
        return t;
      else
        {
          T return_value{};
          internal::all_reduce(MPI_SUM,
                               ArrayView<const T>(&t, 1),
                               mpi_communicator,
                               ArrayView<T>(&return_value, 1));
          return return_value;
        }
    }



    template <typename T, typename U>
    void
    sum(const T &values, const MPI_Comm mpi_communicator, U &sums)
    {
      static_assert(std::is_same_v<std::decay_t<T>, std::decay_t<U>>,
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
        const MPI_Comm            mpi_communicator,
        const ArrayView<T>       &sums)
    {
      if (mpi_communicator == MPI_COMM_SELF)
        for (unsigned int i = 0; i < values.size(); ++i)
          sums[i] = values[i];
      else
        internal::all_reduce(MPI_SUM, values, mpi_communicator, sums);
    }



    template <int rank, int dim, typename Number>
    Tensor<rank, dim, Number>
    sum(const Tensor<rank, dim, Number> &t, const MPI_Comm mpi_communicator)
    {
      // Copy the tensor into a C-style array with which we can then
      // call the other sum() function.
      Number array[Tensor<rank, dim, Number>::n_independent_components];
      for (unsigned int i = 0;
           i < Tensor<rank, dim, Number>::n_independent_components;
           ++i)
        array[i] =
          t[Tensor<rank, dim, Number>::unrolled_to_component_indices(i)];

      sum(array, mpi_communicator, array);

      return Tensor<rank, dim, Number>(make_array_view(array));
    }



    template <int rank, int dim, typename Number>
    SymmetricTensor<rank, dim, Number>
    sum(const SymmetricTensor<rank, dim, Number> &local,
        const MPI_Comm                            mpi_communicator)
    {
      // Copy the tensor into a C-style array with which we can then
      // call the other sum() function.
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
        const MPI_Comm              mpi_communicator,
        SparseMatrix<Number>       &global)
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
    max(const T &t, const MPI_Comm mpi_communicator)
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
    max(const T &values, const MPI_Comm mpi_communicator, U &maxima)
    {
      static_assert(std::is_same_v<std::decay_t<T>, std::decay_t<U>>,
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
        const MPI_Comm            mpi_communicator,
        const ArrayView<T>       &maxima)
    {
      internal::all_reduce(MPI_MAX, values, mpi_communicator, maxima);
    }



    template <typename T>
    T
    min(const T &t, const MPI_Comm mpi_communicator)
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
    min(const T &values, const MPI_Comm mpi_communicator, U &minima)
    {
      static_assert(std::is_same_v<std::decay_t<T>, std::decay_t<U>>,
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
        const MPI_Comm            mpi_communicator,
        const ArrayView<T>       &minima)
    {
      internal::all_reduce(MPI_MIN, values, mpi_communicator, minima);
    }



    template <typename T>
    T
    logical_or(const T &t, const MPI_Comm mpi_communicator)
    {
      static_assert(std::is_integral_v<T>,
                    "The MPI_LOR operation only allows integral data types.");

      T return_value{};
      internal::all_reduce(MPI_LOR,
                           ArrayView<const T>(&t, 1),
                           mpi_communicator,
                           ArrayView<T>(&return_value, 1));
      return return_value;
    }



    template <typename T, typename U>
    void
    logical_or(const T &values, const MPI_Comm mpi_communicator, U &results)
    {
      static_assert(std::is_same_v<std::decay_t<T>, std::decay_t<U>>,
                    "Input and output arguments must have the same type!");

      static_assert(std::is_integral_v<typename T::value_type>,
                    "The MPI_LOR operation only allows integral data types.");

      // Specializations of std containers for the data type bool do not
      // necessarily store its elements as a contiguous array. Thus we will use
      // the make_array_view() function with iterators here, which verifies
      // this.
      const auto array_view_values =
        make_array_view(values.cbegin(), values.cend());
      const auto array_view_results =
        make_array_view(results.begin(), results.end());

      using const_type =
        ArrayView<const typename decltype(array_view_values)::value_type>;
      logical_or(static_cast<const_type>(array_view_values),
                 mpi_communicator,
                 array_view_results);
    }



    template <typename T>
    void
    logical_or(const ArrayView<const T> &values,
               const MPI_Comm            mpi_communicator,
               const ArrayView<T>       &results)
    {
      static_assert(std::is_integral_v<T>,
                    "The MPI_LOR operation only allows integral data types.");

      internal::all_reduce(MPI_LOR, values, mpi_communicator, results);
    }



    template <typename T>
    T
    reduce(const T                                      &vec,
           const MPI_Comm                                comm,
           const std::function<T(const T &, const T &)> &combiner,
           const unsigned int                            root_process)
    {
#ifdef DEAL_II_WITH_MPI
      if (n_mpi_processes(comm) > 1)
        {
          // 1) perform custom reduction
          T result = vec;

          const unsigned int rank  = this_mpi_process(comm);
          const unsigned int nproc = n_mpi_processes(comm);

          for (unsigned int stride = 1; stride < nproc; stride *= 2)
            {
              unsigned int rank_recv =
                (2 * stride) *
                  ((rank + nproc - root_process) % nproc / (2 * stride)) +
                root_process;
              unsigned int rank_send = rank_recv + stride;

              if (rank_send >= nproc + root_process) // nothing to do
                continue;

              rank_recv = rank_recv % nproc;
              rank_send = rank_send % nproc;

              if (rank_recv == rank) // process receives data
                {
                  MPI_Status status;
                  const auto ierr_1 = MPI_Probe(rank_send,
                                                internal::Tags::compute_union,
                                                comm,
                                                &status);
                  AssertThrowMPI(ierr_1);

                  int        amount;
                  const auto ierr_2 = MPI_Get_count(&status, MPI_CHAR, &amount);
                  AssertThrowMPI(ierr_2);

                  std::vector<char> temp(amount);

                  const auto ierr_3 = MPI_Recv(temp.data(),
                                               amount,
                                               MPI_CHAR,
                                               status.MPI_SOURCE,
                                               status.MPI_TAG,
                                               comm,
                                               MPI_STATUS_IGNORE);
                  AssertThrowMPI(ierr_3);

                  result = combiner(result, Utilities::unpack<T>(temp, false));
                }
              else if (rank_send == rank) // process sends data
                {
                  const std::vector<char> temp = Utilities::pack(result, false);

                  const auto ierr = MPI_Send(temp.data(),
                                             temp.size(),
                                             MPI_CHAR,
                                             rank_recv,
                                             internal::Tags::compute_union,
                                             comm);
                  AssertThrowMPI(ierr);
                }
            }

          if (rank == root_process)
            return result;
          else
            return {};
        }
#endif
      (void)comm;
      (void)combiner;
      (void)root_process;
      return vec;
    }



    template <typename T>
    T
    all_reduce(const T                                      &vec,
               const MPI_Comm                                comm,
               const std::function<T(const T &, const T &)> &combiner)
    {
      if (n_mpi_processes(comm) > 1)
        {
          // 1) perform reduction
          const auto result = Utilities::MPI::reduce<T>(vec, comm, combiner);

          // 2) broadcast result
          return Utilities::MPI::broadcast(comm, result);
        }
      else
        return vec;
    }


    template <typename T>
    std::vector<T>
    compute_set_union(const std::vector<T> &vec, const MPI_Comm comm)
    {
      return Utilities::MPI::all_reduce<std::vector<T>>(
        vec, comm, [](const auto &set_1, const auto &set_2) {
          // merge vectors, sort, and eliminate duplicates -> mimic std::set<T>
          auto result = set_1;
          result.insert(result.end(), set_2.begin(), set_2.end());
          std::sort(result.begin(), result.end());
          result.erase(std::unique(result.begin(), result.end()), result.end());
          return result;
        });
    }



    template <typename T>
    std::set<T>
    compute_set_union(const std::set<T> &set_in, const MPI_Comm comm)
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
