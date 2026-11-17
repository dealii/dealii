// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/base/memory_space.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/solver_gmres.h>

#include <Kokkos_Core.hpp>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace SolverGMRESImplementation
  {
    /**
     * Value type for the team-level reduction of a block of @p block_size
     * inner products in do_Tvmult_add_device(), plus the same number of
     * correction terms for the delayed reorthogonalization variant and the
     * norm of the new vector.
     */
    template <int block_size, bool delayed_reorthogonalization>
    struct BlockSums
    {
      double dot[block_size];
      double correct[delayed_reorthogonalization ? block_size : 1];
      double norm_sqr;

      KOKKOS_INLINE_FUNCTION
      BlockSums()
      {
        for (int b = 0; b < block_size; ++b)
          dot[b] = 0.;
        for (int b = 0; b < (delayed_reorthogonalization ? block_size : 1); ++b)
          correct[b] = 0.;
        norm_sqr = 0.;
      }

      KOKKOS_INLINE_FUNCTION BlockSums &
      operator+=(const BlockSums &other)
      {
        for (int b = 0; b < block_size; ++b)
          dot[b] += other.dot[b];
        if constexpr (delayed_reorthogonalization)
          for (int b = 0; b < block_size; ++b)
            correct[b] += other.correct[b];
        norm_sqr += other.norm_sqr;
        return *this;
      }
    };
  } // namespace SolverGMRESImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

namespace Kokkos
{
  template <int block_size, bool delayed_reorthogonalization>
  struct reduction_identity<
    dealii::internal::SolverGMRESImplementation::
      BlockSums<block_size, delayed_reorthogonalization>>
  {
    KOKKOS_FORCEINLINE_FUNCTION static dealii::internal::
      SolverGMRESImplementation::BlockSums<block_size,
                                           delayed_reorthogonalization>
      sum()
    {
      return {};
    }
  };
} // namespace Kokkos

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace SolverGMRESImplementation
  {
    template <bool delayed_reorthogonalization, typename Number>
    void
    do_Tvmult_add_host(const unsigned int             n_vectors,
                       const std::size_t              locally_owned_size,
                       const Number                  *current_vector,
                       OrthogonalizationData<Number> &temp_data,
                       Vector<double>                &h)
    {
      unsigned int j                  = 0;
      const auto  &orthogonal_vectors = temp_data.vector_ptrs;

      if (n_vectors <= 128)
        {
          // optimized path
          static constexpr unsigned int n_lanes =
            VectorizedArray<double>::size();

          VectorizedArray<double> hs[128];
          for (unsigned int i = 0; i < n_vectors; ++i)
            hs[i] = 0.0;
          VectorizedArray<double>
            correct[delayed_reorthogonalization ? 129 : 1];
          if (delayed_reorthogonalization)
            for (unsigned int i = 0; i < n_vectors + 1; ++i)
              correct[i] = 0.0;

          unsigned int c = 0;

          constexpr unsigned int inner_batch_size =
            delayed_reorthogonalization ? 6 : 12;

          const unsigned int loop_length_c =
            locally_owned_size / n_lanes / inner_batch_size;

          // At this point, we would like to run a loop over the variable 'c'
          // from 0 all the way to loop_length_c, and then run a loop over 'i'
          // to work on all vectors we want to compute the inner product
          // against. In other words, the loop layout would be:
          //
          // for (unsigned int c = 0; c < loop_length_c; ++c)
          //   {
          //     // do some work that only depends on the index c or the
          //     // derived index j = c * n_lanes * inner_batch_size
          //     ...
          //     for (unsigned int i = 0; i < n_vectors - 1; ++i)
          //       {
          //         // do the work with orthonormal_vectors[i][j] and
          //         // current_vector[j] to fill the summation variables
          //       }
          //   }
          //
          // However, this access pattern leads to relatively poor memory
          // behavior with low performance because we would access only a few
          // entries of n_vectors 'orthonormal_vectors[i][j]' for an index j
          // at a time, before we pass on to the next vector i. More
          // precisely, we access inner_batch_size * n_lanes many entries
          // before moving on to the next vector. Modern CPUs derive a good
          // deal of performance from so-called hardware prefetching, which is
          // a mechanism that speculatively initiates loads from main memory
          // that the hardware guesses will be accessed soon, in order to
          // reduce the waiting time once the access actually happens. This
          // data is put into fast cache memory in the meantime. Prefetching
          // gets typically initiated when we access the entries of an array
          // consecutively, or also when looping over the data elements of a
          // few vectors at the same time. However, for more than around 10
          // vectors looped over simultaneously, the hardware gets to see too
          // many streams and the capacity of the loop stream detectors gets
          // exceeded. As a result, the hardware will first initiate a load
          // when the entry is actually requested by the respective code with
          // that loop index. As it takes many clock cycles for data to arrive
          // from memory (on the order of 200-500 clock cycles on 2024
          // hardware, whereas caches can deliver data in 4-20 cycles and
          // computations can be done in 3-6 cycles), and since even CPUs with
          // good out-of-order execution capabilities can issue only a limited
          // number of loads, the memory interface will become under-utilized,
          // which is counter-intuitive for a loop we expect to be very
          // memory-bandwidth heavy.
          //
          // The solution is to perform so-called loop blocking, see, e.g.,
          // https://en.wikipedia.org/wiki/Loop_nest_optimization - we split
          // the loop over the variable i into tiles (or blocks) of size 8 to
          // make sure the hardware prefetchers are able to follow all 8
          // streams, using a variable called 'i_block', and then run the loop
          // over 'c' inside. Since we do not want to re-load the entries of
          // 'current_vector' every time we work on blocks of size 8 for the
          // 'i' variable, but want to make sure to obtained it from faster
          // cache memory, we also tile the loop over 'c' into blocks. These
          // blocks have size 64, which is big enough to get most of the
          // effect of prefetchers (64 * inner_block_size * n_lanes is already
          // several kB of data) but small enough for 'current_vector' to
          // still be in fast cache memory for the next round of
          // 'i_block'. The end result are thus three nested loops visible
          // here, over 'c_block', 'i_block', and 'c', and an inner loop over
          // i and the inner batch size for the actual work. Not the prettiest
          // code, but giving adequate performance.
          for (unsigned int c_block = 0; c_block < (loop_length_c + 63) / 64;
               ++c_block)
            for (unsigned int i_block = 0; i_block < (n_vectors + 7) / 8;
                 ++i_block)
              for (c = c_block * 64, j = c * n_lanes * inner_batch_size;
                   c < std::min(loop_length_c, (c_block + 1) * 64);
                   ++c, j += n_lanes * inner_batch_size)
                {
                  VectorizedArray<double> vvec[inner_batch_size];
                  for (unsigned int k = 0; k < inner_batch_size; ++k)
                    vvec[k].load(current_vector + j + k * n_lanes);
                  VectorizedArray<double> prev_vector[inner_batch_size];
                  if (delayed_reorthogonalization || i_block == 0)
                    for (unsigned int k = 0; k < inner_batch_size; ++k)
                      prev_vector[k].load(orthogonal_vectors[n_vectors - 1] +
                                          j + k * n_lanes);

                  if (i_block == 0)
                    {
                      VectorizedArray<double> local_sum_0 =
                        prev_vector[0] * vvec[0];
                      VectorizedArray<double> local_sum_1 =
                        prev_vector[0] * prev_vector[0];
                      VectorizedArray<double> local_sum_2 = vvec[0] * vvec[0];
                      for (unsigned int k = 1; k < inner_batch_size; ++k)
                        {
                          local_sum_0 += prev_vector[k] * vvec[k];
                          if (delayed_reorthogonalization)
                            {
                              local_sum_1 += prev_vector[k] * prev_vector[k];
                              local_sum_2 += vvec[k] * vvec[k];
                            }
                        }
                      hs[n_vectors - 1] += local_sum_0;
                      if (delayed_reorthogonalization)
                        {
                          correct[n_vectors - 1] += local_sum_1;
                          correct[n_vectors] += local_sum_2;
                        }
                    }

                  for (unsigned int i = i_block * 8;
                       i < std::min(n_vectors - 1, (i_block + 1) * 8);
                       ++i)
                    {
                      // break the dependency chain into the field hs[i] for
                      // small sizes i by first accumulating 6 or 12 results
                      // into a local variable
                      VectorizedArray<double> temp;
                      temp.load(orthogonal_vectors[i] + j);
                      VectorizedArray<double> local_sum_0 = temp * vvec[0];
                      VectorizedArray<double> local_sum_1 =
                        delayed_reorthogonalization ? temp * prev_vector[0] :
                                                      0.;
                      for (unsigned int k = 1; k < inner_batch_size; ++k)
                        {
                          temp.load(orthogonal_vectors[i] + j + k * n_lanes);
                          local_sum_0 += temp * vvec[k];
                          if (delayed_reorthogonalization)
                            local_sum_1 += temp * prev_vector[k];
                        }
                      hs[i] += local_sum_0;
                      if (delayed_reorthogonalization)
                        correct[i] += local_sum_1;
                    }
                }

          c *= inner_batch_size;
          for (; c < locally_owned_size / n_lanes; ++c, j += n_lanes)
            {
              VectorizedArray<double> vvec, prev_vector;
              vvec.load(current_vector + j);
              prev_vector.load(orthogonal_vectors[n_vectors - 1] + j);
              hs[n_vectors - 1] += prev_vector * vvec;
              if (delayed_reorthogonalization)
                {
                  correct[n_vectors - 1] += prev_vector * prev_vector;
                  correct[n_vectors] += vvec * vvec;
                }

              for (unsigned int i = 0; i < n_vectors - 1; ++i)
                {
                  VectorizedArray<double> temp;
                  temp.load(orthogonal_vectors[i] + j);
                  hs[i] += temp * vvec;
                  if (delayed_reorthogonalization)
                    correct[i] += temp * prev_vector;
                }
            }

          for (unsigned int i = 0; i < n_vectors; ++i)
            {
              h(i) += hs[i].sum();
              if (delayed_reorthogonalization)
                h(i + n_vectors) += correct[i].sum();
            }
          if (delayed_reorthogonalization)
            h(n_vectors + n_vectors) += correct[n_vectors].sum();
        }

      // remainder loop of optimized path or non-optimized path (if
      // n>128)
      for (; j < locally_owned_size; ++j)
        {
          const double vvec        = current_vector[j];
          const double prev_vector = orthogonal_vectors[n_vectors - 1][j];
          h(n_vectors - 1) += prev_vector * vvec;
          if (delayed_reorthogonalization)
            {
              h(n_vectors + n_vectors - 1) += prev_vector * prev_vector;
              h(n_vectors + n_vectors) += vvec * vvec;
            }
          for (unsigned int i = 0; i < n_vectors - 1; ++i)
            {
              const double temp = orthogonal_vectors[i][j];
              h(i) += temp * vvec;
              if (delayed_reorthogonalization)
                h(n_vectors + i) += temp * prev_vector;
            }
        }
    }



    template <bool delayed_reorthogonalization, typename Number>
    double
    do_subtract_and_norm_host(const unsigned int             n_vectors,
                              const std::size_t              locally_owned_size,
                              OrthogonalizationData<Number> &temp_data,
                              const Vector<double>          &h,
                              Number                        *current_vector)
    {
      const auto &orthogonal_vectors = temp_data.vector_ptrs;
      double      norm_vv_temp       = 0;

      Number *previous_vector =
        const_cast<Number *>(orthogonal_vectors[n_vectors - 1]);
      const double inverse_norm_previous =
        delayed_reorthogonalization ? 1. / h(n_vectors + n_vectors - 1) : 0.;
      const double scaling_factor_vv =
        delayed_reorthogonalization ?
          (h(n_vectors + n_vectors) > 0.0 ?
             inverse_norm_previous / h(n_vectors + n_vectors) :
             inverse_norm_previous / h(n_vectors + n_vectors - 1)) :
          0.;
      const double last_factor = h(n_vectors - 1);

      VectorizedArray<double> norm_vv_temp_vectorized = 0.0;

      static constexpr unsigned int n_lanes = VectorizedArray<double>::size();
      constexpr unsigned int        inner_batch_size =
        delayed_reorthogonalization ? 6 : 12;

      // As for the do_Tvmult_add loop above, we perform loop blocking on both
      // the 'i' and 'c' variable to help hardware prefetchers to perform
      // adequately, and get three nested loops here plus the inner loops. See
      // the extensive comments above for the full rationale.
      unsigned int       j = 0;
      unsigned int       c = 0;
      const unsigned int loop_length_c =
        locally_owned_size / n_lanes / inner_batch_size;
      const unsigned int loop_length_i = (n_vectors + 7) / 8;
      for (unsigned int c_block = 0; c_block < (loop_length_c + 63) / 64;
           ++c_block)
        for (unsigned int i_block = 0; i_block < (n_vectors + 7) / 8; ++i_block)
          for (c = c_block * 64, j = c * n_lanes * inner_batch_size;
               c < std::min(loop_length_c, (c_block + 1) * 64);
               ++c, j += n_lanes * inner_batch_size)
            {
              VectorizedArray<double> temp[inner_batch_size];
              for (unsigned int k = 0; k < inner_batch_size; ++k)
                temp[k].load(current_vector + j + k * n_lanes);
              VectorizedArray<double> prev_vector[inner_batch_size];
              if (delayed_reorthogonalization)
                for (unsigned int k = 0; k < inner_batch_size; ++k)
                  prev_vector[k].load(previous_vector + j + k * n_lanes);

              for (unsigned int i = i_block * 8;
                   i < std::min(n_vectors - 1, (i_block + 1) * 8);
                   ++i)
                {
                  const double factor = h(i);
                  const double correction_factor =
                    (delayed_reorthogonalization ? h(n_vectors + i) : 0.0);
                  for (unsigned int k = 0; k < inner_batch_size; ++k)
                    {
                      VectorizedArray<double> vec;
                      vec.load(orthogonal_vectors[i] + j + k * n_lanes);
                      temp[k] -= factor * vec;
                      if (delayed_reorthogonalization)
                        prev_vector[k] -= correction_factor * vec;
                    }
                }

              if (delayed_reorthogonalization)
                {
                  if (i_block + 1 == loop_length_i)
                    for (unsigned int k = 0; k < inner_batch_size; ++k)
                      {
                        prev_vector[k] = prev_vector[k] * inverse_norm_previous;
                        prev_vector[k].store(previous_vector + j + k * n_lanes);
                        temp[k] -= last_factor * prev_vector[k];
                        temp[k] = temp[k] * scaling_factor_vv;
                        temp[k].store(current_vector + j + k * n_lanes);
                      }
                  else
                    for (unsigned int k = 0; k < inner_batch_size; ++k)
                      {
                        prev_vector[k].store(previous_vector + j + k * n_lanes);
                        temp[k].store(current_vector + j + k * n_lanes);
                      }
                }
              else
                {
                  if (i_block + 1 == loop_length_i)
                    for (unsigned int k = 0; k < inner_batch_size; ++k)
                      {
                        prev_vector[k].load(previous_vector + j + k * n_lanes);
                        temp[k] -= last_factor * prev_vector[k];
                        temp[k].store(current_vector + j + k * n_lanes);
                        norm_vv_temp_vectorized += temp[k] * temp[k];
                      }
                  else
                    for (unsigned int k = 0; k < inner_batch_size; ++k)
                      temp[k].store(current_vector + j + k * n_lanes);
                }
            }

      c *= inner_batch_size;
      for (; c < locally_owned_size / n_lanes; ++c, j += n_lanes)
        {
          VectorizedArray<double> temp, prev_vector;
          temp.load(current_vector + j);
          prev_vector.load(previous_vector + j);
          if (!delayed_reorthogonalization)
            temp -= h(n_vectors - 1) * prev_vector;

          for (unsigned int i = 0; i < n_vectors - 1; ++i)
            {
              VectorizedArray<double> vec;
              vec.load(orthogonal_vectors[i] + j);
              temp -= h(i) * vec;
              if (delayed_reorthogonalization)
                prev_vector -= h(n_vectors + i) * vec;
            }

          if (delayed_reorthogonalization)
            {
              prev_vector = prev_vector * inverse_norm_previous;
              prev_vector.store(previous_vector + j);
              temp -= h(n_vectors - 1) * prev_vector;
              temp = temp * scaling_factor_vv;
              temp.store(current_vector + j);
            }
          else
            {
              temp.store(current_vector + j);
              norm_vv_temp_vectorized += temp * temp;
            }
        }

      if (!delayed_reorthogonalization)
        norm_vv_temp += norm_vv_temp_vectorized.sum();

      for (; j < locally_owned_size; ++j)
        {
          double temp        = current_vector[j];
          double prev_vector = previous_vector[j];
          if (delayed_reorthogonalization)
            {
              for (unsigned int i = 0; i < n_vectors - 1; ++i)
                {
                  const double vec = orthogonal_vectors[i][j];
                  temp -= h(i) * vec;
                  prev_vector -= h(n_vectors + i) * vec;
                }
              prev_vector *= inverse_norm_previous;
              previous_vector[j] = prev_vector;
              temp -= h(n_vectors - 1) * prev_vector;
              temp *= scaling_factor_vv;
            }
          else
            {
              temp -= h(n_vectors - 1) * prev_vector;
              for (unsigned int i = 0; i < n_vectors - 1; ++i)
                temp -= h(i) * orthogonal_vectors[i][j];
              norm_vv_temp += temp * temp;
            }
          current_vector[j] = temp;
        }

      return norm_vv_temp;
    }



    template <typename Number>
    void
    do_add_host(const unsigned int             n_vectors,
                const std::size_t              locally_owned_size,
                OrthogonalizationData<Number> &temp_data,
                const Vector<double>          &h,
                const bool                     zero_out,
                Number                        *output)
    {
      const auto &orthogonal_vectors = temp_data.vector_ptrs;

      static constexpr unsigned int n_lanes = VectorizedArray<double>::size();
      constexpr unsigned int        inner_batch_size = 12;

      unsigned int j = 0;
      unsigned int c = 0;
      for (; c < locally_owned_size / n_lanes / inner_batch_size;
           ++c, j += n_lanes * inner_batch_size)
        {
          VectorizedArray<double> temp[inner_batch_size];
          if (zero_out)
            for (VectorizedArray<double> &a : temp)
              a = {};
          else
            for (unsigned int k = 0; k < inner_batch_size; ++k)
              temp[k].load(output + j + k * n_lanes);

          for (unsigned int i = 0; i < n_vectors; ++i)
            {
              const double h_i = h(i);
              for (unsigned int k = 0; k < inner_batch_size; ++k)
                {
                  VectorizedArray<double> v_ij;
                  v_ij.load(orthogonal_vectors[i] + j + k * n_lanes);
                  temp[k] += v_ij * h_i;
                }
            }

          for (unsigned int k = 0; k < inner_batch_size; ++k)
            temp[k].store(output + j + k * n_lanes);
        }

      c *= inner_batch_size;
      for (; c < locally_owned_size / n_lanes; ++c, j += n_lanes)
        {
          VectorizedArray<double> temp = {};
          if (!zero_out)
            temp.load(output + j);

          for (unsigned int i = 0; i < n_vectors; ++i)
            {
              VectorizedArray<double> v_ij;
              v_ij.load(orthogonal_vectors[i] + j);
              temp += v_ij * h(i);
            }

          temp.store(output + j);
        }

      for (; j < locally_owned_size; ++j)
        {
          double temp = zero_out ? 0.0 : output[j];
          for (unsigned int i = 0; i < n_vectors; ++i)
            temp += orthogonal_vectors[i][j] * h(i);
          output[j] = temp;
        }
    }



    /**
     * Wrapper around a pointer to the data of a vector on the device. This
     * is necessary because a Kokkos::View<Number **> would be interpreted as
     * a two-dimensional array rather than an array of pointers.
     */
    template <typename Number>
    struct DevicePointer
    {
      Number *ptr;
    };



    template <typename Number>
    struct DeviceData
    {
      using ExecutionSpace =
        MemorySpace::Default::kokkos_space::execution_space;
      using KokkosSpace = MemorySpace::Default::kokkos_space;

      /**
       * Number of vector entries processed by one team in the kernel of
       * do_Tvmult_add_device(). The entries of the current vector (and the
       * previous vector in the delayed variant) of a tile are re-read from
       * the L1/L2 caches for every block of vectors, whereas the
       * orthogonal vectors are read from main memory only once.
       */
      static constexpr std::size_t rows_per_team = 2048;

      /**
       * Number of vectors handled at once in do_Tvmult_add_device().
       */
      static constexpr int block_size = 8;

      void
      set_pointers(const ExecutionSpace              &exec,
                   const std::vector<const Number *> &vector_ptrs)
      {
        const std::size_t n_vectors = vector_ptrs.size();
        if (ptrs.extent(0) < n_vectors)
          {
            // grow geometrically to avoid re-allocations as the basis grows
            ptrs = Kokkos::View<DevicePointer<Number> *, KokkosSpace>(
              Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                 "dealii::gmres_vector_pointers"),
              std::max(n_vectors, 2 * ptrs.extent(0)));
            ptrs_host = Kokkos::create_mirror_view(ptrs);
          }
        // the pointers are only used for read access, except for the
        // in-place update of the last vector in the delayed variant
        for (std::size_t i = 0; i < n_vectors; ++i)
          ptrs_host(i).ptr = const_cast<Number *>(vector_ptrs[i]);
        const auto range = std::make_pair(std::size_t(0), n_vectors);
        Kokkos::deep_copy(exec,
                          Kokkos::subview(ptrs, range),
                          Kokkos::subview(ptrs_host, range));
      }

      void
      set_coefficients(const ExecutionSpace &exec,
                       const std::size_t     n_coefficients,
                       const double         *h)
      {
        if (coefficients.extent(0) < n_coefficients)
          {
            coefficients = Kokkos::View<double *, KokkosSpace>(
              Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                 "dealii::gmres_coefficients"),
              std::max(n_coefficients, 2 * coefficients.extent(0)));
            coefficients_host = Kokkos::create_mirror_view(coefficients);
          }
        for (std::size_t i = 0; i < n_coefficients; ++i)
          coefficients_host(i) = h[i];
        const auto range = std::make_pair(std::size_t(0), n_coefficients);
        Kokkos::deep_copy(exec,
                          Kokkos::subview(coefficients, range),
                          Kokkos::subview(coefficients_host, range));
      }

      void
      reserve_results(const std::size_t n_teams, const std::size_t n_results)
      {
        if (partial_sums.extent(0) < n_teams ||
            partial_sums.extent(1) < n_results)
          partial_sums =
            Kokkos::View<double **, Kokkos::LayoutLeft, KokkosSpace>(
              Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                 "dealii::gmres_partial_sums"),
              std::max(n_teams, partial_sums.extent(0)),
              std::max(2 * n_results, partial_sums.extent(1)));
        if (results.extent(0) < n_results)
          {
            results = Kokkos::View<double *, KokkosSpace>(
              Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                 "dealii::gmres_results"),
              2 * n_results);
            results_host = Kokkos::create_mirror_view(results);
          }
      }

      Kokkos::View<DevicePointer<Number> *, KokkosSpace>   ptrs;
      typename Kokkos::View<DevicePointer<Number> *,
                            KokkosSpace>::host_mirror_type ptrs_host;

      Kokkos::View<double *, KokkosSpace> coefficients;
      typename Kokkos::View<double *, KokkosSpace>::host_mirror_type
        coefficients_host;

      // layout (team, result) with LayoutLeft for coalesced access in the
      // final reduction over the teams
      Kokkos::View<double **, Kokkos::LayoutLeft, KokkosSpace> partial_sums;

      Kokkos::View<double *, KokkosSpace> results;
      typename Kokkos::View<double *, KokkosSpace>::host_mirror_type
        results_host;
    };



    template <typename Number>
    DeviceData<Number> &
    get_device_data(OrthogonalizationData<Number> &temp_data)
    {
      if (temp_data.device_data == nullptr)
        temp_data.device_data = std::make_shared<DeviceData<Number>>();
      return *temp_data.device_data;
    }



    template <bool delayed_reorthogonalization, typename Number>
    void
    do_Tvmult_add_device(const unsigned int             n_vectors,
                         const std::size_t              locally_owned_size,
                         const Number                  *current_vector,
                         OrthogonalizationData<Number> &temp_data,
                         Vector<double>                &h)
    {
      using Data                          = DeviceData<Number>;
      constexpr int         block_size    = Data::block_size;
      constexpr std::size_t rows_per_team = Data::rows_per_team;

      Assert(n_vectors > 0, ExcInternalError());
      AssertDimension(temp_data.vector_ptrs.size(), n_vectors);

      // Compute the following quantities, with v_i the orthogonal vectors
      // and w the current vector:
      //   h(i)     += v_i * w        for i < n
      // and for the delayed variant additionally
      //   h(n + i) += v_i * v_{n-1}  for i < n
      //   h(2n)    += w * w
      const unsigned int n_results =
        delayed_reorthogonalization ? 2 * n_vectors + 1 : n_vectors;
      AssertDimension(h.size(), n_results);
      const std::size_t n_teams =
        (locally_owned_size + rows_per_team - 1) / rows_per_team;
      if (n_teams == 0)
        return;

      Data                         &data = get_device_data(temp_data);
      typename Data::ExecutionSpace exec;
      data.set_pointers(exec, temp_data.vector_ptrs);
      data.reserve_results(n_teams, n_results);

      const auto         ptrs         = data.ptrs;
      const auto         partial_sums = data.partial_sums;
      const auto         results      = data.results;
      const unsigned int n_blocks = (n_vectors + block_size - 1) / block_size;
      const std::size_t  size     = locally_owned_size;

      using TeamPolicy = Kokkos::TeamPolicy<typename Data::ExecutionSpace>;
      using Member     = typename TeamPolicy::member_type;

      // Each team works on a tile of the vectors and loops over the
      // orthogonal vectors in blocks of size block_size, which are kept in
      // registers. The current vector is re-read from caches for each block.
      // Partial sums of all teams are written to an intermediate array and
      // summed in a second kernel, which makes the result deterministic in
      // contrast to an approach with atomics.
      Kokkos::parallel_for(
        "dealii::gmres_multi_dot",
        TeamPolicy(exec, n_teams, Kokkos::AUTO),
        KOKKOS_LAMBDA(const Member &team) {
          const std::size_t team_id = team.league_rank();
          const std::size_t begin   = team_id * rows_per_team;
          const std::size_t end     = Kokkos::min(begin + rows_per_team, size);
          [[maybe_unused]] const Number *prev = ptrs(n_vectors - 1).ptr;

          for (unsigned int block = 0; block < n_blocks; ++block)
            {
              const unsigned int i0 = block * block_size;
              const unsigned int nb =
                Kokkos::min<unsigned int>(block_size, n_vectors - i0);

              // pad the last block with a valid pointer (data is in caches,
              // results are discarded) to keep the inner loop unrolled
              const Number *v[block_size];
              for (int b = 0; b < block_size; ++b)
                v[b] = ptrs(i0 + (b < static_cast<int>(nb) ? b : 0)).ptr;

              BlockSums<block_size, delayed_reorthogonalization> sums;
              Kokkos::parallel_reduce(
                Kokkos::TeamThreadRange(team, begin, end),
                [&](const std::size_t                                   j,
                    BlockSums<block_size, delayed_reorthogonalization> &s) {
                  const double w_j    = current_vector[j];
                  double       prev_j = 0.;
                  if constexpr (delayed_reorthogonalization)
                    prev_j = prev[j];
                  for (int b = 0; b < block_size; ++b)
                    {
                      const double v_ij = v[b][j];
                      s.dot[b] += v_ij * w_j;
                      if constexpr (delayed_reorthogonalization)
                        s.correct[b] += v_ij * prev_j;
                    }
                  if (delayed_reorthogonalization && block == 0)
                    s.norm_sqr += w_j * w_j;
                },
                sums);

              Kokkos::single(Kokkos::PerTeam(team), [&]() {
                for (unsigned int b = 0; b < nb; ++b)
                  {
                    partial_sums(team_id, i0 + b) = sums.dot[b];
                    if constexpr (delayed_reorthogonalization)
                      partial_sums(team_id, n_vectors + i0 + b) =
                        sums.correct[b];
                  }
                if (delayed_reorthogonalization && block == 0)
                  partial_sums(team_id, 2 * n_vectors) = sums.norm_sqr;
              });
            }
        });

      Kokkos::parallel_for(
        "dealii::gmres_multi_dot_finalize",
        TeamPolicy(exec, n_results, Kokkos::AUTO),
        KOKKOS_LAMBDA(const Member &team) {
          const unsigned int i   = team.league_rank();
          double             sum = 0.;
          Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team, n_teams),
            [&](const std::size_t t, double &s) { s += partial_sums(t, i); },
            sum);
          Kokkos::single(Kokkos::PerTeam(team), [&]() { results(i) = sum; });
        });

      const auto range = std::make_pair(0u, n_results);
      Kokkos::deep_copy(exec,
                        Kokkos::subview(data.results_host, range),
                        Kokkos::subview(results, range));
      exec.fence();
      for (unsigned int i = 0; i < n_results; ++i)
        h(i) += data.results_host(i);
    }



    template <bool delayed_reorthogonalization, typename Number>
    double
    do_subtract_and_norm_device(const unsigned int n_vectors,
                                const std::size_t  locally_owned_size,
                                OrthogonalizationData<Number> &temp_data,
                                const Vector<double>          &h,
                                Number                        *current_vector)
    {
      using Data   = DeviceData<Number>;
      using Policy = Kokkos::RangePolicy<typename Data::ExecutionSpace,
                                         Kokkos::IndexType<std::size_t>>;

      Assert(n_vectors > 0, ExcInternalError());
      AssertDimension(temp_data.vector_ptrs.size(), n_vectors);
      AssertDimension(h.size(),
                      delayed_reorthogonalization ? 2 * n_vectors + 1 :
                                                    n_vectors);
      if (locally_owned_size == 0)
        return 0.;

      Data                         &data = get_device_data(temp_data);
      typename Data::ExecutionSpace exec;
      data.set_pointers(exec, temp_data.vector_ptrs);
      data.set_coefficients(exec, h.size(), h.data());

      const auto ptrs   = data.ptrs;
      const auto coeffs = data.coefficients;

      if constexpr (delayed_reorthogonalization)
        {
          // Same operations as in do_subtract_and_norm: subtract the
          // projections from the current vector and the correction terms
          // from the previous vector, which gets normalized, all in one sweep
          const double inverse_norm_previous =
            1. / h(n_vectors + n_vectors - 1);
          const double scaling_factor_vv =
            h(n_vectors + n_vectors) > 0.0 ?
              inverse_norm_previous / h(n_vectors + n_vectors) :
              inverse_norm_previous / h(n_vectors + n_vectors - 1);
          const double last_factor = h(n_vectors - 1);

          Kokkos::parallel_for(
            "dealii::gmres_subtract_delayed",
            Policy(exec, 0, locally_owned_size),
            KOKKOS_LAMBDA(const std::size_t j) {
              Number *previous_vector = ptrs(n_vectors - 1).ptr;
              double  temp            = current_vector[j];
              double  prev            = previous_vector[j];
              for (unsigned int i = 0; i + 1 < n_vectors; ++i)
                {
                  const double vec = ptrs(i).ptr[j];
                  temp -= coeffs(i) * vec;
                  prev -= coeffs(n_vectors + i) * vec;
                }
              prev *= inverse_norm_previous;
              previous_vector[j] = prev;
              temp -= last_factor * prev;
              current_vector[j] = temp * scaling_factor_vv;
            });

          // the norm is computed from other quantities in the delayed
          // variant, and no synchronization with the device is needed here
          return 0.;
        }
      else
        {
          double norm_sqr = 0.;
          Kokkos::parallel_reduce(
            "dealii::gmres_subtract_and_norm",
            Policy(exec, 0, locally_owned_size),
            KOKKOS_LAMBDA(const std::size_t j, double &sum) {
              double temp = current_vector[j];
              for (unsigned int i = 0; i < n_vectors; ++i)
                temp -= coeffs(i) * ptrs(i).ptr[j];
              current_vector[j] = temp;
              sum += temp * temp;
            },
            norm_sqr);
          return norm_sqr;
        }
    }



    template <typename Number>
    void
    do_add_device(const unsigned int             n_vectors,
                  const std::size_t              locally_owned_size,
                  OrthogonalizationData<Number> &temp_data,
                  const Vector<double>          &h,
                  const bool                     zero_out,
                  Number                        *output)
    {
      using Data = DeviceData<Number>;

      AssertDimension(temp_data.vector_ptrs.size(), n_vectors);
      AssertIndexRange(n_vectors, h.size() + 1);
      if (locally_owned_size == 0 || n_vectors == 0)
        return;

      Data                         &data = get_device_data(temp_data);
      typename Data::ExecutionSpace exec;
      data.set_pointers(exec, temp_data.vector_ptrs);
      data.set_coefficients(exec, n_vectors, h.data());

      const auto ptrs   = data.ptrs;
      const auto coeffs = data.coefficients;

      Kokkos::parallel_for(
        "dealii::gmres_add",
        Kokkos::RangePolicy<typename Data::ExecutionSpace,
                            Kokkos::IndexType<std::size_t>>(exec,
                                                            0,
                                                            locally_owned_size),
        KOKKOS_LAMBDA(const std::size_t j) {
          double temp = zero_out ? 0. : static_cast<double>(output[j]);
          for (unsigned int i = 0; i < n_vectors; ++i)
            temp += coeffs(i) * ptrs(i).ptr[j];
          output[j] = temp;
        });
    }



    template <bool delayed_reorthogonalization,
              typename MemorySpace,
              typename Number>
    void
    do_Tvmult_add(const unsigned int             n_vectors,
                  const std::size_t              locally_owned_size,
                  const Number                  *current_vector,
                  OrthogonalizationData<Number> &temp_data,
                  Vector<double>                &h)
    {
      if constexpr (std::is_same_v<MemorySpace, dealii::MemorySpace::Default>)
        do_Tvmult_add_device<delayed_reorthogonalization>(
          n_vectors, locally_owned_size, current_vector, temp_data, h);
      else
        do_Tvmult_add_host<delayed_reorthogonalization>(
          n_vectors, locally_owned_size, current_vector, temp_data, h);
    }



    template <bool delayed_reorthogonalization,
              typename MemorySpace,
              typename Number>
    double
    do_subtract_and_norm(const unsigned int             n_vectors,
                         const std::size_t              locally_owned_size,
                         OrthogonalizationData<Number> &temp_data,
                         const Vector<double>          &h,
                         Number                        *current_vector)
    {
      if constexpr (std::is_same_v<MemorySpace, dealii::MemorySpace::Default>)
        return do_subtract_and_norm_device<delayed_reorthogonalization>(
          n_vectors, locally_owned_size, temp_data, h, current_vector);
      else
        return do_subtract_and_norm_host<delayed_reorthogonalization>(
          n_vectors, locally_owned_size, temp_data, h, current_vector);
    }



    template <typename MemorySpace, typename Number>
    void
    do_add(const unsigned int             n_vectors,
           const std::size_t              locally_owned_size,
           OrthogonalizationData<Number> &temp_data,
           const Vector<double>          &h,
           const bool                     zero_out,
           Number                        *output)
    {
      if constexpr (std::is_same_v<MemorySpace, dealii::MemorySpace::Default>)
        do_add_device(
          n_vectors, locally_owned_size, temp_data, h, zero_out, output);
      else
        do_add_host(
          n_vectors, locally_owned_size, temp_data, h, zero_out, output);
    }
  } // namespace SolverGMRESImplementation
} // namespace internal

#include "lac/solver_gmres.inst"

DEAL_II_NAMESPACE_CLOSE
