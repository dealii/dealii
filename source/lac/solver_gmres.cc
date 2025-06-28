// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/vectorization.h>

#include <deal.II/lac/solver_gmres.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace SolverGMRESImplementation
  {
    template <bool delayed_reorthogonalization, typename Number>
    void
    do_Tvmult_add(const unsigned int                 n_vectors,
                  const std::size_t                  locally_owned_size,
                  const Number                      *current_vector,
                  const std::vector<const Number *> &orthogonal_vectors,
                  Vector<double>                    &h)
    {
      unsigned int j = 0;

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
    do_subtract_and_norm(const unsigned int                 n_vectors,
                         const std::size_t                  locally_owned_size,
                         const std::vector<const Number *> &orthogonal_vectors,
                         const Vector<double>              &h,
                         Number                            *current_vector)
    {
      double norm_vv_temp = 0;

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
    do_add(const unsigned int                 n_vectors,
           const std::size_t                  locally_owned_size,
           const std::vector<const Number *> &tmp_vectors,
           const Vector<double>              &h,
           const bool                         zero_out,
           Number                            *output)
    {
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
                  v_ij.load(tmp_vectors[i] + j + k * n_lanes);
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
              v_ij.load(tmp_vectors[i] + j);
              temp += v_ij * h(i);
            }

          temp.store(output + j);
        }

      for (; j < locally_owned_size; ++j)
        {
          double temp = zero_out ? 0.0 : output[j];
          for (unsigned int i = 0; i < n_vectors; ++i)
            temp += tmp_vectors[i][j] * h(i);
          output[j] = temp;
        }
    }
  } // namespace SolverGMRESImplementation
} // namespace internal

#include "lac/solver_gmres.inst"

DEAL_II_NAMESPACE_CLOSE
