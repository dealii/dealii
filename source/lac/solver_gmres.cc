// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
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

#include "solver_gmres.inst"

DEAL_II_NAMESPACE_CLOSE
