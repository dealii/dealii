// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

#include <deal.II/lac/cuda_kernels.templates.h>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
    namespace kernel
    {
      /////////////////////////////////////////////////////////////////////////
      // Explicit instantiation                                              //
      /////////////////////////////////////////////////////////////////////////

      template __global__ void
      vec_scale<float>(float *, const float a, const size_type);
      template __global__ void
      vector_bin_op<float, Binop_Addition>(float *         v1,
                                           const float *   v2,
                                           const size_type N);
      template __global__ void
      vector_bin_op<float, Binop_Subtraction>(float *         v1,
                                              const float *   v2,
                                              const size_type N);
      template __global__ void
      masked_vector_bin_op<float, Binop_Addition>(const unsigned int *mask,
                                                  float *             v1,
                                                  const float *       v2,
                                                  const size_type     N);
      template __global__ void
      masked_vector_bin_op<float, Binop_Subtraction>(const unsigned int *mask,
                                                     float *             v1,
                                                     const float *       v2,
                                                     const size_type     N);
      template struct ElemSum<float>;
      template struct L1Norm<float>;
      template struct LInfty<float>;
      template __global__ void
      reduction<float, ElemSum<float>>(float *         result,
                                       const float *   v,
                                       const size_type N);
      template __global__ void
      reduction<float, L1Norm<float>>(float *         result,
                                      const float *   v,
                                      const size_type N);
      template __global__ void
      reduction<float, LInfty<float>>(float *         result,
                                      const float *   v,
                                      const size_type N);
      template struct DotProduct<float>;
      template __global__ void
      double_vector_reduction<float, DotProduct<float>>(float *         result,
                                                        const float *   v1,
                                                        const float *   v2,
                                                        const size_type N);
      template __global__ void
      vec_add<float>(float *val, const float, const size_type N);
      template __global__ void
      add_aV<float>(float *         val,
                    const float     a,
                    const float *   V_val,
                    const size_type N);
      template __global__ void
      add_aVbW<float>(float *         val,
                      const float     a,
                      const float *   V_val,
                      const float     b,
                      const float *   W_val,
                      const size_type N);
      template __global__ void
      sadd<float>(const float     s,
                  float *         val,
                  const float     a,
                  const float *   V_val,
                  const size_type N);
      template __global__ void
      sadd<float>(const float     s,
                  float *         val,
                  const float     a,
                  const float *   V_val,
                  const float     b,
                  const float *   W_val,
                  const size_type N);
      template __global__ void
      scale<float>(float *val, const float *V_val, const size_type N);
      template __global__ void
      equ<float>(float *         val,
                 const float     a,
                 const float *   V_val,
                 const size_type N);
      template __global__ void
      equ<float>(float *         val,
                 const float     a,
                 const float *   V_val,
                 const float     b,
                 const float *   W_val,
                 const size_type N);
      template __global__ void
      add_and_dot<float>(float *         res,
                         float *         v1,
                         const float *   v2,
                         const float *   v3,
                         const float     a,
                         const size_type N);
      template __global__ void
      set<float>(float *val, const float s, const size_type N);
      template __global__ void
      set_permutated<float, size_type>(const size_type *indices,
                                       float *          val,
                                       const float *    v,
                                       const size_type  N);
      template __global__ void
      gather<float, size_type>(float *          val,
                               const size_type *indices,
                               const float *    v,
                               const size_type  N);
      template __global__ void
      add_permutated<float>(const size_type *indices,
                            float *          val,
                            const float *    v,
                            const size_type  N);



      template __global__ void
      vec_scale<double>(double *, const double a, const size_type);
      template __global__ void
      vector_bin_op<double, Binop_Addition>(double *        v1,
                                            const double *  v2,
                                            const size_type N);
      template __global__ void
      vector_bin_op<double, Binop_Subtraction>(double *        v1,
                                               const double *  v2,
                                               const size_type N);
      template __global__ void
      masked_vector_bin_op<double, Binop_Addition>(const unsigned int *mask,
                                                   double *            v1,
                                                   const double *      v2,
                                                   const size_type     N);
      template __global__ void
      masked_vector_bin_op<double, Binop_Subtraction>(const unsigned int *mask,
                                                      double *            v1,
                                                      const double *      v2,
                                                      const size_type     N);
      template struct ElemSum<double>;
      template struct L1Norm<double>;
      template struct LInfty<double>;
      template __global__ void
      reduction<double, ElemSum<double>>(double *        result,
                                         const double *  v,
                                         const size_type N);
      template __global__ void
      reduction<double, L1Norm<double>>(double *        result,
                                        const double *  v,
                                        const size_type N);
      template __global__ void
      reduction<double, LInfty<double>>(double *        result,
                                        const double *  v,
                                        const size_type N);
      template struct DotProduct<double>;
      template __global__ void
      double_vector_reduction<double, DotProduct<double>>(double *      result,
                                                          const double *v1,
                                                          const double *v2,
                                                          const size_type N);
      template __global__ void
      vec_add<double>(double *val, const double, const size_type N);
      template __global__ void
      add_aV<double>(double *        val,
                     const double    a,
                     const double *  V_val,
                     const size_type N);
      template __global__ void
      add_aVbW<double>(double *        val,
                       const double    a,
                       const double *  V_val,
                       const double    b,
                       const double *  W_val,
                       const size_type N);
      template __global__ void
      sadd<double>(const double    s,
                   double *        val,
                   const double    a,
                   const double *  V_val,
                   const size_type N);
      template __global__ void
      sadd<double>(const double    s,
                   double *        val,
                   const double    a,
                   const double *  V_val,
                   const double    b,
                   const double *  W_val,
                   const size_type N);
      template __global__ void
      scale<double>(double *val, const double *V_val, const size_type N);
      template __global__ void
      equ<double>(double *        val,
                  const double    a,
                  const double *  V_val,
                  const size_type N);
      template __global__ void
      equ<double>(double *        val,
                  const double    a,
                  const double *  V_val,
                  const double    b,
                  const double *  W_val,
                  const size_type N);
      template __global__ void
      add_and_dot<double>(double *        res,
                          double *        v1,
                          const double *  v2,
                          const double *  v3,
                          const double    a,
                          const size_type N);
      template __global__ void
      set<double>(double *val, const double s, const size_type N);
      template __global__ void
      set_permutated<double, size_type>(const size_type *indices,
                                        double *         val,
                                        const double *   v,
                                        const size_type  N);
      template __global__ void
      gather<double, size_type>(double *         val,
                                const size_type *indices,
                                const double *   v,
                                const size_type  N);
      template __global__ void
      add_permutated<double>(const size_type *indices,
                             double *         val,
                             const double *   v,
                             const size_type  N);
    } // namespace kernel
  }   // namespace CUDAWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE
