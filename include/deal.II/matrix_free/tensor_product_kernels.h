// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_tensor_product_kernels_h
#define dealii_matrix_free_tensor_product_kernels_h

#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/utilities.h>

#include <deal.II/matrix_free/shape_info.h>


DEAL_II_NAMESPACE_OPEN



namespace internal
{
  /**
   * In this namespace, the evaluator routines that evaluate the tensor
   * products are implemented.
   */
  enum EvaluatorVariant
  {
    /**
     * Do not use anything more than the tensor product structure of the
     * finite element.
     */
    evaluate_general,
    /**
     * Perform evaluation by exploiting symmetry in the finite element: i.e.,
     * skip some computations by utilizing the symmetry in the shape functions
     * and quadrature points.
     */
    evaluate_symmetric,
    /**
     * Use symmetry to apply the operator to even and odd parts of the input
     * vector separately: see the documentation of the EvaluatorTensorProduct
     * specialization for more information.
     */
    evaluate_evenodd,
    /**
     * Use symmetry in Legendre and similar polynomial spaces where the shape
     * functions with even number are symmetric about the center of the
     * quadrature points (think about even polynomial degrees) and the shape
     * functions with odd number are anti-symmetric about the center of the
     * quadrature points (think about odd polynomial degrees). This allows to
     * use a strategy similar to the even-odd technique but without separate
     * coefficient arrays. See the documentation of the EvaluatorTensorProduct
     * specialization for more information.
     */
    evaluate_symmetric_hierarchical
  };



  /**
   * Determine which quantity should be computed via the tensor product kernels.
   */
  enum class EvaluatorQuantity
  {
    /**
     * Evaluate/integrate by shape functions.
     */
    value,
    /**
     * Evaluate/integrate by gradients of the shape functions.
     */
    gradient,
    /**
     * Evaluate/integrate by hessians of the shape functions.
     */
    hessian
  };



  /**
   * One-dimensional kernel for use by the generic tensor product
   * interpolation as provided by the class EvaluatorTensorProduct,
   * implementing a matrix-vector product along this dimension, controlled by
   * the number of rows and columns and the stride in the input and output
   * arrays, which are embedded into some lexicographic ordering of unknowns
   * in a tensor-product arrangement.
   *
   * Besides this generic function for templated loop lengths, there are
   * several specializations of this class to account for run-time matrix
   * sizes as well as some symmetries that reduce the data access or
   * arithmetic operations. The specializations are technically realized by
   * conditional function overloading with std::enable_if_t based on the first
   * template parameter.
   */
  template <EvaluatorVariant  variant,
            EvaluatorQuantity quantity,
            int               n_rows,
            int               n_columns,
            int               stride_in,
            int               stride_out,
            bool              transpose_matrix,
            bool              add,
            typename Number,
            typename Number2>
  std::enable_if_t<(variant == evaluate_general), void>
  apply_matrix_vector_product(const Number2 *matrix,
                              const Number  *in,
                              Number        *out)
  {
    // We can only statically assert that one argument is non-zero because
    // face evaluation might instantiate some functions, so we need to use the
    // run-time assert to verify that we do not end up involuntarily.
    static_assert(n_rows > 0 || n_columns > 0,
                  "Specialization only for n_rows, n_columns > 0");
    Assert(n_rows > 0 && n_columns > 0,
           ExcInternalError("The evaluation needs n_rows, n_columns > 0, but " +
                            std::to_string(n_rows) + ", " +
                            std::to_string(n_columns) + " was passed!"));
    static_assert(quantity == EvaluatorQuantity::value,
                  "This function should only use EvaluatorQuantity::value");

    constexpr int mm = transpose_matrix ? n_rows : n_columns,
                  nn = transpose_matrix ? n_columns : n_rows;

    std::array<Number, mm> x;
    for (int i = 0; i < mm; ++i)
      x[i] = in[stride_in * i];
    for (int col = 0; col < nn; ++col)
      {
        Number res0;
        if (transpose_matrix == true)
          {
            res0 = matrix[col] * x[0];
            for (int i = 1; i < mm; ++i)
              {
                const Number2 mji = matrix[i * n_columns + col];
                if constexpr (numbers::NumberTraits<Number>::is_complex &&
                              numbers::NumberTraits<Number2>::is_complex)
                  {
                    res0.real(res0.real() + mji.real() * x[i].real() -
                              mji.imag() * x[i].imag());
                    res0.imag(res0.imag() + mji.imag() * x[i].real() +
                              mji.real() * x[i].imag());
                  }
                else
                  res0 += mji * x[i];
              }
          }
        else
          {
            res0 = matrix[col * n_columns] * x[0];
            for (int i = 1; i < mm; ++i)
              {
                const Number2 mij = matrix[col * n_columns + i];
                if constexpr (numbers::NumberTraits<Number>::is_complex &&
                              numbers::NumberTraits<Number2>::is_complex)
                  {
                    res0.real(res0.real() + mij.real() * x[i].real() -
                              mij.imag() * x[i].imag());
                    res0.imag(res0.imag() + mij.imag() * x[i].real() +
                              mij.real() * x[i].imag());
                  }
                else
                  res0 += mij * x[i];
              }
          }
        if (add)
          out[stride_out * col] += res0;
        else
          out[stride_out * col] = res0;
      }
  }



  /**
   * Specialization of the matrix-vector kernel for run-time loop bounds in
   * the generic evaluator.
   */
  template <EvaluatorVariant  variant,
            EvaluatorQuantity quantity,
            bool              transpose_matrix,
            bool              add,
            bool              consider_strides,
            typename Number,
            typename Number2,
            int n_components = 1>
  std::enable_if_t<(variant == evaluate_general), void>
  apply_matrix_vector_product(const Number2 *matrix,
                              const Number  *in,
                              Number        *out,
                              const int      n_rows,
                              const int      n_columns,
                              const int      stride_in_given,
                              const int      stride_out_given)
  {
    const int mm = transpose_matrix ? n_rows : n_columns,
              nn = transpose_matrix ? n_columns : n_rows;
    Assert(n_rows > 0 && n_columns > 0,
           ExcInternalError("Empty evaluation task!"));
    Assert(n_rows > 0 && n_columns > 0,
           ExcInternalError("The evaluation needs n_rows, n_columns > 0, but " +
                            std::to_string(n_rows) + ", " +
                            std::to_string(n_columns) + " was passed!"));

    static_assert(quantity == EvaluatorQuantity::value,
                  "This function should only use EvaluatorQuantity::value");

    Assert(consider_strides || (stride_in_given == 1 && stride_out_given == 1),
           ExcInternalError());
    const int stride_in  = consider_strides ? stride_in_given : 1;
    const int stride_out = consider_strides ? stride_out_given : 1;

    static_assert(n_components > 0 && n_components < 4,
                  "Invalid number of components");

    // specialization for n_rows = 2 that manually unrolls the innermost loop
    // to make the operation perform better (not completely as good as the
    // templated one, but much better than the generic version down below,
    // because the loop over col can be more effectively unrolled by the
    // compiler)
    if (transpose_matrix && n_rows == 2 && n_components == 1)
      {
        const Number2 *matrix_1 = matrix + n_columns;
        const Number   x0 = in[0], x1 = in[stride_in];
        for (int col = 0; col < nn; ++col)
          {
            const Number result = matrix[col] * x0 + matrix_1[col] * x1;
            if (add)
              out[stride_out * col] += result;
            else
              out[stride_out * col] = result;
          }
      }
    else if (transpose_matrix && n_rows == 3 && n_components == 1)
      {
        const Number2 *matrix_1 = matrix + n_columns;
        const Number2 *matrix_2 = matrix_1 + n_columns;
        const Number   x0 = in[0], x1 = in[stride_in], x2 = in[2 * stride_in];
        for (int col = 0; col < nn; ++col)
          {
            const Number result =
              matrix[col] * x0 + matrix_1[col] * x1 + matrix_2[col] * x2;
            if (add)
              out[stride_out * col] += result;
            else
              out[stride_out * col] = result;
          }
      }
    else if (std::abs(in - out) < std::min(stride_out * nn, stride_in * mm) &&
             n_components == 1)
      {
        Assert(mm <= 128,
               ExcNotImplemented("For large sizes, arrays may not overlap"));
        std::array<Number, 129> x;
        for (int i = 0; i < mm; ++i)
          x[i] = in[stride_in * i];

        for (int col = 0; col < nn; ++col)
          {
            Number res0;
            if (transpose_matrix == true)
              {
                res0 = matrix[col] * x[0];
                for (int i = 1; i < mm; ++i)
                  res0 += matrix[i * n_columns + col] * x[i];
              }
            else
              {
                res0 = matrix[col * n_columns] * x[0];
                for (int i = 1; i < mm; ++i)
                  res0 += matrix[col * n_columns + i] * x[i];
              }
            if (add)
              out[stride_out * col] += res0;
            else
              out[stride_out * col] = res0;
          }
      }
    else
      {
        const Number *in0 = in;
        const Number *in1 = n_components > 1 ? in + mm : nullptr;
        const Number *in2 = n_components > 2 ? in + 2 * mm : nullptr;

        Number *out0 = out;
        Number *out1 = n_components > 1 ? out + nn : nullptr;
        Number *out2 = n_components > 2 ? out + 2 * nn : nullptr;

        int nn_regular = (nn / 4) * 4;
        for (int col = 0; col < nn_regular; col += 4)
          {
            Number res[12];
            if (transpose_matrix == true)
              {
                const Number2 *matrix_ptr = matrix + col;
                const Number   a          = in0[0];
                res[0]                    = matrix_ptr[0] * a;
                res[1]                    = matrix_ptr[1] * a;
                res[2]                    = matrix_ptr[2] * a;
                res[3]                    = matrix_ptr[3] * a;

                if (n_components > 1)
                  {
                    const Number b = in1[0];
                    res[4]         = matrix_ptr[0] * b;
                    res[5]         = matrix_ptr[1] * b;
                    res[6]         = matrix_ptr[2] * b;
                    res[7]         = matrix_ptr[3] * b;
                  }

                if (n_components > 2)
                  {
                    const Number c = in2[0];
                    res[8]         = matrix_ptr[0] * c;
                    res[9]         = matrix_ptr[1] * c;
                    res[10]        = matrix_ptr[2] * c;
                    res[11]        = matrix_ptr[3] * c;
                  }

                matrix_ptr += n_columns;
                for (int i = 1; i < mm; ++i, matrix_ptr += n_columns)
                  {
                    const Number a = in0[stride_in * i];
                    res[0] += matrix_ptr[0] * a;
                    res[1] += matrix_ptr[1] * a;
                    res[2] += matrix_ptr[2] * a;
                    res[3] += matrix_ptr[3] * a;

                    if (n_components > 1)
                      {
                        const Number b = in1[stride_in * i];
                        res[4] += matrix_ptr[0] * b;
                        res[5] += matrix_ptr[1] * b;
                        res[6] += matrix_ptr[2] * b;
                        res[7] += matrix_ptr[3] * b;
                      }
                    if (n_components > 2)
                      {
                        const Number c = in2[stride_in * i];
                        res[8] += matrix_ptr[0] * c;
                        res[9] += matrix_ptr[1] * c;
                        res[10] += matrix_ptr[2] * c;
                        res[11] += matrix_ptr[3] * c;
                      }
                  }
              }
            else
              {
                const Number2 *matrix_0 = matrix + col * n_columns;
                const Number2 *matrix_1 = matrix + (col + 1) * n_columns;
                const Number2 *matrix_2 = matrix + (col + 2) * n_columns;
                const Number2 *matrix_3 = matrix + (col + 3) * n_columns;

                const Number a = in0[0];
                res[0]         = matrix_0[0] * a;
                res[1]         = matrix_1[0] * a;
                res[2]         = matrix_2[0] * a;
                res[3]         = matrix_3[0] * a;

                if (n_components > 1)
                  {
                    const Number b = in1[0];
                    res[4]         = matrix_0[0] * b;
                    res[5]         = matrix_1[0] * b;
                    res[6]         = matrix_2[0] * b;
                    res[7]         = matrix_3[0] * b;
                  }

                if (n_components > 2)
                  {
                    const Number c = in2[0];
                    res[8]         = matrix_0[0] * c;
                    res[9]         = matrix_1[0] * c;
                    res[10]        = matrix_2[0] * c;
                    res[11]        = matrix_3[0] * c;
                  }

                for (int i = 1; i < mm; ++i)
                  {
                    const Number a = in0[stride_in * i];
                    res[0] += matrix_0[i] * a;
                    res[1] += matrix_1[i] * a;
                    res[2] += matrix_2[i] * a;
                    res[3] += matrix_3[i] * a;

                    if (n_components > 1)
                      {
                        const Number b = in1[stride_in * i];
                        res[4] += matrix_0[i] * b;
                        res[5] += matrix_1[i] * b;
                        res[6] += matrix_2[i] * b;
                        res[7] += matrix_3[i] * b;
                      }

                    if (n_components > 2)
                      {
                        const Number c = in2[stride_in * i];
                        res[8] += matrix_0[i] * c;
                        res[9] += matrix_1[i] * c;
                        res[10] += matrix_2[i] * c;
                        res[11] += matrix_3[i] * c;
                      }
                  }
              }
            if (add)
              {
                out0[0] += res[0];
                out0[stride_out] += res[1];
                out0[2 * stride_out] += res[2];
                out0[3 * stride_out] += res[3];
                if (n_components > 1)
                  {
                    out1[0] += res[4];
                    out1[stride_out] += res[5];
                    out1[2 * stride_out] += res[6];
                    out1[3 * stride_out] += res[7];
                  }
                if (n_components > 2)
                  {
                    out2[0] += res[8];
                    out2[stride_out] += res[9];
                    out2[2 * stride_out] += res[10];
                    out2[3 * stride_out] += res[11];
                  }
              }
            else
              {
                out0[0]              = res[0];
                out0[stride_out]     = res[1];
                out0[2 * stride_out] = res[2];
                out0[3 * stride_out] = res[3];
                if (n_components > 1)
                  {
                    out1[0]              = res[4];
                    out1[stride_out]     = res[5];
                    out1[2 * stride_out] = res[6];
                    out1[3 * stride_out] = res[7];
                  }
                if (n_components > 2)
                  {
                    out2[0]              = res[8];
                    out2[stride_out]     = res[9];
                    out2[2 * stride_out] = res[10];
                    out2[3 * stride_out] = res[11];
                  }
              }
            out0 += 4 * stride_out;
            if (n_components > 1)
              out1 += 4 * stride_out;
            if (n_components > 2)
              out2 += 4 * stride_out;
          }
        if (nn - nn_regular == 3)
          {
            Number res0, res1, res2, res3, res4, res5, res6, res7, res8;
            if (transpose_matrix == true)
              {
                const Number2 *matrix_ptr = matrix + nn_regular;
                res0                      = matrix_ptr[0] * in0[0];
                res1                      = matrix_ptr[1] * in0[0];
                res2                      = matrix_ptr[2] * in0[0];
                if (n_components > 1)
                  {
                    res3 = matrix_ptr[0] * in1[0];
                    res4 = matrix_ptr[1] * in1[0];
                    res5 = matrix_ptr[2] * in1[0];
                  }
                if (n_components > 2)
                  {
                    res6 = matrix_ptr[0] * in2[0];
                    res7 = matrix_ptr[1] * in2[0];
                    res8 = matrix_ptr[2] * in2[0];
                  }
                matrix_ptr += n_columns;
                for (int i = 1; i < mm; ++i, matrix_ptr += n_columns)
                  {
                    res0 += matrix_ptr[0] * in0[stride_in * i];
                    res1 += matrix_ptr[1] * in0[stride_in * i];
                    res2 += matrix_ptr[2] * in0[stride_in * i];
                    if (n_components > 1)
                      {
                        res3 += matrix_ptr[0] * in1[stride_in * i];
                        res4 += matrix_ptr[1] * in1[stride_in * i];
                        res5 += matrix_ptr[2] * in1[stride_in * i];
                      }
                    if (n_components > 2)
                      {
                        res6 += matrix_ptr[0] * in2[stride_in * i];
                        res7 += matrix_ptr[1] * in2[stride_in * i];
                        res8 += matrix_ptr[2] * in2[stride_in * i];
                      }
                  }
              }
            else
              {
                const Number2 *matrix_0 = matrix + nn_regular * n_columns;
                const Number2 *matrix_1 = matrix + (nn_regular + 1) * n_columns;
                const Number2 *matrix_2 = matrix + (nn_regular + 2) * n_columns;

                res0 = matrix_0[0] * in0[0];
                res1 = matrix_1[0] * in0[0];
                res2 = matrix_2[0] * in0[0];
                if (n_components > 1)
                  {
                    res3 = matrix_0[0] * in1[0];
                    res4 = matrix_1[0] * in1[0];
                    res5 = matrix_2[0] * in1[0];
                  }
                if (n_components > 2)
                  {
                    res6 = matrix_0[0] * in2[0];
                    res7 = matrix_1[0] * in2[0];
                    res8 = matrix_2[0] * in2[0];
                  }
                for (int i = 1; i < mm; ++i)
                  {
                    res0 += matrix_0[i] * in0[stride_in * i];
                    res1 += matrix_1[i] * in0[stride_in * i];
                    res2 += matrix_2[i] * in0[stride_in * i];
                    if (n_components > 1)
                      {
                        res3 += matrix_0[i] * in1[stride_in * i];
                        res4 += matrix_1[i] * in1[stride_in * i];
                        res5 += matrix_2[i] * in1[stride_in * i];
                      }
                    if (n_components > 2)
                      {
                        res6 += matrix_0[i] * in2[stride_in * i];
                        res7 += matrix_1[i] * in2[stride_in * i];
                        res8 += matrix_2[i] * in2[stride_in * i];
                      }
                  }
              }
            if (add)
              {
                out0[0] += res0;
                out0[stride_out] += res1;
                out0[2 * stride_out] += res2;
                if (n_components > 1)
                  {
                    out1[0] += res3;
                    out1[stride_out] += res4;
                    out1[2 * stride_out] += res5;
                  }
                if (n_components > 2)
                  {
                    out2[0] += res6;
                    out2[stride_out] += res7;
                    out2[2 * stride_out] += res8;
                  }
              }
            else
              {
                out0[0]              = res0;
                out0[stride_out]     = res1;
                out0[2 * stride_out] = res2;
                if (n_components > 1)
                  {
                    out1[0]              = res3;
                    out1[stride_out]     = res4;
                    out1[2 * stride_out] = res5;
                  }
                if (n_components > 2)
                  {
                    out2[0]              = res6;
                    out2[stride_out]     = res7;
                    out2[2 * stride_out] = res8;
                  }
              }
          }
        else if (nn - nn_regular == 2)
          {
            Number res0, res1, res2, res3, res4, res5;
            if (transpose_matrix == true)
              {
                const Number2 *matrix_ptr = matrix + nn_regular;
                res0                      = matrix_ptr[0] * in0[0];
                res1                      = matrix_ptr[1] * in0[0];
                if (n_components > 1)
                  {
                    res2 = matrix_ptr[0] * in1[0];
                    res3 = matrix_ptr[1] * in1[0];
                  }
                if (n_components > 2)
                  {
                    res4 = matrix_ptr[0] * in2[0];
                    res5 = matrix_ptr[1] * in2[0];
                  }
                matrix_ptr += n_columns;
                for (int i = 1; i < mm; ++i, matrix_ptr += n_columns)
                  {
                    res0 += matrix_ptr[0] * in0[stride_in * i];
                    res1 += matrix_ptr[1] * in0[stride_in * i];
                    if (n_components > 1)
                      {
                        res2 += matrix_ptr[0] * in1[stride_in * i];
                        res3 += matrix_ptr[1] * in1[stride_in * i];
                      }
                    if (n_components > 2)
                      {
                        res4 += matrix_ptr[0] * in2[stride_in * i];
                        res5 += matrix_ptr[1] * in2[stride_in * i];
                      }
                  }
              }
            else
              {
                const Number2 *matrix_0 = matrix + nn_regular * n_columns;
                const Number2 *matrix_1 = matrix + (nn_regular + 1) * n_columns;

                res0 = matrix_0[0] * in0[0];
                res1 = matrix_1[0] * in0[0];
                if (n_components > 1)
                  {
                    res2 = matrix_0[0] * in1[0];
                    res3 = matrix_1[0] * in1[0];
                  }
                if (n_components > 2)
                  {
                    res4 = matrix_0[0] * in2[0];
                    res5 = matrix_1[0] * in2[0];
                  }
                for (int i = 1; i < mm; ++i)
                  {
                    res0 += matrix_0[i] * in0[stride_in * i];
                    res1 += matrix_1[i] * in0[stride_in * i];
                    if (n_components > 1)
                      {
                        res2 += matrix_0[i] * in1[stride_in * i];
                        res3 += matrix_1[i] * in1[stride_in * i];
                      }
                    if (n_components > 2)
                      {
                        res4 += matrix_0[i] * in2[stride_in * i];
                        res5 += matrix_1[i] * in2[stride_in * i];
                      }
                  }
              }
            if (add)
              {
                out0[0] += res0;
                out0[stride_out] += res1;
                if (n_components > 1)
                  {
                    out1[0] += res2;
                    out1[stride_out] += res3;
                  }
                if (n_components > 2)
                  {
                    out2[0] += res4;
                    out2[stride_out] += res5;
                  }
              }
            else
              {
                out0[0]          = res0;
                out0[stride_out] = res1;
                if (n_components > 1)
                  {
                    out1[0]          = res2;
                    out1[stride_out] = res3;
                  }
                if (n_components > 2)
                  {
                    out2[0]          = res4;
                    out2[stride_out] = res5;
                  }
              }
          }
        else if (nn - nn_regular == 1)
          {
            Number res0, res1, res2;
            if (transpose_matrix == true)
              {
                const Number2 *matrix_ptr = matrix + nn_regular;
                res0                      = matrix_ptr[0] * in0[0];
                if (n_components > 1)
                  res1 = matrix_ptr[0] * in1[0];
                if (n_components > 2)
                  res2 = matrix_ptr[0] * in2[0];
                matrix_ptr += n_columns;
                for (int i = 1; i < mm; ++i, matrix_ptr += n_columns)
                  {
                    res0 += matrix_ptr[0] * in0[stride_in * i];
                    if (n_components > 1)
                      res1 += matrix_ptr[0] * in1[stride_in * i];
                    if (n_components > 2)
                      res2 += matrix_ptr[0] * in2[stride_in * i];
                  }
              }
            else
              {
                const Number2 *matrix_ptr = matrix + nn_regular * n_columns;
                res0                      = matrix_ptr[0] * in0[0];
                if (n_components > 1)
                  res1 = matrix_ptr[0] * in1[0];
                if (n_components > 2)
                  res2 = matrix_ptr[0] * in2[0];
                for (int i = 1; i < mm; ++i)
                  {
                    res0 += matrix_ptr[i] * in0[stride_in * i];
                    if (n_components > 1)
                      res1 += matrix_ptr[i] * in1[stride_in * i];
                    if (n_components > 2)
                      res2 += matrix_ptr[i] * in2[stride_in * i];
                  }
              }
            if (add)
              {
                out0[0] += res0;
                if (n_components > 1)
                  out1[0] += res1;
                if (n_components > 2)
                  out2[0] += res2;
              }
            else
              {
                out0[0] = res0;
                if (n_components > 1)
                  out1[0] = res1;
                if (n_components > 2)
                  out2[0] = res2;
              }
          }
      }
  }



  /**
   * Internal evaluator specialized for "symmetric" finite elements, i.e.,
   * when the shape functions and quadrature points are symmetric about the
   * middle point, making the matrix entries the same when starting to read in
   * the (1,1) entry forward compared to the (N,N) entry backward.
   */
  template <EvaluatorVariant  variant,
            EvaluatorQuantity quantity,
            int               n_rows,
            int               n_columns,
            int               stride_in,
            int               stride_out,
            bool              transpose_matrix,
            bool              add,
            typename Number,
            typename Number2>
  std::enable_if_t<(variant == evaluate_symmetric), void>
  apply_matrix_vector_product(const Number2 *matrix,
                              const Number  *in,
                              Number        *out)
  {
    // We can only statically assert that one argument is non-zero because
    // face evaluation might instantiate some functions, so we need to use the
    // run-time assert to verify that we do not end up involuntarily.
    static_assert(n_rows > 0 || n_columns > 0,
                  "Specialization only for n_rows, n_columns > 0");
    Assert(n_rows > 0 && n_columns > 0,
           ExcInternalError("The evaluation needs n_rows, n_columns > 0, but " +
                            std::to_string(n_rows) + ", " +
                            std::to_string(n_columns) + " was passed!"));

    constexpr int mm     = transpose_matrix ? n_rows : n_columns,
                  nn     = transpose_matrix ? n_columns : n_rows;
    constexpr int n_cols = nn / 2;
    constexpr int mid    = mm / 2;

    std::array<Number, mm> x;
    for (int i = 0; i < mm; ++i)
      x[i] = in[stride_in * i];

    if (quantity == EvaluatorQuantity::value)
      {
        // In this case, the 1d shape values read (sorted lexicographically,
        // rows run over 1d dofs, columns over quadrature points):
        // Q2 --> [ 0.687  0 -0.087 ]
        //        [ 0.4    1  0.4   ]
        //        [-0.087  0  0.687 ]
        // Q3 --> [ 0.66   0.003  0.002  0.049 ]
        //        [ 0.521  1.005 -0.01  -0.230 ]
        //        [-0.230 -0.01   1.005  0.521 ]
        //        [ 0.049  0.002  0.003  0.66  ]
        // Q4 --> [ 0.658  0.022  0 -0.007 -0.032 ]
        //        [ 0.608  1.059  0  0.039  0.176 ]
        //        [-0.409 -0.113  1 -0.113 -0.409 ]
        //        [ 0.176  0.039  0  1.059  0.608 ]
        //        [-0.032 -0.007  0  0.022  0.658 ]
        //
        // In these matrices, we want to use avoid computations involving
        // zeros and ones and use the symmetry in entries starting from (1,1)
        // forward and (N,N) backward, respectively to reduce the number of
        // read operations.
        for (int col = 0; col < n_cols; ++col)
          {
            Number2 val0, val1;
            Number  res0, res1;
            if (transpose_matrix == true)
              {
                val0 = matrix[col];
                val1 = matrix[nn - 1 - col];
              }
            else
              {
                val0 = matrix[col * n_columns];
                val1 = matrix[(col + 1) * n_columns - 1];
              }
            if (mid > 0)
              {
                res0 = val0 * x[0];
                res1 = val1 * x[0];
                res0 += val1 * x[mm - 1];
                res1 += val0 * x[mm - 1];
                for (int ind = 1; ind < mid; ++ind)
                  {
                    if (transpose_matrix == true)
                      {
                        val0 = matrix[ind * n_columns + col];
                        val1 = matrix[ind * n_columns + nn - 1 - col];
                      }
                    else
                      {
                        val0 = matrix[col * n_columns + ind];
                        val1 = matrix[(col + 1) * n_columns - 1 - ind];
                      }
                    res0 += val0 * x[ind];
                    res1 += val1 * x[ind];
                    res0 += val1 * x[mm - 1 - ind];
                    res1 += val0 * x[mm - 1 - ind];
                  }
              }
            else
              res0 = res1 = Number();
            if (transpose_matrix == true)
              {
                if (mm % 2 == 1)
                  {
                    const Number tmp = matrix[mid * n_columns + col] * x[mid];
                    res0 += tmp;
                    res1 += tmp;
                  }
              }
            else
              {
                if (mm % 2 == 1 && nn % 2 == 0)
                  {
                    const Number tmp = matrix[col * n_columns + mid] * x[mid];
                    res0 += tmp;
                    res1 += tmp;
                  }
              }
            if (add)
              {
                out[stride_out * col] += res0;
                out[stride_out * (nn - 1 - col)] += res1;
              }
            else
              {
                out[stride_out * col]            = res0;
                out[stride_out * (nn - 1 - col)] = res1;
              }
          }
        if (transpose_matrix == true && nn % 2 == 1 && mm % 2 == 1)
          {
            if (add)
              out[stride_out * n_cols] += x[mid];
            else
              out[stride_out * n_cols] = x[mid];
          }
        else if (transpose_matrix == true && nn % 2 == 1)
          {
            Number res0;
            if (mid > 0)
              {
                res0 = matrix[n_cols] * (x[0] + x[mm - 1]);
                for (int ind = 1; ind < mid; ++ind)
                  {
                    const Number2 val0 = matrix[ind * n_columns + n_cols];
                    res0 += val0 * (x[ind] + in[mm - 1 - ind]);
                  }
              }
            else
              res0 = Number();
            if (add)
              out[stride_out * n_cols] += res0;
            else
              out[stride_out * n_cols] = res0;
          }
        else if (transpose_matrix == false && nn % 2 == 1)
          {
            Number res0;
            if (mid > 0)
              {
                res0 = matrix[n_cols * n_columns] * (x[0] + x[mm - 1]);
                for (int ind = 1; ind < mid; ++ind)
                  {
                    const Number2 val0 = matrix[n_cols * n_columns + ind];
                    res0 += val0 * (x[ind] + x[mm - 1 - ind]);
                  }
                if (mm % 2)
                  res0 += x[mid];
              }
            else
              res0 = in[0];
            if (add)
              out[stride_out * n_cols] += res0;
            else
              out[stride_out * n_cols] = res0;
          }
      }
    else if (quantity == EvaluatorQuantity::gradient)
      {
        // For the specialized loop used for gradient computations we again
        // exploit symmetries according to the following entries (sorted
        // lexicographically, rows run over 1d dofs, columns over quadrature
        // points):
        // Q2 --> [-2.549 -1  0.549 ]
        //        [ 3.098  0 -3.098 ]
        //        [-0.549  1  2.549 ]
        // Q3 --> [-4.315 -1.03  0.5  -0.44  ]
        //        [ 6.07  -1.44 -2.97  2.196 ]
        //        [-2.196  2.97  1.44 -6.07  ]
        //        [ 0.44  -0.5   1.03  4.315 ]
        // Q4 --> [-6.316 -1.3    0.333 -0.353  0.413 ]
        //        [10.111 -2.76  -2.667  2.066 -2.306 ]
        //        [-5.688  5.773  0     -5.773  5.688 ]
        //        [ 2.306 -2.066  2.667  2.76 -10.111 ]
        //        [-0.413  0.353 -0.333 -0.353  0.413 ]
        for (int col = 0; col < n_cols; ++col)
          {
            Number2 val0, val1;
            Number  res0, res1;
            if (transpose_matrix == true)
              {
                val0 = matrix[col];
                val1 = matrix[nn - 1 - col];
              }
            else
              {
                val0 = matrix[col * n_columns];
                val1 = matrix[(nn - col - 1) * n_columns];
              }
            if (mid > 0)
              {
                res0 = val0 * x[0];
                res1 = val1 * x[0];
                res0 -= val1 * x[mm - 1];
                res1 -= val0 * x[mm - 1];
                for (int ind = 1; ind < mid; ++ind)
                  {
                    if (transpose_matrix == true)
                      {
                        val0 = matrix[ind * n_columns + col];
                        val1 = matrix[ind * n_columns + nn - 1 - col];
                      }
                    else
                      {
                        val0 = matrix[col * n_columns + ind];
                        val1 = matrix[(nn - col - 1) * n_columns + ind];
                      }
                    res0 += val0 * x[ind];
                    res1 += val1 * x[ind];
                    res0 -= val1 * x[mm - 1 - ind];
                    res1 -= val0 * x[mm - 1 - ind];
                  }
              }
            else
              res0 = res1 = Number();
            if (mm % 2 == 1)
              {
                if (transpose_matrix == true)
                  val0 = matrix[mid * n_columns + col];
                else
                  val0 = matrix[col * n_columns + mid];
                const Number tmp = val0 * x[mid];
                res0 += tmp;
                res1 -= tmp;
              }
            if (add)
              {
                out[stride_out * col] += res0;
                out[stride_out * (nn - 1 - col)] += res1;
              }
            else
              {
                out[stride_out * col]            = res0;
                out[stride_out * (nn - 1 - col)] = res1;
              }
          }
        if (nn % 2 == 1)
          {
            Number2 val0;
            Number  res0;
            if (transpose_matrix == true)
              val0 = matrix[n_cols];
            else
              val0 = matrix[n_cols * n_columns];
            res0 = val0 * (x[0] - x[mm - 1]);
            for (int ind = 1; ind < mid; ++ind)
              {
                if (transpose_matrix == true)
                  val0 = matrix[ind * n_columns + n_cols];
                else
                  val0 = matrix[n_cols * n_columns + ind];
                Number in1 = val0 * (x[ind] - x[mm - 1 - ind]);
                res0 += in1;
              }
            if (add)
              out[stride_out * n_cols] += res0;
            else
              out[stride_out * n_cols] = res0;
          }
      }
    else
      {
        // Hessians are almost the same as values, apart from some missing '1'
        // entries
        for (int col = 0; col < n_cols; ++col)
          {
            Number2 val0, val1;
            Number  res0, res1;
            if (transpose_matrix == true)
              {
                val0 = matrix[col];
                val1 = matrix[nn - 1 - col];
              }
            else
              {
                val0 = matrix[col * n_columns];
                val1 = matrix[(col + 1) * n_columns - 1];
              }
            if (mid > 0)
              {
                res0 = val0 * x[0];
                res1 = val1 * x[0];
                res0 += val1 * x[mm - 1];
                res1 += val0 * x[mm - 1];
                for (int ind = 1; ind < mid; ++ind)
                  {
                    if (transpose_matrix == true)
                      {
                        val0 = matrix[ind * n_columns + col];
                        val1 = matrix[ind * n_columns + nn - 1 - col];
                      }
                    else
                      {
                        val0 = matrix[col * n_columns + ind];
                        val1 = matrix[(col + 1) * n_columns - 1 - ind];
                      }
                    res0 += val0 * x[ind];
                    res1 += val1 * x[ind];
                    res0 += val1 * x[mm - 1 - ind];
                    res1 += val0 * x[mm - 1 - ind];
                  }
              }
            else
              res0 = res1 = Number();
            if (mm % 2 == 1)
              {
                if (transpose_matrix == true)
                  val0 = matrix[mid * n_columns + col];
                else
                  val0 = matrix[col * n_columns + mid];
                const Number tmp = val0 * x[mid];
                res0 += tmp;
                res1 += tmp;
              }
            if (add)
              {
                out[stride_out * col] += res0;
                out[stride_out * (nn - 1 - col)] += res1;
              }
            else
              {
                out[stride_out * col]            = res0;
                out[stride_out * (nn - 1 - col)] = res1;
              }
          }
        if (nn % 2 == 1)
          {
            Number2 val0;
            Number  res0;
            if (transpose_matrix == true)
              val0 = matrix[n_cols];
            else
              val0 = matrix[n_cols * n_columns];
            if (mid > 0)
              {
                res0 = val0 * (x[0] + x[mm - 1]);
                for (int ind = 1; ind < mid; ++ind)
                  {
                    if (transpose_matrix == true)
                      val0 = matrix[ind * n_columns + n_cols];
                    else
                      val0 = matrix[n_cols * n_columns + ind];
                    Number in1 = val0 * (x[ind] + x[mm - 1 - ind]);
                    res0 += in1;
                  }
              }
            else
              res0 = Number();
            if (mm % 2 == 1)
              {
                if (transpose_matrix == true)
                  val0 = matrix[mid * n_columns + n_cols];
                else
                  val0 = matrix[n_cols * n_columns + mid];
                res0 += val0 * x[mid];
              }
            if (add)
              out[stride_out * n_cols] += res0;
            else
              out[stride_out * n_cols] = res0;
          }
      }
  }



  /**
   * Internal evaluator specialized for "symmetric" finite elements in the
   * evenodd matrix format.
   *
   * This function implements a different approach to the symmetric case for
   * values, gradients, and Hessians as in the above matrices: It is possible
   * to reduce the cost per dimension from N^2 to N^2/2, where N is the number
   * of 1d dofs (there are only N^2/2 different entries in the shape matrix,
   * so this is plausible). The approach is based on the idea of applying the
   * operator on the even and odd part of the input vectors separately, given
   * that the basis of shape functions evaluated at quadrature points is
   * symmetric. This method is presented e.g. in the book "Implementing
   * Spectral Methods for Partial Differential Equations" by David A. Kopriva,
   * Springer, 2009, section 3.5.3 (Even-Odd-Decomposition). Even though the
   * experiments in the book say that the method is not efficient for N<20, it
   * is more efficient in the context where the loop bounds are compile-time
   * constants (templates).
   */
  template <EvaluatorVariant  variant,
            EvaluatorQuantity quantity,
            int               n_rows_static,
            int               n_columns_static,
            int               stride_in_static,
            int               stride_out_static,
            bool              transpose_matrix,
            bool              add,
            typename Number,
            typename Number2>
#ifndef DEBUG
  inline DEAL_II_ALWAYS_INLINE
#endif
    std::enable_if_t<(variant == evaluate_evenodd), void>
    apply_matrix_vector_product(const Number2 *DEAL_II_RESTRICT matrix,
                                const Number                   *in,
                                Number                         *out,
                                int n_rows_runtime     = 0,
                                int n_columns_runtime  = 0,
                                int stride_in_runtime  = 0,
                                int stride_out_runtime = 0)
  {
    static_assert(n_rows_static >= 0 && n_columns_static >= 0,
                  "Negative loop ranges are not allowed!");

    const int n_rows = n_rows_static == 0 ? n_rows_runtime : n_rows_static;
    const int n_columns =
      n_rows_static == 0 ? n_columns_runtime : n_columns_static;
    const int stride_in =
      stride_in_static == 0 ? stride_in_runtime : stride_in_static;
    const int stride_out =
      stride_out_static == 0 ? stride_out_runtime : stride_out_static;

    Assert(n_rows > 0 && n_columns > 0,
           ExcInternalError("The evaluation needs n_rows, n_columns > 0, but " +
                            std::to_string(n_rows) + ", " +
                            std::to_string(n_columns) + " was passed!"));

    const int mm     = transpose_matrix ? n_rows : n_columns,
              nn     = transpose_matrix ? n_columns : n_rows;
    const int n_half = nn / 2;
    const int m_half = mm / 2;

    constexpr int array_length =
      (n_rows_static == 0) ?
        16 // for non-templated execution
        :
        (1 + (transpose_matrix ? n_rows_static : n_columns_static) / 2);
    const int offset = (n_columns + 1) / 2;

    Assert(m_half <= array_length, ExcNotImplemented());

    std::array<Number, array_length> xp, xm;
    for (int i = 0; i < m_half; ++i)
      {
        if (transpose_matrix == true && quantity == EvaluatorQuantity::gradient)
          {
            xp[i] = in[stride_in * i] - in[stride_in * (mm - 1 - i)];
            xm[i] = in[stride_in * i] + in[stride_in * (mm - 1 - i)];
          }
        else
          {
            xp[i] = in[stride_in * i] + in[stride_in * (mm - 1 - i)];
            xm[i] = in[stride_in * i] - in[stride_in * (mm - 1 - i)];
          }
      }
    Number xmid = in[stride_in * m_half];
    for (int col = 0; col < n_half; ++col)
      {
        Number r0, r1;
        if (m_half > 0)
          {
            if (transpose_matrix == true)
              {
                r0 = matrix[col] * xp[0];
                r1 = matrix[(n_rows - 1) * offset + col] * xm[0];
              }
            else
              {
                r0 = matrix[col * offset] * xp[0];
                r1 = matrix[(n_rows - 1 - col) * offset] * xm[0];
              }
            for (int ind = 1; ind < m_half; ++ind)
              {
                if (transpose_matrix == true)
                  {
                    r0 += matrix[ind * offset + col] * xp[ind];
                    r1 += matrix[(n_rows - 1 - ind) * offset + col] * xm[ind];
                  }
                else
                  {
                    r0 += matrix[col * offset + ind] * xp[ind];
                    r1 += matrix[(n_rows - 1 - col) * offset + ind] * xm[ind];
                  }
              }
          }
        else
          r0 = r1 = Number();
        if (mm % 2 == 1 && transpose_matrix == true)
          {
            if (quantity == EvaluatorQuantity::gradient)
              r1 += matrix[m_half * offset + col] * xmid;
            else
              r0 += matrix[m_half * offset + col] * xmid;
          }
        else if (mm % 2 == 1 &&
                 (nn % 2 == 0 || quantity != EvaluatorQuantity::value ||
                  mm == 3))
          r0 += matrix[col * offset + m_half] * xmid;

        if (add)
          {
            out[stride_out * col] += r0 + r1;
            if (quantity == EvaluatorQuantity::gradient &&
                transpose_matrix == false)
              out[stride_out * (nn - 1 - col)] += r1 - r0;
            else
              out[stride_out * (nn - 1 - col)] += r0 - r1;
          }
        else
          {
            out[stride_out * col] = r0 + r1;
            if (quantity == EvaluatorQuantity::gradient &&
                transpose_matrix == false)
              out[stride_out * (nn - 1 - col)] = r1 - r0;
            else
              out[stride_out * (nn - 1 - col)] = r0 - r1;
          }
      }
    if (quantity == EvaluatorQuantity::value && transpose_matrix == true &&
        nn % 2 == 1 && mm % 2 == 1 && mm > 3)
      {
        if (add)
          out[stride_out * n_half] += matrix[m_half * offset + n_half] * xmid;
        else
          out[stride_out * n_half] = matrix[m_half * offset + n_half] * xmid;
      }
    else if (transpose_matrix == true && nn % 2 == 1)
      {
        Number r0;
        if (m_half > 0)
          {
            r0 = matrix[n_half] * xp[0];
            for (int ind = 1; ind < m_half; ++ind)
              r0 += matrix[ind * offset + n_half] * xp[ind];
          }
        else
          r0 = Number();
        if (quantity != EvaluatorQuantity::gradient && mm % 2 == 1)
          r0 += matrix[m_half * offset + n_half] * xmid;

        if (add)
          out[stride_out * n_half] += r0;
        else
          out[stride_out * n_half] = r0;
      }
    else if (transpose_matrix == false && nn % 2 == 1)
      {
        Number r0;
        if (m_half > 0)
          {
            if (quantity == EvaluatorQuantity::gradient)
              {
                r0 = matrix[n_half * offset] * xm[0];
                for (int ind = 1; ind < m_half; ++ind)
                  r0 += matrix[n_half * offset + ind] * xm[ind];
              }
            else
              {
                r0 = matrix[n_half * offset] * xp[0];
                for (int ind = 1; ind < m_half; ++ind)
                  r0 += matrix[n_half * offset + ind] * xp[ind];
              }
          }
        else
          r0 = Number();

        if (quantity != EvaluatorQuantity::gradient && mm % 2 == 1)
          r0 += matrix[n_half * offset + m_half] * xmid;

        if (add)
          out[stride_out * n_half] += r0;
        else
          out[stride_out * n_half] = r0;
      }
  }



  /**
   * Internal evaluator specialized for "symmetric" finite elements in the
   * evenodd matrix format with run-time bounds.
   */
  template <EvaluatorVariant  variant,
            EvaluatorQuantity quantity,
            bool              transpose_matrix,
            bool              add,
            bool              consider_strides,
            typename Number,
            typename Number2>
  std::enable_if_t<(variant == evaluate_evenodd), void>
  apply_matrix_vector_product(const Number2 *matrix,
                              const Number  *in,
                              Number        *out,
                              int            n_rows,
                              int            n_columns,
                              int            stride_in,
                              int            stride_out)
  {
    apply_matrix_vector_product<evaluate_evenodd,
                                quantity,
                                0,
                                0,
                                consider_strides ? 0 : 1,
                                consider_strides ? 0 : 1,
                                transpose_matrix,
                                add>(
      matrix, in, out, n_rows, n_columns, stride_in, stride_out);
  }



  /**
   * Internal evaluator specialized for "symmetric" finite elements in the
   * symmetric_hierarchical matrix format.
   *
   * This class implements an approach similar to the even-odd decomposition
   * but with a different type of symmetry. In this case, we assume that a
   * single shape function already shows the symmetry over the quadrature
   * points, rather than the complete basis that is considered in the even-odd
   * case. In particular, we assume that the shape functions are ordered as in
   * the Legendre basis, with symmetric shape functions in the even slots
   * (rows of the values array) and point-symmetric in the odd slots. Like the
   * even-odd decomposition, the number of operations are N^2/2 rather than
   * N^2 FMAs (fused multiply-add), where N is the number of 1d dofs. The
   * difference is in the way the input and output quantities are symmetrized.
   */
  template <EvaluatorVariant  variant,
            EvaluatorQuantity quantity,
            int               n_rows,
            int               n_columns,
            int               stride_in,
            int               stride_out,
            bool              transpose_matrix,
            bool              add,
            typename Number,
            typename Number2>
  std::enable_if_t<(variant == evaluate_symmetric_hierarchical), void>
  apply_matrix_vector_product(const Number2 *matrix,
                              const Number  *in,
                              Number        *out)
  {
    static_assert(n_rows > 0 && n_columns > 0,
                  "Specialization requires n_rows, n_columns > 0");

    constexpr bool evaluate_antisymmetric =
      (quantity == EvaluatorQuantity::gradient);

    constexpr int mm     = transpose_matrix ? n_rows : n_columns,
                  nn     = transpose_matrix ? n_columns : n_rows;
    constexpr int n_half = nn / 2;
    constexpr int m_half = mm / 2;

    if (transpose_matrix)
      {
        std::array<Number, mm> x;
        for (unsigned int i = 0; i < mm; ++i)
          x[i] = in[stride_in * i];
        for (unsigned int col = 0; col < n_half; ++col)
          {
            Number r0, r1;
            if (m_half > 0)
              {
                r0 = matrix[col] * x[0];
                r1 = matrix[col + n_columns] * x[1];
                for (unsigned int ind = 1; ind < m_half; ++ind)
                  {
                    r0 += matrix[col + 2 * ind * n_columns] * x[2 * ind];
                    r1 +=
                      matrix[col + (2 * ind + 1) * n_columns] * x[2 * ind + 1];
                  }
              }
            else
              r0 = r1 = Number();
            if (mm % 2 == 1)
              r0 += matrix[col + (mm - 1) * n_columns] * x[mm - 1];
            if (add)
              {
                out[stride_out * col] += r0 + r1;
                if (evaluate_antisymmetric)
                  out[stride_out * (nn - 1 - col)] += r1 - r0;
                else
                  out[stride_out * (nn - 1 - col)] += r0 - r1;
              }
            else
              {
                out[stride_out * col] = r0 + r1;
                if (evaluate_antisymmetric)
                  out[stride_out * (nn - 1 - col)] = r1 - r0;
                else
                  out[stride_out * (nn - 1 - col)] = r0 - r1;
              }
          }
        if (nn % 2 == 1)
          {
            Number             r0;
            const unsigned int shift = evaluate_antisymmetric ? 1 : 0;
            if (m_half > 0)
              {
                r0 = matrix[n_half + shift * n_columns] * x[shift];
                for (unsigned int ind = 1; ind < m_half; ++ind)
                  r0 += matrix[n_half + (2 * ind + shift) * n_columns] *
                        x[2 * ind + shift];
              }
            else
              r0 = 0;
            if (!evaluate_antisymmetric && mm % 2 == 1)
              r0 += matrix[n_half + (mm - 1) * n_columns] * x[mm - 1];
            if (add)
              out[stride_out * n_half] += r0;
            else
              out[stride_out * n_half] = r0;
          }
      }
    else
      {
        std::array<Number, m_half + 1> xp, xm;
        for (int i = 0; i < m_half; ++i)
          if (!evaluate_antisymmetric)
            {
              xp[i] = in[stride_in * i] + in[stride_in * (mm - 1 - i)];
              xm[i] = in[stride_in * i] - in[stride_in * (mm - 1 - i)];
            }
          else
            {
              xp[i] = in[stride_in * i] - in[stride_in * (mm - 1 - i)];
              xm[i] = in[stride_in * i] + in[stride_in * (mm - 1 - i)];
            }
        if (mm % 2 == 1)
          xp[m_half] = in[stride_in * m_half];
        for (unsigned int col = 0; col < n_half; ++col)
          {
            Number r0, r1;
            if (m_half > 0)
              {
                r0 = matrix[2 * col * n_columns] * xp[0];
                r1 = matrix[(2 * col + 1) * n_columns] * xm[0];
                for (unsigned int ind = 1; ind < m_half; ++ind)
                  {
                    r0 += matrix[2 * col * n_columns + ind] * xp[ind];
                    r1 += matrix[(2 * col + 1) * n_columns + ind] * xm[ind];
                  }
              }
            else
              r0 = r1 = Number();
            if (mm % 2 == 1)
              {
                if (evaluate_antisymmetric)
                  r1 += matrix[(2 * col + 1) * n_columns + m_half] * xp[m_half];
                else
                  r0 += matrix[2 * col * n_columns + m_half] * xp[m_half];
              }
            if (add)
              {
                out[stride_out * (2 * col)] += r0;
                out[stride_out * (2 * col + 1)] += r1;
              }
            else
              {
                out[stride_out * (2 * col)]     = r0;
                out[stride_out * (2 * col + 1)] = r1;
              }
          }
        if (nn % 2 == 1)
          {
            Number r0;
            if (m_half > 0)
              {
                r0 = matrix[(nn - 1) * n_columns] * xp[0];
                for (unsigned int ind = 1; ind < m_half; ++ind)
                  r0 += matrix[(nn - 1) * n_columns + ind] * xp[ind];
              }
            else
              r0 = Number();
            if (mm % 2 == 1 && !evaluate_antisymmetric)
              r0 += matrix[(nn - 1) * n_columns + m_half] * xp[m_half];
            if (add)
              out[stride_out * (nn - 1)] += r0;
            else
              out[stride_out * (nn - 1)] = r0;
          }
      }
  }



  /**
   * Generic evaluator framework that valuates the given shape data in general
   * dimensions using the tensor product form. Depending on the particular
   * layout in the matrix entries, this corresponds to a usual matrix-matrix
   * product or a matrix-matrix product including some symmetries. The actual
   * work is implemented by functions of type apply_matrix_vector_product
   * working on a single dimension, controlled by suitable strides, using the
   * kernel specified via variant.
   *
   * @tparam dim Space dimension in which this class is applied
   * @tparam n_rows Number of rows in the transformation matrix, which corresponds
   *                to the number of 1d shape functions in the usual tensor
   *                contraction setting
   * @tparam n_columns Number of columns in the transformation matrix, which
   *                   corresponds to the number of 1d shape functions in the
   *                   usual tensor contraction setting
   * @tparam Number Abstract number type for input and output arrays
   * @tparam Number2 Abstract number type for coefficient arrays (defaults to
   *                 same type as the input/output arrays); must implement
   *                 operator* with Number and produce Number as an output to
   *                 be a valid type
   */
  template <EvaluatorVariant variant,
            int              dim,
            int              n_rows,
            int              n_columns,
            typename Number,
            typename Number2 = Number>
  struct EvaluatorTensorProduct
  {
    static constexpr unsigned int n_rows_of_product =
      Utilities::pow(n_rows, dim);
    static constexpr unsigned int n_columns_of_product =
      Utilities::pow(n_columns, dim);

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other pointer
     */
    EvaluatorTensorProduct()
      : shape_values(nullptr)
      , shape_gradients(nullptr)
      , shape_hessians(nullptr)
    {}

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct(const AlignedVector<Number2> &shape_values,
                           const AlignedVector<Number2> &shape_gradients,
                           const AlignedVector<Number2> &shape_hessians,
                           const unsigned int = 0,
                           const unsigned int = 0)
      : shape_values(shape_values.begin())
      , shape_gradients(shape_gradients.begin())
      , shape_hessians(shape_hessians.begin())
    {
      if (variant == evaluate_evenodd)
        {
          if (!shape_values.empty())
            AssertDimension(shape_values.size(),
                            n_rows * ((n_columns + 1) / 2));
          if (!shape_gradients.empty())
            AssertDimension(shape_gradients.size(),
                            n_rows * ((n_columns + 1) / 2));
          if (!shape_hessians.empty())
            AssertDimension(shape_hessians.size(),
                            n_rows * ((n_columns + 1) / 2));
        }
      else
        {
          Assert(shape_values.empty() ||
                   shape_values.size() == n_rows * n_columns,
                 ExcDimensionMismatch(shape_values.size(), n_rows * n_columns));
          Assert(shape_gradients.empty() ||
                   shape_gradients.size() == n_rows * n_columns,
                 ExcDimensionMismatch(shape_gradients.size(),
                                      n_rows * n_columns));
          Assert(shape_hessians.empty() ||
                   shape_hessians.size() == n_rows * n_columns,
                 ExcDimensionMismatch(shape_hessians.size(),
                                      n_rows * n_columns));
        }
    }

    /**
     * Constructor, taking the data from ShapeInfo via raw pointers
     */
    EvaluatorTensorProduct(const Number2     *shape_values,
                           const Number2     *shape_gradients,
                           const Number2     *shape_hessians,
                           const unsigned int dummy1 = 0,
                           const unsigned int dummy2 = 0)
      : shape_values(shape_values)
      , shape_gradients(shape_gradients)
      , shape_hessians(shape_hessians)
    {
      (void)dummy1;
      (void)dummy2;
    }

    /**
     * This interpolates values with sum factorization, according to the
     * following parameters:
     *
     * @tparam direction The direction along which the one-dimensional
     * operations should be performed, using 0 as the fastest running
     * direction x, 1 for y, and so on.
     * @tparam contract_over_rows Describes whether the interpolation is over
     * the rows or the columns in the underlying interpolation matrix. With
     * the chosen convention in the data fields, `contract_over_rows==true`
     * means interpolation from DoF values to quadrature points, whereas using
     * the argument `false` will perform summation over quadrature points to
     * produce integrals for each test function.
     * @tparam add Specify whether to add into the output array or overwrite
     * the previous content.
     * @tparam stride This parameter can specify an additional stride in the
     * array associated with quadrature points (`in` for `dof_to_quad==false`,
     * otherwise `out`) from consecutive points in the x-direction. This is
     * used to place results from different interpolation steps next to each
     * other in memory.
     *
     * @param in Input array for the operation, needs to be backed up by a
     * sufficiently large memory region.
     * @param out Array holding the result of the sum factorization operation.
     */
    template <int direction, bool contract_over_rows, bool add, int stride = 1>
    void
    values(const Number in[], Number out[]) const
    {
      constexpr EvaluatorQuantity value_type = EvaluatorQuantity::value;
      apply<direction, contract_over_rows, add, false, value_type, stride>(
        shape_values, in, out);
    }

    /**
     * This interpolates gradients with sum factorization, based on the second
     * argument given to the constructor of this class. For the documentation
     * of the template and function parameters, see the values function.
     */
    template <int direction, bool contract_over_rows, bool add, int stride = 1>
    void
    gradients(const Number in[], Number out[]) const
    {
      constexpr EvaluatorQuantity gradient_type =
        (variant == evaluate_general ? EvaluatorQuantity::value :
                                       EvaluatorQuantity::gradient);
      apply<direction, contract_over_rows, add, false, gradient_type, stride>(
        shape_gradients, in, out);
    }

    /**
     * This interpolates hessians with sum factorization, based on the third
     * argument given to the constructor of this class. For the documentation
     * of the template and function parameters, see the values function.
     */
    template <int direction, bool contract_over_rows, bool add>
    void
    hessians(const Number in[], Number out[]) const
    {
      constexpr EvaluatorQuantity hessian_type =
        (((variant == evaluate_general) |
          (variant == evaluate_symmetric_hierarchical)) ?
           EvaluatorQuantity::value :
           EvaluatorQuantity::hessian);
      apply<direction, contract_over_rows, add, false, hessian_type>(
        shape_hessians, in, out);
    }

    /**
     * A variant of interpolation with sum factorization that only applies the
     * operation to a single 1d line, leaving all other entries untouched,
     * rather than expanding the loop over all other directions of a tensor
     * product mesh. For the documentation of the template and function
     * parameters, see the other values function.
     */
    template <int direction, bool contract_over_rows, bool add>
    void
    values_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_values != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, true, EvaluatorQuantity::value>(
        shape_values, in, out);
    }

    /**
     * A variant of interpolation with sum factorization that only applies the
     * gradient operation to a single 1d line, leaving all other entries
     * untouched, rather than expanding the loop over all other directions of
     * a tensor product mesh. For the documentation of the template and
     * function parameters, see the values function.
     */
    template <int direction, bool contract_over_rows, bool add>
    void
    gradients_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_gradients != nullptr, ExcNotInitialized());
      constexpr EvaluatorQuantity gradient_type =
        (variant == evaluate_general ? EvaluatorQuantity::value :
                                       EvaluatorQuantity::gradient);
      apply<direction, contract_over_rows, add, true, gradient_type>(
        shape_gradients, in, out);
    }

    /**
     * A variant of interpolation with sum factorization that only applies the
     * Hessian operation to a single 1d line, leaving all other entries
     * untouched, rather than expanding the loop over all other directions of
     * a tensor product mesh. For the documentation of the template and
     * function parameters, see the values function.
     */
    template <int direction, bool contract_over_rows, bool add>
    void
    hessians_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_hessians != nullptr, ExcNotInitialized());
      constexpr EvaluatorQuantity hessian_type =
        (((variant == evaluate_general) |
          (variant == evaluate_symmetric_hierarchical)) ?
           EvaluatorQuantity::value :
           EvaluatorQuantity::hessian);
      apply<direction, contract_over_rows, add, true, hessian_type>(
        shape_hessians, in, out);
    }

    /**
     * This function applies the tensor product kernel with sum factorization,
     * corresponding to a matrix-vector multiplication of 1d stripes, along
     * the given @p direction of the tensor data in the input array. This
     * function allows the @p in and @p out arrays to alias for the case
     * n_rows == n_columns, i.e., it is safe to perform the contraction in
     * place where @p in and @p out point to the same address. For the case
     * `n_rows != n_columns`, the output is in general not correct.
     *
     * @tparam direction Direction that is evaluated
     * @tparam contract_over_rows If true, the tensor contraction sums
     *                            over the rows in the given @p shape_data
     *                            array, otherwise it sums over the columns
     * @tparam add If true, the result is added to the output vector, else
     *             the computed values overwrite the content in the output
     * @tparam one_line If true, the kernel is only applied along a single 1d
     *                  stripe within a dim-dimensional tensor, not the full
     *                  n_rows^dim points as in the @p false case.
     * @tparam quantity Specify whether values, gradients or Hessians should
     *                  be interpolated, allowing specialized algorithms
     *                  for some class template parameters of `variant` to
     *                  find the right path.
     * @tparam stride This parameter enables to place the result of the
     *                tensor product evaluation in the output array (if
     *                `contract_over_rows == true`) or input array (if
     *                `contract_over_rows == false`) with additional strides
     *                between adjacent points in x direction, which is used
     *                to group all components of a gradient adjacent in
     *                memory. If the stride is one, the data will form a
     *                contiguous range in memory.
     *
     * @param shape_data Transformation matrix with @p n_rows rows and
     *                   @p n_columns columns, stored in row-major format.
     * @param in Pointer to the start of the input data vector.
     * @param out Pointer to the start of the output data vector.
     */
    template <int               direction,
              bool              contract_over_rows,
              bool              add,
              bool              one_line = false,
              EvaluatorQuantity quantity = EvaluatorQuantity::value,
              int               stride   = 1>
    static void
    apply(const Number2 *DEAL_II_RESTRICT shape_data,
          const Number                   *in,
          Number                         *out);

  private:
    const Number2 *shape_values;
    const Number2 *shape_gradients;
    const Number2 *shape_hessians;
  };



  template <EvaluatorVariant variant,
            int              dim,
            int              n_rows,
            int              n_columns,
            typename Number,
            typename Number2>
  template <int               direction,
            bool              contract_over_rows,
            bool              add,
            bool              one_line,
            EvaluatorQuantity quantity,
            int               stride>
  inline void
  EvaluatorTensorProduct<variant, dim, n_rows, n_columns, Number, Number2>::
    apply(const Number2 *DEAL_II_RESTRICT shape_data,
          const Number                   *in,
          Number                         *out)
  {
    static_assert(one_line == false || direction == dim - 1,
                  "Single-line evaluation only works for direction=dim-1.");
    Assert(shape_data != nullptr,
           ExcMessage(
             "The given array shape_data must not be the null pointer!"));
    Assert(dim == direction + 1 || one_line == true || n_rows == n_columns ||
             in != out,
           ExcMessage("In-place operation only supported for "
                      "n_rows==n_columns or single-line interpolation"));
    AssertIndexRange(direction, dim);
    constexpr int mm = contract_over_rows ? n_rows : n_columns,
                  nn = contract_over_rows ? n_columns : n_rows;

    constexpr int stride_operation = Utilities::pow(n_columns, direction);
    constexpr int n_blocks1        = one_line ? 1 : stride_operation;
    constexpr int n_blocks2 =
      Utilities::pow(n_rows, (direction >= dim) ? 0 : (dim - direction - 1));

    constexpr int stride_in  = !contract_over_rows ? stride : 1;
    constexpr int stride_out = contract_over_rows ? stride : 1;
    for (int i2 = 0; i2 < n_blocks2; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks1; ++i1)
          {
            apply_matrix_vector_product<variant,
                                        quantity,
                                        n_rows,
                                        n_columns,
                                        stride_operation * stride_in,
                                        stride_operation * stride_out,
                                        contract_over_rows,
                                        add>(shape_data, in, out);

            if (one_line == false)
              {
                in += stride_in;
                out += stride_out;
              }
          }
        if (one_line == false)
          {
            in += stride_operation * (mm - 1) * stride_in;
            out += stride_operation * (nn - 1) * stride_out;
          }
      }
  }



  /**
   * Internal evaluator for shape function using the tensor product form
   * of the basis functions. The same as the other templated class but
   * without making use of template arguments and variable loop bounds
   * instead.
   *
   * @tparam dim Space dimension in which this class is applied
   * @tparam Number Abstract number type for input and output arrays
   * @tparam Number2 Abstract number type for coefficient arrays (defaults to
   *                 same type as the input/output arrays); must implement
   *                 operator* with Number and produce Number as an output to
   *                 be a valid type
   */
  template <EvaluatorVariant variant,
            int              dim,
            typename Number,
            typename Number2>
  struct EvaluatorTensorProduct<variant, dim, 0, 0, Number, Number2>
  {
    static constexpr unsigned int n_rows_of_product =
      numbers::invalid_unsigned_int;
    static constexpr unsigned int n_columns_of_product =
      numbers::invalid_unsigned_int;

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other constructor
     */
    EvaluatorTensorProduct()
      : shape_values(nullptr)
      , shape_gradients(nullptr)
      , shape_hessians(nullptr)
      , n_rows(numbers::invalid_unsigned_int)
      , n_columns(numbers::invalid_unsigned_int)
    {}

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct(const AlignedVector<Number2> &shape_values,
                           const AlignedVector<Number2> &shape_gradients,
                           const AlignedVector<Number2> &shape_hessians,
                           const unsigned int            n_rows    = 0,
                           const unsigned int            n_columns = 0)
      : shape_values(shape_values.begin())
      , shape_gradients(shape_gradients.begin())
      , shape_hessians(shape_hessians.begin())
      , n_rows(n_rows)
      , n_columns(n_columns)
    {
      if (variant == evaluate_evenodd)
        {
          if (!shape_values.empty())
            AssertDimension(shape_values.size(),
                            n_rows * ((n_columns + 1) / 2));
          if (!shape_gradients.empty())
            AssertDimension(shape_gradients.size(),
                            n_rows * ((n_columns + 1) / 2));
          if (!shape_hessians.empty())
            AssertDimension(shape_hessians.size(),
                            n_rows * ((n_columns + 1) / 2));
        }
      else
        {
          Assert(shape_values.empty() ||
                   shape_values.size() == n_rows * n_columns,
                 ExcDimensionMismatch(shape_values.size(), n_rows * n_columns));
          Assert(shape_gradients.empty() ||
                   shape_gradients.size() == n_rows * n_columns,
                 ExcDimensionMismatch(shape_gradients.size(),
                                      n_rows * n_columns));
          Assert(shape_hessians.empty() ||
                   shape_hessians.size() == n_rows * n_columns,
                 ExcDimensionMismatch(shape_hessians.size(),
                                      n_rows * n_columns));
        }
    }

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct(const Number2     *shape_values,
                           const Number2     *shape_gradients,
                           const Number2     *shape_hessians,
                           const unsigned int n_rows    = 0,
                           const unsigned int n_columns = 0)
      : shape_values(shape_values)
      , shape_gradients(shape_gradients)
      , shape_hessians(shape_hessians)
      , n_rows(n_rows)
      , n_columns(n_columns)
    {}

    template <int direction, bool contract_over_rows, bool add, int stride = 1>
    void
    values(const Number *in, Number *out) const
    {
      constexpr EvaluatorQuantity value_type = EvaluatorQuantity::value;
      apply<direction, contract_over_rows, add, false, value_type, stride>(
        shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add, int stride = 1>
    void
    gradients(const Number *in, Number *out) const
    {
      constexpr EvaluatorQuantity gradient_type =
        (variant != evaluate_evenodd ? EvaluatorQuantity::value :
                                       EvaluatorQuantity::gradient);
      apply<direction, contract_over_rows, add, false, gradient_type, stride>(
        shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians(const Number *in, Number *out) const
    {
      constexpr EvaluatorQuantity hessian_type =
        (variant != evaluate_evenodd ? EvaluatorQuantity::value :
                                       EvaluatorQuantity::hessian);
      apply<direction, contract_over_rows, add, false, hessian_type>(
        shape_hessians, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_values != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, true, EvaluatorQuantity::value>(
        shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_gradients != nullptr, ExcNotInitialized());
      constexpr EvaluatorQuantity gradient_type =
        (variant != evaluate_evenodd ? EvaluatorQuantity::value :
                                       EvaluatorQuantity::gradient);
      apply<direction, contract_over_rows, add, true, gradient_type>(
        shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_hessians != nullptr, ExcNotInitialized());
      constexpr EvaluatorQuantity hessian_type =
        (variant != evaluate_evenodd ? EvaluatorQuantity::value :
                                       EvaluatorQuantity::hessian);
      apply<direction, contract_over_rows, add, true, hessian_type>(
        shape_hessians, in, out);
    }

    template <int               direction,
              bool              contract_over_rows,
              bool              add,
              bool              one_line = false,
              EvaluatorQuantity quantity = EvaluatorQuantity::value,
              int               stride   = 1>
    void
    apply(const Number2 *DEAL_II_RESTRICT shape_data,
          const Number                   *in,
          Number                         *out) const;

    const Number2     *shape_values;
    const Number2     *shape_gradients;
    const Number2     *shape_hessians;
    const unsigned int n_rows;
    const unsigned int n_columns;
  };



  template <EvaluatorVariant variant,
            int              dim,
            typename Number,
            typename Number2>
  template <int               direction,
            bool              contract_over_rows,
            bool              add,
            bool              one_line,
            EvaluatorQuantity quantity,
            int               stride>
  inline void
  EvaluatorTensorProduct<variant, dim, 0, 0, Number, Number2>::apply(
    const Number2 *DEAL_II_RESTRICT shape_data,
    const Number                   *in,
    Number                         *out) const
  {
    static_assert(one_line == false || direction == dim - 1,
                  "Single-line evaluation only works for direction=dim-1.");
    Assert(shape_data != nullptr,
           ExcMessage(
             "The given array shape_data must not be the null pointer!"));
    Assert(dim == direction + 1 || one_line == true || n_rows == n_columns ||
             in != out,
           ExcMessage("In-place operation only supported for "
                      "n_rows==n_columns or single-line interpolation"));
    AssertIndexRange(direction, dim);
    const int mm = contract_over_rows ? n_rows : n_columns,
              nn = contract_over_rows ? n_columns : n_rows;

    const int stride_operation =
      direction == 0 ? 1 : Utilities::fixed_power<direction>(n_columns);
    const int n_blocks1 = one_line ? 1 : stride_operation;
    const int n_blocks2 = direction >= dim - 1 ?
                            1 :
                            Utilities::fixed_power<dim - direction - 1>(n_rows);
    Assert(n_rows <= 128, ExcNotImplemented());

    constexpr int stride_in  = !contract_over_rows ? stride : 1;
    constexpr int stride_out = contract_over_rows ? stride : 1;
    for (int i2 = 0; i2 < n_blocks2; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks1; ++i1)
          {
            // the empty template case can only run the general evaluator or
            // evenodd
            constexpr EvaluatorVariant restricted_variant =
              variant == evaluate_evenodd ? evaluate_evenodd : evaluate_general;
            apply_matrix_vector_product<restricted_variant,
                                        quantity,
                                        contract_over_rows,
                                        add,
                                        (direction != 0 || stride != 1)>(
              shape_data,
              in,
              out,
              n_rows,
              n_columns,
              stride_operation * stride_in,
              stride_operation * stride_out);

            if (one_line == false)
              {
                in += stride_in;
                out += stride_out;
              }
          }
        if (one_line == false)
          {
            in += stride_operation * (mm - 1) * stride_in;
            out += stride_operation * (nn - 1) * stride_out;
          }
      }
  }



  template <int  dim,
            int  fe_degree,
            int  n_q_points_1d,
            bool contract_over_rows,
            bool symmetric_evaluate = true>
  struct EvaluatorTensorProductAnisotropic
  {
    template <int direction,
              int stride       = 1,
              typename Number  = double,
              typename Number2 = double>
    static void
    normal(const MatrixFreeFunctions::UnivariateShapeData<Number2> &data,
           const Number                                            *in,
           Number                                                  *out,
           const bool add_into_result  = false,
           const int  subface_index_1d = 0)
    {
      AssertIndexRange(direction, dim);
      AssertDimension(fe_degree, data.fe_degree);
      AssertDimension(n_q_points_1d, data.n_q_points_1d);
      constexpr int  n_rows    = fe_degree + 1;
      constexpr int  n_columns = n_q_points_1d;
      constexpr int  mm        = contract_over_rows ? n_rows : n_columns;
      constexpr int  nn        = contract_over_rows ? n_columns : n_rows;
      const Number2 *shape_data =
        symmetric_evaluate ?
          data.shape_values_eo.data() :
          data.values_within_subface[subface_index_1d].data();
      Assert(shape_data != nullptr, ExcNotInitialized());
      Assert(contract_over_rows == false || !add_into_result,
             ExcMessage("Cannot add into result if contract_over_rows = true"));

      constexpr int n_blocks1  = Utilities::pow(fe_degree, direction);
      constexpr int n_blocks2  = Utilities::pow(fe_degree, dim - direction - 1);
      constexpr int stride_in  = contract_over_rows ? 1 : stride;
      constexpr int stride_out = contract_over_rows ? stride : 1;
      constexpr EvaluatorVariant variant =
        symmetric_evaluate ? evaluate_evenodd : evaluate_general;

      for (int i2 = 0; i2 < n_blocks2; ++i2)
        {
          for (int i1 = 0; i1 < n_blocks1; ++i1)
            {
              if (contract_over_rows == false && add_into_result)
                apply_matrix_vector_product<variant,
                                            EvaluatorQuantity::value,
                                            n_rows,
                                            n_columns,
                                            n_blocks1 * stride_in,
                                            n_blocks1 * stride_out,
                                            contract_over_rows,
                                            true>(shape_data, in, out);
              else
                apply_matrix_vector_product<variant,
                                            EvaluatorQuantity::value,
                                            n_rows,
                                            n_columns,
                                            n_blocks1 * stride_in,
                                            n_blocks1 * stride_out,
                                            contract_over_rows,
                                            false>(shape_data, in, out);

              in += stride_in;
              out += stride_out;
            }
          in += n_blocks1 * (mm - 1) * stride_in;
          out += n_blocks1 * (nn - 1) * stride_out;
        }
    }

    template <int direction,
              int normal_direction,
              int stride       = 1,
              typename Number  = double,
              typename Number2 = double>
    static void
    tangential(const MatrixFreeFunctions::UnivariateShapeData<Number2> &data,
               const Number                                            *in,
               Number                                                  *out,
               const int subface_index_1d = 0)
    {
      AssertIndexRange(direction, dim);
      AssertDimension(fe_degree - 1, data.fe_degree);
      AssertDimension(n_q_points_1d, data.n_q_points_1d);
      static_assert(direction != normal_direction,
                    "Cannot interpolate tangentially in normal direction");

      constexpr int  n_rows    = std::max(fe_degree, 0);
      constexpr int  n_columns = n_q_points_1d;
      const Number2 *shape_data =
        symmetric_evaluate ?
          data.shape_values_eo.data() :
          data.values_within_subface[subface_index_1d].data();
      Assert(shape_data != nullptr, ExcNotInitialized());

      constexpr int n_blocks1 =
        (direction > normal_direction) ?
          Utilities::pow(n_q_points_1d, direction) :
          (direction > 0 ?
             (Utilities::pow(fe_degree, direction - 1) * n_q_points_1d) :
             1);
      constexpr int n_blocks2 =
        (direction > normal_direction) ?
          Utilities::pow(fe_degree, dim - 1 - direction) :
          ((direction + 1 < dim) ?
             (Utilities::pow(fe_degree, dim - 2 - direction) * n_q_points_1d) :
             1);

      constexpr EvaluatorVariant variant =
        symmetric_evaluate ? evaluate_evenodd : evaluate_general;

      // Since we may perform an in-place interpolation, we must run the step
      // expanding the size of the basis backward ('contract_over_rows' aka
      // 'evaluate' case), so shift the pointers and decrement during the loop
      if (contract_over_rows)
        {
          in += (n_blocks2 - 1) * n_blocks1 * n_rows + n_blocks1 - 1;
          out +=
            stride * ((n_blocks2 - 1) * n_blocks1 * n_columns + n_blocks1 - 1);
          for (int i2 = 0; i2 < n_blocks2; ++i2)
            {
              for (int i1 = 0; i1 < n_blocks1; ++i1)
                {
                  apply_matrix_vector_product<variant,
                                              EvaluatorQuantity::value,
                                              n_rows,
                                              n_columns,
                                              n_blocks1,
                                              n_blocks1 * stride,
                                              true,
                                              false>(shape_data, in, out);

                  --in;
                  out -= stride;
                }
              in -= n_blocks1 * (n_rows - 1);
              out -= n_blocks1 * (n_columns - 1) * stride;
            }
        }
      else
        {
          for (int i2 = 0; i2 < n_blocks2; ++i2)
            {
              for (int i1 = 0; i1 < n_blocks1; ++i1)
                {
                  apply_matrix_vector_product<variant,
                                              EvaluatorQuantity::value,
                                              n_rows,
                                              n_columns,
                                              n_blocks1 * stride,
                                              n_blocks1,
                                              false,
                                              false>(shape_data, in, out);

                  in += stride;
                  ++out;
                }
              in += n_blocks1 * (n_columns - 1) * stride;
              out += n_blocks1 * (n_rows - 1);
            }
        }
    }
  };



  /**
   * This function applies the tensor product operation to produce face values
   * from cell values. The algorithm involved here can be interpreted as the
   * first sweep in sum factorization, reducing the dimensionality of the data
   * set from dim-dimensional cell values to (dim-1)-dimensional face
   * values. This step is always done before we evaluate within the face, as
   * it reduces the length of the loops for the successive steps.
   *
   * @tparam n_rows_template The number of entries within the interpolation,
   *             typically equal to the polynomial degree plus one, if known
   *             at compile time, otherwise n_rows_runtime is used.
   * @tparam stride_template The stride between successive entries in the
   *             one-dimensional operation of sum factorization, if known at
   *             compile time, otherwise stride_runtime is used.
   * @tparam contract_onto_face If true, the input vector is of size n_rows^dim
   *                            and interpolation into n_rows^(dim-1) points
   *                            is performed. This is a typical scenario in
   *                            FEFaceEvaluation::evaluate() calls. If false,
   *                            data from n_rows^(dim-1) points is expanded
   *                            into the n_rows^dim points of the
   *                            higher-dimensional data array. Derivatives in
   * the case contract_onto_face==false are summed together.
   * @tparam add If true, the result is added to the output vector, else
   *             the computed values overwrite the content in the output.
   * @tparam max_derivative Sets the number of derivatives that should be
   *             computed. 0 means only values, 1 means values and first
   *             derivatives, 2 up to second derivatives. Note that all the
   *             derivatives access the data in @p shape_values passed to
   *             the constructor of the class.
   *
   * @param shape_values Address of the interpolation matrix.
   * @param n_blocks Number of interpolation layers used along the up to two
   *             dimensions tangential to the interpolation direction.
   * @param steps Increments in the input array from one step to the next,
   *             varied in conjunction with the @p stride_template variable: We
   *             increment by @p stride_template along the 1d interpolation,
   *             and then increment by @p steps when passing from one line
   *             to the next.
   * @param input Address of the input data vector.
   * @param output Address of the output data vector.
   * @param n_rows_runtime Alternative number of rows to be used if the
   *             variable @p n_rows_template is 0, enabling a run-time path.
   * @param stride_runtime Alternative number for the stride to be used if the
   *             variable @p n_rows_template is 0.
   */
  template <int  n_rows_template,
            int  stride_template,
            bool contract_onto_face,
            bool add,
            int  max_derivative,
            typename Number,
            typename Number2>
  inline std::enable_if_t<contract_onto_face, void>
  interpolate_to_face(const Number2            *shape_values,
                      const std::array<int, 2> &n_blocks,
                      const std::array<int, 2> &steps,
                      const Number             *input,
                      Number *DEAL_II_RESTRICT  output,
                      const int                 n_rows_runtime = 0,
                      const int                 stride_runtime = 1)
  {
    const int n_rows = n_rows_template > 0 ? n_rows_template : n_rows_runtime;
    const int stride = n_rows_template > 0 ? stride_template : stride_runtime;

    Number *output1 = output + n_blocks[0] * n_blocks[1];
    Number *output2 = output1 + n_blocks[0] * n_blocks[1];
    for (int i2 = 0; i2 < n_blocks[1]; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks[0]; ++i1)
          {
            Number res0 = shape_values[0] * input[0];
            Number res1, res2;
            if (max_derivative > 0)
              res1 = shape_values[n_rows] * input[0];
            if (max_derivative > 1)
              res2 = shape_values[2 * n_rows] * input[0];
            for (int ind = 1; ind < n_rows; ++ind)
              {
                res0 += shape_values[ind] * input[stride * ind];
                if (max_derivative > 0)
                  res1 += shape_values[ind + n_rows] * input[stride * ind];
                if (max_derivative > 1)
                  res2 += shape_values[ind + 2 * n_rows] * input[stride * ind];
              }
            if (add)
              {
                output[i1] += res0;
                if (max_derivative > 0)
                  output1[i1] += res1;
                if (max_derivative > 1)
                  output2[i2] += res2;
              }
            else
              {
                output[i1] = res0;
                if (max_derivative > 0)
                  output1[i1] = res1;
                if (max_derivative > 1)
                  output2[i1] = res2;
              }
            input += steps[0];
          }
        output += n_blocks[0];
        if (max_derivative > 0)
          output1 += n_blocks[0];
        if (max_derivative > 1)
          output2 += n_blocks[0];
        input += steps[1];
      }
  }



  /**
   * Helper function to specify whether a transformation to collocation should
   * be used: It should give correct results (first condition), we need to be
   * able to initialize the fields in shape_info.templates.h from the
   * polynomials (second condition), and it should be the most efficient
   * choice in terms of operation counts (third condition).
   */
  constexpr bool
  use_collocation_evaluation(const unsigned int fe_degree,
                             const unsigned int n_q_points_1d)
  {
    return (n_q_points_1d > fe_degree) && (n_q_points_1d < 200) &&
           (n_q_points_1d <= 3 * fe_degree / 2 + 1);
  }



  /**
   * This function performs the opposite operation to the interpolate_to_face
   * function, done as the last step in sum factorization to embed face values
   * and gradients back to values on all degrees of freedom of the cell.
   */
  template <int  n_rows_template,
            int  stride_template,
            bool contract_onto_face,
            bool add,
            int  max_derivative,
            typename Number,
            typename Number2>
  inline std::enable_if_t<!contract_onto_face, void>
  interpolate_to_face(const Number2            *shape_values,
                      const std::array<int, 2> &n_blocks,
                      const std::array<int, 2> &steps,
                      const Number             *input,
                      Number *DEAL_II_RESTRICT  output,
                      const int                 n_rows_runtime = 0,
                      const int                 stride_runtime = 1)
  {
    const int n_rows = n_rows_template > 0 ? n_rows_template : n_rows_runtime;
    const int stride = n_rows_template > 0 ? stride_template : stride_runtime;

    const Number *input1 = input + n_blocks[0] * n_blocks[1];
    const Number *input2 = input1 + n_blocks[0] * n_blocks[1];
    for (int i2 = 0; i2 < n_blocks[1]; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks[0]; ++i1)
          {
            const Number in = input[i1];
            Number       in1, in2;
            if (max_derivative > 0)
              in1 = input1[i1];
            if (max_derivative > 1)
              in2 = input2[i1];
            for (int col = 0; col < n_rows; ++col)
              {
                Number result =
                  add ? (output[col * stride] + shape_values[col] * in) :
                        (shape_values[col] * in);
                if (max_derivative > 0)
                  result += shape_values[col + n_rows] * in1;
                if (max_derivative > 1)
                  result += shape_values[col + 2 * n_rows] * in2;

                output[col * stride] = result;
              }
            output += steps[0];
          }
        input += n_blocks[0];
        if (max_derivative > 0)
          input1 += n_blocks[0];
        if (max_derivative > 1)
          input2 += n_blocks[0];
        output += steps[1];
      }
  }

  template <int dim, int n_points_1d_template, typename Number>
  inline void
  weight_fe_q_dofs_by_entity(const Number      *weights,
                             const unsigned int n_components,
                             const int          n_points_1d_non_template,
                             Number            *data)
  {
    const int n_points_1d = n_points_1d_template != -1 ?
                              n_points_1d_template :
                              n_points_1d_non_template;

    Assert(n_points_1d > 0, ExcNotImplemented());
    Assert(n_points_1d < 100, ExcNotImplemented());

    unsigned int compressed_index[100];
    compressed_index[0] = 0;
    for (int i = 1; i < n_points_1d - 1; ++i)
      compressed_index[i] = 1;
    compressed_index[n_points_1d - 1] = 2;

    for (unsigned int c = 0; c < n_components; ++c)
      for (int k = 0; k < (dim > 2 ? n_points_1d : 1); ++k)
        for (int j = 0; j < (dim > 1 ? n_points_1d : 1); ++j)
          {
            const unsigned int shift =
              9 * compressed_index[k] + 3 * compressed_index[j];
            data[0] *= weights[shift];
            // loop bound as int avoids compiler warnings in case n_points_1d
            // == 1 (polynomial degree 0)
            const Number weight = weights[shift + 1];
            for (int i = 1; i < n_points_1d - 1; ++i)
              data[i] *= weight;
            data[n_points_1d - 1] *= weights[shift + 2];
            data += n_points_1d;
          }
  }


  template <int dim, int n_points_1d_template, typename Number>
  inline void
  weight_fe_q_dofs_by_entity_shifted(const Number      *weights,
                                     const unsigned int n_components,
                                     const int n_points_1d_non_template,
                                     Number   *data)
  {
    const int n_points_1d = n_points_1d_template != -1 ?
                              n_points_1d_template :
                              n_points_1d_non_template;

    Assert((n_points_1d % 2) == 1,
           ExcMessage("The function can only with add number of points"));
    Assert(n_points_1d > 0, ExcNotImplemented());
    Assert(n_points_1d < 100, ExcNotImplemented());

    const unsigned int n_inside_1d = n_points_1d / 2;

    unsigned int compressed_index[100];

    unsigned int c = 0;
    for (int i = 0; i < n_inside_1d; ++i)
      compressed_index[c++] = 0;
    compressed_index[c++] = 1;
    for (int i = 0; i < n_inside_1d; ++i)
      compressed_index[c++] = 2;

    for (unsigned int c = 0; c < n_components; ++c)
      for (int k = 0; k < (dim > 2 ? n_points_1d : 1); ++k)
        for (int j = 0; j < (dim > 1 ? n_points_1d : 1); ++j)
          {
            const unsigned int shift =
              9 * compressed_index[k] + 3 * compressed_index[j];

            unsigned int c       = 0;
            const Number weight1 = weights[shift];
            for (int i = 0; i < n_inside_1d; ++i)
              data[c++] *= weight1;
            data[c++] *= weights[shift + 1];
            const Number weight2 = weights[shift + 2];
            for (int i = 0; i < n_inside_1d; ++i)
              data[c++] *= weight2;
            data += n_points_1d;
          }
  }


  template <int dim, int n_points_1d_template, typename Number>
  inline bool
  compute_weights_fe_q_dofs_by_entity(const Number      *data,
                                      const unsigned int n_components,
                                      const int n_points_1d_non_template,
                                      Number   *weights)
  {
    const int n_points_1d = n_points_1d_template != -1 ?
                              n_points_1d_template :
                              n_points_1d_non_template;

    Assert(n_points_1d > 0, ExcNotImplemented());
    Assert(n_points_1d < 100, ExcNotImplemented());

    unsigned int compressed_index[100];
    compressed_index[0] = 0;
    for (int i = 1; i < n_points_1d - 1; ++i)
      compressed_index[i] = 1;
    compressed_index[n_points_1d - 1] = 2;

    // Insert the number data into a storage position for weight,
    // ensuring that the array has either not been touched before
    // or the previous content is the same. In case the previous
    // content has a different value, we exit this function and
    // signal to outer functions that the compression was not possible.
    const auto check_and_set = [](Number &weight, const Number &data) {
      if (weight == Number(-1.0) || weight == data)
        {
          weight = data;
          return true; // success for the entry
        }

      return false; // failure for the entry
    };

    for (unsigned int c = 0; c < Utilities::pow<unsigned int>(3, dim); ++c)
      weights[c] = Number(-1.0);

    for (unsigned int c = 0; c < n_components; ++c)
      for (int k = 0; k < (dim > 2 ? n_points_1d : 1); ++k)
        for (int j = 0; j < (dim > 1 ? n_points_1d : 1);
             ++j, data += n_points_1d)
          {
            const unsigned int shift =
              9 * compressed_index[k] + 3 * compressed_index[j];

            if (!check_and_set(weights[shift], data[0]))
              return false; // failure

            for (int i = 1; i < n_points_1d - 1; ++i)
              if (!check_and_set(weights[shift + 1], data[i]))
                return false; // failure

            if (!check_and_set(weights[shift + 2], data[n_points_1d - 1]))
              return false; // failure
          }

    return true; // success
  }


  template <int dim, int n_points_1d_template, typename Number>
  inline bool
  compute_weights_fe_q_dofs_by_entity_shifted(
    const Number      *data,
    const unsigned int n_components,
    const int          n_points_1d_non_template,
    Number            *weights)
  {
    const int n_points_1d = n_points_1d_template != -1 ?
                              n_points_1d_template :
                              n_points_1d_non_template;

    Assert((n_points_1d % 2) == 1,
           ExcMessage("The function can only with add number of points"));
    Assert(n_points_1d > 0, ExcNotImplemented());
    Assert(n_points_1d < 100, ExcNotImplemented());

    const unsigned int n_inside_1d = n_points_1d / 2;

    unsigned int compressed_index[100];

    unsigned int c = 0;
    for (int i = 0; i < n_inside_1d; ++i)
      compressed_index[c++] = 0;
    compressed_index[c++] = 1;
    for (int i = 0; i < n_inside_1d; ++i)
      compressed_index[c++] = 2;

    // Insert the number data into a storage position for weight,
    // ensuring that the array has either not been touched before
    // or the previous content is the same. In case the previous
    // content has a different value, we exit this function and
    // signal to outer functions that the compression was not possible.
    const auto check_and_set = [](Number &weight, const Number &data) {
      if (weight == Number(-1.0) || weight == data)
        {
          weight = data;
          return true; // success for the entry
        }

      return false; // failure for the entry
    };

    for (unsigned int c = 0; c < Utilities::pow<unsigned int>(3, dim); ++c)
      weights[c] = Number(-1.0);

    for (unsigned int comp = 0; comp < n_components; ++comp)
      for (int k = 0; k < (dim > 2 ? n_points_1d : 1); ++k)
        for (int j = 0; j < (dim > 1 ? n_points_1d : 1);
             ++j, data += n_points_1d)
          {
            const unsigned int shift =
              9 * compressed_index[k] + 3 * compressed_index[j];

            unsigned int c = 0;

            for (int i = 0; i < n_inside_1d; ++i)
              if (!check_and_set(weights[shift], data[c++]))
                return false; // failure

            if (!check_and_set(weights[shift + 1], data[c++]))
              return false; // failure

            for (int i = 0; i < n_inside_1d; ++i)
              if (!check_and_set(weights[shift + 2], data[c++]))
                return false; // failure
          }

    return true; // success
  }


} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
