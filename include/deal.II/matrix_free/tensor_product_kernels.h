// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2023 by the deal.II authors
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


#ifndef dealii_matrix_free_tensor_product_kernels_h
#define dealii_matrix_free_tensor_product_kernels_h

#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/ndarray.h>
#include <deal.II/base/polynomial.h>
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



  /**
   * Specialization of the matrix-vector kernel for run-time loop bounds in
   * the generic evaluator.
   */
  template <EvaluatorVariant  variant,
            EvaluatorQuantity quantity,
            bool              transpose_matrix,
            bool              add,
            typename Number,
            typename Number2>
  std::enable_if_t<(variant == evaluate_general), void>
  apply_matrix_vector_product(const Number2 *matrix,
                              const Number  *in,
                              Number        *out,
                              const int      n_rows,
                              const int      n_columns,
                              const int      stride_in,
                              const int      stride_out)
  {
    const int mm = transpose_matrix ? n_rows : n_columns,
              nn = transpose_matrix ? n_columns : n_rows;
    Assert(n_rows <= 128, ExcNotImplemented());
    Assert(n_rows > 0 && n_columns > 0,
           ExcInternalError("Empty evaluation task!"));
    Assert(n_rows > 0 && n_columns > 0,
           ExcInternalError("The evaluation needs n_rows, n_columns > 0, but " +
                            std::to_string(n_rows) + ", " +
                            std::to_string(n_columns) + " was passed!"));

    static_assert(quantity == EvaluatorQuantity::value,
                  "This function should only use EvaluatorQuantity::value");

    // specialization for n_rows = 2 that manually unrolls the innermost loop
    // to make the operation perform better (not completely as good as the
    // templated one, but much better than the generic version down below,
    // because the loop over col can be more effectively unrolled by the
    // compiler)
    if (transpose_matrix && n_rows == 2)
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
    else if (transpose_matrix && n_rows == 3)
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
    else
      {
        std::array<Number, 129> x;
        for (int i = 0; i < mm; ++i)
          x[i] = in[stride_in * i];

        Number res0;
        for (int col = 0; col < nn; ++col)
          {
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
                    ;
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
      n_rows_static == 0 ? stride_in_runtime : stride_in_static;
    const int stride_out =
      n_rows_static == 0 ? stride_out_runtime : stride_out_static;

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
                                0,
                                0,
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
                                        add>(shape_data,
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
   *             derivatives, 2 up to second derivates. Note that all the
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



  /**
   * Struct to avoid using Tensor<1, dim, Point<dim2>> in
   * evaluate_tensor_product_value_and_gradient because a Point cannot be used
   * within Tensor. Instead, a specialization of this struct upcasts the point
   * to a Tensor<1,dim>.
   */
  template <typename Number, typename Number2>
  struct ProductTypeNoPoint
  {
    using type = typename ProductType<Number, Number2>::type;
  };

  template <int dim, typename Number, typename Number2>
  struct ProductTypeNoPoint<Point<dim, Number>, Number2>
  {
    using type = typename ProductType<Tensor<1, dim, Number>, Number2>::type;
  };



  /**
   * Computes the values and derivatives of the 1d polynomials @p poly at the
   * specified point @p p and stores it in @p shapes.
   */
  template <int dim, typename Number>
  inline void
  compute_values_of_array(
    dealii::ndarray<Number, 2, dim>                    *shapes,
    const std::vector<Polynomials::Polynomial<double>> &poly,
    const Point<dim, Number>                           &p,
    const unsigned int                                  derivative = 1)
  {
    const int n_shapes = poly.size();

    // Evaluate 1d polynomials and their derivatives
    std::array<Number, dim> point;
    for (unsigned int d = 0; d < dim; ++d)
      point[d] = p[d];
    for (int i = 0; i < n_shapes; ++i)
      poly[i].values_of_array(point, derivative, shapes[i].data());
  }



  /**
   * Specialization of above function for dim = 0. Should not be called.
   */
  template <typename Number>
  inline void
  compute_values_of_array(dealii::ndarray<Number, 2, 0> *,
                          const std::vector<Polynomials::Polynomial<double>> &,
                          const Point<0, Number> &,
                          const unsigned int)
  {
    Assert(false, ExcInternalError());
  }



  /**
   * Interpolate inner dimensions of tensor product shape functions.
   */
  template <int dim,
            int length,
            typename Number2,
            typename Number,
            int  n_values    = 1,
            bool do_renumber = true>
  inline
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
      std::array<typename ProductTypeNoPoint<Number, Number2>::type,
                 2 + n_values>
      do_interpolate_xy(const Number                           *values,
                        const std::vector<unsigned int>        &renumber,
                        const dealii::ndarray<Number2, 2, dim> *shapes,
                        const int n_shapes_runtime,
                        int      &i)
  {
    static_assert(0 <= dim && dim <= 3, "Only dim=0,1,2,3 implemented");
    static_assert(1 <= n_values && n_values <= 2,
                  "Only n_values=1,2 implemented");

    const int n_shapes = length > 0 ? length : n_shapes_runtime;

    // If n_values > 1, we want to interpolate from a second array,
    // placed in the same array immediately after the main data. This
    // is used to interpolate normal derivatives onto faces.
    const Number *values_2 =
      n_values > 1 ?
        values + (length > 0 ? Utilities::pow(length, dim) :
                               Utilities::fixed_power<dim>(n_shapes_runtime)) :
        nullptr;
    using Number3 = typename ProductTypeNoPoint<Number, Number2>::type;
    std::array<Number3, 2 + n_values> result = {};
    for (int i1 = 0; i1 < (dim > 1 ? n_shapes : 1); ++i1)
      {
        // Interpolation + derivative x direction
        std::array<Number3, 1 + n_values> inner_result = {};

        // Distinguish the inner loop based on whether we have a
        // renumbering or not
        if (do_renumber && !renumber.empty())
          for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
            {
              // gradient
              inner_result[0] += shapes[i0][1][0] * values[renumber[i]];
              // values
              inner_result[1] += shapes[i0][0][0] * values[renumber[i]];
              if (n_values > 1)
                inner_result[2] += shapes[i0][0][0] * values_2[renumber[i]];
            }
        else
          for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
            {
              // gradient
              inner_result[0] += shapes[i0][1][0] * values[i];
              // values
              inner_result[1] += shapes[i0][0][0] * values[i];
              if (n_values > 1)
                inner_result[2] += shapes[i0][0][0] * values_2[i];
            }

        if (dim > 1)
          {
            // Interpolation + derivative in y direction
            // gradient
            result[0] += inner_result[0] * shapes[i1][0][1];
            result[1] += inner_result[1] * shapes[i1][1][1];
            // values
            result[2] += inner_result[1] * shapes[i1][0][1];
            if (n_values > 1)
              result[3] += inner_result[2] * shapes[i1][0][1];
          }
        else
          {
            // gradient
            result[0] = inner_result[0];
            // values
            result[1] = inner_result[1];
            if (n_values > 1)
              result[2] = inner_result[2];
          }
      }
    return result;
  }



  /**
   * Interpolates the values and gradients into the points specified in
   * @p compute_values_of_array() with help of the precomputed @p shapes.
   */
  template <int dim,
            typename Number,
            typename Number2,
            int  n_values    = 1,
            bool do_renumber = true>
  inline std::array<typename ProductTypeNoPoint<Number, Number2>::type,
                    dim + n_values>
  evaluate_tensor_product_value_and_gradient_shapes(
    const dealii::ndarray<Number2, 2, dim> *shapes,
    const int                               n_shapes,
    const Number                           *values,
    const std::vector<unsigned int>        &renumber = {})
  {
    static_assert(0 <= dim && dim <= 3, "Only dim=0,1,2,3 implemented");
    static_assert(1 <= n_values && n_values <= 2,
                  "Only n_values=1,2 implemented");

    using Number3 = typename ProductTypeNoPoint<Number, Number2>::type;

    std::array<Number3, dim + n_values> result = {};
    if (dim == 0)
      {
        // We only need the interpolation of the value and normal derivatives on
        // faces of a 1d element. As the interpolation is the value at the
        // point, simply set the result vector accordingly.
        result[0] = values[0];
        if (n_values > 1)
          result[1] = values[1];
        return result;
      }

    // Go through the tensor product of shape functions and interpolate
    // with optimal algorithm
    for (int i2 = 0, i = 0; i2 < (dim > 2 ? n_shapes : 1); ++i2)
      {
        std::array<Number3, 2 + n_values> inner_result;
        // Generate separate code with known loop bounds for the most common
        // cases
        if (n_shapes == 2)
          inner_result =
            do_interpolate_xy<dim, 2, Number2, Number, n_values, do_renumber>(
              values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 3)
          inner_result =
            do_interpolate_xy<dim, 3, Number2, Number, n_values, do_renumber>(
              values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 4)
          inner_result =
            do_interpolate_xy<dim, 4, Number2, Number, n_values, do_renumber>(
              values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 5)
          inner_result =
            do_interpolate_xy<dim, 5, Number2, Number, n_values, do_renumber>(
              values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 6)
          inner_result =
            do_interpolate_xy<dim, 6, Number2, Number, n_values, do_renumber>(
              values, renumber, shapes, n_shapes, i);
        else
          inner_result =
            do_interpolate_xy<dim, -1, Number2, Number, n_values, do_renumber>(
              values, renumber, shapes, n_shapes, i);
        if (dim == 3)
          {
            // derivative + interpolation in z direction
            // gradient
            result[0] += inner_result[0] * shapes[i2][0][2];
            result[1] += inner_result[1] * shapes[i2][0][2];
            result[2] += inner_result[2] * shapes[i2][1][2];
            // values
            result[3] += inner_result[2] * shapes[i2][0][2];
            if (n_values > 1)
              result[4] += inner_result[3] * shapes[i2][0][2];
          }
        else if (dim == 2)
          {
            // gradient
            result[0] = inner_result[0];
            result[1] = inner_result[1];
            // values
            result[2] = inner_result[2];
            if (n_values > 1)
              result[3] = inner_result[3];
          }
        else
          {
            // gradient
            result[0] = inner_result[0];
            // values
            result[1] = inner_result[1];
            if (n_values > 1)
              result[2] = inner_result[2];
          }
      }

    return result;
  }



  /**
   * Specializes @p evaluate_tensor_product_value_and_gradient() for linear
   * polynomials which massively reduces the necessary instructions.
   */
  template <int dim,
            typename Number,
            typename Number2,
            int n_values = 1,
            int stride   = 1>
  inline std::array<typename ProductTypeNoPoint<Number, Number2>::type,
                    dim + n_values>
  evaluate_tensor_product_value_and_gradient_linear(
    const Number              *values,
    const Point<dim, Number2> &p)
  {
    static_assert(0 <= dim && dim <= 3, "Only dim=0,1,2,3 implemented");
    static_assert(1 <= n_values && n_values <= 2,
                  "Only n_values=1,2 implemented");

    using Number3 = typename ProductTypeNoPoint<Number, Number2>::type;

    static_assert(
      n_values == 1 || stride == 1,
      "Either n_values or stride has to be one for correct data access!");
    // If n_values > 1, we want to interpolate from a second array,
    // placed in the same array immediately after the main data. This
    // is used to interpolate normal derivatives onto faces.

    std::array<Number3, dim + n_values> result;
    if (dim == 0)
      {
        // we only need the value on faces of a 1d element
        result[0] = values[0];
        if (n_values > 1)
          result[1] = values[1];
      }
    else if (dim == 1)
      {
        // gradient
        result[0] = Number3(values[stride] - values[0]);
        // values
        result[1] = Number3(values[0]) + p[0] * result[0];
        if (n_values > 1)
          result[2] = Number3(values[2]) + p[0] * (values[3] - values[2]);
      }
    else if (dim == 2)
      {
        const Number3 val10 = Number3(values[stride] - values[0]);
        const Number3 val32 = Number3(values[3 * stride] - values[2 * stride]);
        const Number3 tmp0  = Number3(values[0]) + p[0] * val10;
        const Number3 tmp1  = Number3(values[2 * stride]) + p[0] * val32;

        // gradient
        result[0] = val10 + p[1] * (val32 - val10);
        result[1] = tmp1 - tmp0;

        // values
        result[2] = tmp0 + p[1] * result[1];

        if (n_values > 1)
          {
            const Number3 tmp0_2 =
              Number3(values[4]) + p[0] * (values[5] - values[4]);
            const Number3 tmp1_2 =
              Number3(values[6]) + p[0] * (values[7] - values[6]);
            result[3] = tmp0_2 + p[1] * (tmp1_2 - tmp0_2);
          }
      }
    else if (dim == 3)
      {
        const Number3 val10 = Number3(values[stride] - values[0]);
        const Number3 val32 = Number3(values[3 * stride] - values[2 * stride]);
        const Number3 tmp0  = Number3(values[0]) + p[0] * val10;
        const Number3 tmp1  = Number3(values[2 * stride]) + p[0] * val32;
        const Number3 tmp10 = tmp1 - tmp0;
        const Number3 tmpy0 = tmp0 + p[1] * tmp10;

        const Number3 val54 = Number3(values[5 * stride] - values[4 * stride]);
        const Number3 val76 = Number3(values[7 * stride] - values[6 * stride]);
        const Number3 tmp2  = Number3(values[4 * stride]) + p[0] * val54;
        const Number3 tmp3  = Number3(values[6 * stride]) + p[0] * val76;
        const Number3 tmp32 = tmp3 - tmp2;
        const Number3 tmpy1 = tmp2 + p[1] * tmp32;

        // gradient
        result[2]           = tmpy1 - tmpy0;
        result[1]           = tmp10 + p[2] * (tmp32 - tmp10);
        const Number3 tmpz0 = val10 + p[1] * (val32 - val10);
        result[0] = tmpz0 + p[2] * (val54 + p[1] * (val76 - val54) - tmpz0);

        // value
        result[3] = tmpy0 + p[2] * result[2];
        Assert(n_values == 1, ExcNotImplemented());
      }

    return result;
  }



  /**
   * Compute the polynomial interpolation of a tensor product shape function
   * $\varphi_i$ given a vector of coefficients $u_i$ in the form
   * $u_h(\mathbf{x}) = \sum_{i=1}^{k^d} \varphi_i(\mathbf{x}) u_i$. The shape
   * functions $\varphi_i(\mathbf{x}) =
   * \prod_{d=1}^{\text{dim}}\varphi_{i_d}^\text{1d}(x_d)$ represent a tensor
   * product. The function returns a pair with the value of the interpolation
   * as the first component and the gradient in reference coordinates as the
   * second component. Note that for compound types (e.g. the `values` field
   * begin a Point<spacedim> argument), the components of the gradient are
   * sorted as Tensor<1, dim, Tensor<1, spacedim>> with the derivatives
   * as the first index; this is a consequence of the generic arguments in the
   * function.
   *
   * @param poly The underlying one-dimensional polynomial basis
   * $\{\varphi^{1d}_{i_1}\}$ given as a vector of polynomials.
   *
   * @param values The expansion coefficients $u_i$ of type `Number` in
   * the polynomial interpolation. The coefficients can be simply `double`
   * variables but e.g. also Point<spacedim> in case they define arithmetic
   * operations with the type `Number2`.
   *
   * @param p The position in reference coordinates where the interpolation
   * should be evaluated.
   *
   * @param d_linear Flag to specify whether a d-linear (linear in 1d,
   * bi-linear in 2d, tri-linear in 3d) interpolation should be made, which
   * allows to unroll loops and considerably speed up evaluation.
   *
   * @param renumber Optional parameter to specify a renumbering in the
   * coefficient vector, assuming that `values[renumber[i]]` returns
   * the lexicographic (tensor product) entry of the coefficients. If the
   * vector is entry, the values are assumed to be sorted lexicographically.
   */
  template <int dim, typename Number, typename Number2>
  inline std::pair<
    typename ProductTypeNoPoint<Number, Number2>::type,
    Tensor<1, dim, typename ProductTypeNoPoint<Number, Number2>::type>>
  evaluate_tensor_product_value_and_gradient(
    const std::vector<Polynomials::Polynomial<double>> &poly,
    const std::vector<Number>                          &values,
    const Point<dim, Number2>                          &p,
    const bool                                          d_linear = false,
    const std::vector<unsigned int>                    &renumber = {})
  {
    using Number3 = typename ProductTypeNoPoint<Number, Number2>::type;

    std::array<Number3, dim + 1> result;
    if (d_linear)
      {
        result =
          evaluate_tensor_product_value_and_gradient_linear(values.data(), p);
      }
    else
      {
        AssertIndexRange(poly.size(), 200);
        std::array<dealii::ndarray<Number2, 2, dim>, 200> shapes;
        compute_values_of_array(shapes.data(), poly, p);
        result = evaluate_tensor_product_value_and_gradient_shapes<dim,
                                                                   Number,
                                                                   Number2>(
          shapes.data(), poly.size(), values.data(), renumber);
      }
    return std::make_pair(result[dim],
                          Tensor<1, dim, Number3>(
                            ArrayView<Number3>(result.data(), dim)));
  }



  template <int dim,
            int length,
            typename Number2,
            typename Number,
            bool do_renumber = true>
  inline
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
    typename ProductTypeNoPoint<Number, Number2>::type
    do_interpolate_xy_value(const Number                           *values,
                            const std::vector<unsigned int>        &renumber,
                            const dealii::ndarray<Number2, 2, dim> *shapes,
                            const int n_shapes_runtime,
                            int      &i)
  {
    const int n_shapes = length > 0 ? length : n_shapes_runtime;
    using Number3      = typename ProductTypeNoPoint<Number, Number2>::type;
    Number3 result     = {};
    for (int i1 = 0; i1 < (dim > 1 ? n_shapes : 1); ++i1)
      {
        // Interpolation x direction
        Number3 value = {};

        // Distinguish the inner loop based on whether we have a
        // renumbering or not
        if (do_renumber && !renumber.empty())
          for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
            value += shapes[i0][0][0] * values[renumber[i]];
        else
          for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
            value += shapes[i0][0][0] * values[i];

        if (dim > 1)
          result += value * shapes[i1][0][1];
        else
          result = value;
      }
    return result;
  }



  template <int dim, typename Number, typename Number2, bool do_renumber = true>
  inline typename ProductTypeNoPoint<Number, Number2>::type
  evaluate_tensor_product_value_shapes(
    const dealii::ndarray<Number2, 2, dim> *shapes,
    const int                               n_shapes,
    const Number                           *values,
    const std::vector<unsigned int>        &renumber = {})
  {
    static_assert(dim >= 0 && dim <= 3, "Only dim=0,1,2,3 implemented");

    // we only need the value on faces of a 1d element
    if (dim == 0)
      {
        return values[0];
      }

    using Number3 = typename ProductTypeNoPoint<Number, Number2>::type;

    // Go through the tensor product of shape functions and interpolate
    // with optimal algorithm
    Number3 result = {};
    for (int i2 = 0, i = 0; i2 < (dim > 2 ? n_shapes : 1); ++i2)
      {
        Number3 inner_result;
        // Generate separate code with known loop bounds for the most common
        // cases
        if (n_shapes == 2)
          inner_result =
            do_interpolate_xy_value<dim, 2, Number2, Number, do_renumber>(
              values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 3)
          inner_result =
            do_interpolate_xy_value<dim, 3, Number2, Number, do_renumber>(
              values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 4)
          inner_result =
            do_interpolate_xy_value<dim, 4, Number2, Number, do_renumber>(
              values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 5)
          inner_result =
            do_interpolate_xy_value<dim, 5, Number2, Number, do_renumber>(
              values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 6)
          inner_result =
            do_interpolate_xy_value<dim, 6, Number2, Number, do_renumber>(
              values, renumber, shapes, n_shapes, i);
        else
          inner_result =
            do_interpolate_xy_value<dim, -1, Number2, Number, do_renumber>(
              values, renumber, shapes, n_shapes, i);
        if (dim == 3)
          {
            // Interpolation + derivative in z direction
            result += inner_result * shapes[i2][0][2];
          }
        else
          {
            result = inner_result;
          }
      }

    return result;
  }



  template <int dim, typename Number, typename Number2, int stride = 1>
  inline typename ProductTypeNoPoint<Number, Number2>::type
  evaluate_tensor_product_value_linear(const Number              *values,
                                       const Point<dim, Number2> &p)
  {
    static_assert(dim >= 0 && dim <= 3, "Only dim=0,1,2,3 implemented");

    using Number3 = typename ProductTypeNoPoint<Number, Number2>::type;

    if (dim == 0)
      {
        // we only need the value on faces of a 1d element
        return values[0];
      }
    else if (dim == 1)
      {
        return Number3(values[0]) + p[0] * Number3(values[stride] - values[0]);
      }
    else if (dim == 2)
      {
        const Number3 val10 = Number3(values[stride] - values[0]);
        const Number3 val32 = Number3(values[3 * stride] - values[2 * stride]);
        const Number3 tmp0  = Number3(values[0]) + p[0] * val10;
        const Number3 tmp1  = Number3(values[2 * stride]) + p[0] * val32;
        return tmp0 + p[1] * (tmp1 - tmp0);
      }
    else if (dim == 3)
      {
        const Number3 val10 = Number3(values[stride] - values[0]);
        const Number3 val32 = Number3(values[3 * stride] - values[2 * stride]);
        const Number3 tmp0  = Number3(values[0]) + p[0] * val10;
        const Number3 tmp1  = Number3(values[2 * stride]) + p[0] * val32;
        const Number3 tmpy0 = tmp0 + p[1] * (tmp1 - tmp0);

        const Number3 val54 = Number3(values[5 * stride] - values[4 * stride]);
        const Number3 val76 = Number3(values[7 * stride] - values[6 * stride]);
        const Number3 tmp2  = Number3(values[4 * stride]) + p[0] * val54;
        const Number3 tmp3  = Number3(values[6 * stride]) + p[0] * val76;
        const Number3 tmpy1 = tmp2 + p[1] * (tmp3 - tmp2);

        return tmpy0 + p[2] * (tmpy1 - tmpy0);
      }

    // work around a compile error: missing return statement
    return Number3();
  }



  template <int dim, typename Number, typename Number2>
  inline typename ProductTypeNoPoint<Number, Number2>::type
  evaluate_tensor_product_value(
    const std::vector<Polynomials::Polynomial<double>> &poly,
    const std::vector<Number>                          &values,
    const Point<dim, Number2>                          &p,
    const bool                                          d_linear = false,
    const std::vector<unsigned int>                    &renumber = {})
  {
    typename ProductTypeNoPoint<Number, Number2>::type result;
    if (d_linear)
      {
        result = evaluate_tensor_product_value_linear(values.data(), p);
      }
    else
      {
        AssertIndexRange(poly.size(), 200);
        std::array<dealii::ndarray<Number2, 2, dim>, 200> shapes;
        const int                n_shapes = poly.size();
        std::array<Number2, dim> point;
        for (unsigned int d = 0; d < dim; ++d)
          point[d] = p[d];
        for (int i = 0; i < n_shapes; ++i)
          poly[i].values_of_array(point, 0, shapes[i].data());
        result = evaluate_tensor_product_value_shapes<dim, Number, Number2>(
          shapes.data(), n_shapes, values.data(), renumber);
      }
    return result;
  }



  /**
   * This function computes derivatives of arbitrary orders in 1d, returning a
   * Tensor with the respective derivative
   */
  template <int derivative_order, typename Number, typename Number2>
  inline Tensor<1, 1, typename ProductTypeNoPoint<Number, Number2>::type>
  evaluate_tensor_product_higher_derivatives(
    const std::vector<Polynomials::Polynomial<double>> &poly,
    const std::vector<Number>                          &values,
    const Point<1, Number2>                            &p,
    const std::vector<unsigned int>                    &renumber = {})
  {
    using Number3 = typename ProductTypeNoPoint<Number, Number2>::type;

    const int n_shapes = poly.size();
    AssertDimension(n_shapes, values.size());
    Assert(renumber.empty() || renumber.size() == values.size(),
           ExcDimensionMismatch(renumber.size(), values.size()));

    std::array<Number2, derivative_order + 1> shapes;
    Tensor<1, 1, Number3>                     result;
    if (renumber.empty())
      for (int i = 0; i < n_shapes; ++i)
        {
          poly[i].value(p[0], derivative_order, shapes.data());
          result[0] += shapes[derivative_order] * values[i];
        }
    else
      for (int i = 0; i < n_shapes; ++i)
        {
          poly[i].value(p[0], derivative_order, shapes.data());
          result[0] += shapes[derivative_order] * values[renumber[i]];
        }
    return result;
  }



  /**
   * This function computes derivatives of arbitrary orders in 2d, returning a
   * Tensor with the respective derivatives
   */
  template <int derivative_order, typename Number, typename Number2>
  inline Tensor<1,
                derivative_order + 1,
                typename ProductTypeNoPoint<Number, Number2>::type>
  evaluate_tensor_product_higher_derivatives(
    const std::vector<Polynomials::Polynomial<double>> &poly,
    const std::vector<Number>                          &values,
    const Point<2, Number2>                            &p,
    const std::vector<unsigned int>                    &renumber = {})
  {
    using Number3     = typename ProductTypeNoPoint<Number, Number2>::type;
    constexpr int dim = 2;

    const int n_shapes = poly.size();
    AssertDimension(Utilities::pow(n_shapes, 2), values.size());
    Assert(renumber.empty() || renumber.size() == values.size(),
           ExcDimensionMismatch(renumber.size(), values.size()));

    AssertIndexRange(n_shapes, 100);
    dealii::ndarray<Number2, 100, derivative_order + 1, dim> shapes;
    // Evaluate 1d polynomials and their derivatives
    std::array<Number2, dim> point;
    for (unsigned int d = 0; d < dim; ++d)
      point[d] = p[d];
    for (int i = 0; i < n_shapes; ++i)
      poly[i].values_of_array(point, derivative_order, &shapes[i][0]);

    Tensor<1, derivative_order + 1, Number3> result;
    for (int i1 = 0, i = 0; i1 < n_shapes; ++i1)
      {
        Tensor<1, derivative_order + 1, Number3> result_x;
        if (renumber.empty())
          for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
            for (unsigned int d = 0; d <= derivative_order; ++d)
              result_x[d] += shapes[i0][d][0] * values[i];
        else
          for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
            for (unsigned int d = 0; d <= derivative_order; ++d)
              result_x[d] += shapes[i0][d][0] * values[renumber[i]];

        for (unsigned int d = 0; d <= derivative_order; ++d)
          result[d] += shapes[i1][d][1] * result_x[derivative_order - d];
      }
    return result;
  }



  /**
   * This function computes derivatives of arbitrary orders in 3d, returning a
   * Tensor with the respective derivatives
   */
  template <int derivative_order, typename Number, typename Number2>
  inline Tensor<1,
                ((derivative_order + 1) * (derivative_order + 2)) / 2,
                typename ProductTypeNoPoint<Number, Number2>::type>
  evaluate_tensor_product_higher_derivatives(
    const std::vector<Polynomials::Polynomial<double>> &poly,
    const std::vector<Number>                          &values,
    const Point<3, Number2>                            &p,
    const std::vector<unsigned int>                    &renumber = {})
  {
    using Number3     = typename ProductTypeNoPoint<Number, Number2>::type;
    constexpr int dim = 3;
    constexpr int n_derivatives =
      ((derivative_order + 1) * (derivative_order + 2)) / 2;

    const int n_shapes = poly.size();
    AssertDimension(Utilities::pow(n_shapes, 3), values.size());
    Assert(renumber.empty() || renumber.size() == values.size(),
           ExcDimensionMismatch(renumber.size(), values.size()));

    AssertIndexRange(n_shapes, 100);
    dealii::ndarray<Number2, 100, derivative_order + 1, dim> shapes;
    // Evaluate 1d polynomials and their derivatives
    std::array<Number2, dim> point;
    for (unsigned int d = 0; d < dim; ++d)
      point[d] = p[d];
    for (int i = 0; i < n_shapes; ++i)
      poly[i].values_of_array(point, derivative_order, &shapes[i][0]);

    Tensor<1, n_derivatives, Number3> result;
    for (int i2 = 0, i = 0; i2 < n_shapes; ++i2)
      {
        Tensor<1, n_derivatives, Number3> result_xy;
        for (int i1 = 0; i1 < n_shapes; ++i1)
          {
            // apply x derivatives
            Tensor<1, derivative_order + 1, Number3> result_x;
            if (renumber.empty())
              for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
                for (unsigned int d = 0; d <= derivative_order; ++d)
                  result_x[d] += shapes[i0][d][0] * values[i];
            else
              for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
                for (unsigned int d = 0; d <= derivative_order; ++d)
                  result_x[d] += shapes[i0][d][0] * values[renumber[i]];

            // multiply by y derivatives, sorting them in upper triangular
            // matrix, starting with highest global derivative order,
            // decreasing the combined order of xy derivatives by one in each
            // row, to be combined with z derivatives in the next step
            for (unsigned int d = 0, c = 0; d <= derivative_order; ++d)
              for (unsigned int e = d; e <= derivative_order; ++e, ++c)
                result_xy[c] +=
                  shapes[i1][e - d][1] * result_x[derivative_order - e];
          }

        // multiply by z derivatives, starting with highest x derivative
        for (unsigned int d = 0, c = 0; d <= derivative_order; ++d)
          for (unsigned int e = d; e <= derivative_order; ++e, ++c)
            result[c] += shapes[i2][d][2] * result_xy[c];
      }
    return result;
  }



  template <int dim, typename Number, typename Number2>
  SymmetricTensor<2, dim, typename ProductTypeNoPoint<Number, Number2>::type>
  evaluate_tensor_product_hessian(
    const std::vector<Polynomials::Polynomial<double>> &poly,
    const std::vector<Number>                          &values,
    const Point<dim, Number2>                          &p,
    const std::vector<unsigned int>                    &renumber = {})
  {
    static_assert(dim >= 1 && dim <= 3, "Only dim=1,2,3 implemented");

    const auto hessian =
      evaluate_tensor_product_higher_derivatives<2>(poly, values, p, renumber);

    using Number3 = typename ProductTypeNoPoint<Number, Number2>::type;
    SymmetricTensor<2, dim, Number3> result;
    if (dim == 1)
      result[0][0] = hessian[0];
    else if (dim >= 2)
      {
        // derivatives in Hessian are xx, xy, yy, xz, yz, zz, so must re-order
        // them for 3D
        for (unsigned int d = 0, c = 0; d < 2; ++d)
          for (unsigned int e = d; e < 2; ++e, ++c)
            result[d][e] = hessian[c];
        if (dim == 3)
          {
            for (unsigned int d = 0; d < 2; ++d)
              result[d][2] = hessian[3 + d];
            result[2][2] = hessian[5];
          }
      }

    return result;
  }



  /**
   * Test inner dimensions of tensor product shape functions and accumulate.
   */
  template <int dim,
            int length,
            typename Number2,
            typename Number,
            bool add,
            int  n_values = 1>
  inline
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
    void
    do_apply_test_functions_xy(
      Number2                                 *values,
      const dealii::ndarray<Number, 2, dim>   *shapes,
      const std::array<Number2, 2 + n_values> &test_grads_value,
      const int                                n_shapes_runtime,
      int                                     &i)
  {
    static_assert(0 <= dim && dim <= 3, "Only dim=0,1,2,3 implemented");
    static_assert(1 <= n_values && n_values <= 2,
                  "Only n_values=1,2 implemented");

    // Note that 'add' is a template argument, so the compiler will remove
    // these checks
    if (length > 0)
      {
        constexpr unsigned int         array_size = length > 0 ? length : 1;
        std::array<Number, array_size> shape_values_x;
        std::array<Number, array_size> shape_derivs_x;
        for (unsigned int j = 0; j < array_size; ++j)
          {
            shape_values_x[j] = shapes[j][0][0];
            shape_derivs_x[j] = shapes[j][1][0];
          }
        for (int i1 = 0; i1 < (dim > 1 ? length : 1); ++i1)
          {
            const Number2 test_value_y =
              dim > 1 ? (test_grads_value[2] * shapes[i1][0][1] +
                         test_grads_value[1] * shapes[i1][1][1]) :
                        test_grads_value[2];
            const Number2 test_grad_xy =
              dim > 1 ? test_grads_value[0] * shapes[i1][0][1] :
                        test_grads_value[0];
            Number2 test_value_y_2;
            if (n_values > 1)
              test_value_y_2 = dim > 1 ?
                                 test_grads_value[3] * shapes[i1][0][1] :
                                 test_grads_value[3];

            Number2 *values_ptr = values + i + i1 * length;
            Number2 *values_ptr_2 =
              n_values > 1 ? values_ptr + Utilities::pow(length, dim) : nullptr;
            for (int i0 = 0; i0 < length; ++i0)
              {
                if (add)
                  values_ptr[i0] += shape_values_x[i0] * test_value_y;
                else
                  values_ptr[i0] = shape_values_x[i0] * test_value_y;
                values_ptr[i0] += shape_derivs_x[i0] * test_grad_xy;
                if (n_values > 1)
                  {
                    if (add)
                      values_ptr_2[i0] += shape_values_x[i0] * test_value_y_2;
                    else
                      values_ptr_2[i0] = shape_values_x[i0] * test_value_y_2;
                  }
              }
          }
        i += (dim > 1 ? length * length : length);
      }
    else
      {
        for (int i1 = 0; i1 < (dim > 1 ? n_shapes_runtime : 1); ++i1)
          {
            const Number2 test_value_y =
              dim > 1 ? (test_grads_value[2] * shapes[i1][0][1] +
                         test_grads_value[1] * shapes[i1][1][1]) :
                        test_grads_value[2];
            const Number2 test_grad_xy =
              dim > 1 ? test_grads_value[0] * shapes[i1][0][1] :
                        test_grads_value[0];
            Number2 test_value_y_2;
            if (n_values > 1)
              test_value_y_2 = dim > 1 ?
                                 test_grads_value[3] * shapes[i1][0][1] :
                                 test_grads_value[3];

            Number2 *values_ptr = values + i + i1 * n_shapes_runtime;
            Number2 *values_ptr_2 =
              n_values > 1 ?
                values_ptr + Utilities::fixed_power<dim>(n_shapes_runtime) :
                nullptr;
            for (int i0 = 0; i0 < n_shapes_runtime; ++i0)
              {
                if (add)
                  values_ptr[i0] += shapes[i0][0][0] * test_value_y;
                else
                  values_ptr[i0] = shapes[i0][0][0] * test_value_y;
                values_ptr[i0] += shapes[i0][1][0] * test_grad_xy;
                if (n_values > 1)
                  {
                    if (add)
                      values_ptr_2[i0] += shapes[i0][0][0] * test_value_y_2;
                    else
                      values_ptr_2[i0] = shapes[i0][0][0] * test_value_y_2;
                  }
              }
          }
        i += (dim > 1 ? n_shapes_runtime * n_shapes_runtime : n_shapes_runtime);
      }
  }



  /**
   * Same as evaluate_tensor_product_value_and_gradient_shapes() but for
   * integration.
   */
  template <int dim,
            typename Number,
            typename Number2,
            bool add,
            int  n_values = 1>
  inline void
  integrate_add_tensor_product_value_and_gradient_shapes(
    const dealii::ndarray<Number, 2, dim> *shapes,
    const int                              n_shapes,
    const Number2                         *value,
    const Tensor<1, dim, Number2>         &gradient,
    Number2                               *values)
  {
    static_assert(0 <= dim && dim <= 3, "Only dim=0,1,2,3 implemented");
    static_assert(1 <= n_values && n_values <= 2,
                  "Only n_values=1,2 implemented");

    // Note that 'add' is a template argument, so the compiler will remove
    // these checks
    if (dim == 0)
      {
        if (add)
          values[0] += value[0];
        else
          values[0] = value[0];
        if (n_values > 1)
          {
            if (add)
              values[1] += value[1];
            else
              values[1] = value[1];
          }
        return;
      }

    // Implement the transpose of the function above
    // as in evaluate, use `int` type to produce better code in this context
    std::array<Number2, 2 + n_values> test_grads_value;
    for (int i2 = 0, i = 0; i2 < (dim > 2 ? n_shapes : 1); ++i2)
      {
        // test grad x
        test_grads_value[0] =
          dim > 2 ? gradient[0] * shapes[i2][0][2] : gradient[0];
        // test grad y
        test_grads_value[1] = dim > 2 ? gradient[1] * shapes[i2][0][2] :
                                        (dim > 1 ? gradient[1] : Number2());
        // test value z
        test_grads_value[2] =
          dim > 2 ?
            (value[0] * shapes[i2][0][2] + gradient[2] * shapes[i2][1][2]) :
            value[0];

        if (n_values > 1)
          test_grads_value[3] =
            dim > 2 ? value[1] * shapes[i2][0][2] : value[1];
        // Generate separate code with known loop bounds for the most common
        // cases
        if (n_shapes == 2)
          do_apply_test_functions_xy<dim, 2, Number2, Number, add, n_values>(
            values, shapes, test_grads_value, n_shapes, i);
        else if (n_shapes == 3)
          do_apply_test_functions_xy<dim, 3, Number2, Number, add, n_values>(
            values, shapes, test_grads_value, n_shapes, i);
        else if (n_shapes == 4)
          do_apply_test_functions_xy<dim, 4, Number2, Number, add, n_values>(
            values, shapes, test_grads_value, n_shapes, i);
        else if (n_shapes == 5)
          do_apply_test_functions_xy<dim, 5, Number2, Number, add, n_values>(
            values, shapes, test_grads_value, n_shapes, i);
        else if (n_shapes == 6)
          do_apply_test_functions_xy<dim, 6, Number2, Number, add, n_values>(
            values, shapes, test_grads_value, n_shapes, i);
        else
          do_apply_test_functions_xy<dim, -1, Number2, Number, add, n_values>(
            values, shapes, test_grads_value, n_shapes, i);
      }
  }



  /**
   * Specializes @p integrate_add_tensor_product_value_and_gradient_shapes() for linear
   * polynomials which massively reduces the necessary instructions.
   */
  template <int dim,
            typename Number,
            typename Number2,
            bool add,
            int  n_values = 1>
  inline void
  integrate_add_tensor_product_value_and_gradient_linear(
    const Number2                 *value,
    const Tensor<1, dim, Number2> &gradient,
    Number2                       *values,
    const Point<dim, Number>      &p)
  {
    static_assert(0 <= dim && dim <= 3, "Only dim=0,1,2,3 implemented");
    static_assert(1 <= n_values && n_values <= 2,
                  "Only n_values=1,2 implemented");

    // Note that 'add' is a template argument, so the compiler will remove
    // these checks
    if (dim == 0)
      {
        if (add)
          values[0] += value[0];
        else
          values[0] = value[0];
        if (n_values > 1)
          {
            if (add)
              values[1] += value[1];
            else
              values[1] = value[1];
          }
      }
    else if (dim == 1)
      {
        const Number2 difference = value[0] * p[0] + gradient[0];
        if (add)
          {
            values[0] += value[0] - difference;
            values[1] += difference;
          }
        else
          {
            values[0] = value[0] - difference;
            values[1] = difference;
          }
        if (n_values > 1)
          {
            const Number2 product = value[1] * p[0];
            if (add)
              {
                values[2] += value[1] - product;
                values[3] += product;
              }
            else
              {
                values[2] = value[1] - product;
                values[3] = product;
              }
          }
      }
    else if (dim == 2)
      {
        const Number2 test_value_y1 = value[0] * p[1] + gradient[1];
        const Number2 test_value_y0 = value[0] - test_value_y1;
        const Number2 test_grad_xy1 = gradient[0] * p[1];
        const Number2 test_grad_xy0 = gradient[0] - test_grad_xy1;
        const Number2 value0        = p[0] * test_value_y0 + test_grad_xy0;
        const Number2 value1        = p[0] * test_value_y1 + test_grad_xy1;

        if (add)
          {
            values[0] += test_value_y0 - value0;
            values[1] += value0;
            values[2] += test_value_y1 - value1;
            values[3] += value1;
          }
        else
          {
            values[0] = test_value_y0 - value0;
            values[1] = value0;
            values[2] = test_value_y1 - value1;
            values[3] = value1;
          }

        if (n_values > 1)
          {
            const Number2 test_value_y1_2 = value[1] * p[1];
            const Number2 test_value_y0_2 = value[1] - test_value_y1_2;
            const Number2 value0_2        = p[0] * test_value_y0_2;
            const Number2 value1_2        = p[0] * test_value_y1_2;

            if (add)
              {
                values[4] += test_value_y0_2 - value0_2;
                values[5] += value0_2;
                values[6] += test_value_y1_2 - value1_2;
                values[7] += value1_2;
              }
            else
              {
                values[4] = test_value_y0_2 - value0_2;
                values[5] = value0_2;
                values[6] = test_value_y1_2 - value1_2;
                values[7] = value1_2;
              }
          }
      }
    else if (dim == 3)
      {
        Assert(n_values == 1, ExcNotImplemented());

        const Number2 test_value_z1 = value[0] * p[2] + gradient[2];
        const Number2 test_value_z0 = value[0] - test_value_z1;
        const Number2 test_grad_x1  = gradient[0] * p[2];
        const Number2 test_grad_x0  = gradient[0] - test_grad_x1;
        const Number2 test_grad_y1  = gradient[1] * p[2];
        const Number2 test_grad_y0  = gradient[1] - test_grad_y1;

        const Number2 test_value_y01 = test_value_z0 * p[1] + test_grad_y0;
        const Number2 test_value_y00 = test_value_z0 - test_value_y01;
        const Number2 test_grad_xy01 = test_grad_x0 * p[1];
        const Number2 test_grad_xy00 = test_grad_x0 - test_grad_xy01;
        const Number2 test_value_y11 = test_value_z1 * p[1] + test_grad_y1;
        const Number2 test_value_y10 = test_value_z1 - test_value_y11;
        const Number2 test_grad_xy11 = test_grad_x1 * p[1];
        const Number2 test_grad_xy10 = test_grad_x1 - test_grad_xy11;

        const Number2 value00 = p[0] * test_value_y00 + test_grad_xy00;
        const Number2 value01 = p[0] * test_value_y01 + test_grad_xy01;
        const Number2 value10 = p[0] * test_value_y10 + test_grad_xy10;
        const Number2 value11 = p[0] * test_value_y11 + test_grad_xy11;

        if (add)
          {
            values[0] += test_value_y00 - value00;
            values[1] += value00;
            values[2] += test_value_y01 - value01;
            values[3] += value01;
            values[4] += test_value_y10 - value10;
            values[5] += value10;
            values[6] += test_value_y11 - value11;
            values[7] += value11;
          }
        else
          {
            values[0] = test_value_y00 - value00;
            values[1] = value00;
            values[2] = test_value_y01 - value01;
            values[3] = value01;
            values[4] = test_value_y10 - value10;
            values[5] = value10;
            values[6] = test_value_y11 - value11;
            values[7] = value11;
          }
      }
  }



  /**
   * Calls the correct @p integrate_add_tensor_product_value_and_gradient_...()
   * function depending on if values should be added to or set and if
   * polynomials are linear.
   */
  template <bool is_linear,
            int  dim,
            typename Number,
            typename Number2,
            int n_values = 1>
  inline void
  integrate_tensor_product_value_and_gradient(
    const dealii::ndarray<Number, 2, dim> *shapes,
    const unsigned int                     n_shapes,
    const Number2                         *value,
    const Tensor<1, dim, Number2>         &gradient,
    Number2                               *values,
    const Point<dim, Number>              &p,
    const bool                             do_add)
  {
    if (do_add)
      {
        if (is_linear)
          internal::integrate_add_tensor_product_value_and_gradient_linear<
            dim,
            Number,
            Number2,
            true,
            n_values>(value, gradient, values, p);
        else
          internal::integrate_add_tensor_product_value_and_gradient_shapes<
            dim,
            Number,
            Number2,
            true,
            n_values>(shapes, n_shapes, value, gradient, values);
      }
    else
      {
        if (is_linear)
          internal::integrate_add_tensor_product_value_and_gradient_linear<
            dim,
            Number,
            Number2,
            false,
            n_values>(value, gradient, values, p);
        else
          internal::integrate_add_tensor_product_value_and_gradient_shapes<
            dim,
            Number,
            Number2,
            false,
            n_values>(shapes, n_shapes, value, gradient, values);
      }
  }



  /**
   * Test inner dimensions of tensor product shape functions and accumulate.
   */
  template <int dim, int length, typename Number2, typename Number, bool add>
  inline
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
    void
    do_apply_test_functions_xy_value(
      Number2                               *values,
      const dealii::ndarray<Number, 2, dim> *shapes,
      const Number2                         &test_value,
      const int                              n_shapes_runtime,
      int                                   &i)
  {
    if (length > 0)
      {
        constexpr unsigned int         array_size = length > 0 ? length : 1;
        std::array<Number, array_size> shape_values_x;
        for (unsigned int i1 = 0; i1 < array_size; ++i1)
          shape_values_x[i1] = shapes[i1][0][0];
        for (unsigned int i1 = 0; i1 < (dim > 1 ? length : 1); ++i1)
          {
            const Number2 test_value_y =
              dim > 1 ? test_value * shapes[i1][0][1] : test_value;

            Number2 *values_ptr = values + i + i1 * length;
            for (unsigned int i0 = 0; i0 < length; ++i0)
              {
                if (add)
                  values_ptr[i0] += shape_values_x[i0] * test_value_y;
                else
                  values_ptr[i0] = shape_values_x[i0] * test_value_y;
              }
          }
        i += (dim > 1 ? length * length : length);
      }
    else
      {
        for (int i1 = 0; i1 < (dim > 1 ? n_shapes_runtime : 1); ++i1)
          {
            const Number2 test_value_y =
              dim > 1 ? test_value * shapes[i1][0][1] : test_value;

            Number2 *values_ptr = values + i + i1 * n_shapes_runtime;
            for (int i0 = 0; i0 < n_shapes_runtime; ++i0)
              {
                if (add)
                  values_ptr[i0] += shapes[i0][0][0] * test_value_y;
                else
                  values_ptr[i0] = shapes[i0][0][0] * test_value_y;
              }
          }
        i += (dim > 1 ? n_shapes_runtime * n_shapes_runtime : n_shapes_runtime);
      }
  }



  /**
   * Same as evaluate_tensor_product_value_shapes() but for integration.
   */
  template <int dim, typename Number, typename Number2, bool add>
  inline void
  integrate_add_tensor_product_value_shapes(
    const dealii::ndarray<Number, 2, dim> *shapes,
    const int                              n_shapes,
    const Number2                         &value,
    Number2                               *values)
  {
    static_assert(dim >= 0 && dim <= 3, "Only dim=0,1,2,3 implemented");

    // as in evaluate, use `int` type to produce better code in this context

    if (dim == 0)
      {
        if (add)
          values[0] += value;
        else
          values[0] = value;
        return;
      }

    // Implement the transpose of the function above
    Number2 test_value;
    for (int i2 = 0, i = 0; i2 < (dim > 2 ? n_shapes : 1); ++i2)
      {
        // test value z
        test_value = dim > 2 ? value * shapes[i2][0][2] : value;

        // Generate separate code with known loop bounds for the most common
        // cases
        if (n_shapes == 2)
          do_apply_test_functions_xy_value<dim, 2, Number2, Number, add>(
            values, shapes, test_value, n_shapes, i);
        else if (n_shapes == 3)
          do_apply_test_functions_xy_value<dim, 3, Number2, Number, add>(
            values, shapes, test_value, n_shapes, i);
        else if (n_shapes == 4)
          do_apply_test_functions_xy_value<dim, 4, Number2, Number, add>(
            values, shapes, test_value, n_shapes, i);
        else if (n_shapes == 5)
          do_apply_test_functions_xy_value<dim, 5, Number2, Number, add>(
            values, shapes, test_value, n_shapes, i);
        else if (n_shapes == 6)
          do_apply_test_functions_xy_value<dim, 6, Number2, Number, add>(
            values, shapes, test_value, n_shapes, i);
        else
          do_apply_test_functions_xy_value<dim, -1, Number2, Number, add>(
            values, shapes, test_value, n_shapes, i);
      }
  }



  /**
   * Specializes @p integrate_tensor_product_value_shapes() for linear
   * polynomials which massively reduces the necessary instructions.
   */
  template <int dim, typename Number, typename Number2, bool add>
  inline void
  integrate_add_tensor_product_value_linear(const Number2            &value,
                                            Number2                  *values,
                                            const Point<dim, Number> &p)
  {
    static_assert(dim >= 0 && dim <= 3, "Only dim=0,1,2,3 implemented");

    if (dim == 0)
      {
        if (add)
          values[0] += value;
        else
          values[0] = value;
      }
    else if (dim == 1)
      {
        const auto x0 = 1. - p[0], x1 = p[0];

        if (add)
          {
            values[0] += value * x0;
            values[1] += value * x1;
          }
        else
          {
            values[0] = value * x0;
            values[1] = value * x1;
          }
      }
    else if (dim == 2)
      {
        const auto x0 = 1. - p[0], x1 = p[0], y0 = 1. - p[1], y1 = p[1];

        const auto test_value_y0 = value * y0;
        const auto test_value_y1 = value * y1;

        if (add)
          {
            values[0] += x0 * test_value_y0;
            values[1] += x1 * test_value_y0;
            values[2] += x0 * test_value_y1;
            values[3] += x1 * test_value_y1;
          }
        else
          {
            values[0] = x0 * test_value_y0;
            values[1] = x1 * test_value_y0;
            values[2] = x0 * test_value_y1;
            values[3] = x1 * test_value_y1;
          }
      }
    else if (dim == 3)
      {
        const auto x0 = 1. - p[0], x1 = p[0], y0 = 1. - p[1], y1 = p[1],
                   z0 = 1. - p[2], z1 = p[2];

        const auto test_value_z0 = value * z0;
        const auto test_value_z1 = value * z1;

        const auto test_value_y00 = test_value_z0 * y0;
        const auto test_value_y01 = test_value_z0 * y1;
        const auto test_value_y10 = test_value_z1 * y0;
        const auto test_value_y11 = test_value_z1 * y1;

        if (add)
          {
            values[0] += x0 * test_value_y00;
            values[1] += x1 * test_value_y00;
            values[2] += x0 * test_value_y01;
            values[3] += x1 * test_value_y01;
            values[4] += x0 * test_value_y10;
            values[5] += x1 * test_value_y10;
            values[6] += x0 * test_value_y11;
            values[7] += x1 * test_value_y11;
          }
        else
          {
            values[0] = x0 * test_value_y00;
            values[1] = x1 * test_value_y00;
            values[2] = x0 * test_value_y01;
            values[3] = x1 * test_value_y01;
            values[4] = x0 * test_value_y10;
            values[5] = x1 * test_value_y10;
            values[6] = x0 * test_value_y11;
            values[7] = x1 * test_value_y11;
          }
      }
  }



  /**
   * Calls the correct @p integrate_add_tensor_product_value_...()
   * function depending on if values should be added to or set and if
   * polynomials are linear.
   */
  template <bool is_linear, int dim, typename Number, typename Number2>
  inline void
  integrate_tensor_product_value(const dealii::ndarray<Number, 2, dim> *shapes,
                                 const unsigned int        n_shapes,
                                 const Number2            &value,
                                 Number2                  *values,
                                 const Point<dim, Number> &p,
                                 const bool                do_add)
  {
    if (do_add)
      {
        if (is_linear)
          internal::integrate_add_tensor_product_value_linear<dim,
                                                              Number,
                                                              Number2,
                                                              true>(value,
                                                                    values,
                                                                    p);
        else
          internal::integrate_add_tensor_product_value_shapes<dim,
                                                              Number,
                                                              Number2,
                                                              true>(shapes,
                                                                    n_shapes,
                                                                    value,
                                                                    values);
      }
    else
      {
        if (is_linear)
          internal::integrate_add_tensor_product_value_linear<dim,
                                                              Number,
                                                              Number2,
                                                              false>(value,
                                                                     values,
                                                                     p);
        else
          internal::integrate_add_tensor_product_value_shapes<dim,
                                                              Number,
                                                              Number2,
                                                              false>(shapes,
                                                                     n_shapes,
                                                                     value,
                                                                     values);
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
