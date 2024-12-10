// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_tensor_product_point_kernels_h
#define dealii_matrix_free_tensor_product_point_kernels_h

#include <deal.II/base/config.h>

#include <deal.II/base/ndarray.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/utilities.h>

#include <deal.II/matrix_free/shape_info.h>


DEAL_II_NAMESPACE_OPEN



namespace internal
{
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
   * Compute the values and derivatives of the 1d polynomials @p poly at the
   * specified point @p p and store them in @p shapes.
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
    DEAL_II_ASSERT_UNREACHABLE();
  }



  /**
   * Evaluate the 1d polynomials @p poly at the two specified points @p p0 and
   * @p p1 and store them in @p shapes. This function can be used as a more
   * efficient alternative to the compute_values_of_array() function, because
   * of reduced overhead when querying the polynomials (which usually have
   * loop bounds that are not known at compile time).
   */
  template <int dim, typename Number>
  inline void
  compute_values_of_array_in_pairs(
    dealii::ndarray<Number, 2, dim>                    *shapes,
    const std::vector<Polynomials::Polynomial<double>> &poly,
    const Point<dim, Number>                           &p0,
    const Point<dim, Number>                           &p1)
  {
    // Use 'int' variable here to let the compiler apply additional
    // optimizations, in particular regarding multiplications and additions in
    // loop increments that are known not to overflow/wrap around (as is the
    // case for unsigned int).
    const int n_shapes = poly.size();

    std::array<Number, 2 * dim> point, result;
    for (unsigned int d = 0; d < dim; ++d)
      point[d] = p0[d];
    for (unsigned int d = 0; d < dim; ++d)
      point[dim + d] = p1[d];
    for (int i = 0; i < n_shapes; ++i)
      {
        poly[i].values_of_array(point, 0, &result);
        for (unsigned int j = 0; j < 2; ++j)
          for (unsigned int d = 0; d < dim; ++d)
            shapes[j * n_shapes + i][0][d] = result[j * dim + d];
      }
  }



  /**
   * Interpolate inner dimensions of tensor product shape functions.
   */
  template <int dim,
            int length,
            typename Number2,
            typename Number,
            int  n_values    = 1,
            bool do_renumber = true,
            int  stride      = 1>
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
        values + stride * (length > 0 ?
                             Utilities::pow(length, dim) :
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
              inner_result[0] +=
                shapes[i0][1][0] * values[renumber[i] * stride];
              // values
              inner_result[1] +=
                shapes[i0][0][0] * values[renumber[i] * stride];
              if (n_values > 1)
                inner_result[2] +=
                  shapes[i0][0][0] * values_2[renumber[i] * stride];
            }
        else
          for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
            {
              // gradient
              inner_result[0] += shapes[i0][1][0] * values[i * stride];
              // values
              inner_result[1] += shapes[i0][0][0] * values[i * stride];
              if (n_values > 1)
                inner_result[2] += shapes[i0][0][0] * values_2[i * stride];
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
            bool do_renumber = true,
            int  stride      = 1>
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
          result[1] = values[1 * stride];
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
            do_interpolate_xy<dim,
                              2,
                              Number2,
                              Number,
                              n_values,
                              do_renumber,
                              stride>(values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 3)
          inner_result =
            do_interpolate_xy<dim,
                              3,
                              Number2,
                              Number,
                              n_values,
                              do_renumber,
                              stride>(values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 4)
          inner_result =
            do_interpolate_xy<dim,
                              4,
                              Number2,
                              Number,
                              n_values,
                              do_renumber,
                              stride>(values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 5)
          inner_result =
            do_interpolate_xy<dim,
                              5,
                              Number2,
                              Number,
                              n_values,
                              do_renumber,
                              stride>(values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 6)
          inner_result =
            do_interpolate_xy<dim,
                              6,
                              Number2,
                              Number,
                              n_values,
                              do_renumber,
                              stride>(values, renumber, shapes, n_shapes, i);
        else
          inner_result =
            do_interpolate_xy<dim,
                              -1,
                              Number2,
                              Number,
                              n_values,
                              do_renumber,
                              stride>(values, renumber, shapes, n_shapes, i);
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

    // If n_values > 1, we want to interpolate from a second array,
    // placed in the same array immediately after the main data. This
    // is used to interpolate normal derivatives onto faces.

    std::array<Number3, dim + n_values> result;
    if (dim == 0)
      {
        // we only need the value on faces of a 1d element
        result[0] = values[0];
        if (n_values > 1)
          result[1] = values[1 * stride];
      }
    else if (dim == 1)
      {
        // gradient
        result[0] = Number3(values[stride] - values[0]);
        // values
        result[1] = Number3(values[0]) + p[0] * result[0];
        if (n_values > 1)
          result[2] = Number3(values[2 * stride]) +
                      p[0] * (values[3 * stride] - values[2 * stride]);
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
              Number3(values[4 * stride]) +
              p[0] * (values[5 * stride] - values[4 * stride]);
            const Number3 tmp1_2 =
              Number3(values[6 * stride]) +
              p[0] * (values[7 * stride] - values[6 * stride]);
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
    const ArrayView<const Number>                      &values,
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
            bool do_renumber = true,
            int  stride      = 1>
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
            value += shapes[i0][0][0] * values[renumber[i] * stride];
        else
          for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
            value += shapes[i0][0][0] * values[i * stride];

        if (dim > 1)
          result += value * shapes[i1][0][1];
        else
          result = value;
      }
    return result;
  }



  template <int dim,
            typename Number,
            typename Number2,
            bool do_renumber = true,
            int  stride      = 1>
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
          inner_result = do_interpolate_xy_value<dim,
                                                 2,
                                                 Number2,
                                                 Number,
                                                 do_renumber,
                                                 stride>(
            values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 3)
          inner_result = do_interpolate_xy_value<dim,
                                                 3,
                                                 Number2,
                                                 Number,
                                                 do_renumber,
                                                 stride>(
            values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 4)
          inner_result = do_interpolate_xy_value<dim,
                                                 4,
                                                 Number2,
                                                 Number,
                                                 do_renumber,
                                                 stride>(
            values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 5)
          inner_result = do_interpolate_xy_value<dim,
                                                 5,
                                                 Number2,
                                                 Number,
                                                 do_renumber,
                                                 stride>(
            values, renumber, shapes, n_shapes, i);
        else if (n_shapes == 6)
          inner_result = do_interpolate_xy_value<dim,
                                                 6,
                                                 Number2,
                                                 Number,
                                                 do_renumber,
                                                 stride>(
            values, renumber, shapes, n_shapes, i);
        else
          inner_result = do_interpolate_xy_value<dim,
                                                 -1,
                                                 Number2,
                                                 Number,
                                                 do_renumber,
                                                 stride>(
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
    const ArrayView<const Number>                      &values,
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
    const ArrayView<const Number>                      &values,
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
    const ArrayView<const Number>                      &values,
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
    const ArrayView<const Number>                      &values,
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
    const ArrayView<const Number>                      &values,
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
        for (int i1 = 0; i1 < (dim > 1 ? length : 1); ++i1)
          {
            const Number2 test_value_y =
              dim > 1 ? test_value * shapes[i1][0][1] : test_value;

            Number2 *values_ptr = values + i + i1 * length;
            for (int i0 = 0; i0 < length; ++i0)
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
        const auto x0 = Number(1.) - p[0], x1 = p[0];

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
        const auto x0 = Number(1.) - p[0], x1 = p[0], y0 = Number(1.) - p[1],
                   y1 = p[1];

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
        const auto x0 = Number(1.) - p[0], x1 = p[0], y0 = Number(1.) - p[1],
                   y1 = p[1], z0 = Number(1.) - p[2], z1 = p[2];

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
} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
