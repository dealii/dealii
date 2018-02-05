// ---------------------------------------------------------------------
//
// Copyright (C) 2017-2018 by the deal.II authors
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


#ifndef dealii_matrix_free_tensor_product_kernels_h
#define dealii_matrix_free_tensor_product_kernels_h

#include <deal.II/base/config.h>
#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/utilities.h>


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
   * Generic evaluator framework that valuates the given shape data in general
   * dimensions using the tensor product form. Depending on the particular
   * layout in the matrix entries, this corresponds to a usual matrix-matrix
   * product or a matrix-matrix product including some symmetries.
  *
   * @tparam variant Variant of evaluation used for creating template
   *                 specializations
   * @tparam dim Dimension of the function
   * @tparam n_rows Number of rows in the transformation matrix, which corresponds
   *                to the number of 1d shape functions in the usual tensor
   *                contraction setting
   * @tparam n_columns Number of columns in the transformation matrix, which
   *                   corresponds to the number of 1d shape functions in the
   *                   usual tensor contraction setting
   * @tparam Number Abstract number type for input and output arrays
   * @tparam Number2 Abstract number type for coefficient arrays (defaults to
   *                 same type as the input/output arrays); must implement
   *                 operator* with Number to be valid
   */
  template <EvaluatorVariant variant, int dim, int n_rows, int n_columns,
            typename Number, typename Number2=Number>
  struct EvaluatorTensorProduct
  {};



  /**
   * Internal evaluator for shape function in arbitrary dimension using the
   * tensor product form of the basis functions.
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
  template <int dim, int n_rows, int n_columns, typename Number, typename Number2>
  struct EvaluatorTensorProduct<evaluate_general,dim,n_rows,n_columns,Number,Number2>
  {
    static const unsigned int n_rows_of_product = Utilities::fixed_int_power<n_rows,dim>::value;
    static const unsigned int n_columns_of_product = Utilities::fixed_int_power<n_columns,dim>::value;

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other pointer
     */
    EvaluatorTensorProduct ()
      :
      shape_values (nullptr),
      shape_gradients (nullptr),
      shape_hessians (nullptr)
    {}

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number2> &shape_values,
                            const AlignedVector<Number2> &shape_gradients,
                            const AlignedVector<Number2> &shape_hessians,
                            const unsigned int            dummy1 = 0,
                            const unsigned int            dummy2 = 0)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {
      // We can enter this function either for the apply() path that has
      // n_rows * n_columns entries or for the apply_face() path that only has
      // n_rows * 3 entries in the array. Since we cannot decide about the use
      // we must allow for both here.
      Assert(shape_values.size() == 0 ||
             shape_values.size() == n_rows*n_columns ||
             shape_values.size() == 3*n_rows,
             ExcDimensionMismatch(shape_values.size(), n_rows*n_columns));
      Assert(shape_gradients.size() == 0 ||
             shape_gradients.size() == n_rows*n_columns,
             ExcDimensionMismatch(shape_gradients.size(), n_rows*n_columns));
      Assert(shape_hessians.size() == 0 ||
             shape_hessians.size() == n_rows*n_columns,
             ExcDimensionMismatch(shape_hessians.size(), n_rows*n_columns));
      (void)dummy1;
      (void)dummy2;
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values (const Number in [],
            Number       out[]) const
    {
      apply<direction,contract_over_rows,add>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const
    {
      apply<direction,contract_over_rows,add>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const
    {
      apply<direction,contract_over_rows,add>(shape_hessians, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values_one_line (const Number in [],
                     Number       out[]) const
    {
      Assert (shape_values != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,true>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients_one_line (const Number in [],
                        Number       out[]) const
    {
      Assert (shape_gradients != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,true>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians_one_line (const Number in [],
                       Number       out[]) const
    {
      Assert (shape_hessians != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,true>(shape_hessians, in, out);
    }

    /**
     * This function applies the tensor product kernel, corresponding to a
     * multiplication of 1D stripes, along the given @p direction of the tensor
     * data in the input array. This function allows the @p in and @p out
     * arrays to alias for the case n_rows == n_columns, i.e., it is safe to
     * perform the contraction in place where @p in and @p out point to the
     * same address. For the case n_rows != n_columns, the output is in general
     * not correct.
     *
     * @tparam direction Direction that is evaluated
     * @tparam contract_over_rows If true, the tensor contraction sums
     *                            over the rows in the given @p shape_data
     *                            array, otherwise it sums over the columns
     * @tparam add If true, the result is added to the output vector, else
     *             the computed values overwrite the content in the output
     * @tparam one_line If true, the kernel is only applied along a single 1D
     *                  stripe within a dim-dimensional tensor, not the full
     *                  n_rows^dim points as in the @p false case.
     *
     * @param shape_data Transformation matrix with @p n_rows rows and
     *                   @p n_columns columns, stored in row-major format
     * @param in Pointer to the start of the input data vector
     * @param out Pointer to the start of the output data vector
     */
    template <int direction, bool contract_over_rows, bool add, bool one_line=false>
    static void apply (const Number2 *DEAL_II_RESTRICT shape_data,
                       const Number  *in,
                       Number        *out);

    /**
     * This function applies the tensor product operation to produce face values
     * from cell values. As opposed to the apply method, this method assumes
     * that the directions orthogonal to the face have n_rows degrees of
     * freedom per direction and not n_columns for those directions lower than
     * the one currently applied. In other words, apply_face() must be called
     * before calling any interpolation within the face.
     *
     * @tparam face_direction Direction of the normal vector (0=x, 1=y, etc)
     * @tparam contract_onto_face If true, the input vector is of size n_rows^dim
     *                            and interpolation into n_rows^(dim-1) points
     *                            is performed. This is a typical scenario in
     *                            FEFaceEvaluation::evaluate() calls. If false,
     *                            data from n_rows^(dim-1) points is expanded
     *                            into the n_rows^dim points of the higher-
     *                            dimensional data array. Derivatives in the
     *                            case contract_onto_face==false are summed
     *                            together
     * @tparam add If true, the result is added to the output vector, else
     *             the computed values overwrite the content in the output
     * @tparam max_derivative Sets the number of derivatives that should be
     *             computed. 0 means only values, 1 means values and first
     *             derivatives, 2 second derivates. Note that all the
     *             derivatives access the data in @p shape_values passed to
     *             the constructor of the class
     *
     * @param in address of the input data vector
     * @param out address of the output data vector
     */
    template <int face_direction, bool contract_onto_face, bool add, int max_derivative>
    void apply_face (const Number *DEAL_II_RESTRICT in,
                     Number       *DEAL_II_RESTRICT out) const;

    const Number2 *shape_values;
    const Number2 *shape_gradients;
    const Number2 *shape_hessians;
  };



  template <int dim, int n_rows, int n_columns, typename Number, typename Number2>
  template <int direction, bool contract_over_rows, bool add, bool one_line>
  inline
  void
  EvaluatorTensorProduct<evaluate_general,dim,n_rows,n_columns,Number,Number2>
  ::apply (const Number2 *DEAL_II_RESTRICT shape_data,
           const Number  *in,
           Number        *out)
  {
    static_assert (one_line == false || direction==dim-1,
                   "Single-line evaluation only works for direction=dim-1.");
    Assert(shape_data != nullptr,
           ExcMessage("The given array shape_data must not be the null pointer!"));
    Assert (dim == direction+1 || one_line == true || n_rows == n_columns || in != out,
            ExcMessage("In-place operation only supported for "
                       "n_rows==n_columns or single-line interpolation"));
    AssertIndexRange (direction, dim);
    constexpr int mm     = contract_over_rows ? n_rows : n_columns,
                  nn     = contract_over_rows ? n_columns : n_rows;

    constexpr int stride    = Utilities::fixed_int_power<n_columns,direction>::value;
    constexpr int n_blocks1 = one_line ? 1 : stride;
    constexpr int n_blocks2 = Utilities::fixed_int_power<n_rows,(direction>=dim)?0:(dim-direction-1)>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            Number x[mm];
            for (int i=0; i<mm; ++i)
              x[i] = in[stride*i];
            for (int col=0; col<nn; ++col)
              {
                Number2 val0;
                if (contract_over_rows == true)
                  val0 = shape_data[col];
                else
                  val0 = shape_data[col*n_columns];
                Number res0 = val0 * x[0];
                for (int i=1; i<mm; ++i)
                  {
                    if (contract_over_rows == true)
                      val0 = shape_data[i*n_columns+col];
                    else
                      val0 = shape_data[col*n_columns+i];
                    res0 += val0 * x[i];
                  }
                if (add == false)
                  out[stride*col]  = res0;
                else
                  out[stride*col] += res0;
              }

            if (one_line == false)
              {
                ++in;
                ++out;
              }
          }
        if (one_line == false)
          {
            in += stride*(mm-1);
            out += stride*(nn-1);
          }
      }
  }



  template <int dim, int n_rows, int n_columns, typename Number, typename Number2>
  template <int face_direction, bool contract_onto_face, bool add, int max_derivative>
  inline
  void
  EvaluatorTensorProduct<evaluate_general,dim,n_rows,n_columns,Number,Number2>
  ::apply_face (const Number *DEAL_II_RESTRICT in,
                Number       *DEAL_II_RESTRICT out) const
  {
    static_assert(dim > 0 && dim<4, "Only dim=1,2,3 supported");
    static_assert(max_derivative >= 0 && max_derivative<3,
                  "Only derivative orders 0-2 implemented");
    Assert(shape_values != nullptr,
           ExcMessage("The given array shape_values must not be the null pointer."));

    constexpr int n_blocks1 = dim > 1 ? n_rows : 1;
    constexpr int n_blocks2 = dim > 2 ? n_rows : 1;

    AssertIndexRange (face_direction, dim);
    constexpr int stride = Utilities::fixed_int_power<n_rows,face_direction>::value;
    constexpr int out_stride = Utilities::fixed_int_power<n_rows,dim-1>::value;
    const Number *DEAL_II_RESTRICT shape_values = this->shape_values;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            if (contract_onto_face == true)
              {
                Number res0 = shape_values[0] * in[0];
                Number res1, res2;
                if (max_derivative > 0)
                  res1 = shape_values[n_rows] * in[0];
                if (max_derivative > 1)
                  res2 = shape_values[2*n_rows] * in[0];
                for (int ind=1; ind<n_rows; ++ind)
                  {
                    res0 += shape_values[ind] * in[stride*ind];
                    if (max_derivative > 0)
                      res1 += shape_values[ind+n_rows] * in[stride*ind];
                    if (max_derivative > 1)
                      res2 += shape_values[ind+2*n_rows] * in[stride*ind];
                  }
                if (add == false)
                  {
                    out[0]              = res0;
                    if (max_derivative > 0)
                      out[out_stride]   = res1;
                    if (max_derivative > 1)
                      out[2*out_stride] = res2;
                  }
                else
                  {
                    out[0]              += res0;
                    if (max_derivative > 0)
                      out[out_stride]   += res1;
                    if (max_derivative > 1)
                      out[2*out_stride] += res2;
                  }
              }
            else
              {
                for (int col=0; col<n_rows; ++col)
                  {
                    if (add == false)
                      out[col*stride]  = shape_values[col] * in[0];
                    else
                      out[col*stride] += shape_values[col] * in[0];
                    if (max_derivative > 0)
                      out[col*stride] += shape_values[col+n_rows] * in[out_stride];
                    if (max_derivative > 1)
                      out[col*stride] += shape_values[col+2*n_rows] * in[2*out_stride];
                  }
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need
            // to jump over to the next layer in z-direction
            switch (face_direction)
              {
              case 0:
                in += contract_onto_face ? n_rows : 1;
                out += contract_onto_face ? 1 : n_rows;
                break;
              case 1:
                ++in;
                ++out;
                // faces 2 and 3 in 3D use local coordinate system zx, which
                // is the other way around compared to the tensor
                // product. Need to take that into account.
                if (dim == 3)
                  {
                    if (contract_onto_face)
                      out += n_rows-1;
                    else
                      in += n_rows-1;
                  }
                break;
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (face_direction == 1 && dim == 3)
          {
            // adjust for local coordinate system zx
            if (contract_onto_face)
              {
                in += n_rows*(n_rows-1);
                out -= n_rows*n_rows-1;
              }
            else
              {
                out += n_rows*(n_rows-1);
                in -= n_rows*n_rows-1;
              }
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
  template <int dim, typename Number, typename Number2>
  struct EvaluatorTensorProduct<evaluate_general,dim,0,0,Number,Number2>
  {
    static const unsigned int n_rows_of_product = numbers::invalid_unsigned_int;
    static const unsigned int n_columns_of_product = numbers::invalid_unsigned_int;

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other constructor
     */
    EvaluatorTensorProduct ()
      :
      shape_values (nullptr),
      shape_gradients (nullptr),
      shape_hessians (nullptr),
      n_rows (numbers::invalid_unsigned_int),
      n_columns (numbers::invalid_unsigned_int)
    {}

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number2> &shape_values,
                            const AlignedVector<Number2> &shape_gradients,
                            const AlignedVector<Number2> &shape_hessians,
                            const unsigned int            n_rows,
                            const unsigned int            n_columns)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin()),
      n_rows (n_rows),
      n_columns (n_columns)
    {
      // We can enter this function either for the apply() path that has
      // n_rows * n_columns entries or for the apply_face() path that only has
      // n_rows * 3 entries in the array. Since we cannot decide about the use
      // we must allow for both here.
      Assert(shape_values.size() == 0 ||
             shape_values.size() == n_rows*n_columns ||
             shape_values.size() == n_rows*3,
             ExcDimensionMismatch(shape_values.size(), n_rows*n_columns));
      Assert(shape_gradients.size() == 0 ||
             shape_gradients.size() == n_rows*n_columns,
             ExcDimensionMismatch(shape_gradients.size(), n_rows*n_columns));
      Assert(shape_hessians.size() == 0 ||
             shape_hessians.size() == n_rows*n_columns,
             ExcDimensionMismatch(shape_hessians.size(), n_rows*n_columns));
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values (const Number *in,
            Number       *out) const
    {
      apply<direction,contract_over_rows,add>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients (const Number *in,
               Number       *out) const
    {
      apply<direction,contract_over_rows,add>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians (const Number *in,
              Number       *out) const
    {
      apply<direction,contract_over_rows,add>(shape_hessians, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values_one_line (const Number in [],
                     Number       out[]) const
    {
      Assert (shape_values != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,true>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients_one_line (const Number in [],
                        Number       out[]) const
    {
      Assert (shape_gradients != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,true>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians_one_line (const Number in [],
                       Number       out[]) const
    {
      Assert (shape_hessians != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,true>(shape_hessians, in, out);
    }

    template <int direction, bool contract_over_rows, bool add, bool one_line=false>
    void apply (const Number2 *DEAL_II_RESTRICT shape_data,
                const Number  *in,
                Number        *out) const;

    template <int face_direction, bool contract_onto_face, bool add, int max_derivative>
    void apply_face (const Number *DEAL_II_RESTRICT in,
                     Number       *DEAL_II_RESTRICT out) const;

    const Number2 *shape_values;
    const Number2 *shape_gradients;
    const Number2 *shape_hessians;
    const unsigned int n_rows;
    const unsigned int n_columns;
  };



  template <int dim, typename Number, typename Number2>
  template <int direction, bool contract_over_rows, bool add, bool one_line>
  inline
  void
  EvaluatorTensorProduct<evaluate_general,dim,0,0,Number,Number2>
  ::apply (const Number2 *DEAL_II_RESTRICT shape_data,
           const Number  *in,
           Number        *out) const
  {
    static_assert (one_line == false || direction==dim-1,
                   "Single-line evaluation only works for direction=dim-1.");
    Assert(shape_data != nullptr,
           ExcMessage("The given array shape_data must not be the null pointer!"));
    Assert (dim == direction+1 || one_line == true || n_rows == n_columns || in != out,
            ExcMessage("In-place operation only supported for "
                       "n_rows==n_columns or single-line interpolation"));
    AssertIndexRange (direction, dim);
    const int mm     = contract_over_rows ? n_rows : n_columns,
              nn     = contract_over_rows ? n_columns : n_rows;

    const int stride    = direction==0 ? 1 : Utilities::fixed_power<direction>(n_columns);
    const int n_blocks1 = one_line ? 1 : stride;
    const int n_blocks2 = direction >= dim-1 ? 1 : Utilities::fixed_power<dim-direction-1>(n_rows);
    Assert(n_rows <= 128, ExcNotImplemented());

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            Number x[129];
            for (int i=0; i<mm; ++i)
              x[i] = in[stride*i];
            for (int col=0; col<nn; ++col)
              {
                Number2 val0;
                if (contract_over_rows == true)
                  val0 = shape_data[col];
                else
                  val0 = shape_data[col*n_columns];
                Number res0 = val0 * x[0];
                for (int i=1; i<mm; ++i)
                  {
                    if (contract_over_rows == true)
                      val0 = shape_data[i*n_columns+col];
                    else
                      val0 = shape_data[col*n_columns+i];
                    res0 += val0 * x[i];
                  }
                if (add == false)
                  out[stride*col]  = res0;
                else
                  out[stride*col] += res0;
              }

            if (one_line == false)
              {
                ++in;
                ++out;
              }
          }
        if (one_line == false)
          {
            in += stride*(mm-1);
            out += stride*(nn-1);
          }
      }
  }



  template <int dim, typename Number, typename Number2>
  template <int face_direction, bool contract_onto_face, bool add, int max_derivative>
  inline
  void
  EvaluatorTensorProduct<evaluate_general,dim,0,0,Number,Number2>
  ::apply_face (const Number *DEAL_II_RESTRICT in,
                Number       *DEAL_II_RESTRICT out) const
  {
    Assert(shape_values != nullptr,
           ExcMessage("The given array shape_data must not be the null pointer!"));
    static_assert(dim > 0 && dim<4, "Only dim=1,2,3 supported");
    const int n_blocks1 = dim > 1 ? n_rows : 1;
    const int n_blocks2 = dim > 2 ? n_rows : 1;

    AssertIndexRange (face_direction, dim);
    const int stride = face_direction > 0 ? Utilities::fixed_power<face_direction>(n_rows) : 1;
    const int out_stride = dim > 1 ? Utilities::fixed_power<dim-1>(n_rows) : 1;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            if (contract_onto_face == true)
              {
                Number res0 = shape_values[0] * in[0];
                Number res1, res2;
                if (max_derivative > 0)
                  res1 = shape_values[n_rows] * in[0];
                if (max_derivative > 1)
                  res2 = shape_values[2*n_rows] * in[0];
                for (unsigned int ind=1; ind<n_rows; ++ind)
                  {
                    res0 += shape_values[ind] * in[stride*ind];
                    if (max_derivative > 0)
                      res1 += shape_values[ind+n_rows] * in[stride*ind];
                    if (max_derivative > 1)
                      res2 += shape_values[ind+2*n_rows] * in[stride*ind];
                  }
                if (add == false)
                  {
                    out[0]              = res0;
                    if (max_derivative > 0)
                      out[out_stride]   = res1;
                    if (max_derivative > 1)
                      out[2*out_stride] = res2;
                  }
                else
                  {
                    out[0]              += res0;
                    if (max_derivative > 0)
                      out[out_stride]   += res1;
                    if (max_derivative > 1)
                      out[2*out_stride] += res2;
                  }
              }
            else
              {
                for (unsigned int col=0; col<n_rows; ++col)
                  {
                    if (add == false)
                      out[col*stride]  = shape_values[col] * in[0];
                    else
                      out[col*stride] += shape_values[col] * in[0];
                    if (max_derivative > 0)
                      out[col*stride] += shape_values[col+n_rows] * in[out_stride];
                    if (max_derivative > 1)
                      out[col*stride] += shape_values[col+2*n_rows] * in[2*out_stride];
                  }
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need
            // to jump over to the next layer in z-direction
            switch (face_direction)
              {
              case 0:
                in += contract_onto_face ? n_rows : 1;
                out += contract_onto_face ? 1 : n_rows;
                break;
              case 1:
                ++in;
                ++out;
                // faces 2 and 3 in 3D use local coordinate system zx, which
                // is the other way around compared to the tensor
                // product. Need to take that into account.
                if (dim == 3)
                  {
                    if (contract_onto_face)
                      out += n_rows-1;
                    else
                      in += n_rows-1;
                  }
                break;
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (face_direction == 1 && dim == 3)
          {
            // adjust for local coordinate system zx
            if (contract_onto_face)
              {
                in += n_rows*(n_rows-1);
                out -= n_rows*n_rows-1;
              }
            else
              {
                out += n_rows*(n_rows-1);
                in -= n_rows*n_rows-1;
              }
          }
      }
  }



  /**
   * Internal evaluator for 1d-3d shape function using the tensor product form
   * of the basis functions. This class specializes the general application of
   * tensor-product based elements for "symmetric" finite elements, i.e., when
   * the shape functions are symmetric about 0.5 and the quadrature points
   * are, too.
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
  template <int dim, int n_rows, int n_columns, typename Number, typename Number2>
  struct EvaluatorTensorProduct<evaluate_symmetric,dim,n_rows,n_columns,Number,Number2>
  {
    static const unsigned int n_rows_of_product = Utilities::fixed_int_power<n_rows,dim>::value;
    static const unsigned int n_columns_of_product = Utilities::fixed_int_power<n_columns,dim>::value;

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number2> &shape_values,
                            const AlignedVector<Number2> &shape_gradients,
                            const AlignedVector<Number2> &shape_hessians,
                            const unsigned int            dummy1 = 0,
                            const unsigned int            dummy2 = 0)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {
      Assert(shape_values.size() == 0 ||
             shape_values.size() == n_rows*n_columns,
             ExcDimensionMismatch(shape_values.size(), n_rows*n_columns));
      Assert(shape_gradients.size() == 0 ||
             shape_gradients.size() == n_rows*n_columns,
             ExcDimensionMismatch(shape_gradients.size(), n_rows*n_columns));
      Assert(shape_hessians.size() == 0 ||
             shape_hessians.size() == n_rows*n_columns,
             ExcDimensionMismatch(shape_hessians.size(), n_rows*n_columns));
      (void)dummy1;
      (void)dummy2;
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values (const Number in [],
            Number       out[]) const;

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const;

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const;

    const Number2 *shape_values;
    const Number2 *shape_gradients;
    const Number2 *shape_hessians;
  };



  // In this case, the 1D shape values read (sorted lexicographically, rows
  // run over 1D dofs, columns over quadrature points):
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
  // In these matrices, we want to use avoid computations involving zeros and
  // ones and in addition use the symmetry in entries to reduce the number of
  // read operations.
  template <int dim, int n_rows, int n_columns, typename Number, typename Number2>
  template <int direction, bool contract_over_rows, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric,dim,n_rows,n_columns,Number,Number2>
  ::values (const Number in [],
            Number       out []) const
  {
    Assert (shape_values != nullptr, ExcNotInitialized());
    AssertIndexRange (direction, dim);
    constexpr int mm     = contract_over_rows ? n_rows : n_columns,
                  nn     = contract_over_rows ? n_columns : n_rows;
    constexpr int n_cols = nn / 2;
    constexpr int mid    = mm / 2;

    constexpr int stride    = Utilities::fixed_int_power<n_columns,direction>::value;
    constexpr int n_blocks1 = stride;
    constexpr int n_blocks2 = Utilities::fixed_int_power<n_rows,(direction>=dim)?0:(dim-direction-1)>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<n_cols; ++col)
              {
                Number2 val0, val1;
                Number in0, in1, res0, res1;
                if (contract_over_rows == true)
                  {
                    val0 = shape_values[col];
                    val1 = shape_values[nn-1-col];
                  }
                else
                  {
                    val0 = shape_values[col*n_columns];
                    val1 = shape_values[(col+1)*n_columns-1];
                  }
                if (mid > 0)
                  {
                    in0 = in[0];
                    in1 = in[stride*(mm-1)];
                    res0 = val0 * in0;
                    res1 = val1 * in0;
                    res0 += val1 * in1;
                    res1 += val0 * in1;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (contract_over_rows == true)
                          {
                            val0 = shape_values[ind*n_columns+col];
                            val1 = shape_values[ind*n_columns+nn-1-col];
                          }
                        else
                          {
                            val0 = shape_values[col*n_columns+ind];
                            val1 = shape_values[(col+1)*n_columns-1-ind];
                          }
                        in0 = in[stride*ind];
                        in1 = in[stride*(mm-1-ind)];
                        res0 += val0 * in0;
                        res1 += val1 * in0;
                        res0 += val1 * in1;
                        res1 += val0 * in1;
                      }
                  }
                else
                  res0 = res1 = Number();
                if (contract_over_rows == true)
                  {
                    if (mm % 2 == 1)
                      {
                        val0 = shape_values[mid*n_columns+col];
                        in1 = val0 * in[stride*mid];
                        res0 += in1;
                        res1 += in1;
                      }
                  }
                else
                  {
                    if (mm % 2 == 1 && nn % 2 == 0)
                      {
                        val0 = shape_values[col*n_columns+mid];
                        in1 = val0 * in[stride*mid];
                        res0 += in1;
                        res1 += in1;
                      }
                  }
                if (add == false)
                  {
                    out[stride*col]         = res0;
                    out[stride*(nn-1-col)]  = res1;
                  }
                else
                  {
                    out[stride*col]        += res0;
                    out[stride*(nn-1-col)] += res1;
                  }
              }
            if ( contract_over_rows == true && nn%2==1 && mm%2==1 )
              {
                if (add==false)
                  out[stride*n_cols]  = in[stride*mid];
                else
                  out[stride*n_cols] += in[stride*mid];
              }
            else if (contract_over_rows == true && nn%2==1)
              {
                Number res0;
                Number2 val0  = shape_values[n_cols];
                if (mid > 0)
                  {
                    res0  = val0 * (in[0] + in[stride*(mm-1)]);
                    for (int ind=1; ind<mid; ++ind)
                      {
                        val0  = shape_values[ind*n_columns+n_cols];
                        res0 += val0 * (in[stride*ind] + in[stride*(mm-1-ind)]);
                      }
                  }
                else
                  res0 = Number();
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }
            else if (contract_over_rows == false && nn%2 == 1)
              {
                Number res0;
                if (mid > 0)
                  {
                    Number2 val0 = shape_values[n_cols*n_columns];
                    res0 = val0 * (in[0] + in[stride*(mm-1)]);
                    for (int ind=1; ind<mid; ++ind)
                      {
                        val0  = shape_values[n_cols*n_columns+ind];
                        Number in1 = val0 * (in[stride*ind] + in[stride*(mm-1-ind)]);
                        res0 += in1;
                      }
                    if (mm % 2)
                      res0 += in[stride*mid];
                  }
                else
                  res0 = in[0];
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }

            ++in;
            ++out;
          }
        in += stride*(mm-1);
        out += stride*(nn-1);
      }
  }



  // For the specialized loop used for the gradient computation in
  // here, the 1D shape values read (sorted lexicographically, rows
  // run over 1D dofs, columns over quadrature points):
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
  //
  // In these matrices, we want to use avoid computations involving
  // zeros and ones and in addition use the symmetry in entries to
  // reduce the number of read operations.
  template <int dim, int n_rows, int n_columns, typename Number, typename Number2>
  template <int direction, bool contract_over_rows, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric,dim,n_rows,n_columns,Number,Number2>
  ::gradients (const Number in [],
               Number       out []) const
  {
    Assert (shape_gradients != nullptr, ExcNotInitialized());
    AssertIndexRange (direction, dim);
    constexpr int mm     = contract_over_rows ? n_rows : n_columns,
                  nn     = contract_over_rows ? n_columns : n_rows;
    constexpr int n_cols = nn / 2;
    constexpr int mid    = mm / 2;

    constexpr int stride    = Utilities::fixed_int_power<n_columns,direction>::value;
    constexpr int n_blocks1 = stride;
    constexpr int n_blocks2 = Utilities::fixed_int_power<n_rows,(direction>=dim)?0:(dim-direction-1)>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<n_cols; ++col)
              {
                Number2 val0, val1;
                Number in0, in1, res0, res1;
                if (contract_over_rows == true)
                  {
                    val0 = shape_gradients[col];
                    val1 = shape_gradients[nn-1-col];
                  }
                else
                  {
                    val0 = shape_gradients[col*n_columns];
                    val1 = shape_gradients[(nn-col-1)*n_columns];
                  }
                if (mid > 0)
                  {
                    in0 = in[0];
                    in1 = in[stride*(mm-1)];
                    res0 = val0 * in0;
                    res1 = val1 * in0;
                    res0 -= val1 * in1;
                    res1 -= val0 * in1;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (contract_over_rows == true)
                          {
                            val0 = shape_gradients[ind*n_columns+col];
                            val1 = shape_gradients[ind*n_columns+nn-1-col];
                          }
                        else
                          {
                            val0 = shape_gradients[col*n_columns+ind];
                            val1 = shape_gradients[(nn-col-1)*n_columns+ind];
                          }
                        in0 = in[stride*ind];
                        in1 = in[stride*(mm-1-ind)];
                        res0 += val0 * in0;
                        res1 += val1 * in0;
                        res0 -= val1 * in1;
                        res1 -= val0 * in1;
                      }
                  }
                else
                  res0 = res1 = Number();
                if (mm % 2 == 1)
                  {
                    if (contract_over_rows == true)
                      val0 = shape_gradients[mid*n_columns+col];
                    else
                      val0 = shape_gradients[col*n_columns+mid];
                    in1 = val0 * in[stride*mid];
                    res0 += in1;
                    res1 -= in1;
                  }
                if (add == false)
                  {
                    out[stride*col]         = res0;
                    out[stride*(nn-1-col)]  = res1;
                  }
                else
                  {
                    out[stride*col]        += res0;
                    out[stride*(nn-1-col)] += res1;
                  }
              }
            if ( nn%2 == 1 )
              {
                Number2 val0;
                Number res0;
                if (contract_over_rows == true)
                  val0 = shape_gradients[n_cols];
                else
                  val0 = shape_gradients[n_cols*n_columns];
                res0  = val0 * (in[0] - in[stride*(mm-1)]);
                for (int ind=1; ind<mid; ++ind)
                  {
                    if (contract_over_rows == true)
                      val0 = shape_gradients[ind*n_columns+n_cols];
                    else
                      val0 = shape_gradients[n_cols*n_columns+ind];
                    Number in1  = val0 * (in[stride*ind] - in[stride*(mm-1-ind)]);
                    res0 += in1;
                  }
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }

            ++in;
            ++out;
          }
        in += stride*(mm-1);
        out += stride*(nn-1);
      }
  }



  // evaluates the given shape data in 1d-3d using the tensor product
  // form assuming the symmetries of unit cell shape hessians for
  // finite elements in FEEvaluation
  template <int dim, int n_rows, int n_columns, typename Number, typename Number2>
  template <int direction, bool contract_over_rows, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric,dim,n_rows,n_columns,Number,Number2>
  ::hessians (const Number in [],
              Number       out []) const
  {
    Assert (shape_hessians != nullptr, ExcNotInitialized());
    AssertIndexRange (direction, dim);
    constexpr int mm     = contract_over_rows ? n_rows : n_columns;
    constexpr int nn     = contract_over_rows ? n_columns : n_rows;
    constexpr int n_cols = nn / 2;
    constexpr int mid    = mm / 2;

    constexpr int stride    = Utilities::fixed_int_power<n_columns,direction>::value;
    constexpr int n_blocks1 = stride;
    constexpr int n_blocks2 = Utilities::fixed_int_power<n_rows,(direction>=dim)?0:(dim-direction-1)>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<n_cols; ++col)
              {
                Number2 val0, val1;
                Number in0, in1, res0, res1;
                if (contract_over_rows == true)
                  {
                    val0 = shape_hessians[col];
                    val1 = shape_hessians[nn-1-col];
                  }
                else
                  {
                    val0 = shape_hessians[col*n_columns];
                    val1 = shape_hessians[(col+1)*n_columns-1];
                  }
                if (mid > 0)
                  {
                    in0 = in[0];
                    in1 = in[stride*(mm-1)];
                    res0 = val0 * in0;
                    res1 = val1 * in0;
                    res0 += val1 * in1;
                    res1 += val0 * in1;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (contract_over_rows == true)
                          {
                            val0 = shape_hessians[ind*n_columns+col];
                            val1 = shape_hessians[ind*n_columns+nn-1-col];
                          }
                        else
                          {
                            val0 = shape_hessians[col*n_columns+ind];
                            val1 = shape_hessians[(col+1)*n_columns-1-ind];
                          }
                        in0 = in[stride*ind];
                        in1 = in[stride*(mm-1-ind)];
                        res0 += val0 * in0;
                        res1 += val1 * in0;
                        res0 += val1 * in1;
                        res1 += val0 * in1;
                      }
                  }
                else
                  res0 = res1 = Number();
                if (mm % 2 == 1)
                  {
                    if (contract_over_rows == true)
                      val0 = shape_hessians[mid*n_columns+col];
                    else
                      val0 = shape_hessians[col*n_columns+mid];
                    in1 = val0 * in[stride*mid];
                    res0 += in1;
                    res1 += in1;
                  }
                if (add == false)
                  {
                    out[stride*col]         = res0;
                    out[stride*(nn-1-col)]  = res1;
                  }
                else
                  {
                    out[stride*col]        += res0;
                    out[stride*(nn-1-col)] += res1;
                  }
              }
            if ( nn%2 == 1 )
              {
                Number2 val0;
                Number res0;
                if (contract_over_rows == true)
                  val0 = shape_hessians[n_cols];
                else
                  val0 = shape_hessians[n_cols*n_columns];
                if (mid > 0)
                  {
                    res0  = val0 * (in[0] + in[stride*(mm-1)]);
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (contract_over_rows == true)
                          val0 = shape_hessians[ind*n_columns+n_cols];
                        else
                          val0 = shape_hessians[n_cols*n_columns+ind];
                        Number in1  = val0*(in[stride*ind] + in[stride*(mm-1-ind)]);
                        res0 += in1;
                      }
                  }
                else
                  res0 = Number();
                if (mm % 2 == 1)
                  {
                    if (contract_over_rows == true)
                      val0 = shape_hessians[mid*n_columns+n_cols];
                    else
                      val0 = shape_hessians[n_cols*n_columns+mid];
                    res0 += val0 * in[stride*mid];
                  }
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }

            ++in;
            ++out;
          }
        in += stride*(mm-1);
        out += stride*(nn-1);
      }
  }



  /**
   * Internal evaluator for 1d-3d shape function using the tensor product form
   * of the basis functions.
   *
   * This class implements a different approach to the symmetric case for
   * values, gradients, and Hessians also treated with the above functions: It
   * is possible to reduce the cost per dimension from N^2 to N^2/2, where N
   * is the number of 1D dofs (there are only N^2/2 different entries in the
   * shape matrix, so this is plausible). The approach is based on the idea of
   * applying the operator on the even and odd part of the input vectors
   * separately, given that the shape functions evaluated on quadrature points
   * are symmetric. This method is presented e.g. in the book "Implementing
   * Spectral Methods for Partial Differential Equations" by David A. Kopriva,
   * Springer, 2009, section 3.5.3 (Even-Odd-Decomposition). Even though the
   * experiments in the book say that the method is not efficient for N<20, it
   * is more efficient in the context where the loop bounds are compile-time
   * constants (templates).
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
  template <int dim, int n_rows, int n_columns, typename Number, typename Number2>
  struct EvaluatorTensorProduct<evaluate_evenodd,dim,n_rows,n_columns,Number,Number2>
  {
    static const unsigned int n_rows_of_product = Utilities::fixed_int_power<n_rows,dim>::value;
    static const unsigned int n_columns_of_product = Utilities::fixed_int_power<n_columns,dim>::value;

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other
     * constructor passing in at least an array for the values.
     */
    EvaluatorTensorProduct ()
      :
      shape_values (nullptr),
      shape_gradients (nullptr),
      shape_hessians (nullptr)
    {}

    /**
     * Constructor, taking the data from ShapeInfo (using the even-odd
     * variants stored there)
     */
    EvaluatorTensorProduct (const AlignedVector<Number2> &shape_values)
      :
      shape_values (shape_values.begin()),
      shape_gradients (nullptr),
      shape_hessians (nullptr)
    {
      AssertDimension(shape_values.size(), n_rows*((n_columns+1)/2));
    }

    /**
     * Constructor, taking the data from ShapeInfo (using the even-odd
     * variants stored there)
     */
    EvaluatorTensorProduct (const AlignedVector<Number2> &shape_values,
                            const AlignedVector<Number2> &shape_gradients,
                            const AlignedVector<Number2> &shape_hessians,
                            const unsigned int            dummy1 = 0,
                            const unsigned int            dummy2 = 0)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {
      // In this function, we allow for dummy pointers if some of values,
      // gradients or hessians should not be computed
      if (!shape_values.empty())
        AssertDimension(shape_values.size(), n_rows*((n_columns+1)/2));
      if (!shape_gradients.empty())
        AssertDimension(shape_gradients.size(), n_rows*((n_columns+1)/2));
      if (!shape_hessians.empty())
        AssertDimension(shape_hessians.size(), n_rows*((n_columns+1)/2));
      (void)dummy1;
      (void)dummy2;
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values (const Number in [],
            Number       out[]) const
    {
      Assert (shape_values != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,0>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const
    {
      Assert (shape_gradients != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,1>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const
    {
      Assert (shape_hessians != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,2>(shape_hessians, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values_one_line (const Number in [],
                     Number       out[]) const
    {
      Assert (shape_values != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,0,true>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients_one_line (const Number in [],
                        Number       out[]) const
    {
      Assert (shape_gradients != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,1,true>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians_one_line (const Number in [],
                       Number       out[]) const
    {
      Assert (shape_hessians != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,2,true>(shape_hessians, in, out);
    }

    /**
     * This function applies the tensor product kernel, corresponding to a
     * multiplication of 1D stripes, along the given @p direction of the tensor
     * data in the input array. This function allows the @p in and @p out
     * arrays to alias for the case n_rows == n_columns, i.e., it is safe to
     * perform the contraction in place where @p in and @p out point to the
     * same address. For the case n_rows != n_columns, the output is only
     * correct if @p one_line is set to true.
     *
     * @tparam direction Direction that is evaluated
     * @tparam contract_over_rows If true, the tensor contraction sums
     *                            over the rows in the given @p shape_data
     *                            array, otherwise it sums over the columns
     * @tparam add If true, the result is added to the output vector, else
     *             the computed values overwrite the content in the output
     * @tparam type Determines whether to use the symmetries appearing in
     *              shape values (type=0), shape gradients (type=1) or
     *              second derivatives (type=2, similar to type 0 but
     *              without two additional zero entries)
     * @tparam one_line If true, the kernel is only applied along a single 1D
     *                  stripe within a dim-dimensional tensor, not the full
     *                  n_rows^dim points as in the @p false case.
     *
     * @param shape_data Transformation matrix with @p n_rows rows and
     *                   @p n_columns columns, stored in row-major format
     * @param in Pointer to the start of the input data vector
     * @param out Pointer to the start of the output data vector
     */
    template <int direction, bool contract_over_rows, bool add, int type,
              bool one_line=false>
    static void apply (const Number2 *DEAL_II_RESTRICT shape_data,
                       const Number  *in,
                       Number        *out);

    const Number2 *shape_values;
    const Number2 *shape_gradients;
    const Number2 *shape_hessians;
  };



  template <int dim, int n_rows, int n_columns, typename Number, typename Number2>
  template <int direction, bool contract_over_rows, bool add, int type, bool one_line>
  inline
  void
  EvaluatorTensorProduct<evaluate_evenodd,dim,n_rows,n_columns,Number,Number2>
  ::apply (const Number2 *DEAL_II_RESTRICT shapes,
           const Number  *in,
           Number        *out)
  {
    static_assert (type < 3, "Only three variants type=0,1,2 implemented");
    static_assert (one_line == false || direction==dim-1,
                   "Single-line evaluation only works for direction=dim-1.");
    Assert (dim == direction+1 || one_line == true || n_rows == n_columns || in != out,
            ExcMessage("In-place operation only supported for "
                       "n_rows==n_columns or single-line interpolation"));

    // We cannot statically assert that direction is less than dim, so must do
    // an additional dynamic check
    AssertIndexRange (direction, dim);

    constexpr int nn = contract_over_rows ? n_columns : n_rows;
    constexpr int mm = contract_over_rows ? n_rows : n_columns;
    constexpr int n_cols = nn / 2;
    constexpr int mid    = mm / 2;

    constexpr int stride    = Utilities::fixed_int_power<n_columns,direction>::value;
    constexpr int n_blocks1 = one_line ? 1 : stride;
    constexpr int n_blocks2 = Utilities::fixed_int_power<n_rows,(direction>=dim)?0:(dim-direction-1)>::value;

    constexpr int offset = (n_columns+1)/2;

    // this code may look very inefficient at first sight due to the many
    // different cases with if's at the innermost loop part, but all of the
    // conditionals can be evaluated at compile time because they are
    // templates, so the compiler should optimize everything away
    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            Number xp[mid>0?mid:1], xm[mid>0?mid:1];
            for (int i=0; i<mid; ++i)
              {
                if (contract_over_rows == true && type == 1)
                  {
                    xp[i] = in[stride*i] - in[stride*(mm-1-i)];
                    xm[i] = in[stride*i] + in[stride*(mm-1-i)];
                  }
                else
                  {
                    xp[i] = in[stride*i] + in[stride*(mm-1-i)];
                    xm[i] = in[stride*i] - in[stride*(mm-1-i)];
                  }
              }
            Number xmid = in[stride*mid];
            for (int col=0; col<n_cols; ++col)
              {
                Number r0, r1;
                if (mid > 0)
                  {
                    if (contract_over_rows == true)
                      {
                        r0 = shapes[col]                     * xp[0];
                        r1 = shapes[(n_rows-1)*offset + col] * xm[0];
                      }
                    else
                      {
                        r0 = shapes[col*offset]            * xp[0];
                        r1 = shapes[(n_rows-1-col)*offset] * xm[0];
                      }
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (contract_over_rows == true)
                          {
                            r0 += shapes[ind*offset+col]            * xp[ind];
                            r1 += shapes[(n_rows-1-ind)*offset+col] * xm[ind];
                          }
                        else
                          {
                            r0 += shapes[col*offset+ind]            * xp[ind];
                            r1 += shapes[(n_rows-1-col)*offset+ind] * xm[ind];
                          }
                      }
                  }
                else
                  r0 = r1 = Number();
                if (mm % 2 == 1 && contract_over_rows == true)
                  {
                    if (type == 1)
                      r1 += shapes[mid*offset+col] * xmid;
                    else
                      r0 += shapes[mid*offset+col] * xmid;
                  }
                else if (mm % 2 == 1 && (nn % 2 == 0 || type > 0))
                  r0 += shapes[col*offset+mid] * xmid;

                if (add == false)
                  {
                    out[stride*col]         = r0 + r1;
                    if (type == 1 && contract_over_rows == false)
                      out[stride*(nn-1-col)]  = r1 - r0;
                    else
                      out[stride*(nn-1-col)]  = r0 - r1;
                  }
                else
                  {
                    out[stride*col]        += r0 + r1;
                    if (type == 1 && contract_over_rows == false)
                      out[stride*(nn-1-col)] += r1 - r0;
                    else
                      out[stride*(nn-1-col)] += r0 - r1;
                  }
              }
            if ( type == 0 && contract_over_rows == true && nn%2==1 && mm%2==1 )
              {
                if (add==false)
                  out[stride*n_cols]  = xmid;
                else
                  out[stride*n_cols] += xmid;
              }
            else if (contract_over_rows == true && nn%2==1)
              {
                Number r0;
                if (mid > 0)
                  {
                    r0  = shapes[n_cols] * xp[0];
                    for (int ind=1; ind<mid; ++ind)
                      r0 += shapes[ind*offset+n_cols] * xp[ind];
                  }
                else
                  r0 = Number();
                if (type != 1 && mm % 2 == 1)
                  r0 += shapes[mid*offset+n_cols] * xmid;

                if (add == false)
                  out[stride*n_cols]  = r0;
                else
                  out[stride*n_cols] += r0;
              }
            else if (contract_over_rows == false && nn%2 == 1)
              {
                Number r0;
                if (mid > 0)
                  {
                    if (type == 1)
                      {
                        r0 = shapes[n_cols*offset] * xm[0];
                        for (int ind=1; ind<mid; ++ind)
                          r0 += shapes[n_cols*offset+ind] * xm[ind];
                      }
                    else
                      {
                        r0 = shapes[n_cols*offset] * xp[0];
                        for (int ind=1; ind<mid; ++ind)
                          r0 += shapes[n_cols*offset+ind] * xp[ind];
                      }
                  }
                else
                  r0 = Number();

                if (type == 0 && mm % 2 == 1)
                  r0 += xmid;
                else if (type == 2 && mm % 2 == 1)
                  r0 += shapes[n_cols*offset+mid] * xmid;

                if (add == false)
                  out[stride*n_cols]  = r0;
                else
                  out[stride*n_cols] += r0;
              }
            if (one_line == false)
              {
                in  += 1;
                out += 1;
              }
          }
        if (one_line == false)
          {
            in += stride * (mm-1);
            out += stride * (nn-1);
          }
      }
  }



  /**
   * Internal evaluator for 1d-3d shape function using the tensor product form
   * of the basis functions.
   *
   * This class implements an approach similar to the even-odd decomposition
   * but with a different type of symmetry. In this case, we assume that a
   * single shape function already shows the symmetry over the quadrature
   * points, rather than the complete basis that is considered in the even-odd
   * case. In particular, we assume that the shape functions are ordered as in
   * the Legendre basis, with symmetric shape functions in the even slots
   * (rows of the values array) and point-symmetric in the odd slots. Like the
   * even-odd decomposition, the number of operations are N^2/2 rather than
   * N^2 FMAs (fused multiply-add), where N is the number of 1D dofs. The
   * difference is in the way the input and output quantities are symmetrized.
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
  template <int dim, int n_rows, int n_columns, typename Number, typename Number2>
  struct EvaluatorTensorProduct<evaluate_symmetric_hierarchical,dim,n_rows,n_columns,Number,Number2>
  {
    static const unsigned int n_rows_of_product = Utilities::fixed_int_power<n_rows,dim>::value;
    static const unsigned int n_columns_of_product = Utilities::fixed_int_power<n_columns,dim>::value;

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other
     * constructor passing in at least an array for the values.
     */
    EvaluatorTensorProduct ()
      :
      shape_values (nullptr),
      shape_gradients (nullptr),
      shape_hessians (nullptr)
    {}

    /**
     * Constructor, taking the data from ShapeInfo (using the even-odd
     * variants stored there)
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values)
      :
      shape_values (shape_values.begin()),
      shape_gradients (nullptr),
      shape_hessians (nullptr)
    {}

    /**
     * Constructor, taking the data from ShapeInfo (using the even-odd
     * variants stored there)
     */
    EvaluatorTensorProduct (const AlignedVector<Number2> &shape_values,
                            const AlignedVector<Number2> &shape_gradients,
                            const AlignedVector<Number2> &shape_hessians,
                            const unsigned int            dummy1 = 0,
                            const unsigned int            dummy2 = 0)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {
      (void)dummy1;
      (void)dummy2;
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values (const Number in [],
            Number       out[]) const
    {
      Assert (shape_values != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,0>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const
    {
      Assert (shape_gradients != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,1>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const
    {
      Assert (shape_hessians != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,0>(shape_hessians, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values_one_line (const Number in [],
                     Number       out[]) const
    {
      Assert (shape_values != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,0,true>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients_one_line (const Number in [],
                        Number       out[]) const
    {
      Assert (shape_gradients != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,1,true>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians_one_line (const Number in [],
                       Number       out[]) const
    {
      Assert (shape_hessians != nullptr, ExcNotInitialized());
      apply<direction,contract_over_rows,add,0,true>(shape_hessians, in, out);
    }

    /**
     * This function applies the tensor product kernel, corresponding to a
     * multiplication of 1D stripes, along the given @p direction of the tensor
     * data in the input array. This function allows the @p in and @p out
     * arrays to alias for the case n_rows == n_columns, i.e., it is safe to
     * perform the contraction in place where @p in and @p out point to the
     * same address. For the case n_rows != n_columns, the output is only
     * correct if @p one_line is set to true.
     *
     * @tparam direction Direction that is evaluated
     * @tparam contract_over_rows If true, the tensor contraction sums
     *                            over the rows in the given @p shape_data
     *                            array, otherwise it sums over the columns
     * @tparam add If true, the result is added to the output vector, else
     *             the computed values overwrite the content in the output
     * @tparam type Determines whether the evaluation is symmetric in even
     *              rows (type=0) or odd rows (type=1) of @p shape_data and
     *              skew-symmetric in odd rows (type=0) or even rows (type=1)
     * @tparam one_line If true, the kernel is only applied along a single 1D
     *                  stripe within a dim-dimensional tensor, not the full
     *                  n_rows^dim points as in the @p false case.
     *
     * @param shape_data Transformation matrix with @p n_rows rows and
     *                   @p n_columns columns, stored in row-major format
     * @param in Pointer to the start of the input data vector
     * @param out Pointer to the start of the output data vector
     */
    template <int direction, bool contract_over_rows, bool add, int type,
              bool one_line=false>
    static void apply (const Number2 *DEAL_II_RESTRICT shape_data,
                       const Number  *in,
                       Number        *out);

    const Number2 *shape_values;
    const Number2 *shape_gradients;
    const Number2 *shape_hessians;
  };



  template <int dim, int n_rows, int n_columns, typename Number, typename Number2>
  template <int direction, bool contract_over_rows, bool add, int type, bool one_line>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric_hierarchical,dim,n_rows,n_columns,Number,Number2>
  ::apply (const Number2 *DEAL_II_RESTRICT shapes,
           const Number  *in,
           Number        *out)
  {
    static_assert (one_line == false || direction==dim-1,
                   "Single-line evaluation only works for direction=dim-1.");
    static_assert (type == 0 || type == 1,
                   "Only types 0 and 1 implemented for evaluate_symmetric_hierarchical.");
    Assert (dim == direction+1 || one_line == true || n_rows == n_columns || in != out,
            ExcMessage("In-place operation only supported for "
                       "n_rows==n_columns or single-line interpolation"));

    // We cannot statically assert that direction is less than dim, so must do
    // an additional dynamic check
    AssertIndexRange (direction, dim);

    constexpr int nn = contract_over_rows ? n_columns : n_rows;
    constexpr int mm = contract_over_rows ? n_rows : n_columns;
    constexpr int n_cols = nn / 2;
    constexpr int mid    = mm / 2;

    constexpr int stride    = Utilities::fixed_int_power<n_columns,direction>::value;
    constexpr int n_blocks1 = one_line ? 1 : stride;
    constexpr int n_blocks2 = Utilities::fixed_int_power<n_rows,(direction>=dim)?0:(dim-direction-1)>::value;

    // this code may look very inefficient at first sight due to the many
    // different cases with if's at the innermost loop part, but all of the
    // conditionals can be evaluated at compile time because they are
    // templates, so the compiler should optimize everything away
    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            if (contract_over_rows)
              {
                Number x[mm];
                for (unsigned int i=0; i<mm; ++i)
                  x[i] = in[stride*i];
                for (unsigned int col=0; col<n_cols; ++col)
                  {
                    Number r0, r1;
                    if (mid > 0)
                      {
                        r0 = shapes[col] * x[0];
                        r1 = shapes[col+n_columns] * x[1];
                        for (unsigned int ind=1; ind<mid; ++ind)
                          {
                            r0 += shapes[col+2*ind*n_columns] * x[2*ind];
                            r1 += shapes[col+(2*ind+1)*n_columns] * x[2*ind+1];
                          }
                      }
                    else
                      r0 = r1 = Number();
                    if (mm%2 == 1)
                      r0 += shapes[col+(mm-1)*n_columns] * x[mm-1];
                    if (add == false)
                      {
                        out[stride*col]           = r0 + r1;
                        if (type == 1)
                          out[stride*(nn-1-col)]  = r1 - r0;
                        else
                          out[stride*(nn-1-col)]  = r0 - r1;
                      }
                    else
                      {
                        out[stride*col]          += r0 + r1;
                        if (type == 1)
                          out[stride*(nn-1-col)] += r1 - r0;
                        else
                          out[stride*(nn-1-col)] += r0 - r1;
                      }
                  }
                if (nn%2 == 1)
                  {
                    Number r0;
                    const unsigned int shift = type==1 ? 1 : 0;
                    if (mid>0)
                      {
                        r0 = shapes[n_cols + shift*n_columns] * x[shift];
                        for (unsigned int ind=1; ind<mid; ++ind)
                          r0 += shapes[n_cols + (2*ind+shift)*n_columns] * x[2*ind+shift];
                      }
                    else
                      r0 = 0;
                    if (type != 1 && mm%2 == 1)
                      r0 += shapes[n_cols + (mm-1)*n_columns] * x[mm-1];
                    if (add == false)
                      out[stride*n_cols]  = r0;
                    else
                      out[stride*n_cols] += r0;
                  }
              }
            else
              {
                Number xp[mid+1], xm[mid>0?mid:1];
                for (int i=0; i<mid; ++i)
                  if (type == 0)
                    {
                      xp[i] = in[stride*i] + in[stride*(mm-1-i)];
                      xm[i] = in[stride*i] - in[stride*(mm-1-i)];
                    }
                  else
                    {
                      xp[i] = in[stride*i] - in[stride*(mm-1-i)];
                      xm[i] = in[stride*i] + in[stride*(mm-1-i)];
                    }
                if (mm%2 == 1)
                  xp[mid] = in[stride*mid];
                for (unsigned int col=0; col<n_cols; ++col)
                  {
                    Number r0, r1;
                    if (mid > 0)
                      {
                        r0 = shapes[2*col*n_columns] * xp[0];
                        r1 = shapes[(2*col+1)*n_columns] * xm[0];
                        for (unsigned int ind=1; ind<mid; ++ind)
                          {
                            r0 += shapes[2*col*n_columns+ind] * xp[ind];
                            r1 += shapes[(2*col+1)*n_columns+ind] * xm[ind];
                          }
                      }
                    else
                      r0 = r1 = Number();
                    if (mm%2 == 1)
                      {
                        if (type == 1)
                          r1 += shapes[(2*col+1)*n_columns+mid] * xp[mid];
                        else
                          r0 += shapes[2*col*n_columns+mid] * xp[mid];
                      }
                    if (add == false)
                      {
                        out[stride*(2*col)]    = r0;
                        out[stride*(2*col+1)]  = r1;
                      }
                    else
                      {
                        out[stride*(2*col)]   += r0;
                        out[stride*(2*col+1)] += r1;
                      }
                  }
                if (nn%2 == 1)
                  {
                    Number r0;
                    if (mid > 0)
                      {
                        r0 = shapes[(nn-1)*n_columns] * xp[0];
                        for (unsigned int ind=1; ind<mid; ++ind)
                          r0 += shapes[(nn-1)*n_columns+ind] * xp[ind];
                      }
                    else
                      r0 = Number();
                    if (mm%2 == 1 && type == 0)
                      r0 += shapes[(nn-1)*n_columns+mid] * xp[mid];
                    if (add == false)
                      out[stride*(nn-1)]      = r0;
                    else
                      out[stride*(nn-1)]     += r0;
                  }
              }
            if (one_line == false)
              {
                in  += 1;
                out += 1;
              }
          }
        if (one_line == false)
          {
            in += stride * (mm-1);
            out += stride * (nn-1);
          }
      }
  }

} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
