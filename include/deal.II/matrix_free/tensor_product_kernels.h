// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
    evaluate_evenodd
  };

  /**
   * Generic evaluator framework
   */
  template <EvaluatorVariant variant, int dim, int fe_degree, int n_q_points_1d,
            typename Number>
  struct EvaluatorTensorProduct
  {};

  /**
   * Internal evaluator for 1d-3d shape function using the tensor product form
   * of the basis functions
   */
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct EvaluatorTensorProduct<evaluate_general,dim,fe_degree,n_q_points_1d,Number>
  {
    static const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    static const unsigned int n_q_points = Utilities::fixed_int_power<n_q_points_1d,dim>::value;

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other pointer
     */
    EvaluatorTensorProduct ()
      :
      shape_values (0),
      shape_gradients (0),
      shape_hessians (0)
    {}

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const AlignedVector<Number> &shape_hessians,
                            const unsigned int           dummy1 = 0,
                            const unsigned int           dummy2 = 0)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {
      (void)dummy1;
      (void)dummy2;
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number in [],
            Number       out[]) const
    {
      apply<direction,dof_to_quad,add>(shape_values, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const
    {
      apply<direction,dof_to_quad,add>(shape_gradients, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const
    {
      apply<direction,dof_to_quad,add>(shape_hessians, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    static void apply (const Number *shape_data,
                       const Number in [],
                       Number       out []);

    const Number *shape_values;
    const Number *shape_gradients;
    const Number *shape_hessians;
  };

  // evaluates the given shape data in 1d-3d using the tensor product
  // form. does not use a particular layout of entries in the matrices
  // like the functions below and corresponds to a usual matrix-matrix
  // product
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_general,dim,fe_degree,n_q_points_1d,Number>
  ::apply (const Number *shape_data,
           const Number in [],
           Number       out [])
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<nn; ++col)
              {
                Number val0;
                if (dof_to_quad == true)
                  val0 = shape_data[col];
                else
                  val0 = shape_data[col*n_q_points_1d];
                Number res0 = val0 * in[0];
                for (int ind=1; ind<mm; ++ind)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_data[ind*n_q_points_1d+col];
                    else
                      val0 = shape_data[col*n_q_points_1d+ind];
                    res0 += val0 * in[stride*ind];
                  }
                if (add == false)
                  out[stride*col]  = res0;
                else
                  out[stride*col] += res0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need
            // to jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
      }
  }



  // This method applies the tensor product operation to produce face values
  // out from cell values. As opposed to the apply_tensor_product method, this
  // method assumes that the directions orthogonal to the face have
  // fe_degree+1 degrees of freedom per direction and not n_q_points_1d for
  // those directions lower than the one currently applied
  template <int dim, int fe_degree, typename Number, int face_direction,
            bool dof_to_quad, bool add>
  inline
  void
  apply_tensor_product_face (const Number *shape_data,
                             const Number in [],
                             Number       out [])
  {
    const int n_blocks1 = dim > 1 ? (fe_degree+1) : 1;
    const int n_blocks2 = dim > 2 ? (fe_degree+1) : 1;

    AssertIndexRange (face_direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : 1,
              nn     = dof_to_quad ? 1 : (fe_degree+1);

    const int stride = Utilities::fixed_int_power<fe_degree+1,face_direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            if (dof_to_quad == true)
              {
                Number res0 = shape_data[0] * in[0];
                for (int ind=1; ind<mm; ++ind)
                  res0 += shape_data[ind] * in[stride*ind];
                if (add == false)
                  out[0]  = res0;
                else
                  out[0] += res0;
              }
            else
              {
                for (int col=0; col<nn; ++col)
                  if (add == false)
                    out[col*stride]  = shape_data[col] * in[0];
                  else
                    out[col*stride] += shape_data[col] * in[0];
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need
            // to jump over to the next layer in z-direction
            switch (face_direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
                ++in;
                ++out;
                // faces 2 and 3 in 3D use local coordinate system zx, which
                // is the other way around compared to the tensor
                // product. Need to take that into account.
                if (dim == 3)
                  {
                    if (dof_to_quad)
                      out += fe_degree;
                    else
                      in += fe_degree;
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
            in += mm*(mm-1);
            out += nn*(nn-1);
            // adjust for local coordinate system zx
            if (dof_to_quad)
              out -= (fe_degree+1)*(fe_degree+1)-1;
            else
              in -= (fe_degree+1)*(fe_degree+1)-1;
          }
      }
  }



  /**
   * Internal evaluator for 1d-3d shape function using the tensor product form
   * of the basis functions. The same as above but without making use of
   * template arguments and rather variable loop bounds.
   */
  template <int dim, typename Number>
  struct EvaluatorTensorProduct<evaluate_general,dim,-1,0,Number>
  {
    static const unsigned int dofs_per_cell = numbers::invalid_unsigned_int;
    static const unsigned int n_q_points = numbers::invalid_unsigned_int;

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other constructor
     */
    EvaluatorTensorProduct ()
      :
      shape_values (0),
      shape_gradients (0),
      shape_hessians (0),
      fe_degree (numbers::invalid_unsigned_int),
      n_q_points_1d (numbers::invalid_unsigned_int)
    {}

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const AlignedVector<Number> &shape_hessians,
                            const unsigned int           fe_degree,
                            const unsigned int           n_q_points_1d)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin()),
      fe_degree (fe_degree),
      n_q_points_1d (n_q_points_1d)
    {}

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number *in,
            Number       *out) const
    {
      apply<direction,dof_to_quad,add>(shape_values, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number *in,
               Number       *out) const
    {
      apply<direction,dof_to_quad,add>(shape_gradients, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    hessians (const Number *in,
              Number       *out) const
    {
      apply<direction,dof_to_quad,add>(shape_hessians, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void apply (const Number *shape_data,
                const Number *in,
                Number       *out) const;

    const Number *shape_values;
    const Number *shape_gradients;
    const Number *shape_hessians;
    const unsigned int fe_degree;
    const unsigned int n_q_points_1d;
  };

  // evaluates the given shape data in 1d-3d using the tensor product
  // form. does not use a particular layout of entries in the matrices
  // like the functions below and corresponds to a usual matrix-matrix
  // product
  template <int dim, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_general,dim,-1,0,Number>
  ::apply (const Number *shape_data,
           const Number *in,
           Number       *out) const
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = direction==0 ? 1 : Utilities::fixed_power<direction>(nn);

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<nn; ++col)
              {
                Number val0;
                if (dof_to_quad == true)
                  val0 = shape_data[col];
                else
                  val0 = shape_data[col*n_q_points_1d];
                Number res0 = val0 * in[0];
                for (int ind=1; ind<mm; ++ind)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_data[ind*n_q_points_1d+col];
                    else
                      val0 = shape_data[col*n_q_points_1d+ind];
                    res0 += val0 * in[stride*ind];
                  }
                if (add == false)
                  out[stride*col]  = res0;
                else
                  out[stride*col] += res0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need
            // to jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
      }
  }



  /**
   * Internal evaluator for 1d-3d shape function using the tensor product form
   * of the basis functions. This class specializes the general application of
   * tensor-product based elements for "symmetric" finite elements, i.e., when
   * the shape functions are symmetric about 0.5 and the quadrature points
   * are, too.
   */
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct EvaluatorTensorProduct<evaluate_symmetric,dim,fe_degree,n_q_points_1d,Number>
  {
    static const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    static const unsigned int n_q_points = Utilities::fixed_int_power<n_q_points_1d,dim>::value;

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const AlignedVector<Number> &shape_hessians,
                            const unsigned int           dummy1 = 0,
                            const unsigned int           dummy2 = 0)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {
      (void)dummy1;
      (void)dummy2;
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number in [],
            Number       out[]) const;

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const;

    template <int direction, bool dof_to_quad, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const;

    const Number *shape_values;
    const Number *shape_gradients;
    const Number *shape_hessians;
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
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric,dim,fe_degree,n_q_points_1d,Number>
  ::values (const Number in [],
            Number       out []) const
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
    const int n_cols = nn / 2;
    const int mid    = mm / 2;

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<n_cols; ++col)
              {
                Number val0, val1, in0, in1, res0, res1;
                if (dof_to_quad == true)
                  {
                    val0 = shape_values[col];
                    val1 = shape_values[nn-1-col];
                  }
                else
                  {
                    val0 = shape_values[col*n_q_points_1d];
                    val1 = shape_values[(col+1)*n_q_points_1d-1];
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
                        if (dof_to_quad == true)
                          {
                            val0 = shape_values[ind*n_q_points_1d+col];
                            val1 = shape_values[ind*n_q_points_1d+nn-1-col];
                          }
                        else
                          {
                            val0 = shape_values[col*n_q_points_1d+ind];
                            val1 = shape_values[(col+1)*n_q_points_1d-1-ind];
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
                if (dof_to_quad == true)
                  {
                    if (mm % 2 == 1)
                      {
                        val0 = shape_values[mid*n_q_points_1d+col];
                        val1 = val0 * in[stride*mid];
                        res0 += val1;
                        res1 += val1;
                      }
                  }
                else
                  {
                    if (mm % 2 == 1 && nn % 2 == 0)
                      {
                        val0 = shape_values[col*n_q_points_1d+mid];
                        val1 = val0 * in[stride*mid];
                        res0 += val1;
                        res1 += val1;
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
            if ( dof_to_quad == true && nn%2==1 && mm%2==1 )
              {
                if (add==false)
                  out[stride*n_cols]  = in[stride*mid];
                else
                  out[stride*n_cols] += in[stride*mid];
              }
            else if (dof_to_quad == true && nn%2==1)
              {
                Number res0;
                Number val0  = shape_values[n_cols];
                if (mid > 0)
                  {
                    res0  = in[0] + in[stride*(mm-1)];
                    res0 *= val0;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        val0  = shape_values[ind*n_q_points_1d+n_cols];
                        Number val1  = in[stride*ind] + in[stride*(mm-1-ind)];
                        val1 *= val0;
                        res0 += val1;
                      }
                  }
                else
                  res0 = Number();
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }
            else if (dof_to_quad == false && nn%2 == 1)
              {
                Number res0;
                if (mid > 0)
                  {
                    Number val0 = shape_values[n_cols*n_q_points_1d];
                    res0 = in[0] + in[stride*(mm-1)];
                    res0 *= val0;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        val0  = shape_values[n_cols*n_q_points_1d+ind];
                        Number val1 = in[stride*ind] + in[stride*(mm-1-ind)];
                        val1 *= val0;
                        res0 += val1;
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

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need to
            // jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
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
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric,dim,fe_degree,n_q_points_1d,Number>
  ::gradients (const Number in [],
               Number       out []) const
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
    const int n_cols = nn / 2;
    const int mid    = mm / 2;

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<n_cols; ++col)
              {
                Number val0, val1, in0, in1, res0, res1;
                if (dof_to_quad == true)
                  {
                    val0 = shape_gradients[col];
                    val1 = shape_gradients[nn-1-col];
                  }
                else
                  {
                    val0 = shape_gradients[col*n_q_points_1d];
                    val1 = shape_gradients[(nn-col-1)*n_q_points_1d];
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
                        if (dof_to_quad == true)
                          {
                            val0 = shape_gradients[ind*n_q_points_1d+col];
                            val1 = shape_gradients[ind*n_q_points_1d+nn-1-col];
                          }
                        else
                          {
                            val0 = shape_gradients[col*n_q_points_1d+ind];
                            val1 = shape_gradients[(nn-col-1)*n_q_points_1d+ind];
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
                    if (dof_to_quad == true)
                      val0 = shape_gradients[mid*n_q_points_1d+col];
                    else
                      val0 = shape_gradients[col*n_q_points_1d+mid];
                    val1 = val0 * in[stride*mid];
                    res0 += val1;
                    res1 -= val1;
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
                Number val0, res0;
                if (dof_to_quad == true)
                  val0 = shape_gradients[n_cols];
                else
                  val0 = shape_gradients[n_cols*n_q_points_1d];
                res0  = in[0] - in[stride*(mm-1)];
                res0 *= val0;
                for (int ind=1; ind<mid; ++ind)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_gradients[ind*n_q_points_1d+n_cols];
                    else
                      val0 = shape_gradients[n_cols*n_q_points_1d+ind];
                    Number val1  = in[stride*ind] - in[stride*(mm-1-ind)];
                    val1 *= val0;
                    res0 += val1;
                  }
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. for y-part in 3D and if we are at the end of one
            // chunk in x-dir, need to jump over to the next layer in
            // z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }

        if (direction == 1)
          {
            in  += nn * (mm-1);
            out += nn * (nn-1);
          }
      }
  }



  // evaluates the given shape data in 1d-3d using the tensor product
  // form assuming the symmetries of unit cell shape hessians for
  // finite elements in FEEvaluation
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric,dim,fe_degree,n_q_points_1d,Number>
  ::hessians (const Number in [],
              Number       out []) const
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
    const int n_cols = nn / 2;
    const int mid    = mm / 2;

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<n_cols; ++col)
              {
                Number val0, val1, in0, in1, res0, res1;
                if (dof_to_quad == true)
                  {
                    val0 = shape_hessians[col];
                    val1 = shape_hessians[nn-1-col];
                  }
                else
                  {
                    val0 = shape_hessians[col*n_q_points_1d];
                    val1 = shape_hessians[(col+1)*n_q_points_1d-1];
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
                        if (dof_to_quad == true)
                          {
                            val0 = shape_hessians[ind*n_q_points_1d+col];
                            val1 = shape_hessians[ind*n_q_points_1d+nn-1-col];
                          }
                        else
                          {
                            val0 = shape_hessians[col*n_q_points_1d+ind];
                            val1 = shape_hessians[(col+1)*n_q_points_1d-1-ind];
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
                    if (dof_to_quad == true)
                      val0 = shape_hessians[mid*n_q_points_1d+col];
                    else
                      val0 = shape_hessians[col*n_q_points_1d+mid];
                    val1 = val0 * in[stride*mid];
                    res0 += val1;
                    res1 += val1;
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
                Number val0, res0;
                if (dof_to_quad == true)
                  val0 = shape_hessians[n_cols];
                else
                  val0 = shape_hessians[n_cols*n_q_points_1d];
                if (mid > 0)
                  {
                    res0  = in[0] + in[stride*(mm-1)];
                    res0 *= val0;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (dof_to_quad == true)
                          val0 = shape_hessians[ind*n_q_points_1d+n_cols];
                        else
                          val0 = shape_hessians[n_cols*n_q_points_1d+ind];
                        Number val1  = in[stride*ind] + in[stride*(mm-1-ind)];
                        val1 *= val0;
                        res0 += val1;
                      }
                  }
                else
                  res0 = Number();
                if (mm % 2 == 1)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_hessians[mid*n_q_points_1d+n_cols];
                    else
                      val0 = shape_hessians[n_cols*n_q_points_1d+mid];
                    res0 += val0 * in[stride*mid];
                  }
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need to
            // jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
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
   */
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct EvaluatorTensorProduct<evaluate_evenodd,dim,fe_degree,n_q_points_1d,Number>
  {
    static const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    static const unsigned int n_q_points = Utilities::fixed_int_power<n_q_points_1d,dim>::value;

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other pointer
     */
    EvaluatorTensorProduct ()
      :
      shape_values (0),
      shape_gradients (0),
      shape_hessians (0)
    {}

    /**
     * Constructor, taking the data from ShapeInfo (using the even-odd
     * variants stored there)
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const AlignedVector<Number> &shape_hessians,
                            const unsigned int           dummy1 = 0,
                            const unsigned int           dummy2 = 0)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {
      (void)dummy1;
      (void)dummy2;
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number in [],
            Number       out[]) const
    {
      apply<direction,dof_to_quad,add,0>(shape_values, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const
    {
      apply<direction,dof_to_quad,add,1>(shape_gradients, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const
    {
      apply<direction,dof_to_quad,add,2>(shape_hessians, in, out);
    }

    template <int direction, bool dof_to_quad, bool add, int type>
    static void apply (const Number *shape_data,
                       const Number  in [],
                       Number        out []);

    const Number *shape_values;
    const Number *shape_gradients;
    const Number *shape_hessians;
  };



  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add, int type>
  inline
  void
  EvaluatorTensorProduct<evaluate_evenodd,dim,fe_degree,n_q_points_1d,Number>
  ::apply (const Number *shapes,
           const Number  in [],
           Number        out [])
  {
    AssertIndexRange (type, 3);
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
    const int n_cols = nn / 2;
    const int mid    = mm / 2;

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    const int offset = (n_q_points_1d+1)/2;

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
                if (dof_to_quad == true && type == 1)
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
            for (int col=0; col<n_cols; ++col)
              {
                Number r0, r1;
                if (mid > 0)
                  {
                    if (dof_to_quad == true)
                      {
                        r0 = shapes[col]                    * xp[0];
                        r1 = shapes[fe_degree*offset + col] * xm[0];
                      }
                    else
                      {
                        r0 = shapes[col*offset]             * xp[0];
                        r1 = shapes[(fe_degree-col)*offset] * xm[0];
                      }
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (dof_to_quad == true)
                          {
                            r0 += shapes[ind*offset+col]             * xp[ind];
                            r1 += shapes[(fe_degree-ind)*offset+col] * xm[ind];
                          }
                        else
                          {
                            r0 += shapes[col*offset+ind]             * xp[ind];
                            r1 += shapes[(fe_degree-col)*offset+ind] * xm[ind];
                          }
                      }
                  }
                else
                  r0 = r1 = Number();
                if (mm % 2 == 1 && dof_to_quad == true)
                  {
                    if (type == 1)
                      r1 += shapes[mid*offset+col] * in[stride*mid];
                    else
                      r0 += shapes[mid*offset+col] * in[stride*mid];
                  }
                else if (mm % 2 == 1 && (nn % 2 == 0 || type > 0))
                  r0 += shapes[col*offset+mid] * in[stride*mid];

                if (add == false)
                  {
                    out[stride*col]         = r0 + r1;
                    if (type == 1 && dof_to_quad == false)
                      out[stride*(nn-1-col)]  = r1 - r0;
                    else
                      out[stride*(nn-1-col)]  = r0 - r1;
                  }
                else
                  {
                    out[stride*col]        += r0 + r1;
                    if (type == 1 && dof_to_quad == false)
                      out[stride*(nn-1-col)] += r1 - r0;
                    else
                      out[stride*(nn-1-col)] += r0 - r1;
                  }
              }
            if ( type == 0 && dof_to_quad == true && nn%2==1 && mm%2==1 )
              {
                if (add==false)
                  out[stride*n_cols]  = in[stride*mid];
                else
                  out[stride*n_cols] += in[stride*mid];
              }
            else if (dof_to_quad == true && nn%2==1)
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
                  r0 += shapes[mid*offset+n_cols] * in[stride*mid];

                if (add == false)
                  out[stride*n_cols]  = r0;
                else
                  out[stride*n_cols] += r0;
              }
            else if (dof_to_quad == false && nn%2 == 1)
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
                  r0 += in[stride*mid];
                else if (type == 2 && mm % 2 == 1)
                  r0 += shapes[n_cols*offset+mid] * in[stride*mid];

                if (add == false)
                  out[stride*n_cols]  = r0;
                else
                  out[stride*n_cols] += r0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need to
            // jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
      }
  }

} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
