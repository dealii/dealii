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


#ifndef dealii__matrix_free_evaluation_kernels_h
#define dealii__matrix_free_evaluation_kernels_h

#include <deal.II/base/config.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/base/utilities.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>
#include <deal.II/matrix_free/shape_info.h>


DEAL_II_NAMESPACE_OPEN



namespace internal
{
  // Select evaluator type from element shape function type
  template <MatrixFreeFunctions::ElementType element, bool is_long>
  struct EvaluatorSelector {};

  template <bool is_long>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_general,is_long>
  {
    static const EvaluatorVariant variant = evaluate_general;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric,false>
  {
    static const EvaluatorVariant variant = evaluate_symmetric;
  };

  template <> struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric,true>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };

  template <bool is_long>
  struct EvaluatorSelector<MatrixFreeFunctions::truncated_tensor,is_long>
  {
    static const EvaluatorVariant variant = evaluate_general;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric_plus_dg0,false>
  {
    static const EvaluatorVariant variant = evaluate_general;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric_plus_dg0,true>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };

  template <bool is_long>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_gausslobatto,is_long>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };



  // This struct performs the evaluation of function values, gradients and
  // Hessians for tensor-product finite elements. The operation is used for
  // both the symmetric and non-symmetric case, which use different apply
  // functions 'values', 'gradients' in the individual coordinate
  // directions. The apply functions for values are provided through one of
  // the template classes EvaluatorTensorProduct which in turn are selected
  // from the MatrixFreeFunctions::ElementType template argument.
  //
  // There is a specialization made for Gauss-Lobatto elements further down
  // where the 'values' operation is identity, which allows us to write
  // shorter code.
  template <MatrixFreeFunctions::ElementType type, int dim, int fe_degree,
            int n_q_points_1d, int n_components, typename Number>
  struct FEEvaluationImpl
  {
    static
    void evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                   VectorizedArray<Number> *values_dofs_actual[],
                   VectorizedArray<Number> *values_quad[],
                   VectorizedArray<Number> *gradients_quad[][dim],
                   VectorizedArray<Number> *hessians_quad[][(dim*(dim+1))/2],
                   VectorizedArray<Number> *scratch_data,
                   const bool               evaluate_val,
                   const bool               evaluate_grad,
                   const bool               evaluate_lapl);

    static
    void integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                    VectorizedArray<Number> *values_dofs_actual[],
                    VectorizedArray<Number> *values_quad[],
                    VectorizedArray<Number> *gradients_quad[][dim],
                    VectorizedArray<Number> *scratch_data,
                    const bool               evaluate_val,
                    const bool               evaluate_grad);
  };


  template <MatrixFreeFunctions::ElementType type, int dim, int fe_degree,
            int n_q_points_1d, int n_components, typename Number>
  inline
  void
  FEEvaluationImpl<type,dim,fe_degree,n_q_points_1d,n_components,Number>
  ::evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
              VectorizedArray<Number> *values_dofs_actual[],
              VectorizedArray<Number> *values_quad[],
              VectorizedArray<Number> *gradients_quad[][dim],
              VectorizedArray<Number> *hessians_quad[][(dim*(dim+1))/2],
              VectorizedArray<Number> *scratch_data,
              const bool               evaluate_val,
              const bool               evaluate_grad,
              const bool               evaluate_lapl)
  {
    if (evaluate_val == false && evaluate_grad == false && evaluate_lapl == false)
      return;

    const EvaluatorVariant variant =
      EvaluatorSelector<type,(fe_degree+n_q_points_1d>4)>::variant;
    typedef EvaluatorTensorProduct<variant, dim, fe_degree, n_q_points_1d,
            VectorizedArray<Number> > Eval;
    Eval eval (variant == evaluate_evenodd ? shape_info.shape_val_evenodd :
               shape_info.shape_values,
               variant == evaluate_evenodd ? shape_info.shape_gra_evenodd :
               shape_info.shape_gradients,
               variant == evaluate_evenodd ? shape_info.shape_hes_evenodd :
               shape_info.shape_hessians,
               shape_info.fe_degree,
               shape_info.n_q_points_1d);

    const unsigned int temp_size = Eval::dofs_per_cell == numbers::invalid_unsigned_int ? 0
                                   : (Eval::dofs_per_cell > Eval::n_q_points ?
                                      Eval::dofs_per_cell : Eval::n_q_points);
    const unsigned int max_stack_size = 100;
    VectorizedArray<Number>  temp_data[(temp_size > 0 && temp_size < max_stack_size) ? 2*temp_size : 1];
    VectorizedArray<Number> *temp1;
    VectorizedArray<Number> *temp2;
    if (temp_size == 0)
      {
        temp1 = scratch_data;
        temp2 = temp1 + std::max(Utilities::fixed_power<dim>(shape_info.fe_degree+1),
                                 Utilities::fixed_power<dim>(shape_info.n_q_points_1d));
      }
    else if (temp_size < max_stack_size)
      {
        temp1 = &temp_data[0];
        temp2 = temp1 + temp_size;
      }
    else
      {
        temp1 = scratch_data;
        temp2 = temp1 + temp_size;
      }

    VectorizedArray<Number> **values_dofs = values_dofs_actual;
    VectorizedArray<Number> *expanded_dof_values[n_components];
    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        values_dofs = expanded_dof_values;
        for (unsigned int c=0; c<n_components; ++c)
          expanded_dof_values[c] = scratch_data+2*(std::max(shape_info.dofs_per_cell,
                                                            shape_info.n_q_points)) +
                                   c*Utilities::fixed_power<dim>(shape_info.fe_degree+1);
        const int degree = fe_degree != -1 ? fe_degree : shape_info.fe_degree;
        unsigned int count_p = 0, count_q = 0;
        for (int i=0; i<(dim>2?degree+1:1); ++i)
          {
            for (int j=0; j<(dim>1?degree+1-i:1); ++j)
              {
                for (int k=0; k<degree+1-j-i; ++k, ++count_p, ++count_q)
                  for (unsigned int c=0; c<n_components; ++c)
                    expanded_dof_values[c][count_q] = values_dofs_actual[c][count_p];
                for (int k=degree+1-j-i; k<degree+1; ++k, ++count_q)
                  for (unsigned int c=0; c<n_components; ++c)
                    expanded_dof_values[c][count_q] = VectorizedArray<Number>();
              }
            for (int j=degree+1-i; j<degree+1; ++j)
              for (int k=0; k<degree+1; ++k, ++count_q)
                for (unsigned int c=0; c<n_components; ++c)
                  expanded_dof_values[c][count_q] = VectorizedArray<Number>();
          }
        AssertDimension(count_q, Utilities::fixed_power<dim>(shape_info.fe_degree+1));
      }

    // These avoid compiler errors; they are only used in sensible context but
    // compilers typically cannot detect when we access something like
    // gradients_quad[2] only for dim==3.
    const unsigned int d1 = dim>1?1:0;
    const unsigned int d2 = dim>2?2:0;
    const unsigned int d3 = dim>2?3:0;
    const unsigned int d4 = dim>2?4:0;
    const unsigned int d5 = dim>2?5:0;

    switch (dim)
      {
      case 1:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (evaluate_val == true)
              eval.template values<0,true,false> (values_dofs[c], values_quad[c]);
            if (evaluate_grad == true)
              eval.template gradients<0,true,false>(values_dofs[c], gradients_quad[c][0]);
            if (evaluate_lapl == true)
              eval.template hessians<0,true,false> (values_dofs[c], hessians_quad[c][0]);
          }
        break;

      case 2:
        for (unsigned int c=0; c<n_components; c++)
          {
            // grad x
            if (evaluate_grad == true)
              {
                eval.template gradients<0,true,false> (values_dofs[c], temp1);
                eval.template values<1,true,false> (temp1, gradients_quad[c][0]);
              }
            if (evaluate_lapl == true)
              {
                // grad xy
                if (evaluate_grad == false)
                  eval.template gradients<0,true,false>(values_dofs[c], temp1);
                eval.template gradients<1,true,false>  (temp1, hessians_quad[c][d1+d1]);

                // grad xx
                eval.template hessians<0,true,false>(values_dofs[c], temp1);
                eval.template values<1,true,false>  (temp1, hessians_quad[c][0]);
              }

            // grad y
            eval.template values<0,true,false> (values_dofs[c], temp1);
            if (evaluate_grad == true)
              eval.template gradients<1,true,false> (temp1, gradients_quad[c][d1]);

            // grad yy
            if (evaluate_lapl == true)
              eval.template hessians<1,true,false> (temp1, hessians_quad[c][d1]);

            // val: can use values applied in x
            if (evaluate_val == true)
              eval.template values<1,true,false> (temp1, values_quad[c]);
          }
        break;

      case 3:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (evaluate_grad == true)
              {
                // grad x
                eval.template gradients<0,true,false> (values_dofs[c], temp1);
                eval.template values<1,true,false> (temp1, temp2);
                eval.template values<2,true,false> (temp2, gradients_quad[c][0]);
              }

            if (evaluate_lapl == true)
              {
                // grad xz
                if (evaluate_grad == false)
                  {
                    eval.template gradients<0,true,false> (values_dofs[c], temp1);
                    eval.template values<1,true,false> (temp1, temp2);
                  }
                eval.template gradients<2,true,false> (temp2, hessians_quad[c][d4]);

                // grad xy
                eval.template gradients<1,true,false> (temp1, temp2);
                eval.template values<2,true,false> (temp2, hessians_quad[c][d3]);

                // grad xx
                eval.template hessians<0,true,false>(values_dofs[c], temp1);
                eval.template values<1,true,false>  (temp1, temp2);
                eval.template values<2,true,false>  (temp2, hessians_quad[c][0]);
              }

            // grad y
            eval.template values<0,true,false> (values_dofs[c], temp1);
            if (evaluate_grad == true)
              {
                eval.template gradients<1,true,false>(temp1, temp2);
                eval.template values<2,true,false>   (temp2, gradients_quad[c][d1]);
              }

            if (evaluate_lapl == true)
              {
                // grad yz
                if (evaluate_grad == false)
                  eval.template gradients<1,true,false>(temp1, temp2);
                eval.template gradients<2,true,false>  (temp2, hessians_quad[c][d5]);

                // grad yy
                eval.template hessians<1,true,false> (temp1, temp2);
                eval.template values<2,true,false> (temp2, hessians_quad[c][d1]);
              }

            // grad z: can use the values applied in x direction stored in temp1
            eval.template values<1,true,false> (temp1, temp2);
            if (evaluate_grad == true)
              eval.template gradients<2,true,false> (temp2, gradients_quad[c][d2]);

            // grad zz: can use the values applied in x and y direction stored
            // in temp2
            if (evaluate_lapl == true)
              eval.template hessians<2,true,false>(temp2, hessians_quad[c][d2]);

            // val: can use the values applied in x & y direction stored in temp2
            if (evaluate_val == true)
              eval.template values<2,true,false> (temp2, values_quad[c]);
          }
        break;

      default:
        AssertThrow(false, ExcNotImplemented());
      }

    // case additional dof for FE_Q_DG0: add values; gradients and second
    // derivatives evaluate to zero
    if (type == MatrixFreeFunctions::tensor_symmetric_plus_dg0 && evaluate_val)
      for (unsigned int c=0; c<n_components; ++c)
        for (unsigned int q=0; q<shape_info.n_q_points; ++q)
          values_quad[c][q] += values_dofs[c][shape_info.dofs_per_cell-1];
  }



  template <MatrixFreeFunctions::ElementType type, int dim, int fe_degree,
            int n_q_points_1d, int n_components, typename Number>
  inline
  void
  FEEvaluationImpl<type,dim,fe_degree,n_q_points_1d,n_components,Number>
  ::integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
               VectorizedArray<Number> *values_dofs_actual[],
               VectorizedArray<Number> *values_quad[],
               VectorizedArray<Number> *gradients_quad[][dim],
               VectorizedArray<Number> *scratch_data,
               const bool               integrate_val,
               const bool               integrate_grad)
  {
    const EvaluatorVariant variant =
      EvaluatorSelector<type,(fe_degree+n_q_points_1d>4)>::variant;
    typedef EvaluatorTensorProduct<variant, dim, fe_degree, n_q_points_1d,
            VectorizedArray<Number> > Eval;
    Eval eval (variant == evaluate_evenodd ? shape_info.shape_val_evenodd :
               shape_info.shape_values,
               variant == evaluate_evenodd ? shape_info.shape_gra_evenodd :
               shape_info.shape_gradients,
               variant == evaluate_evenodd ? shape_info.shape_hes_evenodd :
               shape_info.shape_hessians,
               shape_info.fe_degree,
               shape_info.n_q_points_1d);

    const unsigned int temp_size = Eval::dofs_per_cell == numbers::invalid_unsigned_int ? 0
                                   : (Eval::dofs_per_cell > Eval::n_q_points ?
                                      Eval::dofs_per_cell : Eval::n_q_points);
    const unsigned int max_stack_size = 100;
    VectorizedArray<Number>  temp_data[(temp_size > 0 && temp_size < max_stack_size) ? 2*temp_size : 1];
    VectorizedArray<Number> *temp1;
    VectorizedArray<Number> *temp2;
    if (temp_size == 0)
      {
        temp1 = scratch_data;
        temp2 = temp1 + std::max(Utilities::fixed_power<dim>(shape_info.fe_degree+1),
                                 Utilities::fixed_power<dim>(shape_info.n_q_points_1d));
      }
    else if (temp_size < max_stack_size)
      {
        temp1 = &temp_data[0];
        temp2 = temp1 + temp_size;
      }
    else
      {
        temp1 = scratch_data;
        temp2 = temp1 + temp_size;
      }

    // expand dof_values to tensor product for truncated tensor products
    VectorizedArray<Number> **values_dofs = values_dofs_actual;
    VectorizedArray<Number> *expanded_dof_values[n_components];
    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        values_dofs = expanded_dof_values;
        for (unsigned int c=0; c<n_components; ++c)
          expanded_dof_values[c] = scratch_data+2*(std::max(shape_info.dofs_per_cell,
                                                            shape_info.n_q_points)) +
                                   c*Utilities::fixed_power<dim>(shape_info.fe_degree+1);
      }

    // These avoid compiler errors; they are only used in sensible context but
    // compilers typically cannot detect when we access something like
    // gradients_quad[2] only for dim==3.
    const unsigned int d1 = dim>1?1:0;
    const unsigned int d2 = dim>2?2:0;

    switch (dim)
      {
      case 1:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (integrate_val == true)
              eval.template values<0,false,false> (values_quad[c], values_dofs[c]);
            if (integrate_grad == true)
              {
                if (integrate_val == true)
                  eval.template gradients<0,false,true> (gradients_quad[c][0], values_dofs[c]);
                else
                  eval.template gradients<0,false,false> (gradients_quad[c][0], values_dofs[c]);
              }
          }
        break;

      case 2:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (integrate_val == true)
              {
                // val
                eval.template values<0,false,false> (values_quad[c], temp1);
                //grad x
                if (integrate_grad == true)
                  eval.template gradients<0,false,true> (gradients_quad[c][0], temp1);
                eval.template values<1,false,false>(temp1, values_dofs[c]);
              }
            if (integrate_grad == true)
              {
                // grad y
                eval.template values<0,false,false>  (gradients_quad[c][d1], temp1);
                if (integrate_val == false)
                  {
                    eval.template gradients<1,false,false>(temp1, values_dofs[c]);
                    //grad x
                    eval.template gradients<0,false,false> (gradients_quad[c][0], temp1);
                    eval.template values<1,false,true> (temp1, values_dofs[c]);
                  }
                else
                  eval.template gradients<1,false,true>(temp1, values_dofs[c]);
              }
          }
        break;

      case 3:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (integrate_val == true)
              {
                // val
                eval.template values<0,false,false> (values_quad[c], temp1);
                //grad x: can sum to temporary value in temp1
                if (integrate_grad == true)
                  eval.template gradients<0,false,true> (gradients_quad[c][0], temp1);
                eval.template values<1,false,false>(temp1, temp2);
                if (integrate_grad == true)
                  {
                    eval.template values<0,false,false> (gradients_quad[c][d1], temp1);
                    eval.template gradients<1,false,true>(temp1, temp2);
                  }
                eval.template values<2,false,false> (temp2, values_dofs[c]);
              }
            else if (integrate_grad == true)
              {
                eval.template gradients<0,false,false>(gradients_quad[c][0], temp1);
                eval.template values<1,false,false> (temp1, temp2);
                eval.template values<0,false,false> (gradients_quad[c][d1], temp1);
                eval.template gradients<1,false,true>(temp1, temp2);
                eval.template values<2,false,false> (temp2, values_dofs[c]);
              }
            if (integrate_grad == true)
              {
                // grad z: can sum to temporary x and y value in output
                eval.template values<0,false,false> (gradients_quad[c][d2], temp1);
                eval.template values<1,false,false> (temp1, temp2);
                eval.template gradients<2,false,true> (temp2, values_dofs[c]);
              }
          }
        break;

      default:
        AssertThrow(false, ExcNotImplemented());
      }

    // case FE_Q_DG0: add values, gradients and second derivatives are zero
    if (type == MatrixFreeFunctions::tensor_symmetric_plus_dg0)
      {
        if (integrate_val)
          for (unsigned int c=0; c<n_components; ++c)
            {
              values_dofs[c][shape_info.dofs_per_cell-1] = values_quad[c][0];
              for (unsigned int q=1; q<shape_info.n_q_points; ++q)
                values_dofs[c][shape_info.dofs_per_cell-1] += values_quad[c][q];
            }
        else
          for (unsigned int c=0; c<n_components; ++c)
            values_dofs[c][shape_info.dofs_per_cell-1] = VectorizedArray<Number>();
      }

    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        unsigned int count_p = 0, count_q = 0;
        const int degree = fe_degree != -1 ? fe_degree : shape_info.fe_degree;
        for (int i=0; i<(dim>2?degree+1:1); ++i)
          {
            for (int j=0; j<(dim>1?degree+1-i:1); ++j)
              {
                for (int k=0; k<degree+1-j-i; ++k, ++count_p, ++count_q)
                  {
                    for (unsigned int c=0; c<n_components; ++c)
                      values_dofs_actual[c][count_p] = expanded_dof_values[c][count_q];
                  }
                count_q += j+i;
              }
            count_q += i*(degree+1);
          }
        AssertDimension(count_q, Utilities::fixed_power<dim>(shape_info.fe_degree+1));
      }
  }

  // This a specialization for Gauss-Lobatto elements where the 'values'
  // operation is identity, which allows us to write shorter code.
  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  struct FEEvaluationImpl<MatrixFreeFunctions::tensor_gausslobatto, dim,
    fe_degree, n_q_points_1d, n_components, Number>
  {
    static
    void evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                   VectorizedArray<Number> *values_dofs[],
                   VectorizedArray<Number> *values_quad[],
                   VectorizedArray<Number> *gradients_quad[][dim],
                   VectorizedArray<Number> *hessians_quad[][(dim*(dim+1))/2],
                   VectorizedArray<Number> *scratch_data,
                   const bool               evaluate_val,
                   const bool               evaluate_grad,
                   const bool               evaluate_lapl);

    static
    void integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                    VectorizedArray<Number> *values_dofs[],
                    VectorizedArray<Number> *values_quad[],
                    VectorizedArray<Number> *gradients_quad[][dim],
                    VectorizedArray<Number> *scratch_data,
                    const bool               integrate_val,
                    const bool               integrate_grad);
  };

  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  inline
  void
  FEEvaluationImpl<MatrixFreeFunctions::tensor_gausslobatto, dim,
                   fe_degree, n_q_points_1d, n_components, Number>
                   ::evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                               VectorizedArray<Number> *values_dofs[],
                               VectorizedArray<Number> *values_quad[],
                               VectorizedArray<Number> *gradients_quad[][dim],
                               VectorizedArray<Number> *hessians_quad[][(dim*(dim+1))/2],
                               VectorizedArray<Number> *scratch_data,
                               const bool               evaluate_val,
                               const bool               evaluate_grad,
                               const bool               evaluate_lapl)
  {
    typedef EvaluatorTensorProduct<evaluate_evenodd, dim, fe_degree, fe_degree+1,
            VectorizedArray<Number> > Eval;
    Eval eval (shape_info.shape_val_evenodd,
               shape_info.shape_gra_evenodd,
               shape_info.shape_hes_evenodd,
               shape_info.fe_degree,
               shape_info.n_q_points_1d);

    // These avoid compiler errors; they are only used in sensible context but
    // compilers typically cannot detect when we access something like
    // gradients_quad[2] only for dim==3.
    const unsigned int d1 = dim>1?1:0;
    const unsigned int d2 = dim>2?2:0;
    const unsigned int d3 = dim>2?3:0;
    const unsigned int d4 = dim>2?4:0;
    const unsigned int d5 = dim>2?5:0;

    switch (dim)
      {
      case 1:
        if (evaluate_val == true)
          std::memcpy (values_quad[0], values_dofs[0],
                       eval.dofs_per_cell * n_components *
                       sizeof (values_dofs[0][0]));
        for (unsigned int c=0; c<n_components; c++)
          {
            if (evaluate_grad == true)
              eval.template gradients<0,true,false>(values_dofs[c], gradients_quad[c][0]);
            if (evaluate_lapl == true)
              eval.template hessians<0,true,false> (values_dofs[c], hessians_quad[c][0]);
          }
        break;

      case 2:
        if (evaluate_val == true)
          {
            std::memcpy (values_quad[0], values_dofs[0],
                         Eval::dofs_per_cell * n_components *
                         sizeof (values_dofs[0][0]));
          }
        if (evaluate_grad == true)
          for (unsigned int comp=0; comp<n_components; comp++)
            {
              // grad x
              eval.template gradients<0,true,false> (values_dofs[comp],
                                                     gradients_quad[comp][0]);
              // grad y
              eval.template gradients<1,true,false> (values_dofs[comp],
                                                     gradients_quad[comp][d1]);
            }
        if (evaluate_lapl == true)
          for (unsigned int comp=0; comp<n_components; comp++)
            {
              // hess x
              eval.template hessians<0,true,false> (values_dofs[comp],
                                                    hessians_quad[comp][0]);
              // hess y
              eval.template hessians<1,true,false> (values_dofs[comp],
                                                    hessians_quad[comp][d1]);

              // grad x grad y
              eval.template gradients<0,true,false> (values_dofs[comp], scratch_data);
              eval.template gradients<1,true,false> (scratch_data, hessians_quad[comp][d1+d1]);
            }
        break;

      case 3:
        if (evaluate_val == true)
          {
            std::memcpy (values_quad[0], values_dofs[0],
                         Eval::dofs_per_cell * n_components *
                         sizeof (values_dofs[0][0]));
          }
        if (evaluate_grad == true)
          for (unsigned int comp=0; comp<n_components; comp++)
            {
              // grad x
              eval.template gradients<0,true,false> (values_dofs[comp],
                                                     gradients_quad[comp][0]);
              // grad y
              eval.template gradients<1,true,false> (values_dofs[comp],
                                                     gradients_quad[comp][d1]);
              // grad y
              eval.template gradients<2,true,false> (values_dofs[comp],
                                                     gradients_quad[comp][d2]);
            }
        if (evaluate_lapl == true)
          for (unsigned int comp=0; comp<n_components; comp++)
            {
              // grad x
              eval.template hessians<0,true,false> (values_dofs[comp],
                                                    hessians_quad[comp][0]);
              // grad y
              eval.template hessians<1,true,false> (values_dofs[comp],
                                                    hessians_quad[comp][d1]);
              // grad y
              eval.template hessians<2,true,false> (values_dofs[comp],
                                                    hessians_quad[comp][d2]);

              VectorizedArray<Number> *temp1 = scratch_data;
              // grad xy
              eval.template gradients<0,true,false> (values_dofs[comp], temp1);
              eval.template gradients<1,true,false> (temp1, hessians_quad[comp][d3]);
              // grad xz
              eval.template gradients<2,true,false> (temp1, hessians_quad[comp][d4]);
              // grad yz
              eval.template gradients<1,true,false> (values_dofs[comp], temp1);
              eval.template gradients<2,true,false> (temp1, hessians_quad[comp][d5]);
            }
        break;
      default:
        AssertThrow(false, ExcNotImplemented());
      }
  }

  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  inline
  void
  FEEvaluationImpl<MatrixFreeFunctions::tensor_gausslobatto, dim,
                   fe_degree, n_q_points_1d, n_components, Number>
                   ::integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                                VectorizedArray<Number> *values_dofs[],
                                VectorizedArray<Number> *values_quad[],
                                VectorizedArray<Number> *gradients_quad[][dim],
                                VectorizedArray<Number> *,
                                const bool               integrate_val,
                                const bool               integrate_grad)
  {
    typedef EvaluatorTensorProduct<evaluate_evenodd, dim, fe_degree, fe_degree+1,
            VectorizedArray<Number> > Eval;
    Eval eval (shape_info.shape_val_evenodd,
               shape_info.shape_gra_evenodd,
               shape_info.shape_hes_evenodd,
               shape_info.fe_degree,
               shape_info.n_q_points_1d);

    // These avoid compiler errors; they are only used in sensible context but
    // compilers typically cannot detect when we access something like
    // gradients_quad[2] only for dim==3.
    const unsigned int d1 = dim>1?1:0;
    const unsigned int d2 = dim>2?2:0;

    if (integrate_val == true)
      std::memcpy (values_dofs[0], values_quad[0],
                   Eval::dofs_per_cell * n_components *
                   sizeof (values_dofs[0][0]));
    switch (dim)
      {
      case 1:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (integrate_grad == true)
              {
                if (integrate_val == true)
                  eval.template gradients<0,false,true> (gradients_quad[c][0],
                                                         values_dofs[c]);
                else
                  eval.template gradients<0,false,false> (gradients_quad[c][0],
                                                          values_dofs[c]);
              }
          }

        break;
      case 2:
        if (integrate_grad == true)
          for (unsigned int comp=0; comp<n_components; comp++)
            {
              // grad x: If integrate_val == true we have to add to the
              // previous output
              if (integrate_val == true)
                eval.template gradients<0, false, true> (gradients_quad[comp][0],
                                                         values_dofs[comp]);
              else
                eval.template gradients<0, false, false> (gradients_quad[comp][0],
                                                          values_dofs[comp]);

              // grad y
              eval.template gradients<1, false, true> (gradients_quad[comp][d1],
                                                       values_dofs[comp]);
            }
        break;

      case 3:
        if (integrate_grad == true)
          for (unsigned int comp=0; comp<n_components; comp++)
            {
              // grad x: If integrate_val == true we have to add to the
              // previous output
              if (integrate_val == true)
                eval.template gradients<0, false, true> (gradients_quad[comp][0],
                                                         values_dofs[comp]);
              else
                eval.template gradients<0, false, false> (gradients_quad[comp][0],
                                                          values_dofs[comp]);

              // grad y
              eval.template gradients<1, false, true> (gradients_quad[comp][d1],
                                                       values_dofs[comp]);

              // grad z
              eval.template gradients<2, false, true> (gradients_quad[comp][d2],
                                                       values_dofs[comp]);
            }
        break;

      default:
        AssertThrow(false, ExcNotImplemented());
      }
  }

} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
