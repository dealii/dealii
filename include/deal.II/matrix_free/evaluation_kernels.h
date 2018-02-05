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


#ifndef dealii_matrix_free_evaluation_kernels_h
#define dealii_matrix_free_evaluation_kernels_h

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
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric_collocation,is_long>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };



  /**
   * This struct performs the evaluation of function values, gradients and
   * Hessians for tensor-product finite elements. The operation is used for
   * both the symmetric and non-symmetric case, which use different apply
   * functions 'values', 'gradients' in the individual coordinate
   * directions. The apply functions for values are provided through one of
   * the template classes EvaluatorTensorProduct which in turn are selected
   * from the MatrixFreeFunctions::ElementType template argument.
   *
   * There are two specialized implementation classes
   * FEEvaluationImplCollocation (for Gauss-Lobatto elements where the nodal
   * points and the quadrature points coincide and the 'values' operation is
   * identity) and FEEvaluationImplTransformToCollocation (which can be
   * transformed to a collocation space and can then use the identity in these
   * spaces), which both allow for shorter code.
   *
   * @author Katharina Kormann, Martin Kronbichler, 2012, 2014, 2017
   */
  template <MatrixFreeFunctions::ElementType type, int dim, int fe_degree,
            int n_q_points_1d, int n_components, typename Number>
  struct FEEvaluationImpl
  {
    static
    void evaluate (const MatrixFreeFunctions::ShapeInfo<Number>  &shape_info,
                   const Number *values_dofs_actual,
                   Number *values_quad,
                   Number *gradients_quad,
                   Number *hessians_quad,
                   Number *scratch_data,
                   const bool               evaluate_values,
                   const bool               evaluate_gradients,
                   const bool               evaluate_hessians);

    static
    void integrate (const MatrixFreeFunctions::ShapeInfo<Number>  &shape_info,
                    Number *values_dofs_actual,
                    Number *values_quad,
                    Number *gradients_quad,
                    Number *scratch_data,
                    const bool               integrate_values,
                    const bool               integrate_gradients,
                    const bool               add_into_values_array);
  };



  template <MatrixFreeFunctions::ElementType type, int dim, int fe_degree,
            int n_q_points_1d, int n_components, typename Number>
  inline
  void
  FEEvaluationImpl<type,dim,fe_degree,n_q_points_1d,n_components,Number>
  ::evaluate (const MatrixFreeFunctions::ShapeInfo<Number>  &shape_info,
              const Number *values_dofs_actual,
              Number *values_quad,
              Number *gradients_quad,
              Number *hessians_quad,
              Number *scratch_data,
              const bool               evaluate_values,
              const bool               evaluate_gradients,
              const bool               evaluate_hessians)
  {
    if (evaluate_values == false && evaluate_gradients == false && evaluate_hessians == false)
      return;

    const EvaluatorVariant variant =
      EvaluatorSelector<type,(fe_degree+n_q_points_1d>4)>::variant;
    typedef EvaluatorTensorProduct<variant, dim, fe_degree+1, n_q_points_1d,
            Number> Eval;
    Eval eval (variant == evaluate_evenodd ? shape_info.shape_values_eo :
               shape_info.shape_values,
               variant == evaluate_evenodd ? shape_info.shape_gradients_eo :
               shape_info.shape_gradients,
               variant == evaluate_evenodd ? shape_info.shape_hessians_eo :
               shape_info.shape_hessians,
               shape_info.fe_degree+1,
               shape_info.n_q_points_1d);

    const unsigned int temp_size = Eval::n_rows_of_product == numbers::invalid_unsigned_int ? 0
                                   : (Eval::n_rows_of_product > Eval::n_columns_of_product ?
                                      Eval::n_rows_of_product : Eval::n_columns_of_product);
    Number *temp1;
    Number *temp2;
    if (temp_size == 0)
      {
        temp1 = scratch_data;
        temp2 = temp1 + std::max(Utilities::fixed_power<dim>(shape_info.fe_degree+1),
                                 Utilities::fixed_power<dim>(shape_info.n_q_points_1d));
      }
    else
      {
        temp1 = scratch_data;
        temp2 = temp1 + temp_size;
      }

    const unsigned int n_q_points = temp_size == 0 ? shape_info.n_q_points : Eval::n_columns_of_product;
    const unsigned int dofs_per_comp = (type == MatrixFreeFunctions::truncated_tensor) ?
                                       Utilities::fixed_power<dim>(shape_info.fe_degree+1) : shape_info.dofs_per_component_on_cell;
    const Number *values_dofs = values_dofs_actual;
    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        Number *values_dofs_tmp = scratch_data+2*(std::max(shape_info.dofs_per_component_on_cell, shape_info.n_q_points));
        const int degree = fe_degree != -1 ? fe_degree : shape_info.fe_degree;
        unsigned int count_p = 0, count_q = 0;
        for (int i=0; i<(dim>2?degree+1:1); ++i)
          {
            for (int j=0; j<(dim>1?degree+1-i:1); ++j)
              {
                for (int k=0; k<degree+1-j-i; ++k, ++count_p, ++count_q)
                  for (unsigned int c=0; c<n_components; ++c)
                    values_dofs_tmp[c*dofs_per_comp+count_q] = values_dofs_actual[c*shape_info.dofs_per_component_on_cell+count_p];
                for (int k=degree+1-j-i; k<degree+1; ++k, ++count_q)
                  for (unsigned int c=0; c<n_components; ++c)
                    values_dofs_tmp[c*dofs_per_comp+count_q] = Number();
              }
            for (int j=degree+1-i; j<degree+1; ++j)
              for (int k=0; k<degree+1; ++k, ++count_q)
                for (unsigned int c=0; c<n_components; ++c)
                  values_dofs_tmp[c*dofs_per_comp+count_q] = Number();
          }
        AssertDimension(count_q, dofs_per_comp);
        values_dofs = values_dofs_tmp;
      }

    switch (dim)
      {
      case 1:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (evaluate_values == true)
              eval.template values<0,true,false> (values_dofs, values_quad);
            if (evaluate_gradients == true)
              eval.template gradients<0,true,false>(values_dofs, gradients_quad);
            if (evaluate_hessians == true)
              eval.template hessians<0,true,false> (values_dofs, hessians_quad);

            // advance the next component in 1D array
            values_dofs += dofs_per_comp;
            values_quad += n_q_points;
            gradients_quad += n_q_points;
            hessians_quad += n_q_points;
          }
        break;

      case 2:
        for (unsigned int c=0; c<n_components; c++)
          {
            // grad x
            if (evaluate_gradients == true)
              {
                eval.template gradients<0,true,false> (values_dofs, temp1);
                eval.template values<1,true,false> (temp1, gradients_quad);
              }
            if (evaluate_hessians == true)
              {
                // grad xy
                if (evaluate_gradients == false)
                  eval.template gradients<0,true,false>(values_dofs, temp1);
                eval.template gradients<1,true,false>  (temp1, hessians_quad+2*n_q_points);

                // grad xx
                eval.template hessians<0,true,false>(values_dofs, temp1);
                eval.template values<1,true,false>  (temp1, hessians_quad);
              }

            // grad y
            eval.template values<0,true,false> (values_dofs, temp1);
            if (evaluate_gradients == true)
              eval.template gradients<1,true,false> (temp1, gradients_quad+n_q_points);

            // grad yy
            if (evaluate_hessians == true)
              eval.template hessians<1,true,false> (temp1, hessians_quad+n_q_points);

            // val: can use values applied in x
            if (evaluate_values == true)
              eval.template values<1,true,false> (temp1, values_quad);

            // advance to the next component in 1D array
            values_dofs += dofs_per_comp;
            values_quad += n_q_points;
            gradients_quad += 2*n_q_points;
            hessians_quad += 3*n_q_points;
          }
        break;

      case 3:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (evaluate_gradients == true)
              {
                // grad x
                eval.template gradients<0,true,false> (values_dofs, temp1);
                eval.template values<1,true,false> (temp1, temp2);
                eval.template values<2,true,false> (temp2, gradients_quad);
              }

            if (evaluate_hessians == true)
              {
                // grad xz
                if (evaluate_gradients == false)
                  {
                    eval.template gradients<0,true,false> (values_dofs, temp1);
                    eval.template values<1,true,false> (temp1, temp2);
                  }
                eval.template gradients<2,true,false> (temp2, hessians_quad+4*n_q_points);

                // grad xy
                eval.template gradients<1,true,false> (temp1, temp2);
                eval.template values<2,true,false> (temp2, hessians_quad+3*n_q_points);

                // grad xx
                eval.template hessians<0,true,false>(values_dofs, temp1);
                eval.template values<1,true,false>  (temp1, temp2);
                eval.template values<2,true,false>  (temp2, hessians_quad);
              }

            // grad y
            eval.template values<0,true,false> (values_dofs, temp1);
            if (evaluate_gradients == true)
              {
                eval.template gradients<1,true,false>(temp1, temp2);
                eval.template values<2,true,false>   (temp2, gradients_quad+n_q_points);
              }

            if (evaluate_hessians == true)
              {
                // grad yz
                if (evaluate_gradients == false)
                  eval.template gradients<1,true,false>(temp1, temp2);
                eval.template gradients<2,true,false>  (temp2, hessians_quad+5*n_q_points);

                // grad yy
                eval.template hessians<1,true,false> (temp1, temp2);
                eval.template values<2,true,false> (temp2, hessians_quad+n_q_points);
              }

            // grad z: can use the values applied in x direction stored in temp1
            eval.template values<1,true,false> (temp1, temp2);
            if (evaluate_gradients == true)
              eval.template gradients<2,true,false> (temp2, gradients_quad+2*n_q_points);

            // grad zz: can use the values applied in x and y direction stored
            // in temp2
            if (evaluate_hessians == true)
              eval.template hessians<2,true,false>(temp2, hessians_quad+2*n_q_points);

            // val: can use the values applied in x & y direction stored in temp2
            if (evaluate_values == true)
              eval.template values<2,true,false> (temp2, values_quad);

            // advance to the next component in 1D array
            values_dofs += dofs_per_comp;
            values_quad += n_q_points;
            gradients_quad += 3*n_q_points;
            hessians_quad += 6*n_q_points;
          }
        break;

      default:
        AssertThrow(false, ExcNotImplemented());
      }

    // case additional dof for FE_Q_DG0: add values; gradients and second
    // derivatives evaluate to zero
    if (type == MatrixFreeFunctions::tensor_symmetric_plus_dg0 && evaluate_values)
      {
        values_quad -= n_components*n_q_points;
        values_dofs -= n_components*dofs_per_comp;
        for (unsigned int c=0; c<n_components; ++c)
          for (unsigned int q=0; q<shape_info.n_q_points; ++q)
            values_quad[c*shape_info.n_q_points+q] +=
              values_dofs[(c+1)*shape_info.dofs_per_component_on_cell-1];
      }
  }



  template <MatrixFreeFunctions::ElementType type, int dim, int fe_degree,
            int n_q_points_1d, int n_components, typename Number>
  inline
  void
  FEEvaluationImpl<type,dim,fe_degree,n_q_points_1d,n_components,Number>
  ::integrate (const MatrixFreeFunctions::ShapeInfo<Number>  &shape_info,
               Number *values_dofs_actual,
               Number *values_quad,
               Number *gradients_quad,
               Number *scratch_data,
               const bool               integrate_values,
               const bool               integrate_gradients,
               const bool               add_into_values_array)
  {
    const EvaluatorVariant variant =
      EvaluatorSelector<type,(fe_degree+n_q_points_1d>4)>::variant;
    typedef EvaluatorTensorProduct<variant, dim, fe_degree+1, n_q_points_1d,
            Number> Eval;
    Eval eval (variant == evaluate_evenodd ? shape_info.shape_values_eo :
               shape_info.shape_values,
               variant == evaluate_evenodd ? shape_info.shape_gradients_eo :
               shape_info.shape_gradients,
               variant == evaluate_evenodd ? shape_info.shape_hessians_eo :
               shape_info.shape_hessians,
               shape_info.fe_degree+1,
               shape_info.n_q_points_1d);

    const unsigned int temp_size = Eval::n_rows_of_product == numbers::invalid_unsigned_int ? 0
                                   : (Eval::n_rows_of_product > Eval::n_columns_of_product ?
                                      Eval::n_rows_of_product : Eval::n_columns_of_product);
    Number *temp1;
    Number *temp2;
    if (temp_size == 0)
      {
        temp1 = scratch_data;
        temp2 = temp1 + std::max(Utilities::fixed_power<dim>(shape_info.fe_degree+1),
                                 Utilities::fixed_power<dim>(shape_info.n_q_points_1d));
      }
    else
      {
        temp1 = scratch_data;
        temp2 = temp1 + temp_size;
      }

    const unsigned int n_q_points = temp_size == 0 ? shape_info.n_q_points : Eval::n_columns_of_product;
    const unsigned int dofs_per_comp = (type == MatrixFreeFunctions::truncated_tensor) ?
                                       Utilities::fixed_power<dim>(shape_info.fe_degree+1) : shape_info.dofs_per_component_on_cell;
    // expand dof_values to tensor product for truncated tensor products
    Number *values_dofs = (type == MatrixFreeFunctions::truncated_tensor) ?
                          scratch_data+2*(std::max(shape_info.dofs_per_component_on_cell,
                                                   shape_info.n_q_points)) :
                          values_dofs_actual;

    switch (dim)
      {
      case 1:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (integrate_values == true)
              {
                if (add_into_values_array == false)
                  eval.template values<0,false,false> (values_quad, values_dofs);
                else
                  eval.template values<0,false,true> (values_quad, values_dofs);
              }
            if (integrate_gradients == true)
              {
                if (integrate_values == true || add_into_values_array == true)
                  eval.template gradients<0,false,true> (gradients_quad, values_dofs);
                else
                  eval.template gradients<0,false,false> (gradients_quad, values_dofs);
              }

            // advance to the next component in 1D array
            values_dofs += dofs_per_comp;
            values_quad += n_q_points;
            gradients_quad += n_q_points;
          }
        break;

      case 2:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (integrate_values == true &&
                integrate_gradients == false)
              {
                eval.template values<1,false,false> (values_quad, temp1);
                if (add_into_values_array == false)
                  eval.template values<0,false,false>(temp1, values_dofs);
                else
                  eval.template values<0,false,true>(temp1, values_dofs);
              }
            if (integrate_gradients == true)
              {
                eval.template gradients<1,false,false> (gradients_quad+n_q_points, temp1);
                if (integrate_values)
                  eval.template values<1,false,true> (values_quad, temp1);
                if (add_into_values_array == false)
                  eval.template values<0,false,false>(temp1, values_dofs);
                else
                  eval.template values<0,false,true>(temp1, values_dofs);
                eval.template values<1,false,false> (gradients_quad, temp1);
                eval.template gradients<0,false,true> (temp1, values_dofs);
              }

            // advance to the next component in 1D array
            values_dofs += dofs_per_comp;
            values_quad += n_q_points;
            gradients_quad += 2*n_q_points;
          }
        break;

      case 3:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (integrate_values == true &&
                integrate_gradients == false)
              {
                eval.template values<2,false,false>  (values_quad, temp1);
                eval.template values<1,false,false>  (temp1, temp2);
                if (add_into_values_array == false)
                  eval.template values<0,false,false>(temp2, values_dofs);
                else
                  eval.template values<0,false,true> (temp2, values_dofs);
              }
            if (integrate_gradients == true)
              {
                eval.template gradients<2,false,false>(gradients_quad+2*n_q_points, temp1);
                if (integrate_values)
                  eval.template values<2,false,true>  (values_quad, temp1);
                eval.template values<1,false,false>   (temp1, temp2);
                eval.template values<2,false,false>   (gradients_quad+n_q_points, temp1);
                eval.template gradients<1,false,true> (temp1, temp2);
                if (add_into_values_array == false)
                  eval.template values<0,false,false> (temp2, values_dofs);
                else
                  eval.template values<0,false,true>  (temp2, values_dofs);
                eval.template values<2,false,false>   (gradients_quad, temp1);
                eval.template values<1,false,false>   (temp1, temp2);
                eval.template gradients<0,false,true> (temp2, values_dofs);
              }

            // advance to the next component in 1D array
            values_dofs += dofs_per_comp;
            values_quad += n_q_points;
            gradients_quad += 3*n_q_points;
          }
        break;

      default:
        AssertThrow(false, ExcNotImplemented());
      }

    // case FE_Q_DG0: add values, gradients and second derivatives are zero
    if (type == MatrixFreeFunctions::tensor_symmetric_plus_dg0)
      {
        values_dofs -= n_components * dofs_per_comp - shape_info.dofs_per_component_on_cell + 1;
        values_quad -= n_components * n_q_points;
        if (integrate_values)
          for (unsigned int c=0; c<n_components; ++c)
            {
              values_dofs[0] = values_quad[0];
              for (unsigned int q=1; q<shape_info.n_q_points; ++q)
                values_dofs[0] += values_quad[q];
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
            }
        else
          {
            for (unsigned int c=0; c<n_components; ++c)
              values_dofs[c*shape_info.dofs_per_component_on_cell] = Number();
            values_dofs += n_components*shape_info.dofs_per_component_on_cell;
          }
      }

    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        values_dofs -= dofs_per_comp*n_components;
        unsigned int count_p = 0, count_q = 0;
        const int degree = fe_degree != -1 ? fe_degree : shape_info.fe_degree;
        for (int i=0; i<(dim>2?degree+1:1); ++i)
          {
            for (int j=0; j<(dim>1?degree+1-i:1); ++j)
              {
                for (int k=0; k<degree+1-j-i; ++k, ++count_p, ++count_q)
                  {
                    for (unsigned int c=0; c<n_components; ++c)
                      values_dofs_actual[c*shape_info.dofs_per_component_on_cell+count_p] = values_dofs[c*dofs_per_comp+count_q];
                  }
                count_q += j+i;
              }
            count_q += i*(degree+1);
          }
        AssertDimension(count_q, Utilities::fixed_power<dim>(shape_info.fe_degree+1));
      }
  }



  /**
   * This struct implements the change between two different bases. This is an
   * ingredient in the FEEvaluationImplTransformToCollocation class where we
   * first transform to the appropriate basis where we can compute the
   * derivative through collocation techniques.
   *
   * This class allows for dimension-independent application of the operation,
   * implemented by template recursion. It has been tested up to 6D.
   *
   * @author Katharina Kormann, Martin Kronbichler, 2017
   */
  template <EvaluatorVariant variant, int dim, int basis_size_1, int basis_size_2, int n_components,
            typename Number, typename Number2>
  struct FEEvaluationImplBasisChange
  {
    static_assert(basis_size_1 == 0 || basis_size_1 <= basis_size_2,
                  "The second dimension must not be smaller than the first");

    /**
     * This applies the transformation that contracts over the rows of the
     * coefficient array, generating values along the columns of the
     * coefficient array.
     *
     * @param transformation_matrix The coefficient matrix handed in as a
     *                     vector, using @p basis_size_1 rows and @p basis_size_2
     *                     columns if interpreted as a matrix.
     * @param values_in    The array of the input of size basis_size_1^dim. It
     *                     may alias with values_out
     * @param values_out   The array of size basis_size_2^dim where the results
     *                     of the transformation are stored. It may alias with
     *                     the values_in array.
     * @param basis_size_1_variable In case the template argument basis_size_1 is
     *                     zero, the size of the first basis can alternatively be
     *                     passed in as a run time argument. The template
     *                     argument takes precedence in case it is nonzero
     *                     for efficiency reasons.
     * @param basis_size_2_variable In case the template argument basis_size_1 is
     *                     zero, the size of the second basis can alternatively be
     *                     passed in as a run time argument.
     */
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
    static void do_forward (const AlignedVector<Number2> &transformation_matrix,
                            const Number      *values_in,
                            Number            *values_out,
                            const unsigned int basis_size_1_variable = numbers::invalid_unsigned_int,
                            const unsigned int basis_size_2_variable = numbers::invalid_unsigned_int)
    {
      Assert(basis_size_1 != 0 ||
             basis_size_1_variable <= basis_size_2_variable,
             ExcMessage("The second dimension must not be smaller than the first"));

      // we do recursion until dim==1 or dim==2 and we have
      // basis_size_1==basis_size_2. The latter optimization increases
      // optimization possibilities for the compiler but does only work for
      // aliased pointers if the sizes are equal.
      constexpr int next_dim = (dim > 2 || ((basis_size_1 == 0 || basis_size_2>basis_size_1)
                                            && dim>1)) ? dim-1 : dim;

      EvaluatorTensorProduct<variant, dim, basis_size_1, (basis_size_1==0 ? 0 : basis_size_2),
                             Number,Number2> eval_val (transformation_matrix,
                                                       AlignedVector<Number2>(),
                                                       AlignedVector<Number2>(),
                                                       basis_size_1_variable,
                                                       basis_size_2_variable);
      const unsigned int np_1 = basis_size_1 > 0 ? basis_size_1 : basis_size_1_variable;
      const unsigned int np_2 = basis_size_1 > 0 ? basis_size_2 : basis_size_2_variable;
      Assert(np_1 > 0 && np_1 != numbers::invalid_unsigned_int,
             ExcMessage("Cannot transform with 0-point basis"));
      Assert(np_2 > 0 && np_2 != numbers::invalid_unsigned_int,
             ExcMessage("Cannot transform with 0-point basis"));

      // run loop backwards to ensure correctness if values_in aliases with
      // values_out in case with basis_size_1 < basis_size_2
      values_in = values_in + n_components*Utilities::fixed_power<dim>(np_1);
      values_out = values_out + n_components*Utilities::fixed_power<dim>(np_2);
      for (unsigned int c=n_components; c!=0; --c)
        {
          values_in -= Utilities::fixed_power<dim>(np_1);
          values_out -= Utilities::fixed_power<dim>(np_2);
          if (next_dim < dim)
            for (unsigned int q=np_1; q!=0; --q)
              FEEvaluationImplBasisChange<variant,next_dim,basis_size_1,basis_size_2,1,Number,Number2>
              ::do_forward(transformation_matrix,
                           values_in + (q-1)*Utilities::fixed_power<next_dim>(np_1),
                           values_out + (q-1)*Utilities::fixed_power<next_dim>(np_2),
                           basis_size_1_variable,
                           basis_size_2_variable);

          // the recursion stops if dim==1 or if dim==2 and
          // basis_size_1==basis_size_2 (the latter is used because the
          // compiler generates nicer code)
          if (basis_size_1 > 0 && basis_size_2 == basis_size_1 && dim == 2)
            {
              eval_val.template values<0,true,false>(values_in, values_out);
              eval_val.template values<1,true,false>(values_out, values_out);
            }
          else if (dim==1)
            eval_val.template values<dim-1,true,false>(values_in, values_out);
          else
            eval_val.template values<dim-1,true,false>(values_out, values_out);
        }
    }

    /**
     * This applies the transformation that contracts over the columns of the
     * coefficient array, generating values along the rows of the coefficient
     * array.
     *
     * @param transformation_matrix The coefficient matrix handed in as a
     *                     vector, using @p basis_size_1 rows and @p basis_size_2
     *                     columns if interpreted as a matrix.
     * @param add_into_result Define whether the result should be added into the
     *                     array @p values_out (if true) or overwrite the
     *                     previous content. The result is undefined in case
     *                     values_in and values_out point to the same array and
     *                     @p add_into_result is true, in which case an
     *                     exception is thrown.
     * @param values_in    The array of the input of size basis_size_2^dim. It
     *                     may alias with values_out. Note that the previous
     *                     content of @p values_in is overwritten within the
     *                     function.
     * @param values_out   The array of size basis_size_1^dim where the results
     *                     of the transformation are stored. It may alias with
     *                     the @p values_in array.
     * @param basis_size_1_variable In case the template argument basis_size_1 is
     *                     zero, the size of the first basis can alternatively be
     *                     passed in as a run time argument. The template
     *                     argument takes precedence in case it is nonzero
     *                     for efficiency reasons.
     * @param basis_size_2_variable In case the template argument basis_size_1 is
     *                     zero, the size of the second basis can alternatively be
     *                     passed in as a run time argument.
     */
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
    static void do_backward (const AlignedVector<Number2> &transformation_matrix,
                             const bool         add_into_result,
                             Number            *values_in,
                             Number            *values_out,
                             const unsigned int basis_size_1_variable = numbers::invalid_unsigned_int,
                             const unsigned int basis_size_2_variable = numbers::invalid_unsigned_int)
    {
      Assert(basis_size_1 != 0 ||
             basis_size_1_variable <= basis_size_2_variable,
             ExcMessage("The second dimension must not be smaller than the first"));
      Assert(add_into_result == false || values_in != values_out,
             ExcMessage("Input and output cannot alias with each other when "
                        "adding the result of the basis change to existing data"));

      constexpr int next_dim = (dim > 2 || ((basis_size_1 == 0 || basis_size_2>basis_size_1)
                                            && dim>1)) ? dim-1 : dim;
      EvaluatorTensorProduct<variant, dim, basis_size_1, (basis_size_1==0 ? 0 : basis_size_2),
                             Number,Number2> eval_val (transformation_matrix,
                                                       AlignedVector<Number2>(),
                                                       AlignedVector<Number2>(),
                                                       basis_size_1_variable,
                                                       basis_size_2_variable);
      const unsigned int np_1 = basis_size_1 > 0 ? basis_size_1 : basis_size_1_variable;
      const unsigned int np_2 = basis_size_1 > 0 ? basis_size_2 : basis_size_2_variable;
      Assert(np_1 > 0 && np_1 != numbers::invalid_unsigned_int,
             ExcMessage("Cannot transform with 0-point basis"));
      Assert(np_2 > 0 && np_2 != numbers::invalid_unsigned_int,
             ExcMessage("Cannot transform with 0-point basis"));

      for (unsigned int c=0; c<n_components; ++c)
        {
          if (basis_size_1 > 0 && basis_size_2 == basis_size_1 && dim == 2)
            {
              eval_val.template values<1,false,false>(values_in, values_in);
              if (add_into_result)
                eval_val.template values<0,false,true>(values_in, values_out);
              else
                eval_val.template values<0,false,false>(values_in, values_out);
            }
          else
            {
              if (dim==1 && add_into_result)
                eval_val.template values<0,false,true>(values_in, values_out);
              else if (dim==1)
                eval_val.template values<0,false,false>(values_in, values_out);
              else
                eval_val.template values<dim-1,false,false>(values_in, values_in);
            }
          if (next_dim < dim)
            for (unsigned int q=0; q<np_1; ++q)
              FEEvaluationImplBasisChange<variant,next_dim,basis_size_1,basis_size_2,1,Number,Number2>
              ::do_backward(transformation_matrix,
                            add_into_result,
                            values_in + q*Utilities::fixed_power<next_dim>(np_2),
                            values_out + q*Utilities::fixed_power<next_dim>(np_1),
                            basis_size_1_variable, basis_size_2_variable);

          values_in += Utilities::fixed_power<dim>(np_2);
          values_out += Utilities::fixed_power<dim>(np_1);
        }
    }

    /**
     * This operation applies a mass-matrix-like operation, consisting of a
     * do_forward() operation, multiplication by the coefficients in the
     * quadrature points, and the do_backward() operation.
     *
     * @param transformation_matrix The coefficient matrix handed in as a
     *                     vector, using @p basis_size_1 rows and @p basis_size_2
     *                     columns if interpreted as a matrix.
     * @param coefficients The array of coefficients by which the result is
     *                     multiplied. Its length must be either
     *                     basis_size_2^dim or n_components*basis_size_2^dim
     * @param values_in    The array of the input of size basis_size_2^dim. It
     *                     may alias with values_out
     * @param scratch_data Array to hold temporary data during the operation.
     *                     Must be of length basis_size_2^dim
     * @param values_out   The array of size basis_size_1^dim where the results
     *                     of the transformation are stored. It may alias with
     *                     the values_in array.
     */
    static void do_mass (const AlignedVector<Number2> &transformation_matrix,
                         const AlignedVector<Number>  &coefficients,
                         const Number *values_in,
                         Number *scratch_data,
                         Number *values_out)
    {
      constexpr int next_dim = dim > 1 ? dim-1 : dim;
      Number *my_scratch = basis_size_1 != basis_size_2 ? scratch_data : values_out;
      for (unsigned int q=basis_size_1; q!=0; --q)
        FEEvaluationImplBasisChange<variant,next_dim,basis_size_1,basis_size_2,n_components,Number,Number2>
        ::do_forward(transformation_matrix,
                     values_in + (q-1)*Utilities::fixed_int_power<basis_size_1,dim-1>::value,
                     my_scratch + (q-1)*Utilities::fixed_int_power<basis_size_2,dim-1>::value);
      EvaluatorTensorProduct<variant, dim, basis_size_1, basis_size_2,
                             Number,Number2> eval_val (transformation_matrix);
      const unsigned int n_inner_blocks = (dim > 1 && basis_size_2 < 10) ? basis_size_2 : 1;
      const unsigned int n_blocks = Utilities::fixed_int_power<basis_size_2,dim-1>::value;
      for (unsigned int ii=0; ii<n_blocks; ii+=n_inner_blocks)
        for (unsigned int c=0; c<n_components; ++c)
          {
            for (unsigned int i=ii; i<ii+n_inner_blocks; ++i)
              eval_val.template values_one_line<dim-1,true,false> (my_scratch+i, my_scratch+i);
            for (unsigned int q=0; q<basis_size_2; ++q)
              for (unsigned int i=ii; i<ii+n_inner_blocks; ++i)
                my_scratch[i+q*n_blocks] *= coefficients[i+q*n_blocks];
            for (unsigned int i=ii; i<ii+n_inner_blocks; ++i)
              eval_val.template values_one_line<dim-1,false,false>(my_scratch+i, my_scratch+i);
          }
      for (unsigned int q=0; q<basis_size_1; ++q)
        FEEvaluationImplBasisChange<variant,next_dim,basis_size_1,basis_size_2,n_components,Number,Number2>
        ::do_backward(transformation_matrix, false,
                      my_scratch + q*Utilities::fixed_int_power<basis_size_2,dim-1>::value,
                      values_out + q*Utilities::fixed_int_power<basis_size_1,dim-1>::value);
    }
  };



  /**
   * This struct performs the evaluation of function values, gradients and
   * Hessians for tensor-product finite elements. This a specialization for
   * elements where the nodal points coincide with the quadrature points like
   * FE_Q shape functions on Gauss-Lobatto elements integrated with
   * Gauss-Lobatto quadrature. The assumption of this class is that the shape
   * 'values' operation is identity, which allows us to write shorter code.
   *
   * In literature, this form of evaluation is often called spectral
   * evaluation, spectral collocation or simply collocation, meaning the same
   * location for shape functions and evaluation space (quadrature points).
   *
   * @author Katharina Kormann, 2012
  */
  template <int dim, int fe_degree, int n_components, typename Number>
  struct FEEvaluationImplCollocation
  {
    static
    void evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                   const Number *values_dofs,
                   Number       *values_quad,
                   Number       *gradients_quad,
                   Number       *hessians_quad,
                   Number       *scratch_data,
                   const bool    evaluate_values,
                   const bool    evaluate_gradients,
                   const bool    evaluate_hessians);

    static
    void integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                    Number      *values_dofs,
                    Number      *values_quad,
                    Number      *gradients_quad,
                    Number      *scratch_data,
                    const bool   integrate_values,
                    const bool   integrate_gradients,
                    const bool   add_into_values_array);
  };



  template <int dim, int fe_degree, int n_components, typename Number>
  inline
  void
  FEEvaluationImplCollocation<dim, fe_degree, n_components, Number>
  ::evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
              const Number *values_dofs,
              Number       *values_quad,
              Number       *gradients_quad,
              Number       *hessians_quad,
              Number *,
              const bool    evaluate_values,
              const bool    evaluate_gradients,
              const bool    evaluate_hessians)
  {
    AssertDimension(shape_info.shape_gradients_collocation_eo.size(),
                    (fe_degree+2)/2*(fe_degree+1));

    EvaluatorTensorProduct<evaluate_evenodd, dim, fe_degree+1, fe_degree+1, Number>
    eval(AlignedVector<Number>(),
         shape_info.shape_gradients_collocation_eo,
         shape_info.shape_hessians_collocation_eo);
    constexpr unsigned int n_q_points = Utilities::fixed_int_power<fe_degree+1,dim>::value;

    for (unsigned int c=0; c<n_components; c++)
      {
        if (evaluate_values == true)
          for (unsigned int i=0; i<n_q_points; ++i)
            values_quad[i] = values_dofs[i];
        if (evaluate_gradients == true || evaluate_hessians == true)
          {
            eval.template gradients<0,true,false>(values_dofs, gradients_quad);
            if (dim > 1)
              eval.template gradients<1,true,false>(values_dofs, gradients_quad+n_q_points);
            if (dim > 2)
              eval.template gradients<2,true,false>(values_dofs, gradients_quad+2*n_q_points);
          }
        if (evaluate_hessians == true)
          {
            eval.template hessians<0,true,false> (values_dofs, hessians_quad);
            if (dim > 1)
              {
                eval.template gradients<1,true,false> (gradients_quad, hessians_quad+dim*n_q_points);
                eval.template hessians<1,true,false> (values_dofs, hessians_quad+n_q_points);
              }
            if (dim > 2)
              {
                eval.template gradients<2,true,false> (gradients_quad, hessians_quad+4*n_q_points);
                eval.template gradients<2,true,false> (gradients_quad+n_q_points, hessians_quad+5*n_q_points);
                eval.template hessians<2,true,false> (values_dofs, hessians_quad+2*n_q_points);
              }
            hessians_quad += (dim*(dim+1))/2*n_q_points;
          }
        gradients_quad += dim*n_q_points;
        values_quad += n_q_points;
        values_dofs += n_q_points;
      }
  }



  template <int dim, int fe_degree, int n_components, typename Number>
  inline
  void
  FEEvaluationImplCollocation<dim, fe_degree, n_components, Number>
  ::integrate (const MatrixFreeFunctions::ShapeInfo<Number>  &shape_info,
               Number     *values_dofs,
               Number     *values_quad,
               Number     *gradients_quad,
               Number *,
               const bool  integrate_values,
               const bool  integrate_gradients,
               const bool  add_into_values_array)
  {
    AssertDimension(shape_info.shape_gradients_collocation_eo.size(),
                    (fe_degree+2)/2*(fe_degree+1));

    EvaluatorTensorProduct<evaluate_evenodd, dim, fe_degree+1, fe_degree+1, Number>
    eval(AlignedVector<Number>(),
         shape_info.shape_gradients_collocation_eo,
         shape_info.shape_hessians_collocation_eo);
    constexpr unsigned int n_q_points = Utilities::fixed_int_power<fe_degree+1,dim>::value;

    for (unsigned int c=0; c<n_components; c++)
      {
        if (integrate_values == true && add_into_values_array == false)
          for (unsigned int i=0; i<n_q_points; ++i)
            values_dofs[i] = values_quad[i];
        else if (integrate_values == true)
          for (unsigned int i=0; i<n_q_points; ++i)
            values_dofs[i] += values_quad[i];
        if (integrate_gradients == true)
          {
            if (integrate_values == true || add_into_values_array == true)
              eval.template gradients<0,false,true>(gradients_quad, values_dofs);
            else
              eval.template gradients<0,false,false>(gradients_quad, values_dofs);
            if (dim > 1)
              eval.template gradients<1,false,true>(gradients_quad+n_q_points, values_dofs);
            if (dim > 2)
              eval.template gradients<2,false,true>(gradients_quad+2*n_q_points, values_dofs);
          }
        gradients_quad += dim*n_q_points;
        values_quad += n_q_points;
        values_dofs += n_q_points;
      }
  }




  /**
   * This struct performs the evaluation of function values, gradients and
   * Hessians for tensor-product finite elements. This a specialization for
   * symmetric basis functions about the mid point 0.5 of the unit interval
   * with the same number of quadrature points as degrees of freedom. In that
   * case, we can first transform the basis to one that has the nodal points
   * in the quadrature points (i.e., the collocation space) and then perform
   * the evaluation of the first and second derivatives in this transformed
   * space, using the identity operation for the shape values.
   *
   * @author Katharina Kormann, Martin Kronbichler, 2017
   */
  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  struct FEEvaluationImplTransformToCollocation
  {
    static
    void evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                   const Number *values_dofs,
                   Number       *values_quad,
                   Number       *gradients_quad,
                   Number       *hessians_quad,
                   Number       *scratch_data,
                   const bool    evaluate_values,
                   const bool    evaluate_gradients,
                   const bool    evaluate_hessians);

    static
    void integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                    Number      *values_dofs,
                    Number      *values_quad,
                    Number      *gradients_quad,
                    Number      *scratch_data,
                    const bool   integrate_values,
                    const bool   integrate_gradients,
                    const bool   add_into_values_array);
  };



  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  inline
  void
  FEEvaluationImplTransformToCollocation<dim, fe_degree, n_q_points_1d, n_components, Number>
  ::evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
              const Number  *values_dofs,
              Number        *values_quad,
              Number        *gradients_quad,
              Number        *hessians_quad,
              Number *,
              const bool     ,
              const bool     evaluate_gradients,
              const bool     evaluate_hessians)
  {
    Assert(n_q_points_1d > fe_degree,
           ExcMessage("You lose information when going to a collocation space "
                      "of lower degree, so the evaluation results would be "
                      "wrong. Thus, this class does not permit the desired "
                      "operation."));
    constexpr unsigned int n_q_points = Utilities::fixed_int_power<n_q_points_1d,dim>::value;

    for (unsigned int c=0; c<n_components; c++)
      {
        FEEvaluationImplBasisChange<evaluate_evenodd, dim,
                                    (fe_degree>=n_q_points_1d?n_q_points_1d:fe_degree+1),
                                    n_q_points_1d,1,Number,Number>
                                    ::do_forward(shape_info.shape_values_eo,
                                                 values_dofs, values_quad);

        // apply derivatives in the collocation space
        if (evaluate_gradients == true || evaluate_hessians == true)
          FEEvaluationImplCollocation<dim,n_q_points_1d-1,1,Number>::
          evaluate(shape_info, values_quad, nullptr, gradients_quad, hessians_quad,
                   nullptr, false, evaluate_gradients, evaluate_hessians);

        values_dofs += shape_info.dofs_per_component_on_cell;
        values_quad += n_q_points;
        gradients_quad += dim*n_q_points;
        hessians_quad += (dim*(dim+1))/2*n_q_points;
      }
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  inline
  void
  FEEvaluationImplTransformToCollocation<dim, fe_degree, n_q_points_1d, n_components, Number>
  ::integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
               Number      *values_dofs,
               Number      *values_quad,
               Number      *gradients_quad,
               Number *,
               const bool   integrate_values,
               const bool   integrate_gradients,
               const bool   add_into_values_array)
  {
    Assert(n_q_points_1d > fe_degree,
           ExcMessage("You lose information when going to a collocation space "
                      "of lower degree, so the evaluation results would be "
                      "wrong. Thus, this class does not permit the desired "
                      "operation."));
    AssertDimension(shape_info.shape_gradients_collocation_eo.size(),
                    (n_q_points_1d+1)/2*n_q_points_1d);
    constexpr unsigned int n_q_points = Utilities::fixed_int_power<n_q_points_1d,dim>::value;

    for (unsigned int c=0; c<n_components; c++)
      {

        // apply derivatives in collocation space
        if (integrate_gradients == true)
          FEEvaluationImplCollocation<dim,n_q_points_1d-1,1,Number>::
          integrate(shape_info, values_quad, nullptr, gradients_quad, nullptr, false,
                    integrate_gradients,/*add_into_values_array=*/integrate_values);

        // transform back to the original space
        FEEvaluationImplBasisChange<evaluate_evenodd, dim,
                                    (fe_degree>=n_q_points_1d?n_q_points_1d:fe_degree+1),
                                    n_q_points_1d,1,Number,Number>
                                    ::do_backward(shape_info.shape_values_eo,
                                                  add_into_values_array,
                                                  values_quad,
                                                  values_dofs);
        gradients_quad += dim*n_q_points;
        values_quad += n_q_points;
        values_dofs += shape_info.dofs_per_component_on_cell;
      }
  }



  template <bool symmetric_evaluate, int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  struct FEFaceEvaluationImpl
  {
    static
    void evaluate_in_face (const MatrixFreeFunctions::ShapeInfo<Number>  &data,
                           Number             *values_dofs,
                           Number             *values_quad,
                           Number             *gradients_quad,
                           Number             *scratch_data,
                           const bool          evaluate_val,
                           const bool          evaluate_grad,
                           const unsigned int  subface_index)
    {
      const AlignedVector<Number> &val1
        = symmetric_evaluate ? data.shape_values_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
           data.shape_values : data.values_within_subface[subface_index%2]);
      const AlignedVector<Number> &val2
        = symmetric_evaluate ? data.shape_values_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
           data.shape_values : data.values_within_subface[subface_index/2]);

      const AlignedVector<Number> &grad1
        = symmetric_evaluate ? data.shape_gradients_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
           data.shape_gradients : data.gradients_within_subface[subface_index%2]);
      const AlignedVector<Number> &grad2
        = symmetric_evaluate ? data.shape_gradients_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
           data.shape_gradients : data.gradients_within_subface[subface_index/2]);

      typedef internal::EvaluatorTensorProduct
      <symmetric_evaluate ? internal::evaluate_evenodd :internal::evaluate_general,
      dim-1,fe_degree+1,n_q_points_1d,Number> Eval;
      typedef internal::EvaluatorTensorProduct
      <internal::evaluate_general,dim-1,fe_degree+1,n_q_points_1d,
      Number> EvalGeneric;
      Eval eval1(val1,grad1,AlignedVector<Number>(),
                 data.fe_degree+1, data.n_q_points_1d);
      Eval eval2(val2,grad2,AlignedVector<Number>(),
                 data.fe_degree+1, data.n_q_points_1d);

      const unsigned int size_deg = fe_degree > -1 ?
                                    Utilities::fixed_int_power<fe_degree+1,dim-1>::value :
                                    (dim > 1 ? Utilities::fixed_power<dim-1>(data.fe_degree+1) : 1);

      const unsigned int n_q_points = fe_degree > -1 ?
                                      Utilities::fixed_int_power<n_q_points_1d,dim-1>::value : data.n_q_points_face;

      if (evaluate_grad == false)
        for (unsigned int c=0; c<n_components; ++c)
          {
            switch (dim)
              {
              case 3:
                eval1.template values<0,true,false>(values_dofs, values_quad);
                eval2.template values<1,true,false>(values_quad, values_quad);
                break;
              case 2:
                eval1.template values<0,true,false>(values_dofs, values_quad);
                break;
              case 1:
                values_quad[c] = values_dofs[2*c];
                break;
              default:
                Assert(false, ExcNotImplemented());
              }
            values_dofs += 2*size_deg;
            values_quad += n_q_points;
          }
      else
        for (unsigned int c=0; c<n_components; ++c)
          {
            switch (dim)
              {
              case 3:
                if (symmetric_evaluate && n_q_points_1d > fe_degree)
                  {
                    eval1.template values<0,true,false>(values_dofs, values_quad);
                    eval1.template values<1,true,false>(values_quad, values_quad);
                    internal::EvaluatorTensorProduct
                    <internal::evaluate_evenodd,dim-1,n_q_points_1d,n_q_points_1d,Number> eval_grad
                    (AlignedVector<Number>(),
                     data.shape_gradients_collocation_eo,
                     AlignedVector<Number>());
                    eval_grad.template gradients<0,true,false>(values_quad, gradients_quad);
                    eval_grad.template gradients<1,true,false>(values_quad,
                                                               gradients_quad+n_q_points);
                  }
                else
                  {
                    eval1.template gradients<0,true,false>(values_dofs, scratch_data);
                    eval2.template values<1,true,false>(scratch_data, gradients_quad);

                    eval1.template values<0,true,false>(values_dofs, scratch_data);
                    eval2.template gradients<1,true,false>(scratch_data, gradients_quad+n_q_points);

                    if (evaluate_val == true)
                      eval2.template values<1,true,false>(scratch_data, values_quad);
                  }
                eval1.template values<0,true,false>(values_dofs+size_deg, scratch_data);
                eval2.template values<1,true,false>(scratch_data,
                                                    gradients_quad+(dim-1)*n_q_points);

                break;
              case 2:
                eval1.template values<0,true,false>(values_dofs+size_deg,
                                                    gradients_quad+(dim-1)*n_q_points);
                eval1.template gradients<0,true,false>(values_dofs, gradients_quad);
                if (evaluate_val == true)
                  eval1.template values<0,true,false>(values_dofs, values_quad);
                break;
              case 1:
                values_quad[0] = values_dofs[0];
                gradients_quad[0] = values_dofs[1];
                break;
              default:
                AssertThrow(false, ExcNotImplemented());
              }
            values_dofs += 2*size_deg;
            values_quad += n_q_points;
            gradients_quad += dim*n_q_points;
          }
    }

    static
    void integrate_in_face (const MatrixFreeFunctions::ShapeInfo<Number> &data,
                            Number             *values_dofs,
                            Number             *values_quad,
                            Number             *gradients_quad,
                            Number             *scratch_data,
                            const bool          integrate_val,
                            const bool          integrate_grad,
                            const unsigned int  subface_index)
    {
      const AlignedVector<Number> &val1
        = symmetric_evaluate ? data.shape_values_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
           data.shape_values : data.values_within_subface[subface_index%2]);
      const AlignedVector<Number> &val2
        = symmetric_evaluate ? data.shape_values_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
           data.shape_values : data.values_within_subface[subface_index/2]);

      const AlignedVector<Number> &grad1
        = symmetric_evaluate ? data.shape_gradients_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
           data.shape_gradients : data.gradients_within_subface[subface_index%2]);
      const AlignedVector<Number> &grad2
        = symmetric_evaluate ? data.shape_gradients_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
           data.shape_gradients : data.gradients_within_subface[subface_index/2]);

      typedef internal::EvaluatorTensorProduct
      <symmetric_evaluate ? internal::evaluate_evenodd :internal::evaluate_general,
      dim-1,fe_degree+1,n_q_points_1d,Number> Eval;
      typedef internal::EvaluatorTensorProduct
      <internal::evaluate_general,dim-1,fe_degree+1,n_q_points_1d,
      Number> EvalGeneric;
      Eval eval1(val1,grad1,val1,data.fe_degree+1, data.n_q_points_1d);
      Eval eval2(val2,grad2,val1,data.fe_degree+1, data.n_q_points_1d);

      const unsigned int size_deg = fe_degree > -1 ?
                                    Utilities::fixed_int_power<fe_degree+1,dim-1>::value :
                                    (dim > 1 ? Utilities::fixed_power<dim-1>(data.fe_degree+1) : 1);

      const unsigned int n_q_points = fe_degree > -1 ?
                                      Utilities::fixed_int_power<n_q_points_1d,dim-1>::value : data.n_q_points_face;

      if (integrate_grad == false)
        for (unsigned int c=0; c<n_components; ++c)
          {
            switch (dim)
              {
              case 3:
                eval2.template values<1,false,false>(values_quad, values_quad);
                eval1.template values<0,false,false>(values_quad, values_dofs);
                break;
              case 2:
                eval1.template values<0,false,false>(values_quad, values_dofs);
                break;
              case 1:
                values_dofs[2*c] = values_quad[c][0];
                break;
              default:
                Assert(false, ExcNotImplemented());
              }
            values_dofs += 2*size_deg;
            values_quad += n_q_points;
          }
      else
        for (unsigned int c=0; c<n_components; ++c)
          {
            switch (dim)
              {
              case 3:
                eval2.template values<1,false,false> (gradients_quad+2*n_q_points,
                                                      gradients_quad+2*n_q_points);
                eval1.template values<0,false,false> (gradients_quad+2*n_q_points,
                                                      values_dofs+size_deg);
                if (symmetric_evaluate && n_q_points_1d > fe_degree)
                  {
                    internal::EvaluatorTensorProduct <internal::evaluate_evenodd,
                             dim-1,n_q_points_1d,n_q_points_1d,Number> eval_grad
                             (AlignedVector<Number>(),
                              data.shape_gradients_collocation_eo,
                              AlignedVector<Number>());
                    if (integrate_val)
                      eval_grad.template gradients<1,false,true>(gradients_quad+n_q_points,
                                                                 values_quad);
                    else
                      eval_grad.template gradients<1,false,false>(gradients_quad+n_q_points,
                                                                  values_quad);
                    eval_grad.template gradients<0,false,true>(gradients_quad,
                                                               values_quad);
                    eval1.template values<1,false,false>(values_quad, values_quad);
                    eval1.template values<0,false,false>(values_quad, values_dofs);
                  }
                else
                  {
                    if (integrate_val)
                      {
                        eval2.template values<1,false,false> (values_quad, scratch_data);
                        eval2.template gradients<1,false,true> (gradients_quad+n_q_points,
                                                                scratch_data);
                      }
                    else
                      eval2.template gradients<1,false,false> (gradients_quad+n_q_points,
                                                               scratch_data);

                    eval1.template values<0,false,false> (scratch_data, values_dofs);
                    eval2.template values<1,false,false> (gradients_quad, scratch_data);
                    eval1.template gradients<0,false,true> (scratch_data, values_dofs);
                  }
                break;
              case 2:
                eval1.template values<0,false,false>(gradients_quad+n_q_points,
                                                     values_dofs+size_deg);
                eval1.template gradients<0,false,false>(gradients_quad, values_dofs);
                if (integrate_val == true)
                  eval1.template values<0,false,true>(values_quad, values_dofs);
                break;
              case 1:
                values_dofs[0] = values_quad[0];
                values_dofs[1] = gradients_quad[0];
                break;
              default:
                AssertThrow(false, ExcNotImplemented());
              }
            values_dofs += 2*size_deg;
            values_quad += n_q_points;
            gradients_quad += dim*n_q_points;
          }
    }
  };



  template <int dim, int fe_degree, int n_components, typename Number>
  struct FEFaceNormalEvaluationImpl
  {
    template <bool do_evaluate, bool add_into_output>
    static void interpolate(const MatrixFreeFunctions::ShapeInfo<Number> &data,
                            const Number       *input,
                            Number             *output,
                            const bool          do_gradients,
                            const unsigned int  face_no)
    {
      internal::EvaluatorTensorProduct<internal::evaluate_general,dim,
               fe_degree+1,0,Number>
               evalf(data.shape_data_on_face[face_no%2],
                     AlignedVector<Number>(),
                     AlignedVector<Number>(),
                     data.fe_degree+1, 0);

      const unsigned int in_stride = do_evaluate ? data.dofs_per_component_on_cell : 2*data.dofs_per_component_on_face;
      const unsigned int out_stride = do_evaluate ? 2*data.dofs_per_component_on_face : data.dofs_per_component_on_cell;
      const unsigned int face_direction = face_no / 2;
      for (unsigned int c=0; c<n_components; c++)
        {
          if (do_gradients)
            {
              if (face_direction == 0)
                evalf.template apply_face<0,do_evaluate,add_into_output,1>(input, output);
              else if (face_direction == 1)
                evalf.template apply_face<1,do_evaluate,add_into_output,1>(input, output);
              else
                evalf.template apply_face<2,do_evaluate,add_into_output,1>(input, output);
            }
          else
            {
              if (face_direction == 0)
                evalf.template apply_face<0,do_evaluate,add_into_output,0>(input, output);
              else if (face_direction == 1)
                evalf.template apply_face<1,do_evaluate,add_into_output,0>(input, output);
              else
                evalf.template apply_face<2,do_evaluate,add_into_output,0>(input, output);
            }
          input += in_stride;
          output += out_stride;
        }
    }
  };
} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
