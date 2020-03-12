// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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


#ifndef dealii_matrix_free_evaluation_selector_h
#define dealii_matrix_free_evaluation_selector_h

#include <deal.II/base/config.h>

#include <deal.II/matrix_free/evaluation_kernels.h>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
namespace internal
{
  namespace EvaluationSelectorImplementation
  {
    // The following classes serve the purpose of choosing the correct template
    // specialization of the FEEvaluationImpl* classes in case fe_degree
    // and n_q_points_1d are only given as runtime parameters.
    // The logic is the following:
    // 1. Start with fe_degree=0, n_q_points_1d=0 and DEPTH=0.
    // 2. If the current assumption on fe_degree doesn't match the runtime
    //    parameter, increase fe_degree  by one and try again.
    //    If fe_degree==10 use the class Default which serves as a fallback.
    // 3. After fixing the fe_degree, DEPTH is increased (DEPTH=1) and we start
    // with
    //    n_q_points=fe_degree+1.
    // 4. If the current assumption on n_q_points_1d doesn't match the runtime
    //    parameter, increase n_q_points_1d by one and try again.
    //    If n_q_points_1d==degree+3 use the class Default which serves as a
    //    fallback.

    /**
     * This class serves as a fallback in case we don't have the appropriate
     * template specialization for the run time and template parameters given.
     */
    template <int dim, int n_components, typename Number>
    struct Default
    {
      static inline void
      evaluate(
        const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
        Number *   values_dofs_actual,
        Number *   values_quad,
        Number *   gradients_quad,
        Number *   hessians_quad,
        Number *   scratch_data,
        const bool evaluate_values,
        const bool evaluate_gradients,
        const bool evaluate_hessians)
      {
        internal::FEEvaluationImpl<
          internal::MatrixFreeFunctions::tensor_general,
          dim,
          -1,
          0,
          n_components,
          Number>::evaluate(shape_info,
                            values_dofs_actual,
                            values_quad,
                            gradients_quad,
                            hessians_quad,
                            scratch_data,
                            evaluate_values,
                            evaluate_gradients,
                            evaluate_hessians);
      }

      static inline void
      integrate(
        const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
        Number *   values_dofs_actual,
        Number *   values_quad,
        Number *   gradients_quad,
        Number *   scratch_data,
        const bool integrate_values,
        const bool integrate_gradients,
        const bool sum_into_values_array = false)
      {
        internal::FEEvaluationImpl<
          internal::MatrixFreeFunctions::tensor_general,
          dim,
          -1,
          0,
          n_components,
          Number>::integrate(shape_info,
                             values_dofs_actual,
                             values_quad,
                             gradients_quad,
                             scratch_data,
                             integrate_values,
                             integrate_gradients,
                             sum_into_values_array);
      }
    };


    /**
     * This class implements the actual choice of the template specialization.
     */
    template <int dim,
              int n_components,
              typename Number,
              int DEPTH         = 0,
              int degree        = 0,
              int n_q_points_1d = 0,
              class Enable      = void>
    struct Factory : Default<dim, n_components, Number>
    {};

    /**
     * This specialization sets the maximal fe_degree for
     * which we want to determine the correct template parameters based at
     * runtime.
     */
    template <int n_q_points_1d, int dim, int n_components, typename Number>
    struct Factory<dim, n_components, Number, 0, 10, n_q_points_1d>
      : Default<dim, n_components, Number>
    {};

    /**
     * This specialization sets the maximal number of n_q_points_1d for
     * which we want to determine the correct template parameters based at
     * runtime.
     */
    template <int degree,
              int n_q_points_1d,
              int dim,
              int n_components,
              typename Number>
    struct Factory<dim,
                   n_components,
                   Number,
                   1,
                   degree,
                   n_q_points_1d,
                   typename std::enable_if<n_q_points_1d == degree + 3>::type>
      : Default<dim, n_components, Number>
    {};

    /**
     * This class chooses the correct template degree.
     */
    template <int degree,
              int n_q_points_1d,
              int dim,
              int n_components,
              typename Number>
    struct Factory<dim, n_components, Number, 0, degree, n_q_points_1d>
    {
      static inline void
      evaluate(
        const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
        Number *   values_dofs_actual,
        Number *   values_quad,
        Number *   gradients_quad,
        Number *   hessians_quad,
        Number *   scratch_data,
        const bool evaluate_values,
        const bool evaluate_gradients,
        const bool evaluate_hessians)
      {
        const unsigned int runtime_degree = shape_info.data.front().fe_degree;
        constexpr unsigned int start_n_q_points = degree + 1;
        if (runtime_degree == degree)
          Factory<dim, n_components, Number, 1, degree, start_n_q_points>::
            evaluate(shape_info,
                     values_dofs_actual,
                     values_quad,
                     gradients_quad,
                     hessians_quad,
                     scratch_data,
                     evaluate_values,
                     evaluate_gradients,
                     evaluate_hessians);
        else
          Factory<dim, n_components, Number, 0, degree + 1, n_q_points_1d>::
            evaluate(shape_info,
                     values_dofs_actual,
                     values_quad,
                     gradients_quad,
                     hessians_quad,
                     scratch_data,
                     evaluate_values,
                     evaluate_gradients,
                     evaluate_hessians);
      }

      static inline void
      integrate(
        const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
        Number *   values_dofs_actual,
        Number *   values_quad,
        Number *   gradients_quad,
        Number *   scratch_data,
        const bool integrate_values,
        const bool integrate_gradients,
        const bool sum_into_values_array = false)
      {
        const int runtime_degree = shape_info.data.front().fe_degree;
        constexpr unsigned int start_n_q_points = degree + 1;
        if (runtime_degree == degree)
          Factory<dim, n_components, Number, 1, degree, start_n_q_points>::
            integrate(shape_info,
                      values_dofs_actual,
                      values_quad,
                      gradients_quad,
                      scratch_data,
                      integrate_values,
                      integrate_gradients,
                      sum_into_values_array);
        else
          Factory<dim, n_components, Number, 0, degree + 1, n_q_points_1d>::
            integrate(shape_info,
                      values_dofs_actual,
                      values_quad,
                      gradients_quad,
                      scratch_data,
                      integrate_values,
                      integrate_gradients,
                      sum_into_values_array);
      }
    };

    /**
     * This class chooses the correct template n_q_points_1d after degree was
     * chosen.
     */
    template <int degree,
              int n_q_points_1d,
              int dim,
              int n_components,
              typename Number>
    struct Factory<dim,
                   n_components,
                   Number,
                   1,
                   degree,
                   n_q_points_1d,
                   typename std::enable_if<(n_q_points_1d < degree + 3)>::type>
    {
      /**
       * We enable a transformation to collocation for derivatives if it gives
       * correct results (first condition), if it is the most efficient choice
       * in terms of operation counts (second condition) and if we were able
       * to initialize the fields in shape_info.templates.h from the
       * polynomials (third condition).
       */
      static constexpr bool      use_collocation =
        n_q_points_1d > degree &&n_q_points_1d <= 3 * degree / 2 + 1 &&
        n_q_points_1d < 200;

      static inline void
      evaluate(
        const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
        Number *   values_dofs_actual,
        Number *   values_quad,
        Number *   gradients_quad,
        Number *   hessians_quad,
        Number *   scratch_data,
        const bool evaluate_values,
        const bool evaluate_gradients,
        const bool evaluate_hessians)
      {
        const int runtime_n_q_points_1d = shape_info.data.front().n_q_points_1d;
        if (runtime_n_q_points_1d == n_q_points_1d)
          {
            if (n_q_points_1d == degree + 1 &&
                shape_info.element_type ==
                  internal::MatrixFreeFunctions::tensor_symmetric_collocation)
              internal::
                FEEvaluationImplCollocation<dim, degree, n_components, Number>::
                  evaluate(shape_info,
                           values_dofs_actual,
                           values_quad,
                           gradients_quad,
                           hessians_quad,
                           scratch_data,
                           evaluate_values,
                           evaluate_gradients,
                           evaluate_hessians);
            else if (use_collocation)
              internal::FEEvaluationImplTransformToCollocation<
                dim,
                degree,
                n_q_points_1d,
                n_components,
                Number>::evaluate(shape_info,
                                  values_dofs_actual,
                                  values_quad,
                                  gradients_quad,
                                  hessians_quad,
                                  scratch_data,
                                  evaluate_values,
                                  evaluate_gradients,
                                  evaluate_hessians);
            else
              internal::FEEvaluationImpl<
                internal::MatrixFreeFunctions::tensor_symmetric,
                dim,
                degree,
                n_q_points_1d,
                n_components,
                Number>::evaluate(shape_info,
                                  values_dofs_actual,
                                  values_quad,
                                  gradients_quad,
                                  hessians_quad,
                                  scratch_data,
                                  evaluate_values,
                                  evaluate_gradients,
                                  evaluate_hessians);
          }
        else
          Factory<dim, n_components, Number, 1, degree, n_q_points_1d + 1>::
            evaluate(shape_info,
                     values_dofs_actual,
                     values_quad,
                     gradients_quad,
                     hessians_quad,
                     scratch_data,
                     evaluate_values,
                     evaluate_gradients,
                     evaluate_hessians);
      }

      static inline void
      integrate(
        const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
        Number *   values_dofs_actual,
        Number *   values_quad,
        Number *   gradients_quad,
        Number *   scratch_data,
        const bool integrate_values,
        const bool integrate_gradients,
        const bool sum_into_values_array)
      {
        const int runtime_n_q_points_1d = shape_info.data.front().n_q_points_1d;
        if (runtime_n_q_points_1d == n_q_points_1d)
          {
            if (n_q_points_1d == degree + 1 &&
                shape_info.element_type ==
                  internal::MatrixFreeFunctions::tensor_symmetric_collocation)
              internal::
                FEEvaluationImplCollocation<dim, degree, n_components, Number>::
                  integrate(shape_info,
                            values_dofs_actual,
                            values_quad,
                            gradients_quad,
                            scratch_data,
                            integrate_values,
                            integrate_gradients,
                            sum_into_values_array);
            else if (use_collocation)
              internal::FEEvaluationImplTransformToCollocation<
                dim,
                degree,
                n_q_points_1d,
                n_components,
                Number>::integrate(shape_info,
                                   values_dofs_actual,
                                   values_quad,
                                   gradients_quad,
                                   scratch_data,
                                   integrate_values,
                                   integrate_gradients,
                                   sum_into_values_array);
            else
              internal::FEEvaluationImpl<
                internal::MatrixFreeFunctions::tensor_symmetric,
                dim,
                degree,
                n_q_points_1d,
                n_components,
                Number>::integrate(shape_info,
                                   values_dofs_actual,
                                   values_quad,
                                   gradients_quad,
                                   scratch_data,
                                   integrate_values,
                                   integrate_gradients,
                                   sum_into_values_array);
          }
        else
          Factory<dim, n_components, Number, 1, degree, n_q_points_1d + 1>::
            integrate(shape_info,
                      values_dofs_actual,
                      values_quad,
                      gradients_quad,
                      scratch_data,
                      integrate_values,
                      integrate_gradients,
                      sum_into_values_array);
      }
    };



    /**
     * This is the entry point for choosing the correct runtime parameters
     * for the 'evaluate' function.
     */
    template <int dim, int n_components, typename Number>
    void
    symmetric_selector_evaluate(
      const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
      Number *   values_dofs_actual,
      Number *   values_quad,
      Number *   gradients_quad,
      Number *   hessians_quad,
      Number *   scratch_data,
      const bool evaluate_values,
      const bool evaluate_gradients,
      const bool evaluate_hessians)
    {
      Assert(shape_info.element_type <=
               internal::MatrixFreeFunctions::tensor_symmetric,
             ExcInternalError());
      Factory<dim, n_components, Number>::evaluate(shape_info,
                                                   values_dofs_actual,
                                                   values_quad,
                                                   gradients_quad,
                                                   hessians_quad,
                                                   scratch_data,
                                                   evaluate_values,
                                                   evaluate_gradients,
                                                   evaluate_hessians);
    }



    /**
     * This is the entry point for choosing the correct runtime parameters
     * for the 'integrate' function.
     */
    template <int dim, int n_components, typename Number>
    void
    symmetric_selector_integrate(
      const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
      Number *   values_dofs_actual,
      Number *   values_quad,
      Number *   gradients_quad,
      Number *   scratch_data,
      const bool integrate_values,
      const bool integrate_gradients,
      const bool sum_into_values_array = false)
    {
      Assert(shape_info.element_type <=
               internal::MatrixFreeFunctions::tensor_symmetric,
             ExcInternalError());
      Factory<dim, n_components, Number>::integrate(shape_info,
                                                    values_dofs_actual,
                                                    values_quad,
                                                    gradients_quad,
                                                    scratch_data,
                                                    integrate_values,
                                                    integrate_gradients,
                                                    sum_into_values_array);
    }
  } // namespace EvaluationSelectorImplementation
} // namespace internal
#endif


/**
 * This class chooses an appropriate evaluation strategy based on the
 * template parameters and the shape_info variable which contains runtime
 * parameters. In case the template parameters fe_degree and n_q_points_1d
 * contain valid information (i.e. fe_degree>-1 and n_q_points_1d>0), we simply
 * pass these values to the respective template specializations.
 * Otherwise, we perform a runtime matching of the runtime parameters to find
 * the correct specialization. This matching currently supports
 * $0\leq fe\_degree \leq 9$ and $degree+1\leq n\_q\_points\_1d\leq
 * fe\_degree+2$.
 */
template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number>
struct SelectEvaluator
{
  /**
   * We enable a transformation to collocation for derivatives if it gives
   * correct results (first condition), if it is the most efficient choice in
   * terms of operation counts (second condition) and if we were able to
   * initialize the fields in shape_info.templates.h from the polynomials
   * (third condition).
   */
  static constexpr bool         use_collocation =
    n_q_points_1d > fe_degree &&n_q_points_1d <= 3 * fe_degree / 2 + 1 &&
    n_q_points_1d < 200;

  /**
   * Chooses an appropriate evaluation strategy for the evaluate function, i.e.
   * this calls internal::FEEvaluationImpl::evaluate(),
   * internal::FEEvaluationImplCollocation::evaluate() or
   * internal::FEEvaluationImplTransformToCollocation::evaluate() with
   * appropriate template parameters.
   */
  static void
  evaluate(const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
           Number *   values_dofs_actual,
           Number *   values_quad,
           Number *   gradients_quad,
           Number *   hessians_quad,
           Number *   scratch_data,
           const bool evaluate_values,
           const bool evaluate_gradients,
           const bool evaluate_hessians);

  /**
   * Chooses an appropriate evaluation strategy for the integrate function, i.e.
   * this calls internal::FEEvaluationImpl::integrate(),
   * internal::FEEvaluationImplCollocation::integrate() or
   * internal::FEEvaluationImplTransformToCollocation::integrate() with
   * appropriate template parameters.
   */
  static void
  integrate(const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
            Number *   values_dofs_actual,
            Number *   values_quad,
            Number *   gradients_quad,
            Number *   scratch_data,
            const bool integrate_values,
            const bool integrate_gradients,
            const bool sum_into_values_array = false);
};

/**
 * This specialization chooses an appropriate evaluation strategy if we
 * don't know the correct template parameters at compile time. Instead
 * the selection is done based on the shape_info variable which contains
 * the relevant runtime parameters.
 * In case these parameters do not satisfy
 * $0\leq fe\_degree \leq 9$ and
 * $degree+1\leq n\_q\_points\_1d\leq fe\_degree+2$, a non-optimized fallback
 * is used.
 */
template <int dim, int n_q_points_1d, int n_components, typename Number>
struct SelectEvaluator<dim, -1, n_q_points_1d, n_components, Number>
{
  /**
   * Based on the run time parameters stored in @p shape_info this function
   * chooses an appropriate evaluation strategy for the integrate function, i.e.
   * this calls internal::FEEvaluationImpl::evaluate(),
   * internal::FEEvaluationImplCollocation::evaluate() or
   * internal::FEEvaluationImplTransformToCollocation::evaluate() with
   * appropriate template parameters.
   */
  static void
  evaluate(const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
           Number *   values_dofs_actual,
           Number *   values_quad,
           Number *   gradients_quad,
           Number *   hessians_quad,
           Number *   scratch_data,
           const bool evaluate_values,
           const bool evaluate_gradients,
           const bool evaluate_hessians);

  /**
   * Based on the run time parameters stored in @p shape_info this function
   * chooses an appropriate evaluation strategy for the integrate function, i.e.
   * this calls internal::FEEvaluationImpl::integrate(),
   * internal::FEEvaluationImplCollocation::integrate() or
   * internal::FEEvaluationImplTransformToCollocation::integrate() with
   * appropriate template parameters.
   */
  static void
  integrate(const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
            Number *   values_dofs_actual,
            Number *   values_quad,
            Number *   gradients_quad,
            Number *   scratch_data,
            const bool integrate_values,
            const bool integrate_gradients,
            const bool sum_into_values_array = false);
};

//----------------------Implementation for SelectEvaluator---------------------
#ifndef DOXYGEN

template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number>
inline void
SelectEvaluator<dim, fe_degree, n_q_points_1d, n_components, Number>::evaluate(
  const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
  Number *                                                values_dofs_actual,
  Number *                                                values_quad,
  Number *                                                gradients_quad,
  Number *                                                hessians_quad,
  Number *                                                scratch_data,
  const bool                                              evaluate_values,
  const bool                                              evaluate_gradients,
  const bool                                              evaluate_hessians)
{
  Assert(fe_degree >= 0 && n_q_points_1d > 0, ExcInternalError());

  if (fe_degree + 1 == n_q_points_1d &&
      shape_info.element_type ==
        internal::MatrixFreeFunctions::tensor_symmetric_collocation)
    {
      internal::
        FEEvaluationImplCollocation<dim, fe_degree, n_components, Number>::
          evaluate(shape_info,
                   values_dofs_actual,
                   values_quad,
                   gradients_quad,
                   hessians_quad,
                   scratch_data,
                   evaluate_values,
                   evaluate_gradients,
                   evaluate_hessians);
    }
  // '<=' on type means tensor_symmetric or tensor_symmetric_hermite, see
  // shape_info.h for more details
  else if (use_collocation && shape_info.element_type <=
                                internal::MatrixFreeFunctions::tensor_symmetric)
    {
      internal::FEEvaluationImplTransformToCollocation<
        dim,
        fe_degree,
        n_q_points_1d,
        n_components,
        Number>::evaluate(shape_info,
                          values_dofs_actual,
                          values_quad,
                          gradients_quad,
                          hessians_quad,
                          scratch_data,
                          evaluate_values,
                          evaluate_gradients,
                          evaluate_hessians);
    }
  else if (shape_info.element_type <=
           internal::MatrixFreeFunctions::tensor_symmetric)
    {
      internal::FEEvaluationImpl<
        internal::MatrixFreeFunctions::tensor_symmetric,
        dim,
        fe_degree,
        n_q_points_1d,
        n_components,
        Number>::evaluate(shape_info,
                          values_dofs_actual,
                          values_quad,
                          gradients_quad,
                          hessians_quad,
                          scratch_data,
                          evaluate_values,
                          evaluate_gradients,
                          evaluate_hessians);
    }
  else if (shape_info.element_type ==
           internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0)
    {
      internal::FEEvaluationImpl<
        internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0,
        dim,
        fe_degree,
        n_q_points_1d,
        n_components,
        Number>::evaluate(shape_info,
                          values_dofs_actual,
                          values_quad,
                          gradients_quad,
                          hessians_quad,
                          scratch_data,
                          evaluate_values,
                          evaluate_gradients,
                          evaluate_hessians);
    }
  else if (shape_info.element_type ==
           internal::MatrixFreeFunctions::truncated_tensor)
    {
      internal::FEEvaluationImpl<
        internal::MatrixFreeFunctions::truncated_tensor,
        dim,
        fe_degree,
        n_q_points_1d,
        n_components,
        Number>::evaluate(shape_info,
                          values_dofs_actual,
                          values_quad,
                          gradients_quad,
                          hessians_quad,
                          scratch_data,
                          evaluate_values,
                          evaluate_gradients,
                          evaluate_hessians);
    }
  else if (shape_info.element_type ==
           internal::MatrixFreeFunctions::tensor_general)
    {
      internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_general,
                                 dim,
                                 fe_degree,
                                 n_q_points_1d,
                                 n_components,
                                 Number>::evaluate(shape_info,
                                                   values_dofs_actual,
                                                   values_quad,
                                                   gradients_quad,
                                                   hessians_quad,
                                                   scratch_data,
                                                   evaluate_values,
                                                   evaluate_gradients,
                                                   evaluate_hessians);
    }
  else
    AssertThrow(false, ExcNotImplemented());
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number>
inline void
SelectEvaluator<dim, fe_degree, n_q_points_1d, n_components, Number>::integrate(
  const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
  Number *                                                values_dofs_actual,
  Number *                                                values_quad,
  Number *                                                gradients_quad,
  Number *                                                scratch_data,
  const bool                                              integrate_values,
  const bool                                              integrate_gradients,
  const bool                                              sum_into_values_array)
{
  Assert(fe_degree >= 0 && n_q_points_1d > 0, ExcInternalError());

  if (fe_degree + 1 == n_q_points_1d &&
      shape_info.element_type ==
        internal::MatrixFreeFunctions::tensor_symmetric_collocation)
    {
      internal::
        FEEvaluationImplCollocation<dim, fe_degree, n_components, Number>::
          integrate(shape_info,
                    values_dofs_actual,
                    values_quad,
                    gradients_quad,
                    scratch_data,
                    integrate_values,
                    integrate_gradients,
                    sum_into_values_array);
    }
  // '<=' on type means tensor_symmetric or tensor_symmetric_hermite, see
  // shape_info.h for more details
  else if (use_collocation && shape_info.element_type <=
                                internal::MatrixFreeFunctions::tensor_symmetric)
    {
      internal::FEEvaluationImplTransformToCollocation<
        dim,
        fe_degree,
        n_q_points_1d,
        n_components,
        Number>::integrate(shape_info,
                           values_dofs_actual,
                           values_quad,
                           gradients_quad,
                           scratch_data,
                           integrate_values,
                           integrate_gradients,
                           sum_into_values_array);
    }
  else if (shape_info.element_type <=
           internal::MatrixFreeFunctions::tensor_symmetric)
    {
      internal::FEEvaluationImpl<
        internal::MatrixFreeFunctions::tensor_symmetric,
        dim,
        fe_degree,
        n_q_points_1d,
        n_components,
        Number>::integrate(shape_info,
                           values_dofs_actual,
                           values_quad,
                           gradients_quad,
                           scratch_data,
                           integrate_values,
                           integrate_gradients,
                           sum_into_values_array);
    }
  else if (shape_info.element_type ==
           internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0)
    {
      internal::FEEvaluationImpl<
        internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0,
        dim,
        fe_degree,
        n_q_points_1d,
        n_components,
        Number>::integrate(shape_info,
                           values_dofs_actual,
                           values_quad,
                           gradients_quad,
                           scratch_data,
                           integrate_values,
                           integrate_gradients,
                           sum_into_values_array);
    }
  else if (shape_info.element_type ==
           internal::MatrixFreeFunctions::truncated_tensor)
    {
      internal::FEEvaluationImpl<
        internal::MatrixFreeFunctions::truncated_tensor,
        dim,
        fe_degree,
        n_q_points_1d,
        n_components,
        Number>::integrate(shape_info,
                           values_dofs_actual,
                           values_quad,
                           gradients_quad,
                           scratch_data,
                           integrate_values,
                           integrate_gradients,
                           sum_into_values_array);
    }
  else if (shape_info.element_type ==
           internal::MatrixFreeFunctions::tensor_general)
    {
      internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_general,
                                 dim,
                                 fe_degree,
                                 n_q_points_1d,
                                 n_components,
                                 Number>::integrate(shape_info,
                                                    values_dofs_actual,
                                                    values_quad,
                                                    gradients_quad,
                                                    scratch_data,
                                                    integrate_values,
                                                    integrate_gradients,
                                                    sum_into_values_array);
    }
  else
    AssertThrow(false, ExcNotImplemented());
}



template <int dim, int dummy, int n_components, typename Number>
inline void
SelectEvaluator<dim, -1, dummy, n_components, Number>::evaluate(
  const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
  Number *                                                values_dofs_actual,
  Number *                                                values_quad,
  Number *                                                gradients_quad,
  Number *                                                hessians_quad,
  Number *                                                scratch_data,
  const bool                                              evaluate_values,
  const bool                                              evaluate_gradients,
  const bool                                              evaluate_hessians)
{
  if (shape_info.element_type ==
      internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0)
    {
      internal::FEEvaluationImpl<
        internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0,
        dim,
        -1,
        0,
        n_components,
        Number>::evaluate(shape_info,
                          values_dofs_actual,
                          values_quad,
                          gradients_quad,
                          hessians_quad,
                          scratch_data,
                          evaluate_values,
                          evaluate_gradients,
                          evaluate_hessians);
    }
  else if (shape_info.element_type ==
           internal::MatrixFreeFunctions::truncated_tensor)
    {
      internal::FEEvaluationImpl<
        internal::MatrixFreeFunctions::truncated_tensor,
        dim,
        -1,
        0,
        n_components,
        Number>::evaluate(shape_info,
                          values_dofs_actual,
                          values_quad,
                          gradients_quad,
                          hessians_quad,
                          scratch_data,
                          evaluate_values,
                          evaluate_gradients,
                          evaluate_hessians);
    }
  else if (shape_info.element_type ==
           internal::MatrixFreeFunctions::tensor_general)
    internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_general,
                               dim,
                               -1,
                               0,
                               n_components,
                               Number>::evaluate(shape_info,
                                                 values_dofs_actual,
                                                 values_quad,
                                                 gradients_quad,
                                                 hessians_quad,
                                                 scratch_data,
                                                 evaluate_values,
                                                 evaluate_gradients,
                                                 evaluate_hessians);
  else
    internal::EvaluationSelectorImplementation::
      symmetric_selector_evaluate<dim, n_components, Number>(shape_info,
                                                             values_dofs_actual,
                                                             values_quad,
                                                             gradients_quad,
                                                             hessians_quad,
                                                             scratch_data,
                                                             evaluate_values,
                                                             evaluate_gradients,
                                                             evaluate_hessians);
}



template <int dim, int dummy, int n_components, typename Number>
inline void
SelectEvaluator<dim, -1, dummy, n_components, Number>::integrate(
  const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
  Number *                                                values_dofs_actual,
  Number *                                                values_quad,
  Number *                                                gradients_quad,
  Number *                                                scratch_data,
  const bool                                              integrate_values,
  const bool                                              integrate_gradients,
  const bool                                              sum_into_values_array)
{
  if (shape_info.element_type ==
      internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0)
    {
      internal::FEEvaluationImpl<
        internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0,
        dim,
        -1,
        0,
        n_components,
        Number>::integrate(shape_info,
                           values_dofs_actual,
                           values_quad,
                           gradients_quad,
                           scratch_data,
                           integrate_values,
                           integrate_gradients,
                           sum_into_values_array);
    }
  else if (shape_info.element_type ==
           internal::MatrixFreeFunctions::truncated_tensor)
    {
      internal::FEEvaluationImpl<
        internal::MatrixFreeFunctions::truncated_tensor,
        dim,
        -1,
        0,
        n_components,
        Number>::integrate(shape_info,
                           values_dofs_actual,
                           values_quad,
                           gradients_quad,
                           scratch_data,
                           integrate_values,
                           integrate_gradients,
                           sum_into_values_array);
    }
  else if (shape_info.element_type ==
           internal::MatrixFreeFunctions::tensor_general)
    internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_general,
                               dim,
                               -1,
                               0,
                               n_components,
                               Number>::integrate(shape_info,
                                                  values_dofs_actual,
                                                  values_quad,
                                                  gradients_quad,
                                                  scratch_data,
                                                  integrate_values,
                                                  integrate_gradients,
                                                  sum_into_values_array);
  else
    internal::EvaluationSelectorImplementation::
      symmetric_selector_integrate<dim, n_components, Number>(
        shape_info,
        values_dofs_actual,
        values_quad,
        gradients_quad,
        scratch_data,
        integrate_values,
        integrate_gradients,
        sum_into_values_array);
}
#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
