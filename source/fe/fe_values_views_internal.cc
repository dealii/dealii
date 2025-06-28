// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/array_view.h>
#include <deal.II/base/numbers.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/fe/fe_values_views_internal.h>

#ifdef DEAL_II_WITH_ADOLC
#  include <adolc/adouble.h>
#  include <adolc/adtl.h>
#endif

#include <type_traits>

DEAL_II_NAMESPACE_OPEN

namespace FEValuesViews
{
  namespace internal
  {
    namespace
    {
      // Check to see if a DoF value is zero, implying that subsequent
      // operations with the value have no effect.
      template <typename Number, typename T = void>
      struct CheckForZero
      {
        static bool
        value(const Number &value)
        {
          return value == dealii::internal::NumberType<Number>::value(0.0);
        }
      };

      // For auto-differentiable numbers, the fact that a DoF value is zero
      // does not imply that its derivatives are zero as well. So we
      // can't filter by value for these number types.
      // Note that we also want to avoid actually checking the value itself,
      // since some AD numbers are not contextually convertible to booleans.
      template <typename Number>
      struct CheckForZero<
        Number,
        std::enable_if_t<Differentiation::AD::is_ad_number<Number>::value>>
      {
        static bool
        value(const Number & /*value*/)
        {
          return false;
        }
      };
    } // namespace

    template <int dim, int spacedim, typename Number>
    void
    do_function_values(
      const ArrayView<const Number> &dof_values,
      const Table<2, double>        &shape_values,
      const std::vector<typename Scalar<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename ProductType<Number, double>::type> &values)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = values.size();

      std::fill(values.begin(),
                values.end(),
                dealii::internal::NumberType<Number>::value(0.0));

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        if (shape_function_data[shape_function]
              .is_nonzero_shape_function_component)
          {
            const Number &value = dof_values[shape_function];
            // For auto-differentiable numbers, the fact that a DoF value is
            // zero does not imply that its derivatives are zero as well. So we
            // can't filter by value for these number types.
            if (CheckForZero<Number>::value(value) == true)
              continue;

            const double *shape_value_ptr =
              &shape_values(shape_function_data[shape_function].row_index, 0);
            for (unsigned int q_point = 0; q_point < n_quadrature_points;
                 ++q_point, ++shape_value_ptr)
              values[q_point] += value * (*shape_value_ptr);
          }
    }



    template <int order, int dim, int spacedim, typename Number>
    void
    do_function_derivatives(
      const ArrayView<const Number>                   &dof_values,
      const Table<2, dealii::Tensor<order, spacedim>> &shape_derivatives,
      const std::vector<typename Scalar<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<
        typename ProductType<Number, dealii::Tensor<order, spacedim>>::type>
        &derivatives)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = derivatives.size();

      std::fill(
        derivatives.begin(),
        derivatives.end(),
        typename ProductType<Number, dealii::Tensor<order, spacedim>>::type());

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        if (shape_function_data[shape_function]
              .is_nonzero_shape_function_component)
          {
            const Number &value = dof_values[shape_function];
            // For auto-differentiable numbers, the fact that a DoF value is
            // zero does not imply that its derivatives are zero as well. So we
            // can't filter by value for these number types.
            if (CheckForZero<Number>::value(value) == true)
              continue;

            const dealii::Tensor<order, spacedim> *shape_derivative_ptr =
              &shape_derivatives[shape_function_data[shape_function].row_index]
                                [0];
            for (unsigned int q_point = 0; q_point < n_quadrature_points;
                 ++q_point)
              derivatives[q_point] += value * (*shape_derivative_ptr++);
          }
    }



    template <int dim, int spacedim, typename Number>
    void
    do_function_laplacians(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<2, spacedim>> &shape_hessians,
      const std::vector<typename Scalar<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename Scalar<dim, spacedim>::
                    template solution_laplacian_type<Number>> &laplacians)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = laplacians.size();

      std::fill(
        laplacians.begin(),
        laplacians.end(),
        typename Scalar<dim,
                        spacedim>::template solution_laplacian_type<Number>());

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        if (shape_function_data[shape_function]
              .is_nonzero_shape_function_component)
          {
            const Number &value = dof_values[shape_function];
            // For auto-differentiable numbers, the fact that a DoF value is
            // zero does not imply that its derivatives are zero as well. So we
            // can't filter by value for these number types.
            if (CheckForZero<Number>::value(value) == true)
              continue;

            const dealii::Tensor<2, spacedim> *shape_hessian_ptr =
              &shape_hessians[shape_function_data[shape_function].row_index][0];
            for (unsigned int q_point = 0; q_point < n_quadrature_points;
                 ++q_point)
              laplacians[q_point] += value * trace(*shape_hessian_ptr++);
          }
    }



    // ----------------------------- vector part ---------------------------

    template <int dim, int spacedim, typename Number>
    void
    do_function_values(
      const ArrayView<const Number> &dof_values,
      const Table<2, double>        &shape_values,
      const std::vector<typename Vector<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<
        typename ProductType<Number, dealii::Tensor<1, spacedim>>::type>
        &values)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = values.size();

      std::fill(
        values.begin(),
        values.end(),
        typename ProductType<Number, dealii::Tensor<1, spacedim>>::type());

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        {
          const int snc =
            shape_function_data[shape_function].single_nonzero_component;

          if (snc == -2)
            // shape function is zero for the selected components
            continue;

          const Number &value = dof_values[shape_function];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const unsigned int comp = shape_function_data[shape_function]
                                          .single_nonzero_component_index;
              const double *shape_value_ptr = &shape_values(snc, 0);
              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point, ++shape_value_ptr)
                values[q_point][comp] += value * (*shape_value_ptr);
            }
          else
            for (unsigned int d = 0; d < spacedim; ++d)
              if (shape_function_data[shape_function]
                    .is_nonzero_shape_function_component[d])
                {
                  const double *shape_value_ptr = &shape_values(
                    shape_function_data[shape_function].row_index[d], 0);
                  for (unsigned int q_point = 0; q_point < n_quadrature_points;
                       ++q_point, ++shape_value_ptr)
                    values[q_point][d] += value * (*shape_value_ptr);
                }
        }
    }



    template <int order, int dim, int spacedim, typename Number>
    void
    do_function_derivatives(
      const ArrayView<const Number>                   &dof_values,
      const Table<2, dealii::Tensor<order, spacedim>> &shape_derivatives,
      const std::vector<typename Vector<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<
        typename ProductType<Number, dealii::Tensor<order + 1, spacedim>>::type>
        &derivatives)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = derivatives.size();

      std::fill(
        derivatives.begin(),
        derivatives.end(),
        typename ProductType<Number,
                             dealii::Tensor<order + 1, spacedim>>::type());

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        {
          const int snc =
            shape_function_data[shape_function].single_nonzero_component;

          if (snc == -2)
            // shape function is zero for the selected components
            continue;

          const Number &value = dof_values[shape_function];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const unsigned int comp = shape_function_data[shape_function]
                                          .single_nonzero_component_index;
              const dealii::Tensor<order, spacedim> *shape_derivative_ptr =
                &shape_derivatives[snc][0];
              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point)
                derivatives[q_point][comp] += value * (*shape_derivative_ptr++);
            }
          else
            for (unsigned int d = 0; d < spacedim; ++d)
              if (shape_function_data[shape_function]
                    .is_nonzero_shape_function_component[d])
                {
                  const dealii::Tensor<order, spacedim> *shape_derivative_ptr =
                    &shape_derivatives[shape_function_data[shape_function]
                                         .row_index[d]][0];
                  for (unsigned int q_point = 0; q_point < n_quadrature_points;
                       ++q_point)
                    derivatives[q_point][d] +=
                      value * (*shape_derivative_ptr++);
                }
        }
    }



    template <int dim, int spacedim, typename Number>
    void
    do_function_symmetric_gradients(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<1, spacedim>> &shape_gradients,
      const std::vector<typename Vector<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<
        typename ProductType<Number,
                             dealii::SymmetricTensor<2, spacedim>>::type>
        &symmetric_gradients)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = symmetric_gradients.size();

      std::fill(
        symmetric_gradients.begin(),
        symmetric_gradients.end(),
        typename ProductType<Number,
                             dealii::SymmetricTensor<2, spacedim>>::type());

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        {
          const int snc =
            shape_function_data[shape_function].single_nonzero_component;

          if (snc == -2)
            // shape function is zero for the selected components
            continue;

          const Number &value = dof_values[shape_function];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const unsigned int comp = shape_function_data[shape_function]
                                          .single_nonzero_component_index;
              const dealii::Tensor<1, spacedim> *shape_gradient_ptr =
                &shape_gradients[snc][0];
              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point)
                {
                  for (unsigned int d = 0; d < dim; ++d)
                    symmetric_gradients[q_point][comp][d] +=
                      0.5 * value * (*shape_gradient_ptr)[d];
                  symmetric_gradients[q_point][comp][comp] +=
                    0.5 * value * (*shape_gradient_ptr++)[comp];
                }
            }
          else
            for (unsigned int q_point = 0; q_point < n_quadrature_points;
                 ++q_point)
              {
                typename ProductType<Number, dealii::Tensor<2, spacedim>>::type
                  grad;
                for (unsigned int d = 0; d < spacedim; ++d)
                  if (shape_function_data[shape_function]
                        .is_nonzero_shape_function_component[d])
                    grad[d] =
                      value *
                      shape_gradients[shape_function_data[shape_function]
                                        .row_index[d]][q_point];
                symmetric_gradients[q_point] += symmetrize(grad);
              }
        }
    }



    template <int dim, int spacedim, typename Number>
    void
    do_function_divergences(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<1, spacedim>> &shape_gradients,
      const std::vector<typename Vector<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename Vector<dim, spacedim>::
                    template solution_divergence_type<Number>> &divergences)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = divergences.size();

      std::fill(
        divergences.begin(),
        divergences.end(),
        typename Vector<dim,
                        spacedim>::template solution_divergence_type<Number>());

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        {
          const int snc =
            shape_function_data[shape_function].single_nonzero_component;

          if (snc == -2)
            // shape function is zero for the selected components
            continue;

          const Number &value = dof_values[shape_function];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const unsigned int comp = shape_function_data[shape_function]
                                          .single_nonzero_component_index;
              const dealii::Tensor<1, spacedim> *shape_gradient_ptr =
                &shape_gradients[snc][0];
              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point)
                divergences[q_point] += value * (*shape_gradient_ptr++)[comp];
            }
          else
            for (unsigned int d = 0; d < spacedim; ++d)
              if (shape_function_data[shape_function]
                    .is_nonzero_shape_function_component[d])
                {
                  const dealii::Tensor<1, spacedim> *shape_gradient_ptr =
                    &shape_gradients[shape_function_data[shape_function]
                                       .row_index[d]][0];
                  for (unsigned int q_point = 0; q_point < n_quadrature_points;
                       ++q_point)
                    divergences[q_point] += value * (*shape_gradient_ptr++)[d];
                }
        }
    }



    template <int dim, int spacedim, typename Number>
    void
    do_function_curls(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<1, spacedim>> &shape_gradients,
      const std::vector<typename Vector<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename ProductType<
        Number,
        typename dealii::internal::CurlType<spacedim>::type>::type> &curls)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = curls.size();

      std::fill(curls.begin(),
                curls.end(),
                typename ProductType<
                  Number,
                  typename dealii::internal::CurlType<spacedim>::type>::type());

      switch (spacedim)
        {
          case 1:
            {
              Assert(false,
                     ExcMessage(
                       "Computing the curl in 1d is not a useful operation"));
              break;
            }

          case 2:
            {
              for (unsigned int shape_function = 0;
                   shape_function < dofs_per_cell;
                   ++shape_function)
                {
                  const int snc = shape_function_data[shape_function]
                                    .single_nonzero_component;

                  if (snc == -2)
                    // shape function is zero for the selected components
                    continue;

                  const Number &value = dof_values[shape_function];
                  // For auto-differentiable numbers, the fact that a DoF value
                  // is zero does not imply that its derivatives are zero as
                  // well. So we can't filter by value for these number types.
                  if (CheckForZero<Number>::value(value) == true)
                    continue;

                  if (snc != -1)
                    {
                      const dealii::Tensor<1, spacedim> *shape_gradient_ptr =
                        &shape_gradients[snc][0];

                      Assert(shape_function_data[shape_function]
                                 .single_nonzero_component >= 0,
                             ExcInternalError());
                      // we're in 2d, so the formula for the curl is simple:
                      if (shape_function_data[shape_function]
                            .single_nonzero_component_index == 0)
                        for (unsigned int q_point = 0;
                             q_point < n_quadrature_points;
                             ++q_point)
                          curls[q_point][0] -=
                            value * (*shape_gradient_ptr++)[1];
                      else
                        for (unsigned int q_point = 0;
                             q_point < n_quadrature_points;
                             ++q_point)
                          curls[q_point][0] +=
                            value * (*shape_gradient_ptr++)[0];
                    }
                  else
                    // we have multiple non-zero components in the shape
                    // functions. not all of them must necessarily be within the
                    // 2-component window this FEValuesViews::Vector object
                    // considers, however.
                    {
                      if (shape_function_data[shape_function]
                            .is_nonzero_shape_function_component[0])
                        {
                          const dealii::Tensor<1,
                                               spacedim> *shape_gradient_ptr =
                            &shape_gradients[shape_function_data[shape_function]
                                               .row_index[0]][0];

                          for (unsigned int q_point = 0;
                               q_point < n_quadrature_points;
                               ++q_point)
                            curls[q_point][0] -=
                              value * (*shape_gradient_ptr++)[1];
                        }

                      if (shape_function_data[shape_function]
                            .is_nonzero_shape_function_component[1])
                        {
                          const dealii::Tensor<1,
                                               spacedim> *shape_gradient_ptr =
                            &shape_gradients[shape_function_data[shape_function]
                                               .row_index[1]][0];

                          for (unsigned int q_point = 0;
                               q_point < n_quadrature_points;
                               ++q_point)
                            curls[q_point][0] +=
                              value * (*shape_gradient_ptr++)[0];
                        }
                    }
                }
              break;
            }

          case 3:
            {
              for (unsigned int shape_function = 0;
                   shape_function < dofs_per_cell;
                   ++shape_function)
                {
                  const int snc = shape_function_data[shape_function]
                                    .single_nonzero_component;

                  if (snc == -2)
                    // shape function is zero for the selected components
                    continue;

                  const Number &value = dof_values[shape_function];
                  // For auto-differentiable numbers, the fact that a DoF value
                  // is zero does not imply that its derivatives are zero as
                  // well. So we can't filter by value for these number types.
                  if (CheckForZero<Number>::value(value) == true)
                    continue;

                  if (snc != -1)
                    {
                      const dealii::Tensor<1, spacedim> *shape_gradient_ptr =
                        &shape_gradients[snc][0];

                      switch (shape_function_data[shape_function]
                                .single_nonzero_component_index)
                        {
                          case 0:
                            {
                              for (unsigned int q_point = 0;
                                   q_point < n_quadrature_points;
                                   ++q_point)
                                {
                                  curls[q_point][1] +=
                                    value * (*shape_gradient_ptr)[2];
                                  curls[q_point][2] -=
                                    value * (*shape_gradient_ptr++)[1];
                                }

                              break;
                            }

                          case 1:
                            {
                              for (unsigned int q_point = 0;
                                   q_point < n_quadrature_points;
                                   ++q_point)
                                {
                                  curls[q_point][0] -=
                                    value * (*shape_gradient_ptr)[2];
                                  curls[q_point][2] +=
                                    value * (*shape_gradient_ptr++)[0];
                                }

                              break;
                            }

                          case 2:
                            {
                              for (unsigned int q_point = 0;
                                   q_point < n_quadrature_points;
                                   ++q_point)
                                {
                                  curls[q_point][0] +=
                                    value * (*shape_gradient_ptr)[1];
                                  curls[q_point][1] -=
                                    value * (*shape_gradient_ptr++)[0];
                                }
                              break;
                            }

                          default:
                            DEAL_II_ASSERT_UNREACHABLE();
                        }
                    }

                  else
                    // we have multiple non-zero components in the shape
                    // functions. not all of them must necessarily be within the
                    // 3-component window this FEValuesViews::Vector object
                    // considers, however.
                    {
                      if (shape_function_data[shape_function]
                            .is_nonzero_shape_function_component[0])
                        {
                          const dealii::Tensor<1,
                                               spacedim> *shape_gradient_ptr =
                            &shape_gradients[shape_function_data[shape_function]
                                               .row_index[0]][0];

                          for (unsigned int q_point = 0;
                               q_point < n_quadrature_points;
                               ++q_point)
                            {
                              curls[q_point][1] +=
                                value * (*shape_gradient_ptr)[2];
                              curls[q_point][2] -=
                                value * (*shape_gradient_ptr++)[1];
                            }
                        }

                      if (shape_function_data[shape_function]
                            .is_nonzero_shape_function_component[1])
                        {
                          const dealii::Tensor<1,
                                               spacedim> *shape_gradient_ptr =
                            &shape_gradients[shape_function_data[shape_function]
                                               .row_index[1]][0];

                          for (unsigned int q_point = 0;
                               q_point < n_quadrature_points;
                               ++q_point)
                            {
                              curls[q_point][0] -=
                                value * (*shape_gradient_ptr)[2];
                              curls[q_point][2] +=
                                value * (*shape_gradient_ptr++)[0];
                            }
                        }

                      if (shape_function_data[shape_function]
                            .is_nonzero_shape_function_component[2])
                        {
                          const dealii::Tensor<1,
                                               spacedim> *shape_gradient_ptr =
                            &shape_gradients[shape_function_data[shape_function]
                                               .row_index[2]][0];

                          for (unsigned int q_point = 0;
                               q_point < n_quadrature_points;
                               ++q_point)
                            {
                              curls[q_point][0] +=
                                value * (*shape_gradient_ptr)[1];
                              curls[q_point][1] -=
                                value * (*shape_gradient_ptr++)[0];
                            }
                        }
                    }
                }
            }
        }
    }



    template <int dim, int spacedim, typename Number>
    void
    do_function_laplacians(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<2, spacedim>> &shape_hessians,
      const std::vector<typename Vector<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename Vector<dim, spacedim>::
                    template solution_laplacian_type<Number>> &laplacians)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = laplacians.size();

      std::fill(
        laplacians.begin(),
        laplacians.end(),
        typename Vector<dim,
                        spacedim>::template solution_laplacian_type<Number>());

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        {
          const int snc =
            shape_function_data[shape_function].single_nonzero_component;

          if (snc == -2)
            // shape function is zero for the selected components
            continue;

          const Number &value = dof_values[shape_function];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const unsigned int comp = shape_function_data[shape_function]
                                          .single_nonzero_component_index;
              const dealii::Tensor<2, spacedim> *shape_hessian_ptr =
                &shape_hessians[snc][0];
              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point)
                laplacians[q_point][comp] +=
                  value * trace(*shape_hessian_ptr++);
            }
          else
            for (unsigned int d = 0; d < spacedim; ++d)
              if (shape_function_data[shape_function]
                    .is_nonzero_shape_function_component[d])
                {
                  const dealii::Tensor<2, spacedim> *shape_hessian_ptr =
                    &shape_hessians[shape_function_data[shape_function]
                                      .row_index[d]][0];
                  for (unsigned int q_point = 0; q_point < n_quadrature_points;
                       ++q_point)
                    laplacians[q_point][d] +=
                      value * trace(*shape_hessian_ptr++);
                }
        }
    }



    // ---------------------- symmetric tensor part ------------------------

    template <int dim, int spacedim, typename Number>
    void
    do_function_values(
      const ArrayView<const Number>  &dof_values,
      const dealii::Table<2, double> &shape_values,
      const std::vector<
        typename SymmetricTensor<2, dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<
        typename ProductType<Number,
                             dealii::SymmetricTensor<2, spacedim>>::type>
        &values)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = values.size();

      std::fill(
        values.begin(),
        values.end(),
        typename ProductType<Number,
                             dealii::SymmetricTensor<2, spacedim>>::type());

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        {
          const int snc =
            shape_function_data[shape_function].single_nonzero_component;

          if (snc == -2)
            // shape function is zero for the selected components
            continue;

          const Number &value = dof_values[shape_function];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const TableIndices<2> comp = dealii::
                SymmetricTensor<2, spacedim>::unrolled_to_component_indices(
                  shape_function_data[shape_function]
                    .single_nonzero_component_index);
              const double *shape_value_ptr = &shape_values(snc, 0);
              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point, ++shape_value_ptr)
                values[q_point][comp] += value * (*shape_value_ptr);
            }
          else
            for (unsigned int d = 0;
                 d <
                 dealii::SymmetricTensor<2, spacedim>::n_independent_components;
                 ++d)
              if (shape_function_data[shape_function]
                    .is_nonzero_shape_function_component[d])
                {
                  const TableIndices<2> comp =
                    dealii::SymmetricTensor<2, spacedim>::
                      unrolled_to_component_indices(d);
                  const double *shape_value_ptr = &shape_values(
                    shape_function_data[shape_function].row_index[d], 0);
                  for (unsigned int q_point = 0; q_point < n_quadrature_points;
                       ++q_point, ++shape_value_ptr)
                    values[q_point][comp] += value * (*shape_value_ptr);
                }
        }
    }



    template <int dim, int spacedim, typename Number>
    void
    do_function_divergences(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<1, spacedim>> &shape_gradients,
      const std::vector<
        typename SymmetricTensor<2, dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename SymmetricTensor<2, dim, spacedim>::
                    template solution_divergence_type<Number>> &divergences)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = divergences.size();

      std::fill(divergences.begin(),
                divergences.end(),
                typename SymmetricTensor<2, dim, spacedim>::
                  template solution_divergence_type<Number>());

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        {
          const int snc =
            shape_function_data[shape_function].single_nonzero_component;

          if (snc == -2)
            // shape function is zero for the selected components
            continue;

          const Number &value = dof_values[shape_function];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const unsigned int comp = shape_function_data[shape_function]
                                          .single_nonzero_component_index;

              const dealii::Tensor<1, spacedim> *shape_gradient_ptr =
                &shape_gradients[snc][0];

              const unsigned int ii = dealii::SymmetricTensor<2, spacedim>::
                unrolled_to_component_indices(comp)[0];
              const unsigned int jj = dealii::SymmetricTensor<2, spacedim>::
                unrolled_to_component_indices(comp)[1];

              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point, ++shape_gradient_ptr)
                {
                  divergences[q_point][ii] += value * (*shape_gradient_ptr)[jj];

                  if (ii != jj)
                    divergences[q_point][jj] +=
                      value * (*shape_gradient_ptr)[ii];
                }
            }
          else
            {
              for (unsigned int d = 0;
                   d <
                   dealii::SymmetricTensor<2,
                                           spacedim>::n_independent_components;
                   ++d)
                if (shape_function_data[shape_function]
                      .is_nonzero_shape_function_component[d])
                  {
                    DEAL_II_NOT_IMPLEMENTED();

                    // the following implementation needs to be looked over -- I
                    // think it can't be right, because we are in a case where
                    // there is no single nonzero component
                    //
                    // the following is not implemented! we need to consider the
                    // interplay between multiple non-zero entries in shape
                    // function and the representation as a symmetric
                    // second-order tensor
                    const unsigned int comp =
                      shape_function_data[shape_function]
                        .single_nonzero_component_index;

                    const dealii::Tensor<1, spacedim> *shape_gradient_ptr =
                      &shape_gradients[shape_function_data[shape_function]
                                         .row_index[d]][0];
                    for (unsigned int q_point = 0;
                         q_point < n_quadrature_points;
                         ++q_point, ++shape_gradient_ptr)
                      {
                        for (unsigned int j = 0; j < spacedim;
                             ++j, ++shape_gradient_ptr)
                          {
                            const unsigned int vector_component =
                              dealii::SymmetricTensor<2, spacedim>::
                                component_to_unrolled_index(
                                  TableIndices<2>(comp, j));
                            divergences[q_point][vector_component] +=
                              value * (*shape_gradient_ptr)[j];
                          }
                      }
                  }
            }
        }
    }

    // ---------------------- non-symmetric tensor part ------------------------

    template <int dim, int spacedim, typename Number>
    void
    do_function_values(
      const ArrayView<const Number>  &dof_values,
      const dealii::Table<2, double> &shape_values,
      const std::vector<typename Tensor<2, dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<
        typename ProductType<Number, dealii::Tensor<2, spacedim>>::type>
        &values)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = values.size();

      std::fill(
        values.begin(),
        values.end(),
        typename ProductType<Number, dealii::Tensor<2, spacedim>>::type());

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        {
          const int snc =
            shape_function_data[shape_function].single_nonzero_component;

          if (snc == -2)
            // shape function is zero for the selected components
            continue;

          const Number &value = dof_values[shape_function];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const unsigned int comp = shape_function_data[shape_function]
                                          .single_nonzero_component_index;

              const TableIndices<2> indices =
                dealii::Tensor<2, spacedim>::unrolled_to_component_indices(
                  comp);

              const double *shape_value_ptr = &shape_values(snc, 0);
              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point, ++shape_value_ptr)
                values[q_point][indices] += value * (*shape_value_ptr);
            }
          else
            for (unsigned int d = 0; d < dim * dim; ++d)
              if (shape_function_data[shape_function]
                    .is_nonzero_shape_function_component[d])
                {
                  const TableIndices<2> indices =
                    dealii::Tensor<2, spacedim>::unrolled_to_component_indices(
                      d);

                  const double *shape_value_ptr = &shape_values(
                    shape_function_data[shape_function].row_index[d], 0);
                  for (unsigned int q_point = 0; q_point < n_quadrature_points;
                       ++q_point, ++shape_value_ptr)
                    values[q_point][indices] += value * (*shape_value_ptr);
                }
        }
    }



    template <int dim, int spacedim, typename Number>
    void
    do_function_divergences(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<1, spacedim>> &shape_gradients,
      const std::vector<typename Tensor<2, dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename Tensor<2, dim, spacedim>::
                    template solution_divergence_type<Number>> &divergences)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = divergences.size();

      std::fill(
        divergences.begin(),
        divergences.end(),
        typename Tensor<2, dim, spacedim>::template solution_divergence_type<
          Number>());

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        {
          const int snc =
            shape_function_data[shape_function].single_nonzero_component;

          if (snc == -2)
            // shape function is zero for the selected components
            continue;

          const Number &value = dof_values[shape_function];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const unsigned int comp = shape_function_data[shape_function]
                                          .single_nonzero_component_index;

              const dealii::Tensor<1, spacedim> *shape_gradient_ptr =
                &shape_gradients[snc][0];

              const TableIndices<2> indices =
                dealii::Tensor<2, spacedim>::unrolled_to_component_indices(
                  comp);
              const unsigned int ii = indices[0];
              const unsigned int jj = indices[1];

              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point, ++shape_gradient_ptr)
                {
                  divergences[q_point][ii] += value * (*shape_gradient_ptr)[jj];
                }
            }
          else
            {
              for (unsigned int d = 0; d < dim * dim; ++d)
                if (shape_function_data[shape_function]
                      .is_nonzero_shape_function_component[d])
                  {
                    DEAL_II_NOT_IMPLEMENTED();
                  }
            }
        }
    }



    template <int dim, int spacedim, typename Number>
    void
    do_function_gradients(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<1, spacedim>> &shape_gradients,
      const std::vector<typename Tensor<2, dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename Tensor<2, dim, spacedim>::
                    template solution_gradient_type<Number>> &gradients)
    {
      const unsigned int dofs_per_cell       = dof_values.size();
      const unsigned int n_quadrature_points = gradients.size();

      std::fill(
        gradients.begin(),
        gradients.end(),
        typename Tensor<2, dim, spacedim>::template solution_gradient_type<
          Number>());

      for (unsigned int shape_function = 0; shape_function < dofs_per_cell;
           ++shape_function)
        {
          const int snc =
            shape_function_data[shape_function].single_nonzero_component;

          if (snc == -2)
            // shape function is zero for the selected components
            continue;

          const Number &value = dof_values[shape_function];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const unsigned int comp = shape_function_data[shape_function]
                                          .single_nonzero_component_index;

              const dealii::Tensor<1, spacedim> *shape_gradient_ptr =
                &shape_gradients[snc][0];

              const TableIndices<2> indices =
                dealii::Tensor<2, spacedim>::unrolled_to_component_indices(
                  comp);
              const unsigned int ii = indices[0];
              const unsigned int jj = indices[1];

              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point, ++shape_gradient_ptr)
                {
                  gradients[q_point][ii][jj] += value * (*shape_gradient_ptr);
                }
            }
          else
            {
              for (unsigned int d = 0; d < dim * dim; ++d)
                if (shape_function_data[shape_function]
                      .is_nonzero_shape_function_component[d])
                  {
                    DEAL_II_NOT_IMPLEMENTED();
                  }
            }
        }
    }
  } // end of namespace internal
} // namespace FEValuesViews



/*------------------------------- Explicit Instantiations -------------*/

#include "fe/fe_values_views_internal.inst"

DEAL_II_NAMESPACE_CLOSE
