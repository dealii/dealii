// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_portable_fe_evaluation_h
#define dealii_portable_fe_evaluation_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/portable_hanging_nodes_internal.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/portable_matrix_free.templates.h>
#include <deal.II/matrix_free/portable_tensor_product_kernels.h>

#include <Kokkos_Core.hpp>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace for portable capabilities
 */
namespace Portable
{
  /**
   * This class provides all the functions necessary to evaluate functions at
   * quadrature points and cell integrations. In functionality, this class is
   * similar to FEValues<dim>.
   *
   * This class has five template arguments:
   *
   * @tparam dim Dimension in which this class is to be used
   *
   * @tparam fe_degree Degree of the tensor prodict finite element with fe_degree+1
   * degrees of freedom per coordinate direction
   *
   * @tparam n_q_points_1d Number of points in the quadrature formular in 1d,
   * defaults to fe_degree+1
   *
   * @tparam n_components Number of vector components when solving a system of
   * PDEs. If the same operation is applied to several components of a PDE (e.g.
   * a vector Laplace equation), they can be applied simultaneously with one
   * call (and often more efficiently). Defaults to 1
   *
   * @tparam Number Number format, @p double or @p float. Defaults to @p
   * double.
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d = fe_degree + 1,
            int n_components_ = 1,
            typename Number   = double>
  class FEEvaluation
  {
  public:
    /**
     * An alias for scalar quantities.
     */
    using value_type = std::conditional_t<(n_components_ == 1),
                                          Number,
                                          Tensor<1, n_components_, Number>>;

    /**
     * An alias for vectorial quantities.
     */
    using gradient_type = std::conditional_t<
      n_components_ == 1,
      Tensor<1, dim, Number>,
      std::conditional_t<n_components_ == dim,
                         Tensor<2, dim, Number>,
                         Tensor<1, n_components_, Tensor<1, dim, Number>>>>;

    /**
     * An alias to kernel specific information.
     */
    using data_type = typename MatrixFree<dim, Number>::Data;

    /**
     * Dimension.
     */
    static constexpr unsigned int dimension = dim;

    /**
     * Number of components.
     */
    static constexpr unsigned int n_components = n_components_;

    /**
     * Number of quadrature points per cell.
     */
    static constexpr unsigned int n_q_points =
      Utilities::pow(n_q_points_1d, dim);

    /**
     * Number of tensor degrees of freedoms per cell.
     */
    static constexpr unsigned int tensor_dofs_per_cell =
      Utilities::pow(fe_degree + 1, dim);

    /**
     * Constructor.
     */
    DEAL_II_HOST_DEVICE
    FEEvaluation(const data_type *data, SharedData<dim, Number> *shdata);

    /**
     * For the vector @p src, read out the values on the degrees of freedom of
     * the current cell, and store them internally. Similar functionality as
     * the function DoFAccessor::get_interpolated_dof_values when no
     * constraints are present, but it also includes constraints from hanging
     * nodes, so once can see it as a similar function to
     * AffineConstraints::read_dof_values() as well.
     */
    DEAL_II_HOST_DEVICE void
    read_dof_values(const Number *src);

    /**
     * Take the value stored internally on dof values of the current cell and
     * sum them into the vector @p dst. The function also applies constraints
     * during the write operation. The functionality is hence similar to the
     * function AffineConstraints::distribute_local_to_global.
     */
    DEAL_II_HOST_DEVICE void
    distribute_local_to_global(Number *dst) const;

    /**
     * Evaluate the function values and the gradients of the FE function given
     * at the DoF values in the input vector at the quadrature points on the
     * unit cell. The function argument @p evaluate_flag specifies which parts
     * shall actually be computed. This function needs to be called before the
     * functions @p get_value() or @p get_gradient() give useful information.
     */
    DEAL_II_HOST_DEVICE void
    evaluate(const EvaluationFlags::EvaluationFlags evaluate_flag);

    /**
     * Evaluate the function values and the gradients of the FE function given
     * at the DoF values in the input vector at the quadrature points on the
     * unit cell. The function arguments specify which parts shall actually be
     * computed. This function needs to be called before the functions
     * @p get_value() or @p get_gradient() give useful information.
     */
    DEAL_II_DEPRECATED_WITH_COMMENT("Use the version taking EvaluationFlags.")
    DEAL_II_HOST_DEVICE
    void
    evaluate(const bool evaluate_val, const bool evaluate_grad);

    /**
     * This function takes the values and/or gradients that are stored on
     * quadrature points, tests them by all the basis functions/gradients on
     * the cell and performs the cell integration as specified by the
     * @p integration_flag argument.
     */
    DEAL_II_HOST_DEVICE void
    integrate(const EvaluationFlags::EvaluationFlags integration_flag);

    /**
     * This function takes the values and/or gradients that are stored on
     * quadrature points, tests them by all the basis functions/gradients on
     * the cell and performs the cell integration. The two function arguments
     * @p integrate_val and @p integrate_grad are used to enable/disable some
     * of the values or the gradients.
     */
    DEAL_II_DEPRECATED_WITH_COMMENT("Use the version taking EvaluationFlags.")
    DEAL_II_HOST_DEVICE
    void
    integrate(const bool integrate_val, const bool integrate_grad);

    /**
     * Same as above, except that the quadrature point is computed from thread
     * id.
     */
    DEAL_II_HOST_DEVICE value_type
    get_value(int q_point) const;

    /**
     * Same as above, except that the local dof index is computed from the
     * thread id.
     */
    DEAL_II_HOST_DEVICE value_type
    get_dof_value(int q_point) const;

    /**
     * Same as above, except that the quadrature point is computed from the
     * thread id.
     */
    DEAL_II_HOST_DEVICE void
    submit_value(const value_type &val_in, int q_point);

    /**
     * Same as above, except that the local dof index is computed from the
     * thread id.
     */
    DEAL_II_HOST_DEVICE void
    submit_dof_value(const value_type &val_in, int q_point);

    /**
     * Same as above, except that the quadrature point is computed from the
     * thread id.
     */
    DEAL_II_HOST_DEVICE gradient_type
    get_gradient(int q_point) const;

    /**
     * Same as above, except that the quadrature point is computed from the
     * thread id.
     */
    DEAL_II_HOST_DEVICE void
    submit_gradient(const gradient_type &grad_in, int q_point);

    // clang-format off
    /**
     * Same as above, except that the functor @p func only takes a single input
     * argument (fe_eval) and computes the quadrature point from the thread id.
     *
     * @p func needs to define
     * \code
     * DEAL_II_HOST_DEVICE void operator()(
     *   Portable::FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, Number> *fe_eval) const;
     * \endcode
     */
    // clang-format on
    template <typename Functor>
    DEAL_II_HOST_DEVICE void
    apply_for_each_quad_point(const Functor &func);

  private:
    const data_type         *data;
    SharedData<dim, Number> *shared_data;
    int                      cell_id;
  };



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    FEEvaluation(const data_type *data, SharedData<dim, Number> *shdata)
    : data(data)
    , shared_data(shdata)
    , cell_id(shared_data->team_member.league_rank())
  {}



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    read_dof_values(const Number *src)
  {
    // Populate the scratch memory
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(shared_data->team_member, n_q_points),
      [&](const int &i) {
        for (unsigned int c = 0; c < n_components_; ++c)
          shared_data->values(i, c) =
            src[data->local_to_global(cell_id, i + tensor_dofs_per_cell * c)];
      });
    shared_data->team_member.team_barrier();

    for (unsigned int c = 0; c < n_components_; ++c)
      {
        internal::resolve_hanging_nodes<dim, fe_degree, false, Number>(
          shared_data->team_member,
          data->constraint_weights,
          data->constraint_mask(cell_id),
          Kokkos::subview(shared_data->values, Kokkos::ALL, c));
      }
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    distribute_local_to_global(Number *dst) const
  {
    for (unsigned int c = 0; c < n_components_; ++c)
      {
        internal::resolve_hanging_nodes<dim, fe_degree, true, Number>(
          shared_data->team_member,
          data->constraint_weights,
          data->constraint_mask(cell_id),
          Kokkos::subview(shared_data->values, Kokkos::ALL, c));
      }

    if (data->use_coloring)
      {
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(shared_data->team_member, n_q_points),
          [&](const int &i) {
            for (unsigned int c = 0; c < n_components_; ++c)
              dst[data->local_to_global(cell_id,
                                        i + tensor_dofs_per_cell * c)] +=
                shared_data->values(i, c);
          });
      }
    else
      {
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(shared_data->team_member, n_q_points),
          [&](const int &i) {
            for (unsigned int c = 0; c < n_components_; ++c)
              Kokkos::atomic_add(&dst[data->local_to_global(
                                   cell_id, i + tensor_dofs_per_cell * c)],
                                 shared_data->values(i, c));
          });
      }
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::evaluate(
    const EvaluationFlags::EvaluationFlags evaluate_flag)
  {
    // First evaluate the gradients because it requires values that will be
    // changed if evaluate_val is true
    internal::EvaluatorTensorProduct<
      internal::EvaluatorVariant::evaluate_general,
      dim,
      fe_degree,
      n_q_points_1d,
      Number>
      evaluator_tensor_product(shared_data->team_member,
                               data->shape_values,
                               data->shape_gradients,
                               data->co_shape_gradients);

    for (unsigned int c = 0; c < n_components_; ++c)
      {
        if ((evaluate_flag & EvaluationFlags::values) &&
            (evaluate_flag & EvaluationFlags::gradients))
          {
            evaluator_tensor_product.evaluate_values_and_gradients(
              Kokkos::subview(shared_data->values, Kokkos::ALL, c),
              Kokkos::subview(
                shared_data->gradients, Kokkos::ALL, Kokkos::ALL, c));
            shared_data->team_member.team_barrier();
          }
        else if (evaluate_flag & EvaluationFlags::gradients)
          {
            evaluator_tensor_product.evaluate_gradients(
              Kokkos::subview(shared_data->values, Kokkos::ALL, c),
              Kokkos::subview(
                shared_data->gradients, Kokkos::ALL, Kokkos::ALL, c));
            shared_data->team_member.team_barrier();
          }
        else if (evaluate_flag & EvaluationFlags::values)
          {
            evaluator_tensor_product.evaluate_values(
              Kokkos::subview(shared_data->values, Kokkos::ALL, c));
            shared_data->team_member.team_barrier();
          }
      }
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::evaluate(
    const bool evaluate_val,
    const bool evaluate_grad)
  {
    evaluate(
      (evaluate_val ? EvaluationFlags::values : EvaluationFlags::nothing) |
      (evaluate_grad ? EvaluationFlags::gradients : EvaluationFlags::nothing));
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::integrate(
    const EvaluationFlags::EvaluationFlags integration_flag)
  {
    internal::EvaluatorTensorProduct<
      internal::EvaluatorVariant::evaluate_general,
      dim,
      fe_degree,
      n_q_points_1d,
      Number>
      evaluator_tensor_product(shared_data->team_member,
                               data->shape_values,
                               data->shape_gradients,
                               data->co_shape_gradients);


    for (unsigned int c = 0; c < n_components_; ++c)
      {
        if ((integration_flag & EvaluationFlags::values) &&
            (integration_flag & EvaluationFlags::gradients))
          {
            evaluator_tensor_product.integrate_values_and_gradients(
              Kokkos::subview(shared_data->values, Kokkos::ALL, c),
              Kokkos::subview(
                shared_data->gradients, Kokkos::ALL, Kokkos::ALL, c));
          }
        else if (integration_flag & EvaluationFlags::values)
          {
            evaluator_tensor_product.integrate_values(
              Kokkos::subview(shared_data->values, Kokkos::ALL, c));
            shared_data->team_member.team_barrier();
          }
        else if (integration_flag & EvaluationFlags::gradients)
          {
            evaluator_tensor_product.template integrate_gradients<false>(
              Kokkos::subview(shared_data->values, Kokkos::ALL, c),
              Kokkos::subview(
                shared_data->gradients, Kokkos::ALL, Kokkos::ALL, c));
            shared_data->team_member.team_barrier();
          }
      }
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::integrate(
    const bool integrate_val,
    const bool integrate_grad)
  {
    integrate(
      (integrate_val ? EvaluationFlags::values : EvaluationFlags::nothing) |
      (integrate_grad ? EvaluationFlags::gradients : EvaluationFlags::nothing));
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE typename FEEvaluation<dim,
                                            fe_degree,
                                            n_q_points_1d,
                                            n_components_,
                                            Number>::value_type
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::get_value(
    int q_point) const
  {
    if constexpr (n_components_ == 1)
      {
        return shared_data->values(q_point, 0);
      }
    else
      {
        value_type result;
        for (unsigned int c = 0; c < n_components; ++c)
          result[c] = shared_data->values(q_point, c);
        return result;
      }
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE typename FEEvaluation<dim,
                                            fe_degree,
                                            n_q_points_1d,
                                            n_components_,
                                            Number>::value_type
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    get_dof_value(int q_point) const
  {
    if constexpr (n_components_ == 1)
      {
        return shared_data->values(q_point, 0);
      }
    else
      {
        value_type result;
        for (unsigned int c = 0; c < n_components; ++c)
          result[c] = shared_data->values(q_point, c);
        return result;
      }
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    submit_value(const value_type &val_in, int q_point)
  {
    if constexpr (n_components_ == 1)
      {
        shared_data->values(q_point, 0) = val_in * data->JxW(cell_id, q_point);
      }
    else
      {
        for (unsigned int c = 0; c < n_components; ++c)
          shared_data->values(q_point, c) =
            val_in[c] * data->JxW(cell_id, q_point);
      }
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    submit_dof_value(const value_type &val_in, int q_point)
  {
    if constexpr (n_components_ == 1)
      {
        shared_data->values(q_point, 0) = val_in;
      }
    else
      {
        for (unsigned int c = 0; c < n_components; ++c)
          shared_data->values(q_point, c) = val_in[c];
      }
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE typename FEEvaluation<dim,
                                            fe_degree,
                                            n_q_points_1d,
                                            n_components_,
                                            Number>::gradient_type
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    get_gradient(int q_point) const
  {
    gradient_type grad;

    if constexpr (n_components_ == 1)
      {
        for (unsigned int d_1 = 0; d_1 < dim; ++d_1)
          {
            Number tmp = 0.;
            for (unsigned int d_2 = 0; d_2 < dim; ++d_2)
              tmp += data->inv_jacobian(cell_id, q_point, d_2, d_1) *
                     shared_data->gradients(q_point, d_2, 0);
            grad[d_1] = tmp;
          }
      }
    else
      {
        for (unsigned int c = 0; c < n_components; ++c)
          for (unsigned int d_1 = 0; d_1 < dim; ++d_1)
            {
              Number tmp = 0.;
              for (unsigned int d_2 = 0; d_2 < dim; ++d_2)
                tmp += data->inv_jacobian(cell_id, q_point, d_2, d_1) *
                       shared_data->gradients(q_point, d_2, c);
              grad[c][d_1] = tmp;
            }
      }

    return grad;
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    submit_gradient(const gradient_type &grad_in, int q_point)
  {
    if constexpr (n_components_ == 1)
      {
        for (unsigned int d_1 = 0; d_1 < dim; ++d_1)
          {
            Number tmp = 0.;
            for (unsigned int d_2 = 0; d_2 < dim; ++d_2)
              tmp +=
                data->inv_jacobian(cell_id, q_point, d_1, d_2) * grad_in[d_2];
            shared_data->gradients(q_point, d_1, 0) =
              tmp * data->JxW(cell_id, q_point);
          }
      }
    else
      {
        for (unsigned int c = 0; c < n_components; ++c)
          for (unsigned int d_1 = 0; d_1 < dim; ++d_1)
            {
              Number tmp = 0.;
              for (unsigned int d_2 = 0; d_2 < dim; ++d_2)
                tmp += data->inv_jacobian(cell_id, q_point, d_1, d_2) *
                       grad_in[c][d_2];
              shared_data->gradients(q_point, d_1, c) =
                tmp * data->JxW(cell_id, q_point);
            }
      }
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  template <typename Functor>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    apply_for_each_quad_point(const Functor &func)
  {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(shared_data->team_member,
                                                 n_q_points),
                         [&](const int &i) { func(this, i); });
    shared_data->team_member.team_barrier();
  }



#ifndef DOXYGEN
  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  constexpr unsigned int
    FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
      n_q_points;
#endif
} // namespace Portable

DEAL_II_NAMESPACE_CLOSE

#endif
