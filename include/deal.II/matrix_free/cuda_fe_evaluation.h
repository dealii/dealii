// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2023 by the deal.II authors
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

#ifndef dealii_cuda_fe_evaluation_h
#define dealii_cuda_fe_evaluation_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/cuda_atomic.h>
#include <deal.II/lac/cuda_vector.h>

#include <deal.II/matrix_free/cuda_hanging_nodes_internal.h>
#include <deal.II/matrix_free/cuda_matrix_free.h>
#include <deal.II/matrix_free/cuda_matrix_free.templates.h>
#include <deal.II/matrix_free/cuda_tensor_product_kernels.h>

#include <Kokkos_Core.hpp>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace for the CUDA wrappers
 */
namespace CUDAWrappers
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
   *
   * @ingroup CUDAWrappers
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
    using value_type = Number;

    /**
     * An alias for vectorial quantities.
     */
    using gradient_type = Tensor<1, dim, Number>;

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
     * unit cell. The function arguments specify which parts shall actually be
     * computed. This function needs to be called before the functions
     * @p get_value() or @p get_gradient() give useful information.
     */
    DEAL_II_HOST_DEVICE void
    evaluate(const bool evaluate_val, const bool evaluate_grad);

    /**
     * This function takes the values and/or gradients that are stored on
     * quadrature points, tests them by all the basis functions/gradients on
     * the cell and performs the cell integration. The two function arguments
     * @p integrate_val and @p integrate_grad are used to enable/disable some
     * of the values or the gradients.
     */
    DEAL_II_HOST_DEVICE void
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
     *   CUDAWrappers::FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, Number> *fe_eval) const;
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
    static_assert(n_components_ == 1, "This function only supports FE with one \
                  components");
    // Populate the scratch memory
    Kokkos::parallel_for(Kokkos::TeamThreadRange(shared_data->team_member,
                                                 n_q_points),
                         [&](const int &i) {
                           shared_data->values(i) =
                             src[data->local_to_global(cell_id, i)];
                         });
    shared_data->team_member.team_barrier();

    internal::resolve_hanging_nodes<dim, fe_degree, false>(
      shared_data->team_member,
      data->constraint_weights,
      data->constraint_mask(cell_id),
      shared_data->values);
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
    static_assert(n_components_ == 1, "This function only supports FE with one \
                  components");

    internal::resolve_hanging_nodes<dim, fe_degree, true>(
      shared_data->team_member,
      data->constraint_weights,
      data->constraint_mask(cell_id),
      shared_data->values);

    if (data->use_coloring)
      {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(shared_data->team_member,
                                                     n_q_points),
                             [&](const int &i) {
                               dst[data->local_to_global(cell_id, i)] +=
                                 shared_data->values(i);
                             });
      }
    else
      {
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(shared_data->team_member, n_q_points),
          [&](const int &i) {
            Kokkos::atomic_add(&dst[data->local_to_global(cell_id, i)],
                               shared_data->values(i));
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
    const bool evaluate_val,
    const bool evaluate_grad)
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
    if (evaluate_val == true && evaluate_grad == true)
      {
        evaluator_tensor_product.value_and_gradient_at_quad_pts(
          shared_data->values, shared_data->gradients);
        shared_data->team_member.team_barrier();
      }
    else if (evaluate_grad == true)
      {
        evaluator_tensor_product.gradient_at_quad_pts(shared_data->values,
                                                      shared_data->gradients);
        shared_data->team_member.team_barrier();
      }
    else if (evaluate_val == true)
      {
        evaluator_tensor_product.value_at_quad_pts(shared_data->values);
        shared_data->team_member.team_barrier();
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
    if (integrate_val == true && integrate_grad == true)
      {
        evaluator_tensor_product.integrate_value_and_gradient(
          shared_data->values, shared_data->gradients);
      }
    else if (integrate_val == true)
      {
        evaluator_tensor_product.integrate_value(shared_data->values);
        shared_data->team_member.team_barrier();
      }
    else if (integrate_grad == true)
      {
        evaluator_tensor_product.template integrate_gradient<false>(
          shared_data->values, shared_data->gradients);
        shared_data->team_member.team_barrier();
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
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::get_value(
    int q_point) const
  {
    return shared_data->values(q_point);
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
    return shared_data->values(q_point);
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
    shared_data->values(q_point) = val_in * data->JxW(cell_id, q_point);
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
    shared_data->values(q_point) = val_in;
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
    static_assert(n_components_ == 1, "This function only supports FE with one \
                  components");

    gradient_type grad;
    for (unsigned int d_1 = 0; d_1 < dim; ++d_1)
      {
        Number tmp = 0.;
        for (unsigned int d_2 = 0; d_2 < dim; ++d_2)
          tmp += data->inv_jacobian(cell_id, q_point, d_2, d_1) *
                 shared_data->gradients(q_point, d_2);
        grad[d_1] = tmp;
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
    for (unsigned int d_1 = 0; d_1 < dim; ++d_1)
      {
        Number tmp = 0.;
        for (unsigned int d_2 = 0; d_2 < dim; ++d_2)
          tmp += data->inv_jacobian(cell_id, q_point, d_1, d_2) * grad_in[d_2];
        shared_data->gradients(q_point, d_1) =
          tmp * data->JxW(cell_id, q_point);
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
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
