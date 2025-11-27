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

#ifndef dealii_portable_fe_evaluation_h
#define dealii_portable_fe_evaluation_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/portable_evaluation_kernels.h>
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
   * @tparam fe_degree Degree of the tensor product finite element with fe_degree+1
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
     * An alias for the value type. This is @p Number for scalar problems
     * and Tensor<1, n_components> for vector-valued problems.
     */
    using value_type = std::conditional_t<(n_components_ == 1),
                                          Number,
                                          Tensor<1, n_components_, Number>>;

    /**
     * An alias for the gradient type.
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
     * Number of tensor degrees of freedom of a scalar component determined
     * from the given template argument `fe_degree`.
     */
    static constexpr unsigned int tensor_dofs_per_component =
      Utilities::pow(fe_degree + 1, dim);

    /**
     * Number of tensor degrees of freedom of all component determined from
     * the given template argument `fe_degree`. This is the total number of
     * local DoFs in a cell.
     */
    static constexpr unsigned int tensor_dofs_per_cell =
      tensor_dofs_per_component * n_components;

    /**
     * Constructor. You will need to provide a pointer to the
     * Portable::MatrixFree::Data object, which is typically provided to the
     * functor inside the
     * Portable::MatrixFree::cell_loop() and the index @p dof_handler_index
     * of the DoFHandler if more than one was provided when the
     * Portable::MatrixFree object was initialized.
     */
    DEAL_II_HOST_DEVICE
    explicit FEEvaluation(const data_type   *data,
                          const unsigned int dof_handler_index = 0);

    /**
     * Return the index of the current cell.
     */
    DEAL_II_HOST_DEVICE
    int
    get_current_cell_index();

    /**
     * Return a pointer to the MatrixFree<dim, Number>::Data object on device
     * that contains necessary constraint, dof index, and shape function
     * information for evaluation used in the matrix-free kernels.
     */
    DEAL_II_HOST_DEVICE
    const data_type *
    get_matrix_free_data();

    /**
     * For the vector @p src, read out the values on the degrees of freedom of
     * the current cell, and store them internally. Similar functionality as
     * the function DoFAccessor::get_interpolated_dof_values when no
     * constraints are present, but it also includes constraints from hanging
     * nodes, so one can see it as a similar function to
     * AffineConstraints::read_dof_values() as well.
     */
    DEAL_II_HOST_DEVICE void
    read_dof_values(const DeviceVector<Number> &src);

    /**
     * Take the value stored internally on dof values of the current cell and
     * sum them into the vector @p dst. The function also applies constraints
     * during the write operation. The functionality is hence similar to the
     * function AffineConstraints::distribute_local_to_global.
     */
    DEAL_II_HOST_DEVICE void
    distribute_local_to_global(DeviceVector<Number> &dst) const;

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
     * This function takes the values and/or gradients that are stored on
     * quadrature points, tests them by all the basis functions/gradients on
     * the cell and performs the cell integration as specified by the
     * @p integration_flag argument.
     */
    DEAL_II_HOST_DEVICE void
    integrate(const EvaluationFlags::EvaluationFlags integration_flag);

    /**
     * Return the value of the finite element function at the quadrature point
     * with index @p q_point after a call to evaluate() with
     * EvaluationFlags::values set.
     */
    DEAL_II_HOST_DEVICE value_type
    get_value(int q_point) const;

    /**
     * Return the value stored for the local degree of freedom with index
     * @p dof_index. This accesses the data loaded by read_dof_values().
     */
    DEAL_II_HOST_DEVICE value_type
    get_dof_value(int dof_index) const;

    /**
     * Submit the value @p val_in at quadrature point @p q_point for
     * subsequent integration via integrate() with EvaluationFlags::values
     * set.
     */
    DEAL_II_HOST_DEVICE void
    submit_value(const value_type &val_in, int q_point);

    /**
     * Submit the value @p value for the local degree of freedom with index
     * @p dof_index, to be written out by a subsequent call to
     * distribute_local_to_global().
     */
    DEAL_II_HOST_DEVICE void
    submit_dof_value(const value_type &value, int dof_index);

    /**
     * Return the gradient of the finite element function at the quadrature
     * point with index @p q_point after a call to evaluate() with
     * EvaluationFlags::gradients set.
     */
    DEAL_II_HOST_DEVICE gradient_type
    get_gradient(int q_point) const;

    /**
     * Submit the gradient @p gradient at quadrature point @p q_point for
     * subsequent integration via integrate() with EvaluationFlags::gradients
     * set.
     */
    DEAL_II_HOST_DEVICE void
    submit_gradient(const gradient_type &gradient, int q_point);

    /**
     * Return the symmetric gradient of the finite element function at
     * quadrature point @p q_point after a call to evaluate() with
     * EvaluationFlags::gradients set. This function is only available when
     * the number of components equals the dimension (n_components==dim).
     */
    DEAL_II_HOST_DEVICE
    SymmetricTensor<2, dim, Number>
    get_symmetric_gradient(int q_point) const;

    /**
     * Submit the symmetric gradient @p sym_grad at quadrature point @p q_point
     * for subsequent integration via integrate() with
     * EvaluationFlags::gradients set. This function is only available when
     * the number of components equals the dimension (n_components_==dim).
     */
    DEAL_II_HOST_DEVICE void
    submit_symmetric_gradient(const SymmetricTensor<2, dim, Number> &sym_grad,
                              int                                    q_point);

    // clang-format off
    /**
     * Apply the functor @p func to each quadrature point in parallel.
     * The functor is invoked with the FEEvaluation object and the quadrature
     * point index as arguments.
     *
     * @p func needs to define
     * \code
     * DEAL_II_HOST_DEVICE void operator()(
     *   Portable::FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, Number> *fe_eval,
     *   int q_point) const;
     * \endcode
     *
     * @deprecated Use MatrixFree::Data::for_each_quad_point() instead.
     */
    // clang-format on
    template <typename Functor>
    DEAL_II_DEPRECATED DEAL_II_HOST_DEVICE void
    apply_for_each_quad_point(const Functor &func);

  private:
    const unsigned int                                       dof_handler_index;
    const typename MatrixFree<dim, Number>::Data            *data;
    const typename MatrixFree<dim, Number>::PrecomputedData *precomputed_data;
    SharedData<dim, Number>                                 *shared_data;
    int                                                      cell_id;
  };



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    FEEvaluation(const data_type *data, const unsigned int dof_handler_index)
    : dof_handler_index(dof_handler_index)
    , data(data)
    , precomputed_data(&data->precomputed_data[dof_handler_index])
    , shared_data(&data->shared_data[dof_handler_index])
    , cell_id(data->team_member.league_rank())
  {
    AssertIndexRange(dof_handler_index, data->n_dof_handler);

    Assert(
      n_components_ == precomputed_data->n_components,
      ExcMessage(
        "Portable::FEEvaluation initialized with wrong number of components. Should be " +
        Utilities::to_string(precomputed_data->n_components) +
        " but the template argument 4 is set to " +
        Utilities::to_string(n_components_)));

    // TODO: check fe_degree is correct by storing the used FE degree
    // inside PrecomputedData.
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE int
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    get_current_cell_index()
  {
    return cell_id;
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE const typename FEEvaluation<dim,
                                                  fe_degree,
                                                  n_q_points_1d,
                                                  n_components_,
                                                  Number>::data_type *
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    get_matrix_free_data()
  {
    return data;
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    read_dof_values(const DeviceVector<Number> &src)
  {
    // Populate the scratch memory
    Kokkos::parallel_for(Kokkos::TeamThreadRange(data->team_member,
                                                 tensor_dofs_per_component),
                         [&](const int &i) {
                           for (unsigned int c = 0; c < n_components_; ++c)
                             shared_data->values(i, c) =
                               src[precomputed_data->local_to_global(
                                 i + tensor_dofs_per_component * c, cell_id)];
                         });
    data->team_member.team_barrier();

    for (unsigned int c = 0; c < n_components_; ++c)
      {
        if (precomputed_data->constraint_mask(cell_id * n_components + c) !=
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::
              unconstrained)
          internal::resolve_hanging_nodes<dim, fe_degree, false, Number>(
            data->team_member,
            precomputed_data->constraint_weights,
            precomputed_data->constraint_mask(cell_id * n_components + c),
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
    distribute_local_to_global(DeviceVector<Number> &dst) const
  {
    for (unsigned int c = 0; c < n_components_; ++c)
      {
        if (precomputed_data->constraint_mask(cell_id * n_components + c) !=
            dealii::internal::MatrixFreeFunctions::ConstraintKinds::
              unconstrained)
          internal::resolve_hanging_nodes<dim, fe_degree, true, Number>(
            data->team_member,
            precomputed_data->constraint_weights,
            precomputed_data->constraint_mask(cell_id * n_components + c),
            Kokkos::subview(shared_data->values, Kokkos::ALL, c));
      }

    if (precomputed_data->use_coloring)
      {
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(data->team_member, tensor_dofs_per_component),
          [&](const int &i) {
            for (unsigned int c = 0; c < n_components_; ++c)
              dst[precomputed_data->local_to_global(
                i + tensor_dofs_per_component * c, cell_id)] +=
                shared_data->values(i, c);
          });
      }
    else
      {
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(data->team_member, tensor_dofs_per_component),
          [&](const int &i) {
            for (unsigned int c = 0; c < n_components_; ++c)
              Kokkos::atomic_add(&dst[precomputed_data->local_to_global(
                                   i + (tensor_dofs_per_component)*c, cell_id)],
                                 shared_data->values(i, c));
          });
      }
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, Number>::evaluate(
    const EvaluationFlags::EvaluationFlags evaluation_flag)
  {
    using ElementType = ::dealii::internal::MatrixFreeFunctions::ElementType;

    if (fe_degree >= 0 && fe_degree + 1 == n_q_points_1d &&
        precomputed_data->element_type ==
          ElementType::tensor_symmetric_collocation)
      {
        internal::FEEvaluationImplCollocation<dim, fe_degree, Number>::evaluate(
          dof_handler_index, n_components, evaluation_flag, data);
      }
    // '<=' on type means tensor_symmetric or tensor_symmetric_hermite, see
    // shape_info.h for more details
    else if (fe_degree >= 0 &&
             internal::use_collocation_evaluation(fe_degree, n_q_points_1d) &&
             precomputed_data->element_type <= ElementType::tensor_symmetric)
      {
        internal::FEEvaluationImplTransformToCollocation<
          dim,
          fe_degree,
          n_q_points_1d,
          Number>::evaluate(dof_handler_index,
                            n_components,
                            evaluation_flag,
                            data);
      }
    else if (fe_degree >= 0 && precomputed_data->element_type <=
                                 ElementType::tensor_symmetric_no_collocation)
      {
        internal::FEEvaluationImpl<dim, fe_degree, n_q_points_1d, Number>::
          evaluate(dof_handler_index, n_components, evaluation_flag, data);
      }
    else
      {
        Kokkos::abort("The element type is not yet supported by the portable "
                      "matrix-free module.");
      }
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
    using ElementType = ::dealii::internal::MatrixFreeFunctions::ElementType;

    if (fe_degree >= 0 && fe_degree + 1 == n_q_points_1d &&
        precomputed_data->element_type ==
          ElementType::tensor_symmetric_collocation)
      {
        internal::FEEvaluationImplCollocation<dim, fe_degree, Number>::
          integrate(dof_handler_index, n_components, integration_flag, data);
      }
    // '<=' on type means tensor_symmetric or tensor_symmetric_hermite, see
    // shape_info.h for more details
    else if (fe_degree >= 0 &&
             internal::use_collocation_evaluation(fe_degree, n_q_points_1d) &&
             precomputed_data->element_type <= ElementType::tensor_symmetric)
      {
        internal::FEEvaluationImplTransformToCollocation<
          dim,
          fe_degree,
          n_q_points_1d,
          Number>::integrate(dof_handler_index,
                             n_components,
                             integration_flag,
                             data);
      }
    else if (fe_degree >= 0 && precomputed_data->element_type <=
                                 ElementType::tensor_symmetric_no_collocation)
      {
        internal::FEEvaluationImpl<dim, fe_degree, n_q_points_1d, Number>::
          integrate(dof_handler_index, n_components, integration_flag, data);
      }
    else
      {
        Kokkos::abort("The element type is not yet supported by the portable "
                      "matrix-free module.");
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
    Assert(q_point >= 0 && q_point < n_q_points, ExcInternalError());
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
    get_dof_value(int dof_index) const
  {
    Assert(dof_index >= 0 &&
             dof_index < static_cast<int>(tensor_dofs_per_component),
           ExcInternalError());
    if constexpr (n_components_ == 1)
      {
        return shared_data->values(dof_index, 0);
      }
    else
      {
        value_type result;
        for (unsigned int c = 0; c < n_components; ++c)
          result[c] = shared_data->values(dof_index, c);
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
    submit_value(const value_type &value, int q_point)
  {
    Assert(q_point >= 0 && q_point < n_q_points, ExcInternalError());
    if constexpr (n_components_ == 1)
      {
        shared_data->values(q_point, 0) =
          value * precomputed_data->JxW(q_point, cell_id);
      }
    else
      {
        for (unsigned int c = 0; c < n_components; ++c)
          shared_data->values(q_point, c) =
            value[c] * precomputed_data->JxW(q_point, cell_id);
      }
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    submit_dof_value(const value_type &value, int dof_index)
  {
    Assert(dof_index >= 0 && dof_index < tensor_dofs_per_component,
           ExcInternalError());
    if constexpr (n_components_ == 1)
      {
        shared_data->values(dof_index, 0) = value;
      }
    else
      {
        for (unsigned int c = 0; c < n_components; ++c)
          shared_data->values(dof_index, c) = value[c];
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
    Assert(q_point >= 0 && q_point < n_q_points, ExcInternalError());
    gradient_type grad;

    if constexpr (n_components_ == 1)
      {
        for (unsigned int d_1 = 0; d_1 < dim; ++d_1)
          {
            Number tmp = 0.;
            for (unsigned int d_2 = 0; d_2 < dim; ++d_2)
              tmp +=
                precomputed_data->inv_jacobian(q_point, cell_id, d_2, d_1) *
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
                tmp +=
                  precomputed_data->inv_jacobian(q_point, cell_id, d_2, d_1) *
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
    submit_gradient(const gradient_type &gradient, int q_point)
  {
    Assert(q_point >= 0 && q_point < n_q_points, ExcInternalError());
    if constexpr (n_components_ == 1)
      {
        for (unsigned int d_1 = 0; d_1 < dim; ++d_1)
          {
            Number tmp = 0.;
            for (unsigned int d_2 = 0; d_2 < dim; ++d_2)
              tmp +=
                precomputed_data->inv_jacobian(q_point, cell_id, d_1, d_2) *
                gradient[d_2];
            shared_data->gradients(q_point, d_1, 0) =
              tmp * precomputed_data->JxW(q_point, cell_id);
          }
      }
    else
      {
        for (unsigned int c = 0; c < n_components; ++c)
          for (unsigned int d_1 = 0; d_1 < dim; ++d_1)
            {
              Number tmp = 0.;
              for (unsigned int d_2 = 0; d_2 < dim; ++d_2)
                tmp +=
                  precomputed_data->inv_jacobian(q_point, cell_id, d_1, d_2) *
                  gradient[c][d_2];
              shared_data->gradients(q_point, d_1, c) =
                tmp * precomputed_data->JxW(q_point, cell_id);
            }
      }
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE SymmetricTensor<2, dim, Number>
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    get_symmetric_gradient(int q_point) const
  {
    Assert(q_point >= 0 && q_point < n_q_points, ExcInternalError());
    Assert(n_components_ == dim,
           ExcMessage("Function get_symmetric_gradient() only works when the "
                      "number of components and the number of dimensions are "
                      "equal."));

    return symmetrize(get_gradient(q_point));
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number>
  DEAL_II_HOST_DEVICE void
  FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number>::
    submit_symmetric_gradient(const SymmetricTensor<2, dim, Number> &sym_grad,
                              int                                    q_point)
  {
    Assert(q_point >= 0 && q_point < n_q_points, ExcInternalError());
    Assert(n_components_ == dim,
           ExcMessage("Function submit_symmetric_gradient() only works when "
                      "the number of components and the number of dimensions "
                      "are equal."));

    for (unsigned int c = 0; c < dim; ++c)
      for (unsigned int d_1 = 0; d_1 < dim; ++d_1)
        {
          Number tmp = 0.;
          for (unsigned int d_2 = 0; d_2 < dim; ++d_2)
            tmp += precomputed_data->inv_jacobian(q_point, cell_id, d_1, d_2) *
                   sym_grad[c][d_2];
          shared_data->gradients(q_point, d_1, c) =
            tmp * precomputed_data->JxW(q_point, cell_id);
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
    Kokkos::parallel_for(Kokkos::TeamThreadRange(data->team_member, n_q_points),
                         [&](const int &i) { func(this, i); });
    data->team_member.team_barrier();
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
