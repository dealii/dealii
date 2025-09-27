// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_values_views_internal_h
#define dealii_fe_values_views_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_values_views.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace FEValuesViews
{
  /**
   * Internal namespace for the utility functions used to actually compute
   * values by FEValuesViews
   */
  namespace internal
  {
    // Given values of degrees of freedom, evaluate the
    // values/gradients/... at quadrature points

    // ------------------------- scalar functions --------------------------

    /**
     * Compute function values for Scalars.
     */
    template <int dim, int spacedim, typename Number>
    void
    do_function_values(
      const ArrayView<const Number> &dof_values,
      const Table<2, double>        &shape_values,
      const std::vector<typename Scalar<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename ProductType<Number, double>::type> &values);

    /**
     * same code for gradient and Hessian, template argument 'order' to give
     * the order of the derivative (= rank of gradient/Hessian tensor)
     */
    template <int order, int dim, int spacedim, typename Number>
    void
    do_function_derivatives(
      const ArrayView<const Number>                   &dof_values,
      const Table<2, dealii::Tensor<order, spacedim>> &shape_derivatives,
      const std::vector<typename Scalar<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<
        typename ProductType<Number, dealii::Tensor<order, spacedim>>::type>
        &derivatives);

    /**
     * Compute Laplacian values for Scalars.
     */
    template <int dim, int spacedim, typename Number>
    void
    do_function_laplacians(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<2, spacedim>> &shape_hessians,
      const std::vector<typename Scalar<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename Scalar<dim, spacedim>::
                    template solution_laplacian_type<Number>> &laplacians);

    // ----------------------------- vector part ---------------------------

    /**
     * Compute function values for Vectors.
     */
    template <int dim, int spacedim, typename Number>
    void
    do_function_values(
      const ArrayView<const Number> &dof_values,
      const Table<2, double>        &shape_values,
      const std::vector<typename Vector<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<
        typename ProductType<Number, dealii::Tensor<1, spacedim>>::type>
        &values);

    /**
     * same code for gradient and Hessian, template argument 'order' to give
     * the order of the derivative (= rank of gradient/Hessian tensor)
     */
    template <int order, int dim, int spacedim, typename Number>
    void
    do_function_derivatives(
      const ArrayView<const Number>                   &dof_values,
      const Table<2, dealii::Tensor<order, spacedim>> &shape_derivatives,
      const std::vector<typename Vector<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<
        typename ProductType<Number, dealii::Tensor<order + 1, spacedim>>::type>
        &derivatives);

    /**
     * Compute gradient values for Vectors.
     */
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
        &symmetric_gradients);

    /**
     * Compute divergence values for Vectors.
     */
    template <int dim, int spacedim, typename Number>
    void
    do_function_divergences(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<1, spacedim>> &shape_gradients,
      const std::vector<typename Vector<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename Vector<dim, spacedim>::
                    template solution_divergence_type<Number>> &divergences);

    /**
     * Compute curl values for Vectors.
     */
    template <int dim, int spacedim, typename Number>
    void
    do_function_curls(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<1, spacedim>> &shape_gradients,
      const std::vector<typename Vector<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename ProductType<
        Number,
        typename dealii::internal::CurlType<spacedim>::type>::type> &curls);

    /**
     * Compute Laplacian values for Vectors.
     */
    template <int dim, int spacedim, typename Number>
    void
    do_function_laplacians(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<2, spacedim>> &shape_hessians,
      const std::vector<typename Vector<dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename Vector<dim, spacedim>::
                    template solution_laplacian_type<Number>> &laplacians);

    // ---------------------- symmetric tensor part ------------------------

    /**
     * Compute values for symmetric tensors.
     */
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
        &values);

    /**
     * Compute divergence values for symmetric tensors.
     */
    template <int dim, int spacedim, typename Number>
    void
    do_function_divergences(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<1, spacedim>> &shape_gradients,
      const std::vector<
        typename SymmetricTensor<2, dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename SymmetricTensor<2, dim, spacedim>::
                    template solution_divergence_type<Number>> &divergences);

    // ---------------------- non-symmetric tensor part ------------------------

    /**
     * Compute values for nonsymmetric tensors.
     */
    template <int dim, int spacedim, typename Number>
    void
    do_function_values(
      const ArrayView<const Number>  &dof_values,
      const dealii::Table<2, double> &shape_values,
      const std::vector<typename Tensor<2, dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<
        typename ProductType<Number, dealii::Tensor<2, spacedim>>::type>
        &values);

    /**
     * Compute divergence values for nonsymmetric tensors.
     */
    template <int dim, int spacedim, typename Number>
    void
    do_function_divergences(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<1, spacedim>> &shape_gradients,
      const std::vector<typename Tensor<2, dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename Tensor<2, dim, spacedim>::
                    template solution_divergence_type<Number>> &divergences);

    /**
     * Compute gradient values for nonsymmetric tensors.
     */
    template <int dim, int spacedim, typename Number>
    void
    do_function_gradients(
      const ArrayView<const Number>               &dof_values,
      const Table<2, dealii::Tensor<1, spacedim>> &shape_gradients,
      const std::vector<typename Tensor<2, dim, spacedim>::ShapeFunctionData>
        &shape_function_data,
      std::vector<typename Tensor<2, dim, spacedim>::
                    template solution_gradient_type<Number>> &gradients);
  } // end of namespace internal
} // namespace FEValuesViews


DEAL_II_NAMESPACE_CLOSE

#endif
