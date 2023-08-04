// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2023 by the deal.II authors
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

#include <deal.II/base/array_view.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values_base.h>
#include <deal.II/fe/fe_values_views.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <boost/container/small_vector.hpp>

#include <iomanip>
#include <memory>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  template <int dim, int spacedim>
  inline std::vector<unsigned int>
  make_shape_function_to_row_table(const FiniteElement<dim, spacedim> &fe)
  {
    std::vector<unsigned int> shape_function_to_row_table(
      fe.n_dofs_per_cell() * fe.n_components(), numbers::invalid_unsigned_int);
    unsigned int row = 0;
    for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
      {
        // loop over all components that are nonzero for this particular
        // shape function. if a component is zero then we leave the
        // value in the table unchanged (at the invalid value)
        // otherwise it is mapped to the next free entry
        unsigned int nth_nonzero_component = 0;
        for (unsigned int c = 0; c < fe.n_components(); ++c)
          if (fe.get_nonzero_components(i)[c] == true)
            {
              shape_function_to_row_table[i * fe.n_components() + c] =
                row + nth_nonzero_component;
              ++nth_nonzero_component;
            }
        row += fe.n_nonzero_components(i);
      }

    return shape_function_to_row_table;
  }

  namespace
  {
    // Check to see if a DoF value is zero, implying that subsequent operations
    // with the value have no effect.
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
} // namespace internal



namespace FEValuesViews
{
  template <int dim, int spacedim>
  Scalar<dim, spacedim>::Scalar(const FEValuesBase<dim, spacedim> &fe_values,
                                const unsigned int                 component)
    : fe_values(&fe_values)
    , component(component)
    , shape_function_data(this->fe_values->fe->n_dofs_per_cell())
  {
    const FiniteElement<dim, spacedim> &fe = *this->fe_values->fe;
    AssertIndexRange(component, fe.n_components());

    // TODO: we'd like to use the fields with the same name as these
    // variables from FEValuesBase, but they aren't initialized yet
    // at the time we get here, so re-create it all
    const std::vector<unsigned int> shape_function_to_row_table =
      dealii::internal::make_shape_function_to_row_table(fe);

    for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
      {
        const bool is_primitive = fe.is_primitive() || fe.is_primitive(i);

        if (is_primitive == true)
          shape_function_data[i].is_nonzero_shape_function_component =
            (component == fe.system_to_component_index(i).first);
        else
          shape_function_data[i].is_nonzero_shape_function_component =
            (fe.get_nonzero_components(i)[component] == true);

        if (shape_function_data[i].is_nonzero_shape_function_component == true)
          shape_function_data[i].row_index =
            shape_function_to_row_table[i * fe.n_components() + component];
        else
          shape_function_data[i].row_index = numbers::invalid_unsigned_int;
      }
  }



  template <int dim, int spacedim>
  Scalar<dim, spacedim>::Scalar()
    : fe_values(nullptr)
    , component(numbers::invalid_unsigned_int)
  {}



  template <int dim, int spacedim>
  Vector<dim, spacedim>::Vector(const FEValuesBase<dim, spacedim> &fe_values,
                                const unsigned int first_vector_component)
    : fe_values(&fe_values)
    , first_vector_component(first_vector_component)
    , shape_function_data(this->fe_values->fe->n_dofs_per_cell())
  {
    const FiniteElement<dim, spacedim> &fe = *this->fe_values->fe;
    AssertIndexRange(first_vector_component + spacedim - 1, fe.n_components());

    // TODO: we'd like to use the fields with the same name as these
    // variables from FEValuesBase, but they aren't initialized yet
    // at the time we get here, so re-create it all
    const std::vector<unsigned int> shape_function_to_row_table =
      dealii::internal::make_shape_function_to_row_table(fe);

    for (unsigned int d = 0; d < spacedim; ++d)
      {
        const unsigned int component = first_vector_component + d;

        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
          {
            const bool is_primitive = fe.is_primitive() || fe.is_primitive(i);

            if (is_primitive == true)
              shape_function_data[i].is_nonzero_shape_function_component[d] =
                (component == fe.system_to_component_index(i).first);
            else
              shape_function_data[i].is_nonzero_shape_function_component[d] =
                (fe.get_nonzero_components(i)[component] == true);

            if (shape_function_data[i].is_nonzero_shape_function_component[d] ==
                true)
              shape_function_data[i].row_index[d] =
                shape_function_to_row_table[i * fe.n_components() + component];
            else
              shape_function_data[i].row_index[d] =
                numbers::invalid_unsigned_int;
          }
      }

    for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
      {
        unsigned int n_nonzero_components = 0;
        for (unsigned int d = 0; d < spacedim; ++d)
          if (shape_function_data[i].is_nonzero_shape_function_component[d] ==
              true)
            ++n_nonzero_components;

        if (n_nonzero_components == 0)
          shape_function_data[i].single_nonzero_component = -2;
        else if (n_nonzero_components > 1)
          shape_function_data[i].single_nonzero_component = -1;
        else
          {
            for (unsigned int d = 0; d < spacedim; ++d)
              if (shape_function_data[i]
                    .is_nonzero_shape_function_component[d] == true)
                {
                  shape_function_data[i].single_nonzero_component =
                    shape_function_data[i].row_index[d];
                  shape_function_data[i].single_nonzero_component_index = d;
                  break;
                }
          }
      }
  }



  template <int dim, int spacedim>
  Vector<dim, spacedim>::Vector()
    : fe_values(nullptr)
    , first_vector_component(numbers::invalid_unsigned_int)
  {}



  template <int dim, int spacedim>
  SymmetricTensor<2, dim, spacedim>::SymmetricTensor(
    const FEValuesBase<dim, spacedim> &fe_values,
    const unsigned int                 first_tensor_component)
    : fe_values(&fe_values)
    , first_tensor_component(first_tensor_component)
    , shape_function_data(this->fe_values->fe->n_dofs_per_cell())
  {
    const FiniteElement<dim, spacedim> &fe = *this->fe_values->fe;
    Assert(first_tensor_component + (dim * dim + dim) / 2 - 1 <
             fe.n_components(),
           ExcIndexRange(
             first_tensor_component +
               dealii::SymmetricTensor<2, dim>::n_independent_components - 1,
             0,
             fe.n_components()));
    // TODO: we'd like to use the fields with the same name as these
    // variables from FEValuesBase, but they aren't initialized yet
    // at the time we get here, so re-create it all
    const std::vector<unsigned int> shape_function_to_row_table =
      dealii::internal::make_shape_function_to_row_table(fe);

    for (unsigned int d = 0;
         d < dealii::SymmetricTensor<2, dim>::n_independent_components;
         ++d)
      {
        const unsigned int component = first_tensor_component + d;

        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
          {
            const bool is_primitive = fe.is_primitive() || fe.is_primitive(i);

            if (is_primitive == true)
              shape_function_data[i].is_nonzero_shape_function_component[d] =
                (component == fe.system_to_component_index(i).first);
            else
              shape_function_data[i].is_nonzero_shape_function_component[d] =
                (fe.get_nonzero_components(i)[component] == true);

            if (shape_function_data[i].is_nonzero_shape_function_component[d] ==
                true)
              shape_function_data[i].row_index[d] =
                shape_function_to_row_table[i * fe.n_components() + component];
            else
              shape_function_data[i].row_index[d] =
                numbers::invalid_unsigned_int;
          }
      }

    for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
      {
        unsigned int n_nonzero_components = 0;
        for (unsigned int d = 0;
             d < dealii::SymmetricTensor<2, dim>::n_independent_components;
             ++d)
          if (shape_function_data[i].is_nonzero_shape_function_component[d] ==
              true)
            ++n_nonzero_components;

        if (n_nonzero_components == 0)
          shape_function_data[i].single_nonzero_component = -2;
        else if (n_nonzero_components > 1)
          shape_function_data[i].single_nonzero_component = -1;
        else
          {
            for (unsigned int d = 0;
                 d < dealii::SymmetricTensor<2, dim>::n_independent_components;
                 ++d)
              if (shape_function_data[i]
                    .is_nonzero_shape_function_component[d] == true)
                {
                  shape_function_data[i].single_nonzero_component =
                    shape_function_data[i].row_index[d];
                  shape_function_data[i].single_nonzero_component_index = d;
                  break;
                }
          }
      }
  }



  template <int dim, int spacedim>
  SymmetricTensor<2, dim, spacedim>::SymmetricTensor()
    : fe_values(nullptr)
    , first_tensor_component(numbers::invalid_unsigned_int)
  {}



  template <int dim, int spacedim>
  Tensor<2, dim, spacedim>::Tensor(const FEValuesBase<dim, spacedim> &fe_values,
                                   const unsigned int first_tensor_component)
    : fe_values(&fe_values)
    , first_tensor_component(first_tensor_component)
    , shape_function_data(this->fe_values->fe->n_dofs_per_cell())
  {
    const FiniteElement<dim, spacedim> &fe = *this->fe_values->fe;
    AssertIndexRange(first_tensor_component + dim * dim - 1, fe.n_components());
    // TODO: we'd like to use the fields with the same name as these
    // variables from FEValuesBase, but they aren't initialized yet
    // at the time we get here, so re-create it all
    const std::vector<unsigned int> shape_function_to_row_table =
      dealii::internal::make_shape_function_to_row_table(fe);

    for (unsigned int d = 0; d < dim * dim; ++d)
      {
        const unsigned int component = first_tensor_component + d;

        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
          {
            const bool is_primitive = fe.is_primitive() || fe.is_primitive(i);

            if (is_primitive == true)
              shape_function_data[i].is_nonzero_shape_function_component[d] =
                (component == fe.system_to_component_index(i).first);
            else
              shape_function_data[i].is_nonzero_shape_function_component[d] =
                (fe.get_nonzero_components(i)[component] == true);

            if (shape_function_data[i].is_nonzero_shape_function_component[d] ==
                true)
              shape_function_data[i].row_index[d] =
                shape_function_to_row_table[i * fe.n_components() + component];
            else
              shape_function_data[i].row_index[d] =
                numbers::invalid_unsigned_int;
          }
      }

    for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
      {
        unsigned int n_nonzero_components = 0;
        for (unsigned int d = 0; d < dim * dim; ++d)
          if (shape_function_data[i].is_nonzero_shape_function_component[d] ==
              true)
            ++n_nonzero_components;

        if (n_nonzero_components == 0)
          shape_function_data[i].single_nonzero_component = -2;
        else if (n_nonzero_components > 1)
          shape_function_data[i].single_nonzero_component = -1;
        else
          {
            for (unsigned int d = 0; d < dim * dim; ++d)
              if (shape_function_data[i]
                    .is_nonzero_shape_function_component[d] == true)
                {
                  shape_function_data[i].single_nonzero_component =
                    shape_function_data[i].row_index[d];
                  shape_function_data[i].single_nonzero_component_index = d;
                  break;
                }
          }
      }
  }



  template <int dim, int spacedim>
  Tensor<2, dim, spacedim>::Tensor()
    : fe_values(nullptr)
    , first_tensor_component(numbers::invalid_unsigned_int)
  {}



  namespace internal
  {
    // Given values of degrees of freedom, evaluate the
    // values/gradients/... at quadrature points

    // ------------------------- scalar functions --------------------------
    template <int dim, int spacedim, typename Number>
    void
    do_function_values(
      const ArrayView<Number> &dof_values,
      const Table<2, double> & shape_values,
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
            if (dealii::internal::CheckForZero<Number>::value(value) == true)
              continue;

            const double *shape_value_ptr =
              &shape_values(shape_function_data[shape_function].row_index, 0);
            for (unsigned int q_point = 0; q_point < n_quadrature_points;
                 ++q_point)
              values[q_point] += value * (*shape_value_ptr++);
          }
    }



    // same code for gradient and Hessian, template argument 'order' to give
    // the order of the derivative (= rank of gradient/Hessian tensor)
    template <int order, int dim, int spacedim, typename Number>
    void
    do_function_derivatives(
      const ArrayView<Number> &                        dof_values,
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
            if (dealii::internal::CheckForZero<Number>::value(value) == true)
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
      const ArrayView<Number> &                    dof_values,
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
            if (dealii::internal::CheckForZero<Number>::value(value) == true)
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
      const ArrayView<Number> &dof_values,
      const Table<2, double> & shape_values,
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
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const unsigned int comp = shape_function_data[shape_function]
                                          .single_nonzero_component_index;
              const double *shape_value_ptr = &shape_values(snc, 0);
              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point)
                values[q_point][comp] += value * (*shape_value_ptr++);
            }
          else
            for (unsigned int d = 0; d < spacedim; ++d)
              if (shape_function_data[shape_function]
                    .is_nonzero_shape_function_component[d])
                {
                  const double *shape_value_ptr = &shape_values(
                    shape_function_data[shape_function].row_index[d], 0);
                  for (unsigned int q_point = 0; q_point < n_quadrature_points;
                       ++q_point)
                    values[q_point][d] += value * (*shape_value_ptr++);
                }
        }
    }



    template <int order, int dim, int spacedim, typename Number>
    void
    do_function_derivatives(
      const ArrayView<Number> &                        dof_values,
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
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
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
      const ArrayView<Number> &                    dof_values,
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
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const unsigned int comp = shape_function_data[shape_function]
                                          .single_nonzero_component_index;
              const dealii::Tensor<1, spacedim> *shape_gradient_ptr =
                &shape_gradients[snc][0];
              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point)
                symmetric_gradients[q_point] +=
                  value * dealii::SymmetricTensor<2, spacedim>(
                            symmetrize_single_row(comp, *shape_gradient_ptr++));
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
      const ArrayView<Number> &                    dof_values,
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
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
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
      const ArrayView<Number> &                    dof_values,
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
                  if (dealii::internal::CheckForZero<Number>::value(value) ==
                      true)
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
                  if (dealii::internal::CheckForZero<Number>::value(value) ==
                      true)
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
                            Assert(false, ExcInternalError());
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
      const ArrayView<Number> &                    dof_values,
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
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
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
      const ArrayView<Number> &       dof_values,
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
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
            continue;

          if (snc != -1)
            {
              const TableIndices<2> comp = dealii::
                SymmetricTensor<2, spacedim>::unrolled_to_component_indices(
                  shape_function_data[shape_function]
                    .single_nonzero_component_index);
              const double *shape_value_ptr = &shape_values(snc, 0);
              for (unsigned int q_point = 0; q_point < n_quadrature_points;
                   ++q_point)
                values[q_point][comp] += value * (*shape_value_ptr++);
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
                       ++q_point)
                    values[q_point][comp] += value * (*shape_value_ptr++);
                }
        }
    }



    template <int dim, int spacedim, typename Number>
    void
    do_function_divergences(
      const ArrayView<Number> &                    dof_values,
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
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
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
                    Assert(false, ExcNotImplemented());

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
                        for (unsigned int j = 0; j < spacedim; ++j)
                          {
                            const unsigned int vector_component =
                              dealii::SymmetricTensor<2, spacedim>::
                                component_to_unrolled_index(
                                  TableIndices<2>(comp, j));
                            divergences[q_point][vector_component] +=
                              value * (*shape_gradient_ptr++)[j];
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
      const ArrayView<Number> &       dof_values,
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
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
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
                   ++q_point)
                values[q_point][indices] += value * (*shape_value_ptr++);
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
                       ++q_point)
                    values[q_point][indices] += value * (*shape_value_ptr++);
                }
        }
    }



    template <int dim, int spacedim, typename Number>
    void
    do_function_divergences(
      const ArrayView<Number> &                    dof_values,
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
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
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
                    Assert(false, ExcNotImplemented());
                  }
            }
        }
    }



    template <int dim, int spacedim, typename Number>
    void
    do_function_gradients(
      const ArrayView<Number> &                    dof_values,
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
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
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
                    Assert(false, ExcNotImplemented());
                  }
            }
        }
    }

  } // end of namespace internal



  template <int dim, int spacedim>
  template <typename Number>
  void
  Scalar<dim, spacedim>::get_function_values(
    const ReadVector<Number> &                fe_function,
    std::vector<solution_value_type<Number>> &values) const
  {
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));

    // get function values of dofs on this cell and call internal worker
    // function
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_values<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_values,
      shape_function_data,
      values);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_function_values_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_values<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_values,
      shape_function_data,
      values);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Scalar<dim, spacedim>::get_function_gradients(
    const ReadVector<Number> &                   fe_function,
    std::vector<solution_gradient_type<Number>> &gradients) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_derivatives<1, dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      gradients);
  }



  template <int dim, int spacedim>
  template <typename InputVector>
  void
  Scalar<dim, spacedim>::get_function_gradients_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_gradient_type<typename InputVector::value_type>>
      &gradients) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_derivatives<1, dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      gradients);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Scalar<dim, spacedim>::get_function_hessians(
    const ReadVector<Number> &                  fe_function,
    std::vector<solution_hessian_type<Number>> &hessians) const
  {
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_derivatives<2, dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_hessians,
      shape_function_data,
      hessians);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_function_hessians_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_hessian_type<typename InputVector::value_type>>
      &hessians) const
  {
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_derivatives<2, dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_hessians,
      shape_function_data,
      hessians);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Scalar<dim, spacedim>::get_function_laplacians(
    const ReadVector<Number> &                    fe_function,
    std::vector<solution_laplacian_type<Number>> &laplacians) const
  {
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_laplacians<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_hessians,
      shape_function_data,
      laplacians);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_function_laplacians_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_laplacian_type<typename InputVector::value_type>>
      &laplacians) const
  {
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_laplacians<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_hessians,
      shape_function_data,
      laplacians);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Scalar<dim, spacedim>::get_function_third_derivatives(
    const ReadVector<Number> &                           fe_function,
    std::vector<solution_third_derivative_type<Number>> &third_derivatives)
    const
  {
    Assert(fe_values->update_flags & update_3rd_derivatives,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_3rd_derivatives")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_derivatives<3, dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_3rd_derivatives,
      shape_function_data,
      third_derivatives);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Scalar<dim, spacedim>::get_function_third_derivatives_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<
      solution_third_derivative_type<typename InputVector::value_type>>
      &third_derivatives) const
  {
    Assert(fe_values->update_flags & update_3rd_derivatives,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_3rd_derivatives")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_derivatives<3, dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_3rd_derivatives,
      shape_function_data,
      third_derivatives);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_values(
    const ReadVector<Number> &                fe_function,
    std::vector<solution_value_type<Number>> &values) const
  {
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_values<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_values,
      shape_function_data,
      values);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_function_values_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_values<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_values,
      shape_function_data,
      values);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_gradients(
    const ReadVector<Number> &                   fe_function,
    std::vector<solution_gradient_type<Number>> &gradients) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_derivatives<1, dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      gradients);
  }



  template <int dim, int spacedim>
  template <typename InputVector>
  void
  Vector<dim, spacedim>::get_function_gradients_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_gradient_type<typename InputVector::value_type>>
      &gradients) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_derivatives<1, dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      gradients);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_symmetric_gradients(
    const ReadVector<Number> &                             fe_function,
    std::vector<solution_symmetric_gradient_type<Number>> &symmetric_gradients)
    const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_symmetric_gradients<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      symmetric_gradients);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_function_symmetric_gradients_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<
      solution_symmetric_gradient_type<typename InputVector::value_type>>
      &symmetric_gradients) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_symmetric_gradients<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      symmetric_gradients);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_divergences(
    const ReadVector<Number> &                     fe_function,
    std::vector<solution_divergence_type<Number>> &divergences) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs
    // on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_divergences<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      divergences);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_function_divergences_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_divergence_type<typename InputVector::value_type>>
      &divergences) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_divergences<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      divergences);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_curls(
    const ReadVector<Number> &               fe_function,
    std::vector<solution_curl_type<Number>> &curls) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           ExcMessage("FEValues object is not reinited to any cell"));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_curls<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      curls);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_function_curls_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_curl_type<typename InputVector::value_type>> &curls)
    const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           ExcMessage("FEValues object is not reinited to any cell"));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_curls<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      curls);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_hessians(
    const ReadVector<Number> &                  fe_function,
    std::vector<solution_hessian_type<Number>> &hessians) const
  {
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_derivatives<2, dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_hessians,
      shape_function_data,
      hessians);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_function_hessians_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_hessian_type<typename InputVector::value_type>>
      &hessians) const
  {
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_derivatives<2, dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_hessians,
      shape_function_data,
      hessians);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_laplacians(
    const ReadVector<Number> &                fe_function,
    std::vector<solution_value_type<Number>> &laplacians) const
  {
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));
    Assert(laplacians.size() == fe_values->n_quadrature_points,
           ExcDimensionMismatch(laplacians.size(),
                                fe_values->n_quadrature_points));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    Assert(
      fe_function.size() == fe_values->present_cell.n_dofs_for_dof_handler(),
      ExcDimensionMismatch(fe_function.size(),
                           fe_values->present_cell.n_dofs_for_dof_handler()));

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_laplacians<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_hessians,
      shape_function_data,
      laplacians);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_function_laplacians_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_laplacian_type<typename InputVector::value_type>>
      &laplacians) const
  {
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));
    Assert(laplacians.size() == fe_values->n_quadrature_points,
           ExcDimensionMismatch(laplacians.size(),
                                fe_values->n_quadrature_points));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_laplacians<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_hessians,
      shape_function_data,
      laplacians);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_third_derivatives(
    const ReadVector<Number> &                           fe_function,
    std::vector<solution_third_derivative_type<Number>> &third_derivatives)
    const
  {
    Assert(fe_values->update_flags & update_3rd_derivatives,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_3rd_derivatives")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_derivatives<3, dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_3rd_derivatives,
      shape_function_data,
      third_derivatives);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Vector<dim, spacedim>::get_function_third_derivatives_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<
      solution_third_derivative_type<typename InputVector::value_type>>
      &third_derivatives) const
  {
    Assert(fe_values->update_flags & update_3rd_derivatives,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_3rd_derivatives")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_derivatives<3, dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_3rd_derivatives,
      shape_function_data,
      third_derivatives);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  SymmetricTensor<2, dim, spacedim>::get_function_values(
    const ReadVector<Number> &                fe_function,
    std::vector<solution_value_type<Number>> &values) const
  {
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_values<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_values,
      shape_function_data,
      values);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  SymmetricTensor<2, dim, spacedim>::get_function_values_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_values<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_values,
      shape_function_data,
      values);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  SymmetricTensor<2, dim, spacedim>::get_function_divergences(
    const ReadVector<Number> &                     fe_function,
    std::vector<solution_divergence_type<Number>> &divergences) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs
    // on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_divergences<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      divergences);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  SymmetricTensor<2, dim, spacedim>::
    get_function_divergences_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_divergence_type<typename InputVector::value_type>>
        &divergences) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_divergences<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      divergences);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Tensor<2, dim, spacedim>::get_function_values(
    const ReadVector<Number> &                fe_function,
    std::vector<solution_value_type<Number>> &values) const
  {
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_values<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_values,
      shape_function_data,
      values);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Tensor<2, dim, spacedim>::get_function_values_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_values<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_values,
      shape_function_data,
      values);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Tensor<2, dim, spacedim>::get_function_divergences(
    const ReadVector<Number> &                     fe_function,
    std::vector<solution_divergence_type<Number>> &divergences) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs
    // on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_divergences<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      divergences);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Tensor<2, dim, spacedim>::get_function_divergences_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_divergence_type<typename InputVector::value_type>>
      &divergences) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_divergences<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      divergences);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Tensor<2, dim, spacedim>::get_function_gradients(
    const ReadVector<Number> &                   fe_function,
    std::vector<solution_gradient_type<Number>> &gradients) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());

    // get function values of dofs
    // on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_gradients<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      gradients);
  }



  template <int dim, int spacedim>
  template <class InputVector>
  void
  Tensor<2, dim, spacedim>::get_function_gradients_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_gradient_type<typename InputVector::value_type>>
      &gradients) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(dof_values.size(), fe_values->dofs_per_cell);

    internal::do_function_gradients<dim, spacedim>(
      make_array_view(dof_values.begin(), dof_values.end()),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      gradients);
  }
} // namespace FEValuesViews


namespace internal
{
  namespace FEValuesViews
  {
    template <int dim, int spacedim>
    Cache<dim, spacedim>::Cache(const FEValuesBase<dim, spacedim> &fe_values)
    {
      const FiniteElement<dim, spacedim> &fe = fe_values.get_fe();

      const unsigned int n_scalars = fe.n_components();
      scalars.reserve(n_scalars);
      for (unsigned int component = 0; component < n_scalars; ++component)
        scalars.emplace_back(fe_values, component);

      // compute number of vectors that we can fit into this finite element.
      // note that this is based on the dimensionality 'dim' of the manifold,
      // not 'spacedim' of the output vector
      const unsigned int n_vectors =
        (fe.n_components() >= Tensor<1, spacedim>::n_independent_components ?
           fe.n_components() - Tensor<1, spacedim>::n_independent_components +
             1 :
           0);
      vectors.reserve(n_vectors);
      for (unsigned int component = 0; component < n_vectors; ++component)
        vectors.emplace_back(fe_values, component);

      // compute number of symmetric tensors in the same way as above
      const unsigned int n_symmetric_second_order_tensors =
        (fe.n_components() >=
             SymmetricTensor<2, spacedim>::n_independent_components ?
           fe.n_components() -
             SymmetricTensor<2, spacedim>::n_independent_components + 1 :
           0);
      symmetric_second_order_tensors.reserve(n_symmetric_second_order_tensors);
      for (unsigned int component = 0;
           component < n_symmetric_second_order_tensors;
           ++component)
        symmetric_second_order_tensors.emplace_back(fe_values, component);


      // compute number of symmetric tensors in the same way as above
      const unsigned int n_second_order_tensors =
        (fe.n_components() >= Tensor<2, spacedim>::n_independent_components ?
           fe.n_components() - Tensor<2, spacedim>::n_independent_components +
             1 :
           0);
      second_order_tensors.reserve(n_second_order_tensors);
      for (unsigned int component = 0; component < n_second_order_tensors;
           ++component)
        second_order_tensors.emplace_back(fe_values, component);
    }
  } // namespace FEValuesViews
} // namespace internal

/*------------------------------- Explicit Instantiations -------------*/

#include "fe_values_views.inst"

DEAL_II_NAMESPACE_CLOSE
