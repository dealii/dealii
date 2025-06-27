// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
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

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values_base.h>
#include <deal.II/fe/fe_values_views.h>
#include <deal.II/fe/fe_values_views_internal.h>

#include <deal.II/lac/vector.h>

#ifdef DEAL_II_WITH_ADOLC
#  include <adolc/adouble.h>
#  include <adolc/adtl.h>
#endif


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace
  {
    template <int dim, int spacedim>
    inline std::vector<unsigned int>
    make_shape_function_to_row_table(const FiniteElement<dim, spacedim> &fe)
    {
      std::vector<unsigned int> shape_function_to_row_table(
        fe.n_dofs_per_cell() * fe.n_components(),
        numbers::invalid_unsigned_int);
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



  template <int dim, int spacedim>
  template <typename Number>
  void
  Scalar<dim, spacedim>::get_function_values(
    const ReadVector<Number>                 &fe_function,
    std::vector<solution_value_type<Number>> &values) const
  {
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(values.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell and call internal worker
    // function
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_values<dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(values.size(), fe_values->n_quadrature_points);

    internal::do_function_values<dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_values,
      shape_function_data,
      values);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Scalar<dim, spacedim>::get_function_gradients(
    const ReadVector<Number>                    &fe_function,
    std::vector<solution_gradient_type<Number>> &gradients) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());
    AssertDimension(gradients.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_derivatives<1, dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(gradients.size(), fe_values->n_quadrature_points);

    internal::do_function_derivatives<1, dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      gradients);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Scalar<dim, spacedim>::get_function_hessians(
    const ReadVector<Number>                   &fe_function,
    std::vector<solution_hessian_type<Number>> &hessians) const
  {
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());
    AssertDimension(hessians.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_derivatives<2, dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(hessians.size(), fe_values->n_quadrature_points);

    internal::do_function_derivatives<2, dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_hessians,
      shape_function_data,
      hessians);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Scalar<dim, spacedim>::get_function_laplacians(
    const ReadVector<Number>                     &fe_function,
    std::vector<solution_laplacian_type<Number>> &laplacians) const
  {
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());
    AssertDimension(laplacians.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_laplacians<dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(laplacians.size(), fe_values->n_quadrature_points);

    internal::do_function_laplacians<dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_hessians,
      shape_function_data,
      laplacians);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Scalar<dim, spacedim>::get_function_third_derivatives(
    const ReadVector<Number>                            &fe_function,
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
    AssertDimension(third_derivatives.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_derivatives<3, dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(third_derivatives.size(), fe_values->n_quadrature_points);

    internal::do_function_derivatives<3, dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_3rd_derivatives,
      shape_function_data,
      third_derivatives);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_values(
    const ReadVector<Number>                 &fe_function,
    std::vector<solution_value_type<Number>> &values) const
  {
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(values.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_values<dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(values.size(), fe_values->n_quadrature_points);

    internal::do_function_values<dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_values,
      shape_function_data,
      values);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_gradients(
    const ReadVector<Number>                    &fe_function,
    std::vector<solution_gradient_type<Number>> &gradients) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());
    AssertDimension(gradients.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_derivatives<1, dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(gradients.size(), fe_values->n_quadrature_points);

    internal::do_function_derivatives<1, dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      gradients);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_symmetric_gradients(
    const ReadVector<Number>                              &fe_function,
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
    AssertDimension(symmetric_gradients.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_symmetric_gradients<dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(symmetric_gradients.size(), fe_values->n_quadrature_points);

    internal::do_function_symmetric_gradients<dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      symmetric_gradients);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_divergences(
    const ReadVector<Number>                      &fe_function,
    std::vector<solution_divergence_type<Number>> &divergences) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());
    AssertDimension(divergences.size(), fe_values->n_quadrature_points);

    // get function values of dofs
    // on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_divergences<dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(divergences.size(), fe_values->n_quadrature_points);

    internal::do_function_divergences<dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      divergences);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_curls(
    const ReadVector<Number>                &fe_function,
    std::vector<solution_curl_type<Number>> &curls) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           ExcMessage("FEValues object is not reinited to any cell"));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());
    AssertDimension(curls.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_curls<dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(curls.size(), fe_values->n_quadrature_points);

    internal::do_function_curls<dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      curls);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_hessians(
    const ReadVector<Number>                   &fe_function,
    std::vector<solution_hessian_type<Number>> &hessians) const
  {
    Assert(fe_values->update_flags & update_hessians,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_hessians")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());
    AssertDimension(hessians.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_derivatives<2, dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(hessians.size(), fe_values->n_quadrature_points);

    internal::do_function_derivatives<2, dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_hessians,
      shape_function_data,
      hessians);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_laplacians(
    const ReadVector<Number>                 &fe_function,
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
    AssertDimension(laplacians.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_laplacians<dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(laplacians.size(), fe_values->n_quadrature_points);

    internal::do_function_laplacians<dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_hessians,
      shape_function_data,
      laplacians);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Vector<dim, spacedim>::get_function_third_derivatives(
    const ReadVector<Number>                            &fe_function,
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
    AssertDimension(third_derivatives.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_derivatives<3, dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(third_derivatives.size(), fe_values->n_quadrature_points);

    internal::do_function_derivatives<3, dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_3rd_derivatives,
      shape_function_data,
      third_derivatives);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  SymmetricTensor<2, dim, spacedim>::get_function_values(
    const ReadVector<Number>                 &fe_function,
    std::vector<solution_value_type<Number>> &values) const
  {
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(values.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_values<dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(values.size(), fe_values->n_quadrature_points);

    internal::do_function_values<dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_values,
      shape_function_data,
      values);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  SymmetricTensor<2, dim, spacedim>::get_function_divergences(
    const ReadVector<Number>                      &fe_function,
    std::vector<solution_divergence_type<Number>> &divergences) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());
    AssertDimension(divergences.size(), fe_values->n_quadrature_points);

    // get function values of dofs
    // on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_divergences<dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(divergences.size(), fe_values->n_quadrature_points);

    internal::do_function_divergences<dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      divergences);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Tensor<2, dim, spacedim>::get_function_values(
    const ReadVector<Number>                 &fe_function,
    std::vector<solution_value_type<Number>> &values) const
  {
    Assert(fe_values->update_flags & update_values,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_values")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(values.size(), fe_values->n_quadrature_points);

    // get function values of dofs on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_values<dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(values.size(), fe_values->n_quadrature_points);

    internal::do_function_values<dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_values,
      shape_function_data,
      values);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Tensor<2, dim, spacedim>::get_function_divergences(
    const ReadVector<Number>                      &fe_function,
    std::vector<solution_divergence_type<Number>> &divergences) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());
    AssertDimension(divergences.size(), fe_values->n_quadrature_points);

    // get function values of dofs
    // on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_divergences<dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(divergences.size(), fe_values->n_quadrature_points);

    internal::do_function_divergences<dim, spacedim>(
      make_const_array_view(dof_values),
      fe_values->finite_element_output.shape_gradients,
      shape_function_data,
      divergences);
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Tensor<2, dim, spacedim>::get_function_gradients(
    const ReadVector<Number>                    &fe_function,
    std::vector<solution_gradient_type<Number>> &gradients) const
  {
    Assert(fe_values->update_flags & update_gradients,
           (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
             "update_gradients")));
    Assert(fe_values->present_cell.is_initialized(),
           (typename FEValuesBase<dim, spacedim>::ExcNotReinited()));
    AssertDimension(fe_function.size(),
                    fe_values->present_cell.n_dofs_for_dof_handler());
    AssertDimension(gradients.size(), fe_values->n_quadrature_points);

    // get function values of dofs
    // on this cell
    dealii::Vector<Number> dof_values(fe_values->dofs_per_cell);
    fe_values->present_cell.get_interpolated_dof_values(fe_function,
                                                        dof_values);
    internal::do_function_gradients<dim, spacedim>(
      make_const_array_view(dof_values),
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
    AssertDimension(gradients.size(), fe_values->n_quadrature_points);

    internal::do_function_gradients<dim, spacedim>(
      make_const_array_view(dof_values),
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
      scalars.resize(n_scalars);

      // compute number of vectors that we can fit into this finite element.
      // note that this is based on the dimensionality 'dim' of the manifold,
      // not 'spacedim' of the output vector
      const unsigned int n_vectors =
        (fe.n_components() >= Tensor<1, spacedim>::n_independent_components ?
           fe.n_components() - Tensor<1, spacedim>::n_independent_components +
             1 :
           0);
      vectors.resize(n_vectors);

      // compute number of symmetric tensors in the same way as above
      const unsigned int n_symmetric_second_order_tensors =
        (fe.n_components() >=
             SymmetricTensor<2, spacedim>::n_independent_components ?
           fe.n_components() -
             SymmetricTensor<2, spacedim>::n_independent_components + 1 :
           0);
      symmetric_second_order_tensors.resize(n_symmetric_second_order_tensors);

      // compute number of symmetric tensors in the same way as above
      const unsigned int n_second_order_tensors =
        (fe.n_components() >= Tensor<2, spacedim>::n_independent_components ?
           fe.n_components() - Tensor<2, spacedim>::n_independent_components +
             1 :
           0);
      second_order_tensors.resize(n_second_order_tensors);
    }
  } // namespace FEValuesViews
} // namespace internal

/*------------------------------- Explicit Instantiations -------------*/

#include "fe/fe_values_views.inst"

DEAL_II_NAMESPACE_CLOSE
