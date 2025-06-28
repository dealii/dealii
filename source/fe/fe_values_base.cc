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
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
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

/* ------------ FEValuesBase<dim,spacedim>::CellIteratorWrapper ----------- */


template <int dim, int spacedim>
FEValuesBase<dim, spacedim>::CellIteratorWrapper::CellIteratorWrapper(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell)
  : cell(cell)
{}



template <int dim, int spacedim>
FEValuesBase<dim, spacedim>::CellIteratorWrapper::CellIteratorWrapper(
  const typename DoFHandler<dim, spacedim>::cell_iterator &cell)
  : cell(cell)
{}



template <int dim, int spacedim>
FEValuesBase<dim, spacedim>::CellIteratorWrapper::CellIteratorWrapper(
  const typename DoFHandler<dim, spacedim>::level_cell_iterator &cell)
  : cell(cell)
{}



template <int dim, int spacedim>
bool
FEValuesBase<dim, spacedim>::CellIteratorWrapper::is_initialized() const
{
  return cell.has_value();
}



template <int dim, int spacedim>
FEValuesBase<dim, spacedim>::CellIteratorWrapper::
operator typename Triangulation<dim, spacedim>::cell_iterator() const
{
  Assert(is_initialized(), ExcNotReinited());

  // We can always convert to a tria iterator, regardless of which of
  // the three types of cell we store.
  return std::visit(
    [](auto &cell_iterator) ->
    typename Triangulation<dim, spacedim>::cell_iterator {
      return cell_iterator;
    },
    cell.value());
}



template <int dim, int spacedim>
types::global_dof_index
FEValuesBase<dim, spacedim>::CellIteratorWrapper::n_dofs_for_dof_handler() const
{
  Assert(is_initialized(), ExcNotReinited());

  switch (cell.value().index())
    {
      case 1:
        return std::get<1>(cell.value())->get_dof_handler().n_dofs();
      case 2:
        return std::get<2>(cell.value())->get_dof_handler().n_dofs();
      default:
        Assert(false, ExcNeedsDoFHandler());
        return numbers::invalid_dof_index;
    }
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::CellIteratorWrapper::get_interpolated_dof_values(
  const ReadVector<Number> &in,
  Vector<Number>           &out) const
{
  Assert(is_initialized(), ExcNotReinited());

  switch (cell.value().index())
    {
      case 1:
        std::get<1>(cell.value())->get_interpolated_dof_values(in, out);
        break;

      case 2:
        std::get<2>(cell.value())->get_interpolated_dof_values(in, out);
        break;

      default:
        Assert(false, ExcNeedsDoFHandler());
        break;
    }
}



/*------------------------------- FEValuesBase ---------------------------*/


template <int dim, int spacedim>
FEValuesBase<dim, spacedim>::FEValuesBase(
  const unsigned int                  n_q_points,
  const unsigned int                  dofs_per_cell,
  const UpdateFlags                   flags,
  const Mapping<dim, spacedim>       &mapping,
  const FiniteElement<dim, spacedim> &fe)
  : n_quadrature_points(n_q_points)
  , max_n_quadrature_points(n_q_points)
  , dofs_per_cell(dofs_per_cell)
  , mapping(&mapping, typeid(*this).name())
  , fe(&fe, typeid(*this).name())
  , cell_similarity(CellSimilarity::Similarity::none)
  , fe_values_views_cache(*this)
  , check_for_cell_similarity_allowed(MultithreadInfo::n_threads() == 1)
{
  Assert(n_q_points > 0,
         ExcMessage("There is nothing useful you can do with an FEValues "
                    "object when using a quadrature formula with zero "
                    "quadrature points!"));
  this->update_flags = flags;
}



template <int dim, int spacedim>
FEValuesBase<dim, spacedim>::~FEValuesBase()
{
  tria_listener_refinement.disconnect();
  tria_listener_mesh_transform.disconnect();
}



template <int dim, int spacedim>
void
FEValuesBase<dim, spacedim>::always_allow_check_for_cell_similarity(
  const bool allow)
{
  check_for_cell_similarity_allowed = allow;
}



namespace internal
{
  // put shape function part of get_function_xxx methods into separate
  // internal functions. this allows us to reuse the same code for several
  // functions (e.g. both the versions with and without indices) as well as
  // the same code for gradients and Hessians. Moreover, this speeds up
  // compilation and reduces the size of the final file since all the
  // different global vectors get channeled through the same code.

  template <typename Number, typename Number2>
  void
  do_function_values(const ArrayView<Number2>       &dof_values,
                     const dealii::Table<2, double> &shape_values,
                     std::vector<Number>            &values)
  {
    // scalar finite elements, so shape_values.size() == dofs_per_cell
    const unsigned int dofs_per_cell       = shape_values.n_rows();
    const unsigned int n_quadrature_points = values.size();

    // initialize with zero
    std::fill_n(values.begin(),
                n_quadrature_points,
                dealii::internal::NumberType<Number>::value(0.0));

    // add up contributions of trial functions. note that here we deal with
    // scalar finite elements, so no need to check for non-primitivity of
    // shape functions. in order to increase the speed of this function, we
    // directly access the data in the shape_values array, and increment
    // pointers for accessing the data. this saves some lookup time and
    // indexing. moreover, the order of the loops is such that we can access
    // the shape_values data stored contiguously
    for (unsigned int shape_func = 0; shape_func < dofs_per_cell; ++shape_func)
      {
        const Number2 value = dof_values[shape_func];
        // For auto-differentiable numbers, the fact that a DoF value is zero
        // does not imply that its derivatives are zero as well. So we
        // can't filter by value for these number types.
        if (!Differentiation::AD::is_ad_number<Number2>::value)
          if (value == dealii::internal::NumberType<Number2>::value(0.0))
            continue;

        const double *shape_value_ptr = &shape_values(shape_func, 0);
        for (unsigned int point = 0; point < n_quadrature_points; ++point)
          values[point] += value * (*shape_value_ptr++);
      }
  }



  template <int dim, int spacedim, typename VectorType>
  void
  do_function_values(
    const ArrayView<typename VectorType::value_type> &dof_values,
    const dealii::Table<2, double>                   &shape_values,
    const FiniteElement<dim, spacedim>               &fe,
    const std::vector<unsigned int> &shape_function_to_row_table,
    const ArrayView<VectorType>     &values,
    const bool                       quadrature_points_fastest = false,
    const unsigned int               component_multiple        = 1)
  {
    using Number = typename VectorType::value_type;
    // initialize with zero
    for (unsigned int i = 0; i < values.size(); ++i)
      std::fill_n(values[i].begin(),
                  values[i].size(),
                  typename VectorType::value_type());

    // see if there the current cell has DoFs at all, and if not
    // then there is nothing else to do.
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    if (dofs_per_cell == 0)
      return;

    const unsigned int n_quadrature_points =
      quadrature_points_fastest ? values[0].size() : values.size();
    const unsigned int n_components = fe.n_components();

    // Assert that we can write all components into the result vectors
    const unsigned result_components = n_components * component_multiple;
    (void)result_components;
    if (quadrature_points_fastest)
      {
        AssertDimension(values.size(), result_components);
        for (unsigned int i = 0; i < values.size(); ++i)
          AssertDimension(values[i].size(), n_quadrature_points);
      }
    else
      {
        AssertDimension(values.size(), n_quadrature_points);
        for (unsigned int i = 0; i < values.size(); ++i)
          AssertDimension(values[i].size(), result_components);
      }

    // add up contributions of trial functions.  now check whether the shape
    // function is primitive or not. if it is, then set its only non-zero
    // component, otherwise loop over components
    for (unsigned int mc = 0; mc < component_multiple; ++mc)
      for (unsigned int shape_func = 0; shape_func < dofs_per_cell;
           ++shape_func)
        {
          const Number &value = dof_values[shape_func + mc * dofs_per_cell];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
            continue;

          if (fe.is_primitive(shape_func))
            {
              const unsigned int comp =
                fe.system_to_component_index(shape_func).first +
                mc * n_components;
              const unsigned int row =
                shape_function_to_row_table[shape_func * n_components + comp];

              const double *shape_value_ptr = &shape_values(row, 0);

              if (quadrature_points_fastest)
                {
                  VectorType &values_comp = values[comp];
                  for (unsigned int point = 0; point < n_quadrature_points;
                       ++point)
                    values_comp[point] += value * (*shape_value_ptr++);
                }
              else
                for (unsigned int point = 0; point < n_quadrature_points;
                     ++point)
                  values[point][comp] += value * (*shape_value_ptr++);
            }
          else
            for (unsigned int c = 0; c < n_components; ++c)
              {
                if (fe.get_nonzero_components(shape_func)[c] == false)
                  continue;

                const unsigned int row =
                  shape_function_to_row_table[shape_func * n_components + c];

                const double      *shape_value_ptr = &shape_values(row, 0);
                const unsigned int comp            = c + mc * n_components;

                if (quadrature_points_fastest)
                  {
                    VectorType &values_comp = values[comp];
                    for (unsigned int point = 0; point < n_quadrature_points;
                         ++point)
                      values_comp[point] += value * (*shape_value_ptr++);
                  }
                else
                  for (unsigned int point = 0; point < n_quadrature_points;
                       ++point)
                    values[point][comp] += value * (*shape_value_ptr++);
              }
        }
  }



  // use the same implementation for gradients and Hessians, distinguish them
  // by the rank of the tensors
  template <int order, int spacedim, typename Number>
  void
  do_function_derivatives(
    const ArrayView<Number>                         &dof_values,
    const dealii::Table<2, Tensor<order, spacedim>> &shape_derivatives,
    std::vector<Tensor<order, spacedim, Number>>    &derivatives)
  {
    const unsigned int dofs_per_cell       = shape_derivatives.size()[0];
    const unsigned int n_quadrature_points = derivatives.size();

    // initialize with zero
    std::fill_n(derivatives.begin(),
                n_quadrature_points,
                Tensor<order, spacedim, Number>());

    // add up contributions of trial functions. note that here we deal with
    // scalar finite elements, so no need to check for non-primitivity of
    // shape functions. in order to increase the speed of this function, we
    // directly access the data in the shape_gradients/hessians array, and
    // increment pointers for accessing the data. this saves some lookup time
    // and indexing. moreover, the order of the loops is such that we can
    // access the shape_gradients/hessians data stored contiguously
    for (unsigned int shape_func = 0; shape_func < dofs_per_cell; ++shape_func)
      {
        const Number &value = dof_values[shape_func];
        // For auto-differentiable numbers, the fact that a DoF value is zero
        // does not imply that its derivatives are zero as well. So we
        // can't filter by value for these number types.
        if (dealii::internal::CheckForZero<Number>::value(value) == true)
          continue;

        const Tensor<order, spacedim> *shape_derivative_ptr =
          &shape_derivatives[shape_func][0];
        for (unsigned int point = 0; point < n_quadrature_points; ++point)
          derivatives[point] += value * (*shape_derivative_ptr++);
      }
  }



  template <int order, int dim, int spacedim, typename Number>
  void
  do_function_derivatives(
    const ArrayView<Number>                         &dof_values,
    const dealii::Table<2, Tensor<order, spacedim>> &shape_derivatives,
    const FiniteElement<dim, spacedim>              &fe,
    const std::vector<unsigned int> &shape_function_to_row_table,
    const ArrayView<std::vector<Tensor<order, spacedim, Number>>> &derivatives,
    const bool         quadrature_points_fastest = false,
    const unsigned int component_multiple        = 1)
  {
    // initialize with zero
    for (unsigned int i = 0; i < derivatives.size(); ++i)
      std::fill_n(derivatives[i].begin(),
                  derivatives[i].size(),
                  Tensor<order, spacedim, Number>());

    // see if there the current cell has DoFs at all, and if not
    // then there is nothing else to do.
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    if (dofs_per_cell == 0)
      return;


    const unsigned int n_quadrature_points =
      quadrature_points_fastest ? derivatives[0].size() : derivatives.size();
    const unsigned int n_components = fe.n_components();

    // Assert that we can write all components into the result vectors
    const unsigned result_components = n_components * component_multiple;
    (void)result_components;
    if (quadrature_points_fastest)
      {
        AssertDimension(derivatives.size(), result_components);
        for (unsigned int i = 0; i < derivatives.size(); ++i)
          AssertDimension(derivatives[i].size(), n_quadrature_points);
      }
    else
      {
        AssertDimension(derivatives.size(), n_quadrature_points);
        for (unsigned int i = 0; i < derivatives.size(); ++i)
          AssertDimension(derivatives[i].size(), result_components);
      }

    // add up contributions of trial functions.  now check whether the shape
    // function is primitive or not. if it is, then set its only non-zero
    // component, otherwise loop over components
    for (unsigned int mc = 0; mc < component_multiple; ++mc)
      for (unsigned int shape_func = 0; shape_func < dofs_per_cell;
           ++shape_func)
        {
          const Number &value = dof_values[shape_func + mc * dofs_per_cell];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
            continue;

          if (fe.is_primitive(shape_func))
            {
              const unsigned int comp =
                fe.system_to_component_index(shape_func).first +
                mc * n_components;
              const unsigned int row =
                shape_function_to_row_table[shape_func * n_components + comp];

              const Tensor<order, spacedim> *shape_derivative_ptr =
                &shape_derivatives[row][0];

              if (quadrature_points_fastest)
                for (unsigned int point = 0; point < n_quadrature_points;
                     ++point)
                  derivatives[comp][point] += value * (*shape_derivative_ptr++);
              else
                for (unsigned int point = 0; point < n_quadrature_points;
                     ++point)
                  derivatives[point][comp] += value * (*shape_derivative_ptr++);
            }
          else
            for (unsigned int c = 0; c < n_components; ++c)
              {
                if (fe.get_nonzero_components(shape_func)[c] == false)
                  continue;

                const unsigned int row =
                  shape_function_to_row_table[shape_func * n_components + c];

                const Tensor<order, spacedim> *shape_derivative_ptr =
                  &shape_derivatives[row][0];
                const unsigned int comp = c + mc * n_components;

                if (quadrature_points_fastest)
                  for (unsigned int point = 0; point < n_quadrature_points;
                       ++point)
                    derivatives[comp][point] +=
                      value * (*shape_derivative_ptr++);
                else
                  for (unsigned int point = 0; point < n_quadrature_points;
                       ++point)
                    derivatives[point][comp] +=
                      value * (*shape_derivative_ptr++);
              }
        }
  }



  template <int spacedim, typename Number, typename Number2>
  void
  do_function_laplacians(
    const ArrayView<Number2>                    &dof_values,
    const dealii::Table<2, Tensor<2, spacedim>> &shape_hessians,
    std::vector<Number>                         &laplacians)
  {
    const unsigned int dofs_per_cell       = shape_hessians.size()[0];
    const unsigned int n_quadrature_points = laplacians.size();

    // initialize with zero
    std::fill_n(laplacians.begin(),
                n_quadrature_points,
                dealii::internal::NumberType<Number>::value(0.0));

    // add up contributions of trial functions. note that here we deal with
    // scalar finite elements and also note that the Laplacian is
    // the trace of the Hessian.
    for (unsigned int shape_func = 0; shape_func < dofs_per_cell; ++shape_func)
      {
        const Number2 value = dof_values[shape_func];
        // For auto-differentiable numbers, the fact that a DoF value is zero
        // does not imply that its derivatives are zero as well. So we
        // can't filter by value for these number types.
        if (!Differentiation::AD::is_ad_number<Number2>::value)
          if (value == dealii::internal::NumberType<Number2>::value(0.0))
            continue;

        const Tensor<2, spacedim> *shape_hessian_ptr =
          &shape_hessians[shape_func][0];
        for (unsigned int point = 0; point < n_quadrature_points; ++point)
          laplacians[point] += value * trace(*shape_hessian_ptr++);
      }
  }



  template <int dim, int spacedim, typename VectorType, typename Number>
  void
  do_function_laplacians(
    const ArrayView<Number>                     &dof_values,
    const dealii::Table<2, Tensor<2, spacedim>> &shape_hessians,
    const FiniteElement<dim, spacedim>          &fe,
    const std::vector<unsigned int>             &shape_function_to_row_table,
    std::vector<VectorType>                     &laplacians,
    const bool         quadrature_points_fastest = false,
    const unsigned int component_multiple        = 1)
  {
    // initialize with zero
    for (unsigned int i = 0; i < laplacians.size(); ++i)
      std::fill_n(laplacians[i].begin(),
                  laplacians[i].size(),
                  typename VectorType::value_type());

    // see if there the current cell has DoFs at all, and if not
    // then there is nothing else to do.
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    if (dofs_per_cell == 0)
      return;


    const unsigned int n_quadrature_points = laplacians.size();
    const unsigned int n_components        = fe.n_components();

    // Assert that we can write all components into the result vectors
    const unsigned result_components = n_components * component_multiple;
    (void)result_components;
    if (quadrature_points_fastest)
      {
        AssertDimension(laplacians.size(), result_components);
        for (unsigned int i = 0; i < laplacians.size(); ++i)
          AssertDimension(laplacians[i].size(), n_quadrature_points);
      }
    else
      {
        AssertDimension(laplacians.size(), n_quadrature_points);
        for (unsigned int i = 0; i < laplacians.size(); ++i)
          AssertDimension(laplacians[i].size(), result_components);
      }

    // add up contributions of trial functions.  now check whether the shape
    // function is primitive or not. if it is, then set its only non-zero
    // component, otherwise loop over components
    for (unsigned int mc = 0; mc < component_multiple; ++mc)
      for (unsigned int shape_func = 0; shape_func < dofs_per_cell;
           ++shape_func)
        {
          const Number &value = dof_values[shape_func + mc * dofs_per_cell];
          // For auto-differentiable numbers, the fact that a DoF value is zero
          // does not imply that its derivatives are zero as well. So we
          // can't filter by value for these number types.
          if (dealii::internal::CheckForZero<Number>::value(value) == true)
            continue;

          if (fe.is_primitive(shape_func))
            {
              const unsigned int comp =
                fe.system_to_component_index(shape_func).first +
                mc * n_components;
              const unsigned int row =
                shape_function_to_row_table[shape_func * n_components + comp];

              const Tensor<2, spacedim> *shape_hessian_ptr =
                &shape_hessians[row][0];
              if (quadrature_points_fastest)
                {
                  VectorType &laplacians_comp = laplacians[comp];
                  for (unsigned int point = 0; point < n_quadrature_points;
                       ++point)
                    laplacians_comp[point] +=
                      value * trace(*shape_hessian_ptr++);
                }
              else
                for (unsigned int point = 0; point < n_quadrature_points;
                     ++point)
                  laplacians[point][comp] +=
                    value * trace(*shape_hessian_ptr++);
            }
          else
            for (unsigned int c = 0; c < n_components; ++c)
              {
                if (fe.get_nonzero_components(shape_func)[c] == false)
                  continue;

                const unsigned int row =
                  shape_function_to_row_table[shape_func * n_components + c];

                const Tensor<2, spacedim> *shape_hessian_ptr =
                  &shape_hessians[row][0];
                const unsigned int comp = c + mc * n_components;

                if (quadrature_points_fastest)
                  {
                    VectorType &laplacians_comp = laplacians[comp];
                    for (unsigned int point = 0; point < n_quadrature_points;
                         ++point)
                      laplacians_comp[point] +=
                        value * trace(*shape_hessian_ptr++);
                  }
                else
                  for (unsigned int point = 0; point < n_quadrature_points;
                       ++point)
                    laplacians[point][comp] +=
                      value * trace(*shape_hessian_ptr++);
              }
        }
  }
} // namespace internal



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_values(
  const ReadVector<Number> &fe_function,
  std::vector<Number>      &values) const
{
  Assert(this->update_flags & update_values,
         ExcAccessToUninitializedField("update_values"));
  AssertDimension(fe->n_components(), 1);
  Assert(present_cell.is_initialized(), ExcNotReinited());
  AssertDimension(fe_function.size(), present_cell.n_dofs_for_dof_handler());

  // get function values of dofs on this cell
  Vector<Number> dof_values(dofs_per_cell);
  present_cell.get_interpolated_dof_values(fe_function, dof_values);
  internal::do_function_values(make_array_view(dof_values.begin(),
                                               dof_values.end()),
                               this->finite_element_output.shape_values,
                               values);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_values(
  const ReadVector<Number>                       &fe_function,
  const ArrayView<const types::global_dof_index> &indices,
  std::vector<Number>                            &values) const
{
  Assert(this->update_flags & update_values,
         ExcAccessToUninitializedField("update_values"));
  AssertDimension(fe->n_components(), 1);
  AssertDimension(indices.size(), dofs_per_cell);

  boost::container::small_vector<Number, 200> dof_values(dofs_per_cell);
  auto view = make_array_view(dof_values.begin(), dof_values.end());
  fe_function.extract_subvector_to(indices, view);
  internal::do_function_values(view,
                               this->finite_element_output.shape_values,
                               values);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_values(
  const ReadVector<Number>    &fe_function,
  std::vector<Vector<Number>> &values) const
{
  Assert(present_cell.is_initialized(), ExcNotReinited());

  Assert(this->update_flags & update_values,
         ExcAccessToUninitializedField("update_values"));
  AssertDimension(fe_function.size(), present_cell.n_dofs_for_dof_handler());

  // get function values of dofs on this cell
  Vector<Number> dof_values(dofs_per_cell);
  present_cell.get_interpolated_dof_values(fe_function, dof_values);
  internal::do_function_values(
    make_array_view(dof_values.begin(), dof_values.end()),
    this->finite_element_output.shape_values,
    *fe,
    this->finite_element_output.shape_function_to_row_table,
    make_array_view(values.begin(), values.end()));
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_values(
  const ReadVector<Number>                       &fe_function,
  const ArrayView<const types::global_dof_index> &indices,
  std::vector<Vector<Number>>                    &values) const
{
  // Size of indices must be a multiple of dofs_per_cell such that an integer
  // number of function values is generated in each point.
  Assert(indices.size() % dofs_per_cell == 0,
         ExcNotMultiple(indices.size(), dofs_per_cell));
  Assert(this->update_flags & update_values,
         ExcAccessToUninitializedField("update_values"));

  boost::container::small_vector<Number, 200> dof_values(dofs_per_cell);
  auto view = make_array_view(dof_values.begin(), dof_values.end());
  fe_function.extract_subvector_to(indices, view);
  internal::do_function_values(
    view,
    this->finite_element_output.shape_values,
    *fe,
    this->finite_element_output.shape_function_to_row_table,
    make_array_view(values.begin(), values.end()),
    false,
    indices.size() / dofs_per_cell);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_values(
  const ReadVector<Number>                       &fe_function,
  const ArrayView<const types::global_dof_index> &indices,
  ArrayView<std::vector<Number>>                  values,
  const bool quadrature_points_fastest) const
{
  Assert(this->update_flags & update_values,
         ExcAccessToUninitializedField("update_values"));

  // Size of indices must be a multiple of dofs_per_cell such that an integer
  // number of function values is generated in each point.
  Assert(indices.size() % dofs_per_cell == 0,
         ExcNotMultiple(indices.size(), dofs_per_cell));

  boost::container::small_vector<Number, 200> dof_values(dofs_per_cell);
  auto view = make_array_view(dof_values.begin(), dof_values.end());
  fe_function.extract_subvector_to(indices, view);
  internal::do_function_values(
    view,
    this->finite_element_output.shape_values,
    *fe,
    this->finite_element_output.shape_function_to_row_table,
    make_array_view(values.begin(), values.end()),
    quadrature_points_fastest,
    indices.size() / dofs_per_cell);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_gradients(
  const ReadVector<Number>                 &fe_function,
  std::vector<Tensor<1, spacedim, Number>> &gradients) const
{
  Assert(this->update_flags & update_gradients,
         ExcAccessToUninitializedField("update_gradients"));
  AssertDimension(fe->n_components(), 1);
  Assert(present_cell.is_initialized(), ExcNotReinited());
  AssertDimension(fe_function.size(), present_cell.n_dofs_for_dof_handler());

  // get function values of dofs on this cell
  Vector<Number> dof_values(dofs_per_cell);
  present_cell.get_interpolated_dof_values(fe_function, dof_values);
  internal::do_function_derivatives(make_array_view(dof_values.begin(),
                                                    dof_values.end()),
                                    this->finite_element_output.shape_gradients,
                                    gradients);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_gradients(
  const ReadVector<Number>                       &fe_function,
  const ArrayView<const types::global_dof_index> &indices,
  std::vector<Tensor<1, spacedim, Number>>       &gradients) const
{
  Assert(this->update_flags & update_gradients,
         ExcAccessToUninitializedField("update_gradients"));
  AssertDimension(fe->n_components(), 1);
  AssertDimension(indices.size(), dofs_per_cell);

  boost::container::small_vector<Number, 200> dof_values(dofs_per_cell);
  auto view = make_array_view(dof_values.begin(), dof_values.end());
  fe_function.extract_subvector_to(indices, view);
  internal::do_function_derivatives(view,
                                    this->finite_element_output.shape_gradients,
                                    gradients);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_gradients(
  const ReadVector<Number>                              &fe_function,
  std::vector<std::vector<Tensor<1, spacedim, Number>>> &gradients) const
{
  Assert(this->update_flags & update_gradients,
         ExcAccessToUninitializedField("update_gradients"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  AssertDimension(fe_function.size(), present_cell.n_dofs_for_dof_handler());

  // get function values of dofs on this cell
  Vector<Number> dof_values(dofs_per_cell);
  present_cell.get_interpolated_dof_values(fe_function, dof_values);
  internal::do_function_derivatives(
    make_array_view(dof_values.begin(), dof_values.end()),
    this->finite_element_output.shape_gradients,
    *fe,
    this->finite_element_output.shape_function_to_row_table,
    make_array_view(gradients.begin(), gradients.end()));
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_gradients(
  const ReadVector<Number>                           &fe_function,
  const ArrayView<const types::global_dof_index>     &indices,
  ArrayView<std::vector<Tensor<1, spacedim, Number>>> gradients,
  const bool quadrature_points_fastest) const
{
  // Size of indices must be a multiple of dofs_per_cell such that an integer
  // number of function values is generated in each point.
  Assert(indices.size() % dofs_per_cell == 0,
         ExcNotMultiple(indices.size(), dofs_per_cell));
  Assert(this->update_flags & update_gradients,
         ExcAccessToUninitializedField("update_gradients"));

  boost::container::small_vector<Number, 200> dof_values(dofs_per_cell);
  auto view = make_array_view(dof_values.begin(), dof_values.end());
  fe_function.extract_subvector_to(indices, view);
  internal::do_function_derivatives(
    view,
    this->finite_element_output.shape_gradients,
    *fe,
    this->finite_element_output.shape_function_to_row_table,
    make_array_view(gradients.begin(), gradients.end()),
    quadrature_points_fastest,
    indices.size() / dofs_per_cell);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_hessians(
  const ReadVector<Number>                 &fe_function,
  std::vector<Tensor<2, spacedim, Number>> &hessians) const
{
  AssertDimension(fe->n_components(), 1);
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  AssertDimension(fe_function.size(), present_cell.n_dofs_for_dof_handler());

  // get function values of dofs on this cell
  Vector<Number> dof_values(dofs_per_cell);
  present_cell.get_interpolated_dof_values(fe_function, dof_values);
  internal::do_function_derivatives(make_array_view(dof_values.begin(),
                                                    dof_values.end()),
                                    this->finite_element_output.shape_hessians,
                                    hessians);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_hessians(
  const ReadVector<Number>                       &fe_function,
  const ArrayView<const types::global_dof_index> &indices,
  std::vector<Tensor<2, spacedim, Number>>       &hessians) const
{
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));
  AssertDimension(fe_function.size(), present_cell.n_dofs_for_dof_handler());
  AssertDimension(indices.size(), dofs_per_cell);

  boost::container::small_vector<Number, 200> dof_values(dofs_per_cell);
  auto view = make_array_view(dof_values.begin(), dof_values.end());
  fe_function.extract_subvector_to(indices, view);
  internal::do_function_derivatives(view,
                                    this->finite_element_output.shape_hessians,
                                    hessians);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_hessians(
  const ReadVector<Number>                              &fe_function,
  std::vector<std::vector<Tensor<2, spacedim, Number>>> &hessians,
  const bool quadrature_points_fastest) const
{
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  AssertDimension(fe_function.size(), present_cell.n_dofs_for_dof_handler());

  // get function values of dofs on this cell
  Vector<Number> dof_values(dofs_per_cell);
  present_cell.get_interpolated_dof_values(fe_function, dof_values);
  internal::do_function_derivatives(
    make_array_view(dof_values.begin(), dof_values.end()),
    this->finite_element_output.shape_hessians,
    *fe,
    this->finite_element_output.shape_function_to_row_table,
    make_array_view(hessians.begin(), hessians.end()),
    quadrature_points_fastest);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_hessians(
  const ReadVector<Number>                           &fe_function,
  const ArrayView<const types::global_dof_index>     &indices,
  ArrayView<std::vector<Tensor<2, spacedim, Number>>> hessians,
  const bool quadrature_points_fastest) const
{
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));
  Assert(indices.size() % dofs_per_cell == 0,
         ExcNotMultiple(indices.size(), dofs_per_cell));

  boost::container::small_vector<Number, 200> dof_values(indices.size());
  auto view = make_array_view(dof_values.begin(), dof_values.end());
  fe_function.extract_subvector_to(indices, view);
  internal::do_function_derivatives(
    view,
    this->finite_element_output.shape_hessians,
    *fe,
    this->finite_element_output.shape_function_to_row_table,
    make_array_view(hessians.begin(), hessians.end()),
    quadrature_points_fastest,
    indices.size() / dofs_per_cell);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_laplacians(
  const ReadVector<Number> &fe_function,
  std::vector<Number>      &laplacians) const
{
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));
  AssertDimension(fe->n_components(), 1);
  Assert(present_cell.is_initialized(), ExcNotReinited());
  AssertDimension(fe_function.size(), present_cell.n_dofs_for_dof_handler());

  // get function values of dofs on this cell
  Vector<Number> dof_values(dofs_per_cell);
  present_cell.get_interpolated_dof_values(fe_function, dof_values);
  internal::do_function_laplacians(make_array_view(dof_values.begin(),
                                                   dof_values.end()),
                                   this->finite_element_output.shape_hessians,
                                   laplacians);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_laplacians(
  const ReadVector<Number>                       &fe_function,
  const ArrayView<const types::global_dof_index> &indices,
  std::vector<Number>                            &laplacians) const
{
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));
  AssertDimension(fe->n_components(), 1);
  AssertDimension(indices.size(), dofs_per_cell);

  boost::container::small_vector<Number, 200> dof_values(dofs_per_cell);
  auto view = make_array_view(dof_values.begin(), dof_values.end());
  fe_function.extract_subvector_to(indices, view);
  internal::do_function_laplacians(view,
                                   this->finite_element_output.shape_hessians,
                                   laplacians);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_laplacians(
  const ReadVector<Number>    &fe_function,
  std::vector<Vector<Number>> &laplacians) const
{
  Assert(present_cell.is_initialized(), ExcNotReinited());
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));
  AssertDimension(fe_function.size(), present_cell.n_dofs_for_dof_handler());

  // get function values of dofs on this cell
  Vector<Number> dof_values(dofs_per_cell);
  present_cell.get_interpolated_dof_values(fe_function, dof_values);
  internal::do_function_laplacians(
    make_array_view(dof_values.begin(), dof_values.end()),
    this->finite_element_output.shape_hessians,
    *fe,
    this->finite_element_output.shape_function_to_row_table,
    laplacians);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_laplacians(
  const ReadVector<Number>                       &fe_function,
  const ArrayView<const types::global_dof_index> &indices,
  std::vector<Vector<Number>>                    &laplacians) const
{
  // Size of indices must be a multiple of dofs_per_cell such that an integer
  // number of function values is generated in each point.
  Assert(indices.size() % dofs_per_cell == 0,
         ExcNotMultiple(indices.size(), dofs_per_cell));
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));

  boost::container::small_vector<Number, 200> dof_values(indices.size());
  auto view = make_array_view(dof_values.begin(), dof_values.end());
  fe_function.extract_subvector_to(indices, view);
  internal::do_function_laplacians(
    view,
    this->finite_element_output.shape_hessians,
    *fe,
    this->finite_element_output.shape_function_to_row_table,
    laplacians,
    false,
    indices.size() / dofs_per_cell);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_laplacians(
  const ReadVector<Number>                       &fe_function,
  const ArrayView<const types::global_dof_index> &indices,
  std::vector<std::vector<Number>>               &laplacians,
  const bool quadrature_points_fastest) const
{
  Assert(indices.size() % dofs_per_cell == 0,
         ExcNotMultiple(indices.size(), dofs_per_cell));
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));

  boost::container::small_vector<Number, 200> dof_values(indices.size());
  auto view = make_array_view(dof_values.begin(), dof_values.end());
  fe_function.extract_subvector_to(indices, view);
  internal::do_function_laplacians(
    view,
    this->finite_element_output.shape_hessians,
    *fe,
    this->finite_element_output.shape_function_to_row_table,
    laplacians,
    quadrature_points_fastest,
    indices.size() / dofs_per_cell);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_third_derivatives(
  const ReadVector<Number>                 &fe_function,
  std::vector<Tensor<3, spacedim, Number>> &third_derivatives) const
{
  AssertDimension(fe->n_components(), 1);
  Assert(this->update_flags & update_3rd_derivatives,
         ExcAccessToUninitializedField("update_3rd_derivatives"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  AssertDimension(fe_function.size(), present_cell.n_dofs_for_dof_handler());

  // get function values of dofs on this cell
  Vector<Number> dof_values(dofs_per_cell);
  present_cell.get_interpolated_dof_values(fe_function, dof_values);
  internal::do_function_derivatives(
    make_array_view(dof_values.begin(), dof_values.end()),
    this->finite_element_output.shape_3rd_derivatives,
    third_derivatives);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_third_derivatives(
  const ReadVector<Number>                       &fe_function,
  const ArrayView<const types::global_dof_index> &indices,
  std::vector<Tensor<3, spacedim, Number>>       &third_derivatives) const
{
  Assert(this->update_flags & update_3rd_derivatives,
         ExcAccessToUninitializedField("update_3rd_derivatives"));
  AssertDimension(fe_function.size(), present_cell.n_dofs_for_dof_handler());
  AssertDimension(indices.size(), dofs_per_cell);

  boost::container::small_vector<Number, 200> dof_values(dofs_per_cell);
  auto view = make_array_view(dof_values.begin(), dof_values.end());
  fe_function.extract_subvector_to(indices, view);
  internal::do_function_derivatives(
    view, this->finite_element_output.shape_3rd_derivatives, third_derivatives);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_third_derivatives(
  const ReadVector<Number>                              &fe_function,
  std::vector<std::vector<Tensor<3, spacedim, Number>>> &third_derivatives,
  const bool quadrature_points_fastest) const
{
  Assert(this->update_flags & update_3rd_derivatives,
         ExcAccessToUninitializedField("update_3rd_derivatives"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  AssertDimension(fe_function.size(), present_cell.n_dofs_for_dof_handler());

  // get function values of dofs on this cell
  Vector<Number> dof_values(dofs_per_cell);
  present_cell.get_interpolated_dof_values(fe_function, dof_values);
  internal::do_function_derivatives(
    make_array_view(dof_values.begin(), dof_values.end()),
    this->finite_element_output.shape_3rd_derivatives,
    *fe,
    this->finite_element_output.shape_function_to_row_table,
    make_array_view(third_derivatives.begin(), third_derivatives.end()),
    quadrature_points_fastest);
}



template <int dim, int spacedim>
template <typename Number>
void
FEValuesBase<dim, spacedim>::get_function_third_derivatives(
  const ReadVector<Number>                           &fe_function,
  const ArrayView<const types::global_dof_index>     &indices,
  ArrayView<std::vector<Tensor<3, spacedim, Number>>> third_derivatives,
  const bool quadrature_points_fastest) const
{
  Assert(this->update_flags & update_3rd_derivatives,
         ExcAccessToUninitializedField("update_3rd_derivatives"));
  Assert(indices.size() % dofs_per_cell == 0,
         ExcNotMultiple(indices.size(), dofs_per_cell));

  boost::container::small_vector<Number, 200> dof_values(indices.size());
  auto view = make_array_view(dof_values.begin(), dof_values.end());
  fe_function.extract_subvector_to(indices, view);
  internal::do_function_derivatives(
    view,
    this->finite_element_output.shape_3rd_derivatives,
    *fe,
    this->finite_element_output.shape_function_to_row_table,
    make_array_view(third_derivatives.begin(), third_derivatives.end()),
    quadrature_points_fastest,
    indices.size() / dofs_per_cell);
}



template <int dim, int spacedim>
typename Triangulation<dim, spacedim>::cell_iterator
FEValuesBase<dim, spacedim>::get_cell() const
{
  return present_cell;
}



template <int dim, int spacedim>
const std::vector<Tensor<1, spacedim>> &
FEValuesBase<dim, spacedim>::get_normal_vectors() const
{
  Assert(this->update_flags & update_normal_vectors,
         (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
           "update_normal_vectors")));

  return this->mapping_output.normal_vectors;
}



template <int dim, int spacedim>
std::size_t
FEValuesBase<dim, spacedim>::memory_consumption() const
{
  return (sizeof(this->update_flags) +
          MemoryConsumption::memory_consumption(n_quadrature_points) +
          MemoryConsumption::memory_consumption(max_n_quadrature_points) +
          sizeof(cell_similarity) +
          MemoryConsumption::memory_consumption(dofs_per_cell) +
          MemoryConsumption::memory_consumption(mapping) +
          MemoryConsumption::memory_consumption(mapping_data) +
          MemoryConsumption::memory_consumption(*mapping_data) +
          MemoryConsumption::memory_consumption(mapping_output) +
          MemoryConsumption::memory_consumption(fe) +
          MemoryConsumption::memory_consumption(fe_data) +
          MemoryConsumption::memory_consumption(*fe_data) +
          MemoryConsumption::memory_consumption(finite_element_output));
}



template <int dim, int spacedim>
UpdateFlags
FEValuesBase<dim, spacedim>::compute_update_flags(
  const UpdateFlags update_flags) const
{
  // first find out which objects need to be recomputed on each
  // cell we visit. this we have to ask the finite element and mapping.
  // elements are first since they might require update in mapping
  //
  // there is no need to iterate since mappings will never require
  // the finite element to compute something for them
  UpdateFlags flags = update_flags | fe->requires_update_flags(update_flags);
  flags |= mapping->requires_update_flags(flags);

  return flags;
}



template <int dim, int spacedim>
void
FEValuesBase<dim, spacedim>::invalidate_present_cell()
{
  // if there is no present cell, then we shouldn't be
  // connected via a signal to a triangulation
  Assert(present_cell.is_initialized(), ExcInternalError());

  // so delete the present cell and
  // disconnect from the signal we have with
  // it
  tria_listener_refinement.disconnect();
  tria_listener_mesh_transform.disconnect();
  present_cell = {};
}



template <int dim, int spacedim>
void
FEValuesBase<dim, spacedim>::maybe_invalidate_previous_present_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell)
{
  if (present_cell.is_initialized())
    {
      if (&cell->get_triangulation() !=
          &present_cell
             .
             operator typename Triangulation<dim, spacedim>::cell_iterator()
             ->get_triangulation())
        {
          // the triangulations for the previous cell and the current cell
          // do not match. disconnect from the previous triangulation and
          // connect to the current one; also invalidate the previous
          // cell because we shouldn't be comparing cells from different
          // triangulations
          invalidate_present_cell();
          tria_listener_refinement =
            cell->get_triangulation().signals.any_change.connect(
              [this]() { this->invalidate_present_cell(); });
          tria_listener_mesh_transform =
            cell->get_triangulation().signals.mesh_movement.connect(
              [this]() { this->invalidate_present_cell(); });
        }
    }
  else
    {
      // if this FEValues has never been set to any cell at all, then
      // at least subscribe to the triangulation to get notified of
      // changes
      tria_listener_refinement =
        cell->get_triangulation().signals.post_refinement.connect(
          [this]() { this->invalidate_present_cell(); });
      tria_listener_mesh_transform =
        cell->get_triangulation().signals.mesh_movement.connect(
          [this]() { this->invalidate_present_cell(); });
    }
}



template <int dim, int spacedim>
inline void
FEValuesBase<dim, spacedim>::check_cell_similarity(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell)
{
  if (check_for_cell_similarity_allowed == false)
    {
      cell_similarity = CellSimilarity::none;
      return;
    }

  // case that there has not been any cell before
  if (this->present_cell.is_initialized() == false)
    cell_similarity = CellSimilarity::none;
  else
    // in MappingQ, data can have been modified during the last call. Then, we
    // can't use that data on the new cell.
    if (cell_similarity == CellSimilarity::invalid_next_cell)
      cell_similarity = CellSimilarity::none;
    else
      cell_similarity =
        (cell->is_translation_of(
           static_cast<
             const typename Triangulation<dim, spacedim>::cell_iterator &>(
             this->present_cell)) ?
           CellSimilarity::translation :
           CellSimilarity::none);

  if ((dim == spacedim - 1) && (cell_similarity == CellSimilarity::translation))
    {
      if (static_cast<const typename Triangulation<dim, spacedim>::cell_iterator
                        &>(this->present_cell)
            ->direction_flag() != cell->direction_flag())
        cell_similarity = CellSimilarity::inverted_translation;
    }
  // TODO: here, one could implement other checks for similarity, e.g. for
  // children of a parallelogram.
}



template <int dim, int spacedim>
CellSimilarity::Similarity
FEValuesBase<dim, spacedim>::get_cell_similarity() const
{
  return cell_similarity;
}



template <int dim, int spacedim>
const unsigned int FEValuesBase<dim, spacedim>::dimension;



template <int dim, int spacedim>
const unsigned int FEValuesBase<dim, spacedim>::space_dimension;

/*-------------------------- Explicit Instantiations -------------------------*/


#include "fe/fe_values_base.inst"

DEAL_II_NAMESPACE_CLOSE
