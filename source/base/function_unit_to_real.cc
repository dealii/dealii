// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2021 by the deal.II authors
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

#include <deal.II/base/function_unit_to_real.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  template <int dim>
  FunctionUnitToReal<dim>::FunctionUnitToReal(const Function<dim> &function)
    : Function<dim>(function.n_components)
    , function(&function)
  {
    cell_index = numbers::invalid_unsigned_int;
    cell_level = numbers::invalid_unsigned_int;
  }



  template <int dim>
  void
  FunctionUnitToReal<dim>::set_active_cell(
    const typename Triangulation<dim>::active_cell_iterator &cell)
  {
    cell_level    = cell->level();
    cell_index    = cell->index();
    triangulation = &(cell->get_triangulation());
  }



  template <int dim>
  double
  FunctionUnitToReal<dim>::value(const Point<dim> & point,
                                 const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);

    typename Triangulation<dim>::active_cell_iterator cell(triangulation,
                                                           cell_level,
                                                           cell_index);

    FEValues<dim> fe_values(mapping,
                            dummy_element,
                            Quadrature<dim>(point),
                            update_quadrature_points);
    fe_values.reinit(cell);

    const Point<dim> &point_in_real_space = fe_values.quadrature_point(0);

    return function->value(point_in_real_space, component);
  }



  template <int dim>
  Tensor<1, dim>
  FunctionUnitToReal<dim>::gradient(const Point<dim> & point,
                                    const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);

    typename Triangulation<dim>::active_cell_iterator cell(triangulation,
                                                           cell_level,
                                                           cell_index);

    const UpdateFlags update_flags =
      update_jacobians | update_quadrature_points;

    FEValues<dim> fe_values(mapping,
                            dummy_element,
                            Quadrature<dim>(point),
                            update_flags);
    fe_values.reinit(cell);

    const Point<dim> &point_in_real_space = fe_values.quadrature_point(0);

    const Tensor<2, dim> &jacobian = fe_values.jacobian(0);

    const Tensor<1, dim> real_gradient =
      function->gradient(point_in_real_space, component);

    return contract<0, 0>(jacobian, real_gradient);
  }



  template <int dim>
  SymmetricTensor<2, dim>
  FunctionUnitToReal<dim>::hessian(const Point<dim> & point,
                                   const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);

    typename Triangulation<dim>::active_cell_iterator cell(triangulation,
                                                           cell_level,
                                                           cell_index);

    const UpdateFlags update_flags =
      update_jacobians | update_jacobian_grads | update_quadrature_points;
    FEValues<dim> fe_values(mapping,
                            dummy_element,
                            Quadrature<dim>(point),
                            update_flags);
    fe_values.reinit(cell);

    const Point<dim> &point_in_real_space = fe_values.quadrature_point(0);

    const Tensor<2, dim> &jacobian      = fe_values.jacobian(0);
    const Tensor<3, dim> &jacobian_grad = fe_values.jacobian_grad(0);

    const Tensor<1, dim> real_gradient =
      function->gradient(point_in_real_space, component);
    const SymmetricTensor<2, dim> real_hessian =
      function->hessian(point_in_real_space, component);

    Tensor<2, dim> ref_space_hessian =
      contract<0, 0>(jacobian, real_hessian * jacobian);
    ref_space_hessian += contract<1, 0>(jacobian_grad, real_gradient);

    // The result should already be symmetric up to floating point.
    return symmetrize(ref_space_hessian);
  }

#include "function_unit_to_real.inst"

} // namespace Functions

DEAL_II_NAMESPACE_CLOSE
