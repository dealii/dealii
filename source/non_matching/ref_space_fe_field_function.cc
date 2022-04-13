// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/non_matching/ref_space_fe_field_function.h>


DEAL_II_NAMESPACE_OPEN

namespace NonMatching
{
  DeclExceptionMsg(
    ExcCellNotSet,
    "The set_active_cell function has to be called before calling this function.");



  template <int dim, class VectorType>
  RefSpaceFEFieldFunction<dim, VectorType>::RefSpaceFEFieldFunction(
    const DoFHandler<dim> &dof_handler,
    const VectorType &     dof_values)
    : dof_handler(&dof_handler)
    , global_dof_values(&dof_values)
  {
    Assert(dof_handler.n_dofs() == dof_values.size(),
           ExcDimensionMismatch(dof_handler.n_dofs(), dof_values.size()));
  }



  template <int dim, class VectorType>
  void
  RefSpaceFEFieldFunction<dim, VectorType>::set_active_cell(
    const typename Triangulation<dim>::active_cell_iterator &cell)
  {
    Assert(
      &cell->get_triangulation() == &dof_handler->get_triangulation(),
      ExcMessage(
        "The incoming cell must belong to the triangulation associated with "
        "the DoFHandler passed to the constructor."));

    const typename DoFHandler<dim>::active_cell_iterator dof_handler_cell(
      &dof_handler->get_triangulation(),
      cell->level(),
      cell->index(),
      dof_handler);

    // Save the element and the local dof values, since this is what we need
    // to evaluate the function.

    // Check if we can use the fast path. In case we have a different
    // element from the one used before we need to set up the data
    // structures again.
    if (element != &dof_handler_cell->get_fe())
      {
        poly.clear();
        element = &dof_handler_cell->get_fe();

        if (element->n_base_elements() == 1 &&
            dealii::internal::FEPointEvaluation::is_fast_path_supported(
              *element, 0))
          {
            dealii::internal::MatrixFreeFunctions::ShapeInfo<double> shape_info;

            shape_info.reinit(QMidpoint<1>(), *element, 0);
            renumber = shape_info.lexicographic_numbering;
            poly = dealii::internal::FEPointEvaluation::get_polynomial_space(
              element->base_element(0));

            polynomials_are_hat_functions =
              (poly.size() == 2 && poly[0].value(0.) == 1. &&
               poly[0].value(1.) == 0. && poly[1].value(0.) == 0. &&
               poly[1].value(1.) == 1.);
          }
      }
    else
      element = &dof_handler_cell->get_fe();

    local_dof_indices.resize(element->dofs_per_cell);
    dof_handler_cell->get_dof_indices(local_dof_indices);

    local_dof_values.resize(element->dofs_per_cell);

    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      local_dof_values[i] =
        dealii::internal::ElementAccess<VectorType>::get(*global_dof_values,
                                                         local_dof_indices[i]);
  }



  template <int dim, class VectorType>
  bool
  RefSpaceFEFieldFunction<dim, VectorType>::cell_is_set() const
  {
    // If set cell hasn't been called the size of local_dof_values will be
    // zero.
    return local_dof_values.size() > 0;
  }



  template <int dim, class VectorType>
  double
  RefSpaceFEFieldFunction<dim, VectorType>::value(
    const Point<dim> & point,
    const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    Assert(cell_is_set(), ExcCellNotSet());

    if (!poly.empty() && component == 0)
      {
        // TODO: this could be extended to a component that is not zero
        return dealii::internal::evaluate_tensor_product_value_and_gradient(
                 poly,
                 local_dof_values,
                 point,
                 polynomials_are_hat_functions,
                 renumber)
          .first;
      }
    else
      {
        double value = 0;
        for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
          value += local_dof_values[i] *
                   element->shape_value_component(i, point, component);

        return value;
      }
  }



  template <int dim, class VectorType>
  Tensor<1, dim>
  RefSpaceFEFieldFunction<dim, VectorType>::gradient(
    const Point<dim> & point,
    const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    Assert(cell_is_set(), ExcCellNotSet());

    if (!poly.empty() && component == 0)
      {
        // TODO: this could be extended to a component that is not zero
        return dealii::internal::evaluate_tensor_product_value_and_gradient(
                 poly,
                 local_dof_values,
                 point,
                 polynomials_are_hat_functions,
                 renumber)
          .second;
      }
    else
      {
        Tensor<1, dim> gradient;
        for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
          gradient += local_dof_values[i] *
                      element->shape_grad_component(i, point, component);

        return gradient;
      }
  }



  template <int dim, class VectorType>
  SymmetricTensor<2, dim>
  RefSpaceFEFieldFunction<dim, VectorType>::hessian(
    const Point<dim> & point,
    const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    Assert(cell_is_set(), ExcCellNotSet());

    if (!poly.empty() && component == 0)
      {
        // TODO: this could be extended to a component that is not zero
        return dealii::internal::evaluate_tensor_product_hessian(
          poly, local_dof_values, point, renumber);
      }
    else
      {
        Tensor<2, dim> hessian;
        for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
          hessian += local_dof_values[i] *
                     element->shape_grad_grad_component(i, point, component);

        return symmetrize(hessian);
      }
  }

} // namespace NonMatching
#include "ref_space_fe_field_function.inst"

DEAL_II_NAMESPACE_CLOSE
