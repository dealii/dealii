// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/non_matching/closest_surface_point.h>

DEAL_II_NAMESPACE_OPEN
namespace NonMatching
{
  template <int dim, class Number>
  ClosestSurfacePoint<dim, Number>::ClosestSurfacePoint(
    const ReadVector<Number> &level_set,
    const DoFHandler<dim>    &dof_handler,
    Mapping<dim>             &mapping,
    const AdditionalData     &data)
    : data(data)
    , dof_handler(&dof_handler)
    , level_set(&level_set)
    , mapping(&mapping)
  {
    if (data.level != numbers::invalid_unsigned_int)
      {
        AssertThrow(data.level <
                      dof_handler.get_triangulation().n_global_levels(),
                    ExcMessage("Level is larger than number of levels in the "
                               "Triangulation"));
      }
    // The only search algorithm that is implemented is for
    // MappingCartesian, so we assert that the mapping is of this type.
    AssertThrow(
      dynamic_cast<const MappingCartesian<dim> *>(&mapping) != nullptr,
      ExcMessage("This class is only implemented with MappingCartesian."));
  }



  template <int dim, class Number>
  std::pair<std::vector<Point<dim>>, std::vector<Point<dim>>>
  ClosestSurfacePoint<dim, Number>::compute_closest_surface_points(
    const typename Triangulation<dim>::cell_iterator &search_cell,
    const typename Triangulation<dim>::cell_iterator &reference_cell,
    const std::vector<Point<dim>>                    &quadrature_points) const
  {
    std::vector<Point<dim>> closest_unit_search_points(quadrature_points);
    for (unsigned int q = 0; q < quadrature_points.size(); ++q)
      closest_unit_search_points[q] =
        mapping->transform_real_to_unit_cell(search_cell, quadrature_points[q]);

    std::vector<Number> dof_values_level_set(
      dof_handler->get_fe().dofs_per_cell);
    std::vector<types::global_dof_index> level_set_dof_indices(
      dof_handler->get_fe().dofs_per_cell);

    typename DoFHandler<dim>::cell_iterator dof_cell(
      &search_cell->get_triangulation(),
      search_cell->level(),
      search_cell->index(),
      dof_handler.get());

    if (data.level != numbers::invalid_unsigned_int)
      dof_cell->get_mg_dof_indices(level_set_dof_indices);
    else
      dof_cell->get_dof_indices(level_set_dof_indices);

    level_set->extract_subvector_to(level_set_dof_indices,
                                    dof_values_level_set);

    for (size_t i = 0; i < closest_unit_search_points.size(); ++i)
      {
        newton_monolithic(closest_unit_search_points[i],
                          dof_handler->get_fe(),
                          dof_values_level_set,
                          closest_unit_search_points[i]);
      }
    std::vector<Point<dim>> closest_real_points(quadrature_points.size());
    std::vector<Point<dim>> closest_unit_reference_points(
      quadrature_points.size());
    // back to absolute coordinates
    for (unsigned int q = 0; q < quadrature_points.size(); ++q)
      closest_real_points[q] =
        mapping->transform_unit_to_real_cell(search_cell,
                                             closest_unit_search_points[q]);

    for (unsigned int q = 0; q < quadrature_points.size(); ++q)
      closest_unit_reference_points[q] =
        mapping->transform_real_to_unit_cell(reference_cell,
                                             closest_real_points[q]);

    return {closest_real_points, closest_unit_reference_points};
  }



  template <int dim, class Number>
  void
  ClosestSurfacePoint<dim, Number>::newton_monolithic(
    const Point<dim>          &point,
    const FiniteElement<dim>  &fe,
    const std::vector<Number> &dof_values,
    Point<dim>                &closest_point) const
  {
    AssertDimension(dof_values.size(), fe.dofs_per_cell);

    Assert(
      fe.degree > 1,
      ExcMessage(
        "The Newton iteration to find closest surface points requires hessians "
        "that are not available when the finite element degree is 1."));

    // X, Y, Z, lambda
    Vector<double> current_solution(dim + 1);
    Vector<double> solution_update(dim + 1);

    Vector<double>     residual(dim + 1);
    FullMatrix<double> hessian(dim + 1, dim + 1);

    for (unsigned int i = 0; i < dim; ++i)
      current_solution[i] = closest_point[i];


    for (unsigned int newton_iter = 0; newton_iter < data.n_iterations;
         ++newton_iter)
      {
        hessian             = 0.0;
        residual            = 0.0;
        const double lambda = current_solution[dim];
        for (unsigned int k = 0; k < dof_values.size(); ++k)
          {
            const auto value_k = fe.shape_value(k, closest_point);
            const auto grad_k  = fe.shape_grad(k, closest_point);
            const auto hess_k  = fe.shape_grad_grad(k, closest_point);
            for (unsigned int i = 0; i < dim; ++i)
              {
                for (unsigned int j = 0; j < dim; ++j)
                  hessian(i, j) += lambda * dof_values[k] * hess_k[i][j];

                hessian(i, dim) += dof_values[k] * grad_k[i];
                hessian(dim, i) += dof_values[k] * grad_k[i];
              }

            for (unsigned int i = 0; i < dim; ++i)
              residual[i] -= lambda * dof_values[k] * grad_k[i];

            residual[dim] -= dof_values[k] * value_k;
          }


        for (unsigned int i = 0; i < dim; ++i)
          {
            residual[i] -= current_solution[i] - point[i];
            hessian[i][i] += 1.0;
          }

        if (residual.l2_norm() < data.tolerance)
          break;


        hessian.gauss_jordan();
        hessian.vmult(solution_update, residual);
        current_solution += solution_update;

        for (unsigned int i = 0; i < dim; ++i)
          closest_point[i] = current_solution[i];
      }

    // Check if the Newton iteration converged
    Assert(residual.l2_norm() < data.tolerance,
           ExcMessage("Newton iteration did not converge"));
  }

#ifndef DOXYGEN
#  include "non_matching/closest_surface_point.inst"
#endif
} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE
