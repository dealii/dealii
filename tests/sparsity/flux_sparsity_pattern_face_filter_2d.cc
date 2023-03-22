/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 - 2022 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */
// @sect3{A 2D test about the flux sparsity pattern for a mesh with hanging
// nodes as in sparsity/face_spartisy_pattern_2d.cc. Here, some of faces are
// excluded from flux sparsity pattern using a face filter. The
// make_flux_sparsity_pattern() with trivial dof_masks is computed but with a
// non-trivial face_has_flux_coupling predicate.}

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "deal.II/../../tests/tests.h"

using namespace dealii;

template <int dim>
bool
is_face_on_OX(const typename Triangulation<dim>::active_cell_iterator &cell,
              const unsigned int face_index)
{
  return std::abs(cell->face(face_index)->center()(1)) < 0.01;
}

template <int dim>
bool
is_face_x_between_zero_one(
  const typename Triangulation<dim>::active_cell_iterator &cell,
  const unsigned int                                       face_index)
{
  return (cell->face(face_index)->center()(0) > 0.0) &&
         (cell->face(face_index)->center()(0) < 1.0);
}

template <int dim>
void
make_grid_with_hanging_nodes(Triangulation<dim> &triangulation)
{
  GridGenerator::hyper_cube(triangulation, -2.0, 2.0);
  triangulation.refine_global(1);

  const auto refinement_subdomain_predicate = [&](const auto &cell) {
    return (cell->center()(0) > 0.0 && cell->center()(1) > 0.0);
  };

  for (auto &cell :
       triangulation.active_cell_iterators() | refinement_subdomain_predicate)
    cell->set_refine_flag();

  triangulation.execute_coarsening_and_refinement();
}

template <int dim>
void
create_and_output_flux_pattern_filtering_faces(
  DoFHandler<dim> &       dof_handler,
  DynamicSparsityPattern &dsp,
  const std::function<
    bool(const typename Triangulation<dim>::active_cell_iterator,
         const unsigned int)> &face_has_flux_coupling)
{
  Table<2, DoFTools::Coupling> coupling(1, 1);
  coupling.fill(DoFTools::always);

  DoFTools::make_flux_sparsity_pattern(dof_handler,
                                       dsp,
                                       AffineConstraints<double>(),
                                       true,
                                       coupling,
                                       coupling,
                                       numbers::invalid_subdomain_id,
                                       face_has_flux_coupling);
  dsp.print(deallog.get_file_stream());
}

template <int dim>
void
check()
{
  Triangulation<dim> triangulation;
  make_grid_with_hanging_nodes(triangulation);

  DoFHandler<dim> dof_handler(triangulation);
  const FE_Q<dim> finite_element(1);
  dof_handler.distribute_dofs(finite_element);

  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  const auto             face_has_flux_coupling = [&](const auto &       cell,
                                          const unsigned int face_index) {
    return (!is_face_on_OX<dim>(cell, face_index) ||
            is_face_x_between_zero_one<dim>(cell, face_index));
  };
  create_and_output_flux_pattern_filtering_faces(dof_handler,
                                                 dsp,
                                                 face_has_flux_coupling);
}

int
main()
{
  initlog();
  deallog.push("2d");
  check<2>();
  deallog.pop();
}
