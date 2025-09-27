/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2023 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// This test considers two elements which share a common, regular face.
// The logic of make_flux_sparsity_pattern() loops over all faces and filters
// some of them according to a predicate face_has_flux_coupling(). We check
// that the only internal face in the current setup is visited exactly once and
// that the predicate face_has_flux_coupling() is evaluated once for this face.

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


template <int dim>
Triangulation<dim>
make_two_elements()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -2.0, 2.0);
  triangulation.begin_active()->set_refine_flag(
    RefinementCase<dim>::cut_axis(0));
  triangulation.execute_coarsening_and_refinement();
  return triangulation;
}


template <int dim>
bool
is_face_on_OY(const typename DoFHandler<dim>::active_cell_iterator &cell,
              const unsigned int                                    face_index)
{
  deallog
    << "This sentence should appear once when the corresponding face is visited only once on cell "
    << cell->index() << std::endl;
  return (std::abs(cell->face(face_index)->center()[0]) < 0.01);
}

template <int dim>
void
create_and_output_flux_sparsity_with_filter(
  DoFHandler<dim>                        &dof_handler,
  std::function<bool(const typename DoFHandler<dim>::active_cell_iterator,
                     const unsigned int)> filter)
{
  DynamicSparsityPattern       dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  Table<2, DoFTools::Coupling> coupling(1, 1);
  coupling.fill(DoFTools::always);

  DoFTools::make_flux_sparsity_pattern(dof_handler,
                                       dsp,
                                       AffineConstraints<double>(),
                                       true,
                                       coupling,
                                       coupling,
                                       numbers::invalid_subdomain_id,
                                       [&](const auto        &cell,
                                           const unsigned int face_index) {
                                         return filter(cell, face_index);
                                       });
  dsp.print(deallog.get_file_stream());
}

template <int dim>
void
check()
{
  const Triangulation<dim> tria = make_two_elements<dim>();

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_Q<dim>(1));

  create_and_output_flux_sparsity_with_filter<dim>(dof_handler,
                                                   is_face_on_OY<dim>);
}

int
main()
{
  initlog();
  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}
