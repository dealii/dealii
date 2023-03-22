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
void
make_anisotropic_grid(Triangulation<dim> &triangulation)
{
  GridGenerator::hyper_cube(triangulation, -2.0, 2.0);
  triangulation.begin_active()->set_refine_flag(
    RefinementCase<dim>::cut_axis(0));
  triangulation.execute_coarsening_and_refinement();

  for (auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->center()(0) > 0.0)
        cell->set_refine_flag(RefinementCase<dim>::cut_axis(0));
      else
        cell->set_refine_flag(RefinementCase<dim>::cut_axis(1));
    }
  triangulation.execute_coarsening_and_refinement();
}

template <int dim>
void
create_and_output_flux_pattern(DoFHandler<dim> &       dof_handler,
                               DynamicSparsityPattern &dsp)
{
  Table<2, DoFTools::Coupling> coupling(1, 1);
  coupling.fill(DoFTools::always);

  DoFTools::make_flux_sparsity_pattern(dof_handler,
                                       dsp,
                                       AffineConstraints<double>(),
                                       true,
                                       coupling,
                                       coupling,
                                       numbers::invalid_subdomain_id);
  dsp.print(deallog.get_file_stream());
}


template <int dim>
void
check()
{
  Triangulation<dim> triangulation;
  make_anisotropic_grid(triangulation);

  DoFHandler<dim> dof_handler(triangulation);
  const FE_Q<dim> finite_element(1);
  dof_handler.distribute_dofs(finite_element);

  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  create_and_output_flux_pattern(dof_handler, dsp);
}

int
main()
{
  initlog();
  deallog.push("2d");
  check<2>();
  deallog.pop();
}
