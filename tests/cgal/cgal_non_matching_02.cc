// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2022 by the deal.II authors
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

// Compute the "exact" sparsity pattern for two non-matching overlapping grids.


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/non_matching/coupling.h>

#include "../tests.h"
using namespace dealii;


template <int dim, int spacedim>
void
test()
{
  constexpr int                     degree = 3;
  Triangulation<spacedim, spacedim> space_tria;
  Triangulation<dim, spacedim>      embedded_tria;

  GridGenerator::hyper_cube(space_tria, -1., 1.);
  GridGenerator::hyper_cube(embedded_tria, -.35, .25);
  space_tria.refine_global(1);
  embedded_tria.refine_global(1);
  GridTools::rotate(numbers::PI_4, 0, embedded_tria);

  DoFHandler<spacedim>      space_dh(space_tria);
  DoFHandler<dim, spacedim> embedded_dh(embedded_tria);

  FE_Q<spacedim>      fe_space(1);
  FE_Q<dim, spacedim> fe_embedded(1);

  space_dh.distribute_dofs(fe_space);
  embedded_dh.distribute_dofs(fe_embedded);


  auto space_cache =
    std::make_unique<GridTools::Cache<spacedim>>(space_tria); // Q1 mapping
  auto embedded_cache = std::make_unique<GridTools::Cache<dim, spacedim>>(
    embedded_tria); // Q1 mapping

  // Compute Quadrature formulas on the intersections of the two
  const double tol             = 1e-10;
  const auto   cells_and_quads = NonMatching::compute_intersection(
    *space_cache, *embedded_cache, degree, tol);


  std::ofstream output_test_space("space_test.vtk");
  std::ofstream output_test_embedded("embedded_test.vtk");
  GridOut().write_vtk(space_tria, output_test_space);
  GridOut().write_vtk(embedded_tria, output_test_embedded);
  remove("space_test.vtk");
  remove("embedded_test.vtk");
  // Print cells ids and points are the printed

  AffineConstraints<double> constraints;
  AffineConstraints<double> embedded_constraints;
  DynamicSparsityPattern    dsp(space_dh.n_dofs(), embedded_dh.n_dofs());
  NonMatching::create_coupling_sparsity_pattern_with_exact_intersections(
    cells_and_quads,
    space_dh,
    embedded_dh,
    dsp,
    constraints,
    ComponentMask(),
    ComponentMask(),
    embedded_constraints);

  dsp.print(deallogfile);
}

int
main()
{
  initlog();
  test<3, 3>();
}
