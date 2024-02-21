// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// there used to be a bug in MappingCartesian, where an assert was triggered
// for valid cells if the cells where too small compared to the dimension
// of the triangulation. The bug would not occur for coarse grids,
// but would reliably be triggered with higher refinement level.
// In this test we move the origin of the box to trigger the bug
// easier, but it would also occur with a zero origin (just requiring higher
// refinement levels).

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
check(const Triangulation<dim> &tria)
{
  MappingCartesian<dim> mapping;
  FE_Q<dim>             fe(1);
  DoFHandler<dim>       dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  QGauss<dim> quadrature(1);

  UpdateFlags update_flags = update_quadrature_points | update_values;

  FEValues<dim> fe_values(mapping, fe, quadrature, update_flags);

  for (auto cell : dof_handler.active_cell_iterators())
    fe_values.reinit(cell);

  deallog << std::endl;

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  {
    Triangulation<2>          coarse_grid;
    std::vector<unsigned int> rep_vec({9, 2});
    Point<2>                  box_origin(3000000., 660000.);
    Point<2>                  extents(3000000., 660000.);
    GridGenerator::subdivided_hyper_rectangle(
      coarse_grid, rep_vec, box_origin, box_origin + extents, true);
    coarse_grid.refine_global(2);
    check(coarse_grid);
  }

  {
    Triangulation<3>          coarse_grid;
    std::vector<unsigned int> rep_vec({9, 9, 2});
    Point<3>                  box_origin(3000000., 3000000., 660000.);
    Point<3>                  extents(3000000., 3000000., 660000.);
    GridGenerator::subdivided_hyper_rectangle(
      coarse_grid, rep_vec, box_origin, box_origin + extents, true);
    coarse_grid.refine_global(2);
    check(coarse_grid);
  }
}
