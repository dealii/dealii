// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test that MappingQEulerian knows how to compute bounding boxes of cells

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <vector>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);

  FESystem<dim>   fe(FE_Q<dim>(2), dim);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> displacements(dof_handler.n_dofs());
  displacements = 1.0;

  MappingQEulerian<dim> euler(2, dof_handler, displacements);

  // now the actual test
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      const auto p1 = cell->bounding_box().get_boundary_points();
      const auto p2 = euler.get_bounding_box(cell).get_boundary_points();
      deallog << "BBox: [" << p1.first << ", " << p1.second
              << "], with mapping [" << p2.first << ", " << p2.second << ']'
              << std::endl;
    }
}



int
main()
{
  initlog(true);

  test<1>();
  test<2>();
  test<3>();
}
