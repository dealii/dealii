// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test that GridTools::Cache::get_cell_bounding_boxes_rtree()
// works in conjunction with MappingQEulerian

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

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

  displacements = 2.0;

  MappingQEulerian<dim> euler(2, dof_handler, displacements);

  GridTools::Cache<dim> cache0(triangulation);
  GridTools::Cache<dim> cache1(triangulation, euler);

  auto &tree0 = cache0.get_cell_bounding_boxes_rtree();
  auto &tree1 = cache1.get_cell_bounding_boxes_rtree();

  const auto box0 = tree0.begin()->first;
  const auto box1 = tree1.begin()->first;

  deallog << "No Mapping      : " << box0.get_boundary_points().first << ", "
          << box0.get_boundary_points().second << std::endl
          << "MappingQEulerian: " << box1.get_boundary_points().first << ", "
          << box1.get_boundary_points().second << std::endl;
}



int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
