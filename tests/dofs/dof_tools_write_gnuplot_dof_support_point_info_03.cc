// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check DoFTools::write_gnuplot_dof_support_point_info with a mapping

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

template <int dim>
void
test()
{
  std::ostream &out = deallog.get_file_stream();

  Point<dim> center;

  Triangulation<dim> triangulation;
  GridGenerator::quarter_hyper_shell(triangulation, center, 0.5, 1.0);

  static const PolarManifold<dim, dim> manifold(center);
  triangulation.set_manifold(0, manifold);
  triangulation.set_all_manifold_ids_on_boundary(0);
  triangulation.set_all_manifold_ids(0);

  FE_Q<dim> fe(2);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  MappingQ<dim> mapping(3);

  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

  out << "set view equal xy" << std::endl
      << "plot '-' using 1:2 with lines, '-' with labels point pt 2 offset 1,1"
      << std::endl;
  GridOut().write_gnuplot(triangulation, deallog.get_file_stream());
  out << 'e' << std::endl;

  DoFTools::write_gnuplot_dof_support_point_info(deallog.get_file_stream(),
                                                 support_points);
  out << 'e' << std::endl;
}



int
main()
{
  initlog();
  test<2>();
}
