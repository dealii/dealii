// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
//

// Test DoFTools::extract_dofs_with_support_on()

#include <deal.II/base/point.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <list>
#include <set>
#include <sstream>

#include "../tests.h"


template <int dim>
bool
pred_d(const typename DoFHandler<dim>::active_cell_iterator &cell)
{
  return (cell->center()[0] < 0.49 && cell->center()[1] < 0.49);
}

template <int dim>
void
test(const unsigned int flag)
{
  // Setup system
  Triangulation<dim> triangulation;

  GridGenerator::hyper_rectangle(triangulation,
                                 Point<dim>(0, 0),
                                 Point<dim>(1, 1));

  if (flag == 0)
    triangulation.refine_global(2);
  else
    triangulation.refine_global(1);

  DoFHandler<dim> dh(triangulation);

  // Extra refinement to generate hanging nodes
  for (typename DoFHandler<dim>::active_cell_iterator cell = dh.begin_active();
       cell != dh.end();
       ++cell)
    if ((flag == 1 && pred_d<dim>(cell)) || (flag == 2 && !pred_d<dim>(cell)))
      cell->set_refine_flag();

  triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement();

  FE_Q<dim> fe(2);
  dh.distribute_dofs(fe);

  AffineConstraints<double> cm;
  DoFTools::make_hanging_node_constraints(dh, cm);

  IndexSet support = DoFTools::extract_dofs_with_support_contained_within(
    dh,
    std::function<bool(const typename DoFHandler<dim>::active_cell_iterator &)>(
      &pred_d<dim>),
    cm);
  support.print(deallog);

  // print grid and DoFs for visual inspection
  if (false)
    {
      std::cout << "-------------------- " << Utilities::int_to_string(dim)
                << Utilities::int_to_string(flag) << std::endl;
      cm.print(std::cout);

      std::map<types::global_dof_index, Point<dim>> support_points;
      MappingQ1<dim>                                mapping;
      DoFTools::map_dofs_to_support_points(mapping, dh, support_points);

      const std::string filename = "grid" + Utilities::int_to_string(dim) +
                                   Utilities::int_to_string(flag) + ".gp";
      std::ofstream f(filename);

      f << "set terminal png size 400,410 enhanced font \"Helvetica,8\""
        << std::endl
        << "set output \"grid" << Utilities::int_to_string(dim)
        << Utilities::int_to_string(flag) << ".png\"" << std::endl
        << "set size square" << std::endl
        << "set view equal xy" << std::endl
        << "unset xtics" << std::endl
        << "unset ytics" << std::endl
        << "plot '-' using 1:2 with lines notitle, '-' with labels point pt 2 offset 1,1 notitle"
        << std::endl;
      GridOut().write_gnuplot(triangulation, f);
      f << 'e' << std::endl;

      DoFTools::write_gnuplot_dof_support_point_info(f, support_points);
      f << 'e' << std::endl;
    }

  dh.clear();
}

int
main(int argc, char **argv)
{
  initlog();

  test<2>(0);
  test<2>(1);
  test<2>(2);

  return 0;
}
