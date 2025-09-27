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


// Check
//   DoFTools::map_dofs_to_support_points (const Mapping<dim> &,
//                   const DoFHandler<dim> &,
//                   std::vector<Point<dim> > &)
// with an hp DoFHandler where one set of cells has no DoFs at all.
//
// The test is adapted from one posted on the mailing list by
// dav.gyulamiryan@gmail.com.


#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> triangulation;

  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(6 - dim);

  FESystem<dim> fe(FE_Q<dim>(1), dim, FE_Q<dim>(1), 1);
  FESystem<dim> fe_nothing(FE_Nothing<dim>(), dim, FE_Nothing<dim>(), 1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(fe);
  fe_collection.push_back(fe_nothing);

  DoFHandler<dim> dof_handler(triangulation);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          if (cell->center()[0] < 0.5)
            cell->set_active_fe_index(0);
          else
            cell->set_active_fe_index(1);
        }
    }

  dof_handler.distribute_dofs(fe_collection);

  const MappingQ1<dim>                                mapping;
  const std::map<types::global_dof_index, Point<dim>> support_points =
    DoFTools::map_dofs_to_support_points(mapping, dof_handler);

  for (const auto p : support_points)
    deallog << p.first << "->" << p.second << '\n';
  deallog << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();
  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
