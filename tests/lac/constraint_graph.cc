// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// take hp/random and generate the graph of constraints


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Q<dim>(1));
  fe_collection.push_back(FE_Q<dim>(2));
  fe_collection.push_back(FE_Q<dim>(QIterated<1>(QTrapezoid<1>(), 3)));
  fe_collection.push_back(FE_Q<dim>(QIterated<1>(QTrapezoid<1>(), 4)));

  DoFHandler<dim> dof_handler(tria);

  for (typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    cell->set_active_fe_index(Testing::rand() % fe_collection.size());

  dof_handler.distribute_dofs(fe_collection);

  AffineConstraints<double> cm;
  DoFTools::make_hanging_node_constraints(dof_handler, cm);
  cm.write_dot(deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallogfile << std::setprecision(8);

  test<3>();

  deallog << "OK" << std::endl;
}
