// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test that FE_Nothing can be called with interpolate_boundary_values
// with vector elements in the hp-context with each element active
// only in a subdomain


#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FESystem<dim>(FE_Q<dim>(1), 2, FE_Nothing<dim>(), 2));
  fe_collection.push_back(FESystem<dim>(FE_Nothing<dim>(), 2, FE_Q<dim>(1), 2));

  DoFHandler<dim> dof_handler(triangulation);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();

  for (; cell != endc; ++cell)
    if (cell->center()[0] > 0)
      cell->set_active_fe_index(1);
    else
      cell->set_active_fe_index(0);

  dof_handler.distribute_dofs(fe_collection);
  deallog << dof_handler.n_dofs() << " dofs" << std::endl;

  std::map<types::global_dof_index, double> bv;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim>(4),
                                           bv);
  for (std::map<types::global_dof_index, double>::iterator p = bv.begin();
       p != bv.end();
       ++p)
    deallog << p->first << ' ' << p->second << std::endl;
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
