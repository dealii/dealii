// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// VectorTools::interpolate_boundary_values still had bugs in 1d after
// switching to a scheme where we can assign boundary indicators also in 1d

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;

  Triangulation<dim> triangulation;
  FE_Q<dim>          fe(1);
  DoFHandler<dim>    dof_handler(triangulation);

  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.begin_active()->face(0)->set_boundary_id(10);
  triangulation.begin_active()->face(1)->set_boundary_id(20);
  triangulation.refine_global(1);

  dof_handler.distribute_dofs(fe);

  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           10,
                                           Functions::SquareFunction<dim>(),
                                           boundary_values);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           20,
                                           Functions::SquareFunction<dim>(),
                                           boundary_values);
  deallog << boundary_values.size() << std::endl;
  for (std::map<types::global_dof_index, double>::const_iterator p =
         boundary_values.begin();
       p != boundary_values.end();
       ++p)
    deallog << p->first << ' ' << p->second << std::endl;
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
