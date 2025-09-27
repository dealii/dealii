// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// interpolate() can not deal with FE_Nothing in an hp-setting


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

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

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -0.5, 0.5);
  triangulation.refine_global(4);

  hp::FECollection<dim> fe_collection;

  fe_collection.push_back(
    FESystem<dim>(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1, FE_Nothing<dim>(), dim));

  fe_collection.push_back(
    FESystem<dim>(FE_Nothing<dim>(dim + 1), 1, FE_Q<dim>(2), dim));

  DoFHandler<dim> dof_handler(triangulation);



  dof_handler.distribute_dofs(fe_collection);

  deallog << "   Number of active cells:       "
          << triangulation.n_active_cells() << std::endl
          << "   Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;


  Vector<double> solution(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler,
                           Functions::ZeroFunction<dim>(2 * dim + 1),
                           solution);

  deallog << "l2_norm = " << solution.l2_norm() << std::endl;
}



int
main()
{
  initlog();

  // test<1> ();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
