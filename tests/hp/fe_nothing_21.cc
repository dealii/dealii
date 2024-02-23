// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// interpolate() can not deal with FE_Nothing, simplified version of
// fe_nothing_20.cc


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

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -0.5, 0.5);
  triangulation.refine_global(4);

  FESystem<dim> fe(FE_Nothing<dim>(), 1, FE_Q<dim>(1), 1);

  deallog << "n support points: " << fe.get_unit_support_points().size()
          << std::endl;


  DoFHandler<dim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);
  deallog << "   Number of active cells:       "
          << triangulation.n_active_cells() << std::endl
          << "   Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;


  Vector<double> solution(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler,
                           Functions::ZeroFunction<dim>(2),
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
