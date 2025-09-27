// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// document a hang in make_hanging_node_constraints with an
// FE_System with 0 components.

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <sstream>
#include <string>

#include "../tests.h"

#define PRECISION 5



template <int dim>
void
test()
{
  // 0 components is not okay
  FESystem<dim>      fe(FE_Q<dim>(1), 1, FE_Q<dim>(2), 0);
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  DoFHandler<dim> dofh(tria);

  dofh.distribute_dofs(fe);

  AffineConstraints<double> cm;

  DoFTools::make_hanging_node_constraints(dofh, cm);
  cm.close();

  std::ostringstream ss;
  cm.print(ss);

  deallog << ss.str() << std::endl;



  deallog << "ok" << std::endl;
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);

  test<2>();

  return 0;
}
