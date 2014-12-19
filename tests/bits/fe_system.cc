// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// document a hang in make_hanging_node_constraints with an
// FE_System with 0 components.

#include "../tests.h"
#include <sstream>
#include <deal.II/base/logstream.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>



#include <fstream>
#include <string>

#define PRECISION 5




template<int dim>
void test()
{
  // 0 components is not okay
  FESystem<dim> fe(FE_Q<dim>(1), 1, FE_Q<dim>(2), 0);
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  DoFHandler<dim> dofh(tria);

  dofh.distribute_dofs(fe);

  ConstraintMatrix cm;

  DoFTools::make_hanging_node_constraints (dofh, cm);
  cm.close ();

  std::ostringstream ss;
  cm.print(ss);

  deallog << ss.str() << std::endl;





  deallog << "ok" << std::endl;
}


int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2>();

  return 0;
}



