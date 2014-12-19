// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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



// test that SolutionTransfer works with FE_Nothing. This was broken
// before r24829, resulting in an exception.

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <fstream>
#include <iostream>
#include <vector>


template <int dim>
void transfer(std::ostream &out)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(5-dim);

  FESystem<dim> fe(FE_Q<dim> (1), 1,
                   FE_Nothing<dim>(), 1);
  DoFHandler<dim> dof_handler(tria);

  Vector<double> solution;
  ConstraintMatrix cm;
  cm.close();

  dof_handler.distribute_dofs (fe);
  solution.reinit(dof_handler.n_dofs());

  for (unsigned int i=0; i<solution.size(); ++i)
    solution(i) = i;

  SolutionTransfer<dim> soltrans(dof_handler);

  // test a): pure refinement
  typename Triangulation<dim>::active_cell_iterator cell=tria.begin_active(),
                                                    endc=tria.end();
  ++cell;
  ++cell;
  for (; cell!=endc; ++cell)
    cell->set_refine_flag();

  tria.prepare_coarsening_and_refinement();
  soltrans.prepare_for_pure_refinement();
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs (fe);

  Vector<double> tmp_q (dof_handler.n_dofs());
  soltrans.refine_interpolate(solution, tmp_q);
  solution.reinit (dof_handler.n_dofs());
  solution = tmp_q;

  // test b): with coarsening
  soltrans.clear();

  cell=tria.begin_active(tria.n_levels()-1);
  endc=tria.end(tria.n_levels()-1);
  cell->set_refine_flag();
  ++cell;
  for (; cell!=endc; ++cell)
    cell->set_coarsen_flag();
  Vector<double> old_solution=solution;
  tria.prepare_coarsening_and_refinement();
  soltrans.prepare_for_coarsening_and_refinement(old_solution);
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs (fe);
  solution.reinit(dof_handler.n_dofs());
  soltrans.interpolate(old_solution, solution);

  deallog << "OK" << std::endl;
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << "   1D solution transfer" << std::endl;
  transfer<1>(logfile);

  deallog << "   2D solution transfer" << std::endl;
  transfer<2>(logfile);

  deallog << "   3D solution transfer" << std::endl;
  transfer<3>(logfile);
}



