// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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



// SolutionTransfer wanted to compute interpolation matrices between
// all pairs of elements used on a mesh in the hp case. unfortunately,
// not all pairs are actually supported, e.g. between FE_Nothing and
// FE_Q, but that shouldn't matter as long as these combinations are
// never exercised on actual cells

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
#include <fstream>
#include <iostream>
#include <vector>


template <int dim>
void transfer(std::ostream &out)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  hp::FECollection<dim> fe;
  fe.push_back (FE_Q<dim> (1));
  fe.push_back (FE_Nothing<dim>());

  // create a DoFHandler on which we
  // have both cells with FE_Q as
  // well as FE_Nothing
  hp::DoFHandler<dim> dof_handler(tria);
  dof_handler.begin(0)->child(0)->set_active_fe_index(1);

  Vector<double> solution;
  ConstraintMatrix cm;
  cm.close();

  dof_handler.distribute_dofs (fe);
  solution.reinit(dof_handler.n_dofs());

  for (unsigned int i=0; i<solution.size(); ++i)
    solution(i) = i;

  SolutionTransfer<dim,Vector<double>,hp::DoFHandler<dim> > soltrans(dof_handler);

  typename Triangulation<dim>::active_cell_iterator cell=tria.begin_active(),
                                                    endc=tria.end();
  ++cell;
  ++cell;
  for (; cell!=endc; ++cell)
    cell->set_refine_flag();

  Vector<double> old_solution=solution;
  // the following line triggered an
  // exception prior to r25670
  tria.prepare_coarsening_and_refinement();
  soltrans.prepare_for_coarsening_and_refinement(old_solution);
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs (fe);
  solution.reinit(dof_handler.n_dofs());
  soltrans.interpolate(old_solution, solution);
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



