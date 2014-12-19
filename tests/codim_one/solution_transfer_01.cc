// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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



/*
  Bug for solution transfer using DG in surfaces


  set solution = 1
  refine
  do solution transfer
  And the magically, solution becomes 0

  This surprisingly turned out to be a problem with the prolongation
  matrices of DGQ in codim-1. (See the fe_dgq_prolongation_01 test.)
*/

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>

#include <fstream>

using namespace dealii;


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  const unsigned int spacedim = 2;
  const unsigned int dim = spacedim-1;

  Triangulation<dim,spacedim> boundary_mesh;

  // create a line in 2d
  GridGenerator::hyper_cube (boundary_mesh);

  // attach a piecewise constant
  // element to this one cell
  FE_DGQ<dim,spacedim> fe(0);
  DoFHandler<dim,spacedim>  dh (boundary_mesh);
  dh.distribute_dofs(fe);

  Vector<double> solution(dh.n_dofs());
  solution = 1.0;

  deallog << "Old values:" << std::endl;
  for (unsigned int i=0; i<solution.size(); i++)
    deallog << solution(i) << std::endl;


  // Do some refinement
  boundary_mesh.begin_active()->set_refine_flag ();

  SolutionTransfer<dim, Vector<double>, DoFHandler<dim,spacedim> >
  soltrans(dh);

  boundary_mesh.prepare_coarsening_and_refinement();

  soltrans.prepare_for_coarsening_and_refinement(solution);
  boundary_mesh.execute_coarsening_and_refinement ();

  dh.distribute_dofs(fe);

  // get the interpolated solution
  // back
  Vector<double> tmp(dh.n_dofs());
  tmp = 2;
  soltrans.interpolate(solution, tmp);

  deallog << "New values:" << std::endl;
  for (unsigned int i=0; i<tmp.size(); i++)
    deallog << tmp(i) << std::endl;

  return 0;
}
