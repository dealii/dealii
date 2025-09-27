// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test to check if the actual fe_index is used in SolutionTransfer

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/solution_transfer.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <int dim>
void
transfer(std::ostream &out)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(1));
  fe.push_back(FE_Q<dim>(2));

  DoFHandler<dim> dof_handler(tria);
  dof_handler.begin(0)->child(0)->set_active_fe_index(1);

  Vector<double> solution;

  dof_handler.distribute_dofs(fe);
  solution.reinit(dof_handler.n_dofs());

  for (unsigned int i = 0; i < solution.size(); ++i)
    solution(i) = i;

  SolutionTransfer<dim, Vector<double>> soltrans(dof_handler);


  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  ++cell;
  ++cell;
  for (; cell != endc; ++cell)
    cell->set_refine_flag();

  Vector<double> old_solution = solution;
  tria.prepare_coarsening_and_refinement();
  soltrans.prepare_for_coarsening_and_refinement(solution);
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs(fe);
  solution.reinit(dof_handler.n_dofs());
  soltrans.interpolate(solution);
}


int
main()
{
  initlog();

  deallog << "   1D solution transfer" << std::endl;
  transfer<1>(deallog.get_file_stream());

  deallog << "   2D solution transfer" << std::endl;
  transfer<2>(deallog.get_file_stream());

  deallog << "   3D solution transfer" << std::endl;
  transfer<3>(deallog.get_file_stream());
}
