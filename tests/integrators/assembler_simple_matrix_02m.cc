// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/**
 * @file Test initialization of Assembler::MatrixSimple and
 * DoFInfo including assigning of local block sizes with multiple matrices
 */

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/local_results.h>
#include <deal.II/meshworker/simple.h>

#include "../tests.h"


template <int dim>
void
test(FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  tr.begin_active()->neighbor(1)->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  dof.initialize_local_block_info();

  typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  typename DoFHandler<dim>::face_iterator        face = cell->face(1);

  std::vector<FullMatrix<double>>                         matrices(2);
  MeshWorker::Assembler::MatrixSimple<FullMatrix<double>> ass;
  ass.initialize(matrices);

  MeshWorker::DoFInfo<dim> info(dof.block_info());
  ass.initialize_info(info, false);

  deallog.push("cell");
  info.reinit(cell);
  info.print_debug(deallog);
  deallog.pop();

  deallog.push("face1");
  info.reinit(cell, face, 1);
  info.print_debug(deallog);
  deallog.pop();

  deallog.push("face2");
  ass.initialize_info(info, true);
  info.reinit(cell, face, 1);
  info.print_debug(deallog);
  deallog.pop();
}

int
main()
{
  const std::string logname = "output";
  std::ofstream     logfile(logname);
  deallog.attach(logfile);

  FE_DGP<2>           p0(0);
  FE_DGP<2>           p1(1);
  FE_DGP<2>           p2(2);
  FE_RaviartThomas<2> rt0(0);

  FESystem<2> sys1(p0, 1, p1, 1);
  FESystem<2> sys2(p2, 2, p0, 3, p1, 1);
  FESystem<2> sys3(p0, 2, rt0, 1);

  test(sys1);
  test(sys2);
  test(sys3);
}
