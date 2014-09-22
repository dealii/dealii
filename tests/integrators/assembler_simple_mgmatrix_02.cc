// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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


/**
 * @file Test initialization of Assembler::MatrixSimple and
 * DoFInfo including assigning of local block sizes
 */

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/meshworker/local_results.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/simple.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

using namespace dealii;

template <int dim>
void test(FiniteElement<dim> &fe)
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
  typename DoFHandler<dim>::face_iterator face = cell->face(1);

  MeshWorker::Assembler::MGMatrixSimple<FullMatrix<double> > ass;
  ass.initialize_local_blocks(dof.block_info().local());
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

int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console (0);

  FE_DGP<2> p0(0);
  FE_DGP<2> p1(1);
  FE_DGP<2> p2(2);
  FE_RaviartThomas<2> rt0(0);

  FESystem<2> sys1(p0,1, p1, 1);
  FESystem<2> sys2(p2,2,p0,3,p1,1);
  FESystem<2> sys3(p0, 2, rt0, 1);

  test(sys1);
  test(sys2);
  test(sys3);
}
