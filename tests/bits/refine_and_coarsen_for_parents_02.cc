// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that, if we take an locally refined mesh, refine it globally once,
// then coarsen it globally again, the parent relation holds

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



void
do_refine(Triangulation<1> &tria)
{
  tria.refine_global(2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
}


void
do_refine(Triangulation<2> &tria)
{
  const int dim = 2;

  tria.refine_global(2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_x);
  tria.execute_coarsening_and_refinement();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_y);
  tria.execute_coarsening_and_refinement();
}


void
do_refine(Triangulation<3> &tria)
{
  const int dim = 3;

  tria.refine_global(2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_x);
  tria.execute_coarsening_and_refinement();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_y);
  tria.execute_coarsening_and_refinement();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_z);
  tria.execute_coarsening_and_refinement();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_xy);
  tria.execute_coarsening_and_refinement();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_xz);
  tria.execute_coarsening_and_refinement();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_yz);
  tria.execute_coarsening_and_refinement();
}


template <int dim>
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  do_refine(tria);
  // refine the mesh globally and
  // verify that the parent relation
  // holds
  tria.refine_global(1);

  DoFHandler<dim> dof_handler(tria);

  for (typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin();
       cell != dof_handler.end();
       ++cell)
    for (unsigned int child = 0; child < cell->n_children(); ++child)
      AssertThrow(cell->child(child)->parent() == cell, ExcInternalError());

  // coarsen the mesh globally and
  // verify that the parent relation
  // holds
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell)
    cell->set_coarsen_flag();

  tria.execute_coarsening_and_refinement();

  for (typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin();
       cell != dof_handler.end();
       ++cell)
    for (unsigned int child = 0; child < cell->n_children(); ++child)
      AssertThrow(cell->child(child)->parent() == cell, ExcInternalError());

  deallog << "OK for " << dim << 'd' << std::endl;
}


int
main()
{
  initlog();

  check<1>();
  check<2>();
  check<3>();
}
