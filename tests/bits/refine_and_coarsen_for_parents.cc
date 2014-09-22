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



// check that, if we take an locally refined mesh, refine it globally once,
// then coarsen it globally again, the parent relation holds

#include "../tests.h"

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <fstream>


void do_refine (Triangulation<1> &tria)
{
  const int dim = 1;

  tria.refine_global (2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement ();
}


void do_refine (Triangulation<2> &tria)
{
  const int dim = 2;

  tria.refine_global (2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement ();
  tria.begin_active ()->set_refine_flag (RefinementPossibilities<dim>::cut_x);
  tria.execute_coarsening_and_refinement ();
  tria.begin_active ()->set_refine_flag (RefinementPossibilities<dim>::cut_y);
  tria.execute_coarsening_and_refinement ();
}


void do_refine (Triangulation<3> &tria)
{
  const int dim = 3;

  tria.refine_global (2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement ();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_x);
  tria.execute_coarsening_and_refinement ();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_y);
  tria.execute_coarsening_and_refinement ();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_z);
  tria.execute_coarsening_and_refinement ();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_xy);
  tria.execute_coarsening_and_refinement ();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_xz);
  tria.execute_coarsening_and_refinement ();
  tria.begin_active()->set_refine_flag(RefinementPossibilities<dim>::cut_yz);
  tria.execute_coarsening_and_refinement ();
}


template <int dim>
void check ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  do_refine (tria);
  // refine the mesh globally and
  // verify that the parent relation
  // holds
  tria.refine_global (1);

  for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
       cell != tria.end (); ++cell)
    for (unsigned int child = 0; child < cell->n_children (); ++child)
      Assert (cell->child (child)->parent () == cell,
              ExcInternalError ());

  // coarsen the mesh globally and
  // verify that the parent relation
  // holds
  for (typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active ();
       cell != tria.end (); ++cell)
    cell->set_coarsen_flag ();

  tria.execute_coarsening_and_refinement ();

  for (typename Triangulation<dim>::cell_iterator cell = tria.begin ();
       cell != tria.end(); ++cell)
    for (unsigned int child = 0; child < cell->n_children (); ++child)
      Assert (cell->child (child)->parent () == cell,
              ExcInternalError());

  deallog << "OK for " << dim << "d" << std::endl;
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
}



