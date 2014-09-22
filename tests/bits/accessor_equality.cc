// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


// test for equality of accessor objects

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>



template <int dim>
void test ()
{
  deallog << dim << 'd' << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria,0,1);
  tria.refine_global (1);

  deallog << tria.begin_active() << std::endl;

  // test a few comparison operators
  {
    const typename Triangulation<dim>::active_cell_iterator
    cell = tria.begin_active();
    Assert (cell == cell, ExcInternalError());
    Assert (! (cell != cell), ExcInternalError());
    Assert (! (cell < cell), ExcInternalError());
  }

  // same with non-active iterators
  {
    const typename Triangulation<dim>::cell_iterator
    cell = tria.begin();
    Assert (cell == cell, ExcInternalError());
    Assert (! (cell != cell), ExcInternalError());
    Assert (! (cell < cell), ExcInternalError());
  }

  FE_Q<dim> fe (1);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  // now do the same tests with the
  // DoFHandler iterators
  {
    const typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active();
    Assert (cell == cell, ExcInternalError());
    Assert (! (cell != cell), ExcInternalError());
    Assert (! (cell < cell), ExcInternalError());
  }
  {
    const typename DoFHandler<dim>::cell_iterator
    cell = dof_handler.begin();
    Assert (cell == cell, ExcInternalError());
    Assert (! (cell != cell), ExcInternalError());
    Assert (! (cell < cell), ExcInternalError());
  }

  // finally check that two iterators
  // pointing into different DoFHandlers
  // can't be the same
  //
  // this test accidentally failed before a
  // fix was checked in on 2005-08-08 that
  // completely forbids comparing iterators
  // into different DoFHandlers at all
  {
    DoFHandler<dim> dof_handler2 (tria);
    dof_handler2.distribute_dofs (fe);
    try
      {
        (dof_handler.begin_active() != dof_handler2.begin_active());
      }
    catch (ExceptionBase &e)
      {
        deallog << e.get_exc_name() << std::endl;
      }
  }
}


int main ()
{
  deal_II_exceptions::disable_abort_on_exception();

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();
  test<3> ();
}
