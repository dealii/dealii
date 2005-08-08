//----------------------------  accessor_equality.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  accessor_equality.cc  ---------------------------

// test for equality of accessor objects

#include "../tests.h"
#include <iostream>
#include <fstream>

#include <fe/fe_q.h>
#include <dofs/dof_handler.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>



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
    Assert (dof_handler.begin_active() != dof_handler2.begin_active(),
            ExcInternalError());
  }
}


int main ()
{
  std::ofstream logfile("accessor_equality.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();
  test<3> ();
}
