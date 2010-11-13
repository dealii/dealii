//----------------------------  dof_accessor_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_accessor_01.cc  ---------------------------

// verify that iterators can be copied. when you just copy the
// iterator itself, everything is alright: the copy operator of the
// iterator calls Accessor::copy_from.
//
// the test was originally written because at that time it was
// possible to do something like *it1 = *it2 for DoF iterators --
// likely not what the author intended since it does not copy the
// cell2 into the cell1, only the accessor. furthermore, the operator=
// was wrongly implemented, so it was removed in the process


#include "../tests.h"
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>

#include <fe/fe_system.h>
#include <fe/fe_q.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>

#include <fstream>


template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);

  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active();
  ++cell;

  typename DoFHandler<dim>::face_iterator
    face = dof_handler.begin_active()->face(0);
  face = cell->face(0);
  deallog << cell->face(0) << ' ' << face << std::endl;
  Assert (cell->face(0) == face, ExcInternalError());
  Assert (!(cell->face(0) != face), ExcInternalError());
}



int main ()
{
  std::ofstream logfile("dof_accessor_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();
  test<3> ();

  return 0;
}
