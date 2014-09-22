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
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

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
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();
  test<3> ();

  return 0;
}
