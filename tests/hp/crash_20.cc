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



// if the mesh is generated after the hp::DoFHandler is attached to the
// triangulation object, then we can't set active fe indices -- which is
// somewhat tragic since we have to assign active fe indices before we can
// call distribute_dofs
//
// originally, this problem was avoided because the hp::DoFHandler listens to
// the refinement listener signal to rebuild its data structures; so if you
// create a triangulation object, attach the hp::DoFHandler, create a coarse
// mesh, then refine the mesh, everything is ok again. the solution is to also
// listen to the creation of triangulations.


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_dgq.h>

#include <fstream>



template <int dim>
void test ()
{
  Triangulation<dim> tria;
  hp::DoFHandler<dim> dof_handler(tria);

  GridGenerator::hyper_cube(tria);

  for (typename hp::DoFHandler<dim>::active_cell_iterator
       cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    cell->set_active_fe_index (0);
}


int main ()
{
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
