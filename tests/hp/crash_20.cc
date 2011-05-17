//----------------------------  crash_20.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  crash_20.cc  ---------------------------


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
  std::ofstream logfile("crash_20/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
