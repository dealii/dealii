//----------------------------  cylinder_shell_02.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  cylinder_shell_02.cc  ---------------------------

// turns out that cylinder_shell_01 wasn't enough: measure just takes the
// positive value it computes. we need to ask a FEValues object for the JxW
// values to get the sign of the Jacobian


#include "../tests.h"
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fstream>

    

int main () 
{
  std::ofstream logfile("cylinder_shell_02.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  logfile.precision (2);

                                   // generate a hyperball in 3d
  Triangulation<3> tria;
  GridGenerator::cylinder_shell (tria, 1, .8, 1);

  FE_Q<3> fe(1);
  DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  QMidpoint<3> q;
  FEValues<3> fe_values (fe, q, update_JxW_values);
  
                                   // make sure that all cells have positive
                                   // volume
  for (DoFHandler<3>::active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    {
      fe_values.reinit (cell);
      deallog << cell << ' ' << fe_values.JxW(0) << std::endl;

      Assert (fe_values.JxW(0) > 0, ExcInternalError());
    }
}
