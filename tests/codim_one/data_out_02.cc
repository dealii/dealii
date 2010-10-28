
//----------------------------  template.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors and Sebastian Pauletti
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// DataOut::build_patches appeared to have a problem with outputting
// lines in 2d where nodes were numbered differently when writing data
// vectors as opposed to writing node locations. in the end this
// turned out to be a feature: the mesh was a circle of lines, so
// there are equally many cells as their were nodes, and consequently
// DataOut assumed that it had cell_data, rather than
// dof_data. passing the correct argument fixed the problem, but it
// won't hurt to have this test anyway.

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>

#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <grid/grid_in.h>
#include <lac/vector.h>
#include <numerics/data_out.h>

std::ofstream logfile("data_out_02/output");


int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  const unsigned int dim = 1;

  Triangulation<dim,dim+1>    triangulation;
  FE_Q<dim,dim+1>             fe(1);
  DoFHandler<dim,dim+1>       dof_handler(triangulation);
  Vector<double> soln;

  GridIn<dim,dim+1> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream fname("data_out_02/square.msh");
  grid_in.read_msh (fname);

  dof_handler.distribute_dofs (fe);
  soln.reinit (dof_handler.n_dofs());
  soln = 0;
  for (unsigned int i=0; i<soln.size(); ++i)
    soln(i) = i;
  DataOut<dim, DoFHandler<dim, dim+1> > data_out;
  data_out.attach_dof_handler (dof_handler);

  data_out.add_data_vector (soln, "scalar_data",
			    DataOut<dim,DoFHandler<dim, dim+1> >::type_dof_data);
  data_out.build_patches ();
  data_out.write_vtk (deallog.get_file_stream());

  return 0;
}
