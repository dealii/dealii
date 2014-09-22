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
#include <deal.II/base/logstream.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>

std::ofstream logfile("output");


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
  std::ifstream fname(SOURCE_DIR "/grids/square.msh");
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
