// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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



// calculates the surface of a sphere

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

// all include files needed for the program

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>

#include <fstream>
#include <string>


std::ofstream logfile("output");


template <int dim, int spacedim>
void test(std::string filename)
{

  Triangulation<dim, spacedim> triangulation;
  GridIn<dim, spacedim> gi;

  gi.attach_triangulation (triangulation);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);

  FE_Q<dim,spacedim>     fe (1);
  DoFHandler<dim,spacedim> dof_handler (triangulation);

  dof_handler.distribute_dofs (fe);

  // Now we interpolate the constant function on the mesh, and check
  // that this is consistent with what we expect.
  Vector<double> interpolated_one(dof_handler.n_dofs());

  Functions::CosineFunction<spacedim> cosine;
  VectorTools::interpolate(dof_handler, cosine, interpolated_one);

  DataOut<dim, DoFHandler<dim,spacedim> > dataout;
  dataout.add_data_vector(dof_handler, interpolated_one, "numbering");
  dataout.build_patches();
  dataout.write_vtk(logfile);
}



int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog << "Test<1,2>" << std::endl;
  test<1,2>(SOURCE_DIR "/grids/circle_2.inp");

  deallog << "Test<1,2>" << std::endl;
  test<2,3>(SOURCE_DIR "/grids/sphere_2.inp");

  return 0;
}
