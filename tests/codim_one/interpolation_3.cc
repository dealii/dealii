//----------------------------  template.cc  ---------------------------
//    $Id: testsuite.html 13373 2006-07-13 13:12:08Z manigrasso $
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2008, 2009 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// calculates the surface of a sphere

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>

// all include files needed for the program

#include <grid/tria.h>
#include <grid/grid_in.h>
#include <grid/grid_out.h>
#include <fe/mapping.h>
#include <fe/mapping_q1.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <numerics/vectors.h>
#include <numerics/data_out.h>
#include <base/function.h>
#include <base/function_parser.h>

#include <fstream>
#include <string>


std::ofstream logfile("interpolation_3/output");

// Test interpolation on system of finite elements. 

template <int dim, int spacedim>
void test(std::string filename) {

  Triangulation<dim, spacedim> triangulation;
  GridIn<dim, spacedim> gi;

  gi.attach_triangulation (triangulation);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);

  FE_Q<dim,spacedim>     fe_base (1);
  FESystem<dim, spacedim> fe(fe_base, spacedim);
  DoFHandler<dim,spacedim> dof_handler (triangulation);

  dof_handler.distribute_dofs (fe);
  
  // Now we interpolate the constant function on the mesh, and check
  // that this is consistent with what we expect.
  Vector<double> interpolated_one(dof_handler.n_dofs());
  
  FunctionParser<spacedim> func(spacedim);
  std::map<std::string, double> maps;
  if(spacedim == 2)
      func.initialize("x,y", "x^2; y^2", maps);
  else
      func.initialize("x,y,z", "x^2; y^2; z^2", maps);
  
  VectorTools::interpolate(dof_handler, func, interpolated_one);
  
  DataOut<dim, DoFHandler<dim,spacedim> > dataout;
  dataout.attach_dof_handler(dof_handler);
  dataout.add_data_vector(interpolated_one, "test");
  dataout.build_patches();
  dataout.write_vtk(logfile);
}



int main () 
{
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  deallog << "Test<1,2>" << std::endl;
  test<1,2>("grids/circle_2.inp");

  deallog << "Test<1,2>" << std::endl;
  test<2,3>("grids/sphere_2.inp");

  return 0;
}
