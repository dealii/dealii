
//----------------------------  template.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2008 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// output the vertex numbering in a vtk file

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
#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/base/function.h>

#include <fstream>
#include <string>

std::ofstream logfile("data_out/output");


template <int dim, int spacedim>
void test(std::string filename) {

  Triangulation<dim, spacedim> triangulation;
  GridIn<dim, spacedim> gi;

  gi.attach_triangulation (triangulation);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);

  FE_Q<dim,spacedim>     fe (1);
  DoFHandler<dim,spacedim> dof_handler (triangulation);

  dof_handler.distribute_dofs (fe);
  
  // Output the vertex numbering
  Vector<double> numbering(dof_handler.n_dofs());
  for(unsigned int i=0; i <numbering.size(); ++i) 
      numbering(i) = i;
  
  DataOut<dim, DoFHandler<dim,spacedim> > dataout;
  dataout.attach_dof_handler(dof_handler);
  dataout.add_data_vector(numbering, "numbering");
  dataout.build_patches();
  dataout.write_vtk(logfile);
}



int main () 
{
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog << "Test<1,2>" << std::endl;
  test<1,2>("grids/circle_2.inp");

  deallog << "Test<2,3>" << std::endl;
  test<2,3>("grids/sphere_2.inp");

  return 0;
}
