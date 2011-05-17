//----------------------------  template.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// interpolation of constant function on the surface of a hyperxosphere

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

std::ofstream logfile("interpolation/output");


template <int dim, int spacedim>
void test(std::string filename, unsigned int n) {

  Triangulation<dim, spacedim> triangulation;
  GridIn<dim, spacedim> gi;

  gi.attach_triangulation (triangulation);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);

  FE_Q<dim,spacedim>     fe (n);
  DoFHandler<dim,spacedim> dof_handler (triangulation);

  dof_handler.distribute_dofs (fe);
  
  // Now we interpolate the constant function on the mesh, and check
  // that this is consistent with what we expect.
  Vector<double> real_one(dof_handler.n_dofs());
  Vector<double> interpolated_one(dof_handler.n_dofs());
  for(unsigned int i=0; i <real_one.size(); ++i) 
      real_one(i) = 1.;
  
  ConstantFunction<spacedim> constant(1.);
  VectorTools::interpolate(dof_handler, constant, interpolated_one);
  
  real_one.add(-1, interpolated_one);
  
  deallog << "L2 Norm of difference: " 
	  << real_one.l2_norm() << std::endl
	  <<  "Linfty Norm of difference: " 
	  << real_one.linfty_norm() << std::endl;
}



int main () 
{
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  for(unsigned int n=1; n<8; ++n) 
    {
      deallog << "Test<1,2>, finite element q_" << n << std::endl;
      test<1,2>("grids/circle_2.inp",n);
      
      deallog << "Test<1,2>, finite element q_" << n << std::endl;
      test<2,3>("grids/sphere_2.inp",n);
    }
  return 0;
}
