//----------------------------  function_manifold_chart ---------------------------
//    Copyright (C) 2011, 2013, 2014 by the mathLab team.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------- function_manifold_chart ---------------------------


// Test the identity Manifold.

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

// Helper function
template <int dim, int spacedim>
void test(unsigned int ref=1)
{
  deallog << "Testing dim " << dim 
	  << ", spacedim " << spacedim << std::endl;
  
  // Here the only allowed axis is z. In cylinder the default is x.
  std::string push_forward_expression;
  std::string pull_back_expression;

  switch(spacedim) {
  case 2:
    push_forward_expression = "x^2; y^2";
    pull_back_expression = "sqrt(x); sqrt(y)";
    break;
  case 3:
    push_forward_expression = "x^2; y^2; z^2";
    pull_back_expression = "sqrt(x); sqrt(y); sqrt(z)";
    break;
  default:
    Assert(false, ExcInternalError());
  }
  
  FunctionManifold<dim,spacedim,spacedim> manifold(push_forward_expression, 
							pull_back_expression);

  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_cube (tria, 0, 1);

  for(typename Triangulation<dim,spacedim>::active_cell_iterator cell = tria.begin_active(); cell != tria.end(); ++cell) {
    cell->set_all_manifold_ids(1);
  }
  
  tria.set_manifold(1, manifold);
  tria.refine_global(2);
  
  GridOut gridout;
  gridout.write_msh(tria, deallog.get_file_stream());
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  
  test<2,2>();
  test<3,3>();
  
  return 0;
}

