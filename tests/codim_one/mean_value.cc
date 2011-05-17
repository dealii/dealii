//----------------------------  template.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2008, 2011 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// continous projection of a function on the surface of a hypersphere

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>

#include <deal.II/base/quadrature_lib.h>

#include <fstream>
#include <string>


std::ofstream logfile("mean_value/output");


template <int dim, int spacedim>
void test()
{
  Triangulation<dim, spacedim> triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);

  FE_Q<dim,spacedim>     fe (1);
  DoFHandler<dim,spacedim> dof_handler (triangulation);

  dof_handler.distribute_dofs (fe);
  Functions::CosineFunction<spacedim> cosine;
  Vector<double> x(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, cosine, x);

  const double mean = VectorTools::compute_mean_value (dof_handler,
						       QGauss<dim>(2),
						       x, 0);
  // we have a symmetric domain and a symmetric function. the result
  // should be close to zero
  Assert (std::fabs(mean) < 1e-15, ExcInternalError());
  deallog << "OK" << std::endl;
}



int main () 
{
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  test<1,2>();
  test<2,3>();
}
