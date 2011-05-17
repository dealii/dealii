
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
#include <deal.II/base/logstream.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/data_out.h>

std::ofstream logfile("fe_values_extractor_01/output");


int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

   const unsigned int dim = 1;

   Triangulation<dim,dim+1>    triangulation;
   FESystem <dim,dim+1> fe (FE_Q<dim,dim+1> (1), dim+1);
   DoFHandler<dim,dim+1>       dof_handler(triangulation);

   Vector<double> soln;

   GridIn<dim,dim+1> grid_in;
   grid_in.attach_triangulation (triangulation);
   std::ifstream fname("data_out_02/square.msh");
   grid_in.read_msh (fname);

   dof_handler.distribute_dofs (fe);
   soln.reinit (dof_handler.n_dofs());
   soln = 1;

   std::vector<Tensor<1,dim+1> > local_velocity_values (1);
   const FEValuesExtractors::Vector velocities (0);
   QGauss<dim>  quadrature_formula(1);
   DoFHandler<dim,dim+1>::active_cell_iterator
     cell     = dof_handler.begin_active(),
     endc     = dof_handler.end();
   FEValues<dim,dim+1> fe_v (fe, quadrature_formula, update_values);

   for (; cell!=endc; ++cell) {
     fe_v.reinit (cell);

     fe_v[velocities].get_function_values (soln, local_velocity_values);

     deallog << local_velocity_values.at(0) << std::endl;
   }
}
