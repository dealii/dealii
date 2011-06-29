//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2008, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// no normal flux constraints on a hyper cube for all faces this caused
// ExcMessage (\"Cycle in constraints detected!\")" in 3d with a higher order
// mapping.  to make things even weirder, mappings of order <4 work.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vectors.h>

#include <fstream>



template<int dim>
void test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  
  FESystem<dim> fe (FE_Q<dim>(2), dim);

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  deallog << "FE=" << fe.get_name()
	  << std::endl;

  std::set<unsigned char> boundary_ids;
  boundary_ids.insert (0);

  ConstraintMatrix cm;
  const MappingQ<dim> mapping(4);
  VectorTools::compute_no_normal_flux_constraints (dof, 0,
						   boundary_ids, cm,
						   mapping);
  cm.close();
      

  cm.print (deallog.get_file_stream ());
}


int main()
{
  std::ofstream logfile ("no_flux_06/output");
  deallog << std::setprecision (2);
  deallog << std::fixed;  
  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-12);

  test_hyper_cube<2>();
  test_hyper_cube<3>();
}
