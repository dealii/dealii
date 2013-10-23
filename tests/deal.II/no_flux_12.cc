// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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



// add test for ExcNotImplemented() error in compute_no_normal_flux_constraints.
// reported by Keith Galvin, mailing list, 2013/10/13. Simplified.
// note that the cylinder boundary is not even the problem here!

/*
5: An error occurred in line <4586> of file </scratch/deal-trunk/deal.II/include/deal.II/numerics/vector_tools.templates.h> in function
5:     void dealii::VectorTools::compute_no_normal_flux_constraints(const DH<dim, spacedim>&, unsigned int, const std::set<unsigned char>&, dealii::ConstraintMatrix&, const dealii::Mapping<dim, spacedim>&) [with int dim = 3; DH = dealii::DoFHandler; int spacedim = 3]
5: The violated condition was: 
5:     contribution->second.size() == dim-1
5: The name and call sequence of the exception was:
5:     ExcNotImplemented()
  
 */

#include "../tests.h"

#include <grid/tria.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_out.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <fe/fe_system.h>
#include <fe/mapping_q.h>
#include <numerics/vector_tools.h>
#include <numerics/data_out.h>
#include <base/exceptions.h>
#include <base/function.h>


template <int dim>
void run()
{
  Triangulation<dim> tria;

  // indicator 6 = cylinder
  GridGenerator::hyper_cube_with_cylindrical_hole (tria, 0.25, 0.5, 0.5, 1, true);
  
  /*  std::string filename = "Mesh.eps";
  std::ofstream output (filename.c_str());
  GridOut grid_out;
  grid_out.write_eps (tria, output);
  */
  
  FESystem<dim> fe(FE_Q<dim>(1), dim);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  ConstraintMatrix constraints;
  std::set<unsigned char> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert (0); // x=0
  no_normal_flux_boundaries.insert (5); // z=1
  
  VectorTools::compute_no_normal_flux_constraints
    (dof_handler, 0,
     no_normal_flux_boundaries,
     constraints);

  constraints.print(deallog.get_file_stream());

  deallog.get_file_stream() << std::flush;
  constraints.close();

  deallog << "OK" << std::endl;
}


int main ()
{
  std::ofstream logfile ("output");
  logfile.precision (7);
  logfile.setf(std::ios::fixed);
  deallog.attach(logfile);
  deallog.depth_console (0);

  run<3> ();
}
