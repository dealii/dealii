// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2014 by the deal.II authors
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


// check the creation of no-flux boundary conditions for a finite
// element that consists of more than dim components and where
// therefore we have to pick the vector components from somewhere in
// the middle


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
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>



template<int dim>
void test (const Triangulation<dim> &tr,
           const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
    {
      deallog << "FE=" << fe.get_name()
              << ", case=" << i
              << std::endl;

      std::set<types::boundary_id> boundary_ids;
      for (unsigned int j=0; j<=i; ++j)
        boundary_ids.insert (j);

      ConstraintMatrix cm;
      VectorTools::compute_no_normal_flux_constraints (dof, 1, boundary_ids, cm);

      cm.print (deallog.get_file_stream ());
    }
}


template<int dim>
void test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
    tr.begin_active()->face(i)->set_boundary_indicator (i);

  tr.refine_global(2);

  for (unsigned int degree=1; degree<4; ++degree)
    {
      FESystem<dim> fe (FE_Q<dim>(degree+1), 1,
                        FE_Q<dim>(degree), dim,
                        FE_Q<dim>(degree+1), 1);
      test(tr, fe);
    }
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-12);

  test_hyper_cube<2>();
  test_hyper_cube<3>();
}
