// ---------------------------------------------------------------------
// $Id: no_flux_11.cc 31349 2013-10-20 19:07:06Z maier $
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

// bug report from mailing list from 11/15/2013 (simplified). no_normal_flux throws
// an exception:
/*
--------------------------------------------------------
An error occurred in line <4314> of file </scratch/deal-trunk/deal.II/include/deal.II/numerics/vector_tools.templates.h> in function
    void dealii::VectorTools::compute_no_normal_flux_constraints(const DH<dim, spacedim>&, unsigned int, const std::set<unsigned char>&, dealii::ConstraintMatrix&, const dealii::Mapping<dim, spacedim>&) [with int dim = 3; DH = dealii::DoFHandler; int spacedim = 3]
The violated condition was: 
    vector_dofs.dof_indices[d] < dof_handler.n_dofs()
The name and call sequence of the exception was:
    ExcInternalError()
Additional Information: 
(none)
*/  

#include "../tests.h"

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

template <int dim>
void test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation,-1.0,1.0);
  triangulation.begin_active()->face(1)->set_all_boundary_indicators(1);

  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  
  ConstraintMatrix constraints;
  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert (1);
  VectorTools::compute_no_normal_flux_constraints (dof_handler,
						   0,
						   no_normal_flux_boundaries,
						   constraints);

  constraints.close();
  constraints.print(deallog.get_file_stream ());
}


int main ()
{
  initlog();
  deallog.depth_console (0);

  test<3> ();
}
