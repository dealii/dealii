// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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



// The function VectorTools::compute_no_normal_flux_constraints had a
// bug that led to an exception whenever we were computing constraints
// for vector fields located on edges shared between two faces of a 3d
// cell if those faces were not perpendicular: in a case like that, we
// computed the unconstrained tangent field as the cross product of
// the two face normals, but the resulting tangent did not have unit
// length
//
// we can test this using a 3d cube that we shear and by computing
// constraints for a Q2^3 field for which there are DoFs on edge
// midpoints.

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

  std::set<types::boundary_id> boundary_ids;
  boundary_ids.insert (0);

  ConstraintMatrix cm;
  VectorTools::compute_no_normal_flux_constraints (dof, 0, boundary_ids, cm);

  cm.print (deallog.get_file_stream ());
}


template<int dim>
void test_hyper_cube()
{
  Triangulation<dim> tr;

  // create a hypercube, then shear
  // its top surface to make the
  // result non-square
  GridGenerator::hyper_cube(tr);
  Point<dim> shift;
  for (unsigned int i=GeometryInfo<dim>::vertices_per_cell/2;
       i < GeometryInfo<dim>::vertices_per_cell; ++i)
    tr.begin_active()->vertex(i)[0] += 0.5;

  FESystem<dim> fe (FE_Q<dim>(2), dim);
  test(tr, fe);
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-12);

  test_hyper_cube<3>();
}
