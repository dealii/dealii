// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
test(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  const std::set<types::boundary_id> boundary_ids = {0};

  AffineConstraints<double> cm;
  VectorTools::compute_no_normal_flux_constraints(dof, 0, boundary_ids, cm);

  cm.print(deallog.get_file_stream());
}


template <int dim>
void
test_hyper_cube()
{
  Triangulation<dim> tr;

  // create a hypercube, then shear
  // its top surface to make the
  // result non-square
  GridGenerator::hyper_cube(tr);
  Point<dim> shift;
  for (unsigned int i = GeometryInfo<dim>::vertices_per_cell / 2;
       i < GeometryInfo<dim>::vertices_per_cell;
       ++i)
    tr.begin_active()->vertex(i)[0] += 0.5;

  FESystem<dim> fe(FE_Q<dim>(2), dim);
  test(tr, fe);
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;

  test_hyper_cube<3>();
}
