// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This function tests the correctness of the matrix-free implementation
// of the FE_RaviartThomasNodal element by evaluating values + gradients
// as well as the divergence on faces and comparing the result with
// FEFaceValues which is considered the reference. The mesh is a hypercube
// mesh with no hanging nodes and no other constraints.

#include "../tests.h"

#include "matrix_vector_rt_face_common.h"


template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  FE_RaviartThomasNodal<dim> fe(fe_degree - 1);

  DoFHandler<dim> dof(tria);

  dof.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  constraints.close();
  deallog << "Using " << dof.get_fe().get_name() << std::endl;
  deallog << "Number of cells: " << dof.get_triangulation().n_active_cells()
          << std::endl;
  deallog << "Number of degrees of freedom: " << dof.n_dofs() << std::endl
          << std::endl;
  do_test<dim, fe_degree, double>(dof, constraints, TestType::values_gradients);
  do_test<dim, fe_degree, double>(dof, constraints, TestType::divergence);
}
