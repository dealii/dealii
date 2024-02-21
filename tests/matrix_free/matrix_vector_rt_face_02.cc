// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This test it the same as matrix_vector_rt_face_01.cc however with
// non-Cartesian (but still affine) cells.

#include "../tests.h"

#include "matrix_vector_rt_face_common.h"


template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  const unsigned int n_subdivisions = 2;
  Point<dim>         corners[dim];
  corners[0] = (dim == 2) ? Point<dim>(1, 0) : Point<dim>(1, 0, 0);
  corners[1] = (dim == 2) ? Point<dim>(0.5, 0.5) : Point<dim>(0.5, 1, 0.25);
  if (dim == 3)
    corners[2] = Point<dim>(0.5, 0, 1);
  GridGenerator::subdivided_parallelepiped(tria, n_subdivisions, corners);

  FE_RaviartThomasNodal<dim> fe(fe_degree - 1);
  DoFHandler<dim>            dof(tria);
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
