// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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
