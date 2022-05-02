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

// This function tests the correctness of the matrix-free implementation
// of the FE_RaviartThomasNodal element by evaluating values + gradients
// as well as the divergence on faces and comparing the result with
// FEFaceVaules which is considered the reference. The mesh is a hypercube
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
