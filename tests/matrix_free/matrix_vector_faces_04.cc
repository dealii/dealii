// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// tests matrix-free face evaluation, matrix-vector products as compared to
// the same implementation with MeshWorker. This example uses a hyperball mesh
// with both hanging nodes and curved boundaries

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_dgq.h>

#include "../tests.h"

#include "matrix_vector_faces_common.h"

template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);
  tria.refine_global(4 - dim);
  tria.begin_active()->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  // test with threads enabled as well
  do_test<dim, fe_degree, fe_degree + 1, double>(dof, constraints, true);
}
