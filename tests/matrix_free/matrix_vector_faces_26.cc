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
// the same implementation with MeshWorker. This example uses a hypercube mesh
// without any hanging nodes

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_dgq.h>

#include "../tests.h"

#include "matrix_vector_faces_common.h"

template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(5 - dim);

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  do_test<dim, fe_degree, fe_degree + 1, double, VectorizedArray<double, 1>>(
    dof, constraints, true);

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
  do_test<dim, fe_degree, fe_degree + 1, double, VectorizedArray<double, 2>>(
    dof, constraints, true);
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
  do_test<dim, fe_degree, fe_degree + 1, double, VectorizedArray<double, 4>>(
    dof, constraints, true);
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
  do_test<dim, fe_degree, fe_degree + 1, double, VectorizedArray<double, 8>>(
    dof, constraints, true);
#endif
}
