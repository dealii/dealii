// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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



// tests matrix-free face evaluation, matrix-vector products as compared to
// the same implementation with MeshWorker. This example uses a hypercube mesh
// without any hanging nodes

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_dgq.h>

#include "../tests.h"

std::ofstream logfile("output");

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

#if (DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 1 && defined(__SSE2__)) || \
  (DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 1 && defined(__ALTIVEC__))
  do_test<dim, fe_degree, fe_degree + 1, double, VectorizedArray<double, 2>>(
    dof, constraints, true);
#endif

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 2 && defined(__AVX__)
  do_test<dim, fe_degree, fe_degree + 1, double, VectorizedArray<double, 4>>(
    dof, constraints, true);
#endif

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 3 && defined(__AVX512F__)
  do_test<dim, fe_degree, fe_degree + 1, double, VectorizedArray<double, 8>>(
    dof, constraints, true);
#endif
}
