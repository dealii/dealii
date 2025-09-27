// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// this function tests the correctness of the implementation of matrix free
// matrix-vector products by comparing with the result of deal.II sparse
// matrix. The mesh uses a hypershell mesh without hanging nodes (only cell
// type: 2)

#include "../tests.h"

#include "matrix_vector_device_common.h"


template <int dim, int fe_degree, typename Number>
void
test()
{
  const SphericalManifold<dim> manifold;
  Triangulation<dim>           tria;
  GridGenerator::hyper_shell(tria, Point<dim>(), 0.5, 1., 96, true);
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, manifold);
  if (dim == 2)
    tria.refine_global(2);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<Number> constraints;
  constraints.close();

  do_test<dim,
          fe_degree,
          Number,
          LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>,
          fe_degree + 1>(dof, constraints, tria.n_active_cells());
}
