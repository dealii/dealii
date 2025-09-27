// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this test is similar to matrix_vector_06, but implements the operations for
// FE_QDG0 instead of FE_Q (where there is an additional discontinuous degree
// of freedom and different routines need to be used). The data is still not
// very useful because the matrix does not include face terms actually present
// in an approximation of the Laplacian. It only contains cell terms.

#include <deal.II/fe/fe_q_dg0.h>

#include "../tests.h"

#include "matrix_vector_common.h"



template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(1);
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->center().norm() < 1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q_DG0<dim>   fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  do_test<dim, fe_degree, double, fe_degree + 1>(dof, constraints);
}
