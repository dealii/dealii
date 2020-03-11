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
// with hanging nodes

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
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(2);
  else
    tria.refine_global(1);
  typename Triangulation<dim>::active_cell_iterator cell, endc;
  cell = tria.begin_active(), endc = tria.end();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 0.2)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active();
  for (unsigned int i = 0; i < 10 - 3 * dim; ++i)
    {
      cell                 = tria.begin_active();
      unsigned int counter = 0;
      for (; cell != endc; ++cell, ++counter)
        if (counter % (7 - i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  // test with threads enabled as well
  do_test<dim, fe_degree, fe_degree + 1, double>(dof, constraints, true);

  // and finally a float test
  do_test<dim, fe_degree, fe_degree + 1, float>(dof, constraints);
}
