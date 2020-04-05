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



// same as matrix_vector_faces_03 but using a Lagrange polynomial based on
// Gauss points rather than Gauss-Lobatto points where the mass matrix is
// diagonal

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

#include "create_mesh.h"
#include "matrix_vector_faces_common.h"


template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  create_mesh(tria);

  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  typename Triangulation<dim>::active_cell_iterator cell, endc;
  cell = tria.begin_active();
  endc = tria.end();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 0.5)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  if (fe_degree == 1)
    tria.refine_global(1);
  cell = tria.begin_active();
  for (unsigned int i = 0; i < 9 - 3 * dim; ++i)
    {
      cell                 = tria.begin_active();
      endc                 = tria.end();
      unsigned int counter = 0;
      for (; cell != endc; ++cell, ++counter)
        if (counter % (7 - i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_DGQArbitraryNodes<dim> fe(QGauss<1>(fe_degree + 1));
  DoFHandler<dim>           dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  // also test with threads
  do_test<dim, fe_degree, fe_degree + 1, double>(dof, constraints, true);
}
