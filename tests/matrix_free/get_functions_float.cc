// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
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
// operations in getting the function values, the function gradients, and the
// function Laplacians on a mesh with all different cell types (same as in
// matrix_vector_06) with hanging nodes, boundary conditions for float values.

#include <deal.II/base/function.h>

#include "../tests.h"

#include "create_mesh.h"
#include "get_functions_common.h"


template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  create_mesh(tria);
  tria.refine_global(4 - dim);

  // refine a few cells
  for (unsigned int i = 0; i < 10 - 3 * dim; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator cell =
                                                          tria.begin_active(),
                                                        endc = tria.end();
      unsigned int counter                                   = 0;
      for (; cell != endc; ++cell, ++counter)
        if (counter % (7 - i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           1,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  do_test<dim, fe_degree, float>(dof, constraints);
}
