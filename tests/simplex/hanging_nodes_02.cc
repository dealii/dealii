// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// verify hanging node constraints on locally p-refined simplex mesh
//
// dofs will be enumerated as follows
//  scenario 1:    scenario 2:
//   6-------4      2---4---3
//   |\      |      |\      |
//   |  \    |      |  \    |
//   3   2   |      |   5   6
//   |    \  |      |    \  |
//   |      \|      |      \|
//   0---1---5      0-------1


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"

#include "hanging_nodes.h"


int
main()
{
  initlog();

  deallog.push("2d");
  {
    const unsigned int dim = 2;

    const auto subdivided_hyper_cube_with_simplices =
      [](Triangulation<dim> &tria) {
        GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);
      };

    test<dim>({0, 0},
              {0, 1},
              hp::FECollection<dim>(FE_SimplexP<dim>(2), FE_SimplexP<dim>(1)),
              subdivided_hyper_cube_with_simplices);
    test<dim>({0, 0},
              {0, 1},
              hp::FECollection<dim>(FE_SimplexP<dim>(1), FE_SimplexP<dim>(2)),
              subdivided_hyper_cube_with_simplices);
  }
  deallog.pop();
}
