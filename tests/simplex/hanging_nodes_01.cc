// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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



// verify hanging node constraints on locally h-refined simplex mesh
//
// dofs will be enumerated as follows
//  scenario 1:
//   1-------0
//   |\      |
//   |  \    |
//   5---6   |
//   |\  |\  |
//   |  \|  \|
//   3---4---2


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

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

    test<dim>({1, 0},
              {0, 0},
              hp::FECollection<dim>(FE_SimplexP<dim>(1)),
              subdivided_hyper_cube_with_simplices);
    test<dim>({1, 0},
              {0, 0},
              hp::FECollection<dim>(FE_SimplexP<dim>(2)),
              subdivided_hyper_cube_with_simplices);
  }
  deallog.pop();
}
