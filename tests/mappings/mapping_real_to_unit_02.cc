// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/*

  FEFieldFunction could not deal with ExcTransformationFailed from
  transform_real_to_unit_cell()
  */

#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/fe_field_function.h>

#include "../tests.h"



template <int dim>
void
test2()
{
  const SphericalManifold<dim> boundary_description;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_ball(triangulation);
  triangulation.set_manifold(0, boundary_description);
  triangulation.refine_global(1);
  MappingQ<dim> mapping(1);


  Point<dim> p(-0.27999999999999992, -0.62999999999999989);

  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double>                solution(dof_handler.n_dofs());
  Functions::FEFieldFunction<2> fe_function(dof_handler, solution);
  fe_function.value(p); // this works <<<<<<<<<<<

  std::vector<Point<dim>> points(19 * 19);
  std::vector<double>     m(19 * 19);

  if (1) // works if changed to "if (0)"   <<<<<<<<<
    for (unsigned int i = 0; i < 19; ++i)
      for (unsigned int j = 0; j < 19; ++j)
        {
          /// all points are inside
          points[19 * i + j][0] = -0.7 + (i + 1) * .07;
          points[19 * i + j][1] = -0.7 + (j + 1) * .07;
        }
  points[95] = p;
  fe_function.value_list(points, m); // <<<< this fails at point[95] but only if
                                     // the other points are filled in?!

  triangulation.reset_manifold(0);

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test2<2>();
}
