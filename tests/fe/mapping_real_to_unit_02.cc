// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


/*

  FEFieldFunction could not deal with ExcTransformationFailed from transform_real_to_unit_cell()
  */

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/fe_field_function.h>





template<int dim>
void test2()
{
  const HyperBallBoundary<dim> boundary_description;

  Triangulation<dim>   triangulation;
  GridGenerator::hyper_ball (triangulation);
  triangulation.set_boundary (0, boundary_description);
  triangulation.refine_global (1);
  MappingQ1< dim > mapping;


  Point<dim> p(-0.27999999999999992, -0.62999999999999989);

  FE_Q<dim> fe(2);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> solution(dof_handler.n_dofs());
  Functions::FEFieldFunction<2> fe_function (dof_handler, solution);
  fe_function.value (p); //this works <<<<<<<<<<<

  std::vector<Point<dim> > points(19*19);
  std::vector<double> m (19*19);

  if (1) //works if changed to "if (0)"   <<<<<<<<<
    for (unsigned int i = 0; i < 19; i++)
      for (unsigned int j = 0; j < 19; j++)
        {
          /// all points are inside
          points[19*i+j] (0) = -0.7 + (i + 1) * .07;
          points[19*i+j] (1) = -0.7 + (j + 1) * .07;
        }
  points[95]=p;
  fe_function.value_list (points, m); // <<<< this fails at point[95] but only if the other points are filled in?!

  triangulation.set_boundary (0);

  deallog << "OK" << std::endl;
}


int
main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test2<2>();

  return 0;
}



