// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

// test Volume of a Ball

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <sstream>

#include "../tests.h"


void
test(const double R)
{
  const unsigned int dim                          = 3;
  const unsigned int global_mesh_refinement_steps = 4;
  const unsigned int fe_degree                    = 2;
  const unsigned int n_q_points_1d                = 3;

  // derived
  Point<dim> center;
  for (unsigned int d = 0; d < dim; d++)
    center[d] = d;

  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler(triangulation);
  FE_Q<dim>          fe(fe_degree);
  QGauss<dim>        quadrature_formula(n_q_points_1d);

  GridGenerator::hyper_ball(triangulation, center, R);
  triangulation.set_all_manifold_ids_on_boundary(0);
  static SphericalManifold<dim> surface_description(center);
  triangulation.set_manifold(0, surface_description);
  triangulation.refine_global(global_mesh_refinement_steps);
  dof_handler.distribute_dofs(fe);
  MappingQ<dim> mapping(fe_degree);

  FEValues<dim> fe_values(mapping, fe, quadrature_formula, update_JxW_values);

  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                        endc = dof_handler.end();
  const unsigned int n_q_points              = quadrature_formula.size();

  double volume = 0.;
  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      for (unsigned int q = 0; q < n_q_points; ++q)
        volume += fe_values.JxW(q);
    }

  deallog << "Volume:       " << volume << std::endl
          << "Exact volume: " << 4.0 * numbers::PI * std::pow(R, 3.0) / 3.
          << std::endl;

  dof_handler.clear();
}



int
main(int argc, char *argv[])
{
  initlog();
  test(15);
  return 0;
}
