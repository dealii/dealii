// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// FEFieldFunction needed changes to work on simplex meshes. Test by
// Sam Scheuerman.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/fe_field_function.h>

#include "../tests.h"


int
main()
{
  initlog();

  const unsigned int dim = 3;

  // Do all the usual setup
  Triangulation<dim>  triangulation;
  DoFHandler<dim>     dof_handler(triangulation);
  const FESystem<dim> fe(FE_SimplexP<dim>(1), 3);

  Vector<double> solution;

  GridGenerator::subdivided_hyper_cube_with_simplices(triangulation, 1);
  dof_handler.distribute_dofs(fe);
  solution.reinit(dof_handler.n_dofs());

  for (unsigned int i = 0; i < solution.size(); ++i)
    solution[i] = static_cast<double>(i);

  // Data storage for the interpolated function
  std::vector<Point<dim>>     input_points;
  std::vector<Vector<double>> output_displacements;

  double top  = 1.;
  double step = top / 20;

  // Fill input points with 21 regular samples along the line
  // from (.5, .5, 0) to (.5, .5, 1)
  for (unsigned int i = 0; i < 21; ++i)
    {
      const double z_value = step * static_cast<double>(i);
      input_points.emplace_back(0.5, 0.5, z_value);
    }

  output_displacements.resize(input_points.size());

  // Create the Field function object
  const MappingFE<dim>            mapping(FE_SimplexP<dim>(1));
  Functions::FEFieldFunction<dim> displacement_function(dof_handler,
                                                        solution,
                                                        mapping);

  // Attempt to get the vector output of the displacement
  // function at each of the input points. This is where the
  // code fails.
  displacement_function.vector_value_list(input_points, output_displacements);

  for (const auto &v : output_displacements)
    deallog << v << ' ';
  deallog << std::endl;
}
