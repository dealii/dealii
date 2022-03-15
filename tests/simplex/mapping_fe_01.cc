// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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


// Distribute FE_WedgeP on a DoFHandler.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

using namespace dealii;

void
test(const unsigned int mapping_degree)
{
  const int dim = 2;

  Triangulation<dim> tria, tria_temp;
  GridGenerator::hyper_shell(tria_temp, Point<dim>(), 1.0, 2.0);
  GridGenerator::convert_hypercube_to_simplex_mesh(tria_temp, tria);
  for (const auto i : tria_temp.get_manifold_ids())
    if (i != numbers::flat_manifold_id)
      tria.set_manifold(i, tria_temp.get_manifold(i));

  FE_SimplexP<dim> fe(2);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  MappingFE<dim> mapping(FE_SimplexP<dim>{mapping_degree});

  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);

    Vector<double> solution(dof_handler.n_dofs());
    data_out.add_data_vector(solution, "solution");

    data_out.build_patches(mapping, 2);

#if false
    std::ofstream output("test." + std::to_string(mapping_degree) +  ".vtk");
    data_out.write_vtk(output);
#else
    data_out.write_vtk(deallog.get_file_stream());
#endif
  }
}

int
main()
{
  initlog();

  test(1);
  test(2);
}
