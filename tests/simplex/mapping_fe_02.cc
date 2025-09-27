// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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


void
test(const unsigned int mapping_degree)
{
  const int dim = 2;

  Triangulation<dim> tria, tria_temp;
  GridGenerator::hyper_ball(tria_temp, Point<dim>(), 1.0, 2.0);
  GridGenerator::convert_hypercube_to_simplex_mesh(tria_temp, tria);
  for (const auto i : tria_temp.get_manifold_ids())
    if (i != numbers::flat_manifold_id)
      tria.set_manifold(i, tria_temp.get_manifold(i));

  FE_SimplexP<dim> fe(2);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  MappingFE<dim> mapping(FE_SimplexP<dim>{mapping_degree});

  {
    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;

    DataOut<dim> data_out;
    data_out.set_flags(flags);

    data_out.attach_dof_handler(dof_handler);

    Vector<double> solution(dof_handler.n_dofs());
    data_out.add_data_vector(solution, "solution");

    data_out.build_patches(mapping, 3);

#if false
    std::ofstream output("test." + std::to_string(mapping_degree) + ".vtk");
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
  test(3);
}
