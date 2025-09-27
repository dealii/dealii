// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// investigate how the order in which we process hanging nodes on faces
// affect hp-discretizations in 3D
//
// Case 0:       Case 1:
//  +---+         +-+-+
//  |   |         |2|2|
//  | 2 |         +-+-+
//  |   |         |2|2|
//  +---+-+-+     +-+-+---+
//  |   |2|2|     |   |   |
//  | 3 +-+-+     | 3 | 2 |
//  |   |2|2|     |   |   |
//  +---+---+     +---+-+-+


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


void
test(const int skip_cell)
{
  constexpr int dim = 3;

  // set up grid
  Triangulation<dim> triangulation;
  {
    const std::vector<unsigned int> repetitions({2, 2, 1});
    const Point<dim>                bottom_left(-1, -1, 0);
    const Point<dim>                top_right(1, 1, 1);
    const std::vector<int>          n_cells_to_remove({-1, -1, 0});

    GridGenerator::subdivided_hyper_L(
      triangulation, repetitions, bottom_left, top_right, n_cells_to_remove);
  }

  // set up hp-refinement
  DoFHandler<dim> dof_handler(triangulation);
  {
    auto cell = dof_handler.begin_active();
    cell->set_active_fe_index(1);
    if (skip_cell)
      ++cell;
    (++cell)->set_refine_flag();

    triangulation.execute_coarsening_and_refinement();
  }

  // enumerate dofs
  hp::FECollection<dim> fe_collection;
  {
    for (unsigned int degree = 2; degree <= 3; ++degree)
      fe_collection.push_back(FE_Q<dim>(degree));
    dof_handler.distribute_dofs(fe_collection);
  }

  // make hanging node constraints
  AffineConstraints<double> constraints;
  {
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    constraints.close();
  }

#if false
  // output vtu
  DataOut<dim> data_out;
  {
    Vector<float> fe_degrees(triangulation.n_active_cells());
    for (const auto &cell : dof_handler.active_cell_iterators())
      fe_degrees(cell->active_cell_index()) = cell->get_fe().degree;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(fe_degrees, "fe_degrees");
    data_out.build_patches();

    std::ofstream output("output-" + Utilities::to_string(skip_cell) + ".vtu");
    data_out.write_vtu(output);
  }
#endif

  constraints.print(deallog.get_file_stream());
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  deallog << "Case 0:" << std::endl;
  test(0);
  deallog << "Case 1:" << std::endl;
  test(1);
}
