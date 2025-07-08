// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this test failed at some stage:
// DoFTools::make_hanging_node_constraints() has not worked with an
// hp::DoFHandler on ghost cells that neighbor artificial cells.
//
// Geometry in 2D divided on two subdomains:
// Subdomain0:          Subdomain1:
// +-------+-------+    +-------+-------+
// |       |       |    |       |       |
// |   g   |   g   |    |   o   |   o   |
// |       |       |    |       |       |
// +---+---+---+!!!+    +---+---+---+---+
// | o | o | g | a |    | g | g | o | o |
// +---+---+---+---+    +---+---+---+---+
// | o | o | g | a |    | a | g | o | o |
// +---+---+---+---+    +---+---+---+---+


#include <deal.II/base/index_set.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"


template <int dim>
void
test()
{
  // setup triangulation
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  // refine half of all cells
  unsigned int i = 0;
  for (const auto &cell : tria.active_cell_iterators())
    if (i++ < .5 * tria.n_active_cells())
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  // setup finite elements
  hp::FECollection<dim> fes;
  fes.push_back(FE_Q<dim>(1));

  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(fes);

  // make constraints
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dh);

  AffineConstraints<double> constraints;
  constraints.reinit(dh.locally_owned_dofs(), locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dh, constraints);

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  MPI_Barrier(MPI_COMM_WORLD);
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
