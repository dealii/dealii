// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test distributed solution_transfer with averaging (without averaging the
// program fails)

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/solution_transfer.h>

#include "../tests.h"


int
main(int argc, char **argv)
{
  constexpr unsigned int dim    = 2;
  constexpr unsigned int degree = 1;
  using Number                  = double;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  MPILogInitAll                    all;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            {2, 1},
                                            {0., 0.0},
                                            {2.0, 1.0});
  if (tria.create_cell_iterator(CellId("0_0:"))->is_locally_owned())
    tria.create_cell_iterator(CellId("0_0:"))->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_Q<dim>(degree));

  LinearAlgebra::distributed::Vector<Number> v;
  v.reinit(dof_handler.locally_owned_dofs(),
           DoFTools::extract_locally_relevant_dofs(dof_handler),
           MPI_COMM_WORLD);

  std::map<CellId, Vector<double>> values;
  values[CellId("1_0:")]  = Vector<double>{0.5, 0.5, 0.5, 0.5};
  values[CellId("0_1:0")] = Vector<double>{0.5, 0.5, 0.5, 0.5};
  values[CellId("0_1:1")] = Vector<double>{0.5, 0.5, 0.5, 0.1};
  values[CellId("0_1:2")] = Vector<double>{0.5, 0.5, 0.5, 0.5};
  values[CellId("0_1:3")] = Vector<double>{0.5, 0.1, 0.5, 0.5};

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_artificial() == false)
      cell->set_dof_values(values[cell->id()], v);

  v.print(deallog.get_file_stream());

  SolutionTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
    solution_trans(dof_handler, true /*enabling averaging*/);

  v.update_ghost_values();
  solution_trans.prepare_for_coarsening_and_refinement(v);

  if (tria.create_cell_iterator(CellId("1_0:"))->is_locally_owned())
    tria.create_cell_iterator(CellId("1_0:"))->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(FE_Q<dim>(degree));
  v.reinit(dof_handler.locally_owned_dofs(),
           DoFTools::extract_locally_relevant_dofs(dof_handler),
           MPI_COMM_WORLD);
  solution_trans.interpolate(v);

  v.print(deallog.get_file_stream());
}
