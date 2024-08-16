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


// Test AffineConstraints::make_consistent_in_parallel().

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"


/**
 * mpirun -np 4 ./constraints_01
 *
 * 6-----7-----8
 * | 0/2 | 1/3 |
 * 2-----3-----5
 * | 0/0 | 0/1 |  ... with dof-index, material-id/rank
 * 0-----1-----4
 *
 * Expected output: {3}, {3, 5}, {3, 7}, {3, 5, 7}
 * Actual output:   {},  {3, 5}, {3, 7}, {3, 5, 7}
 *
 * ... code to reproduce https://github.com/dealii/dealii/issues/11725
 */
template <typename Number>
IndexSet
collect_lines(const AffineConstraints<Number> &constraints,
              const unsigned int               size)
{
  IndexSet lines_local(size);
  for (const auto &line : constraints.get_lines())
    lines_local.add_index(line.index);
  return lines_local;
}

template <int dim, int spacedim>
void
test(const DoFHandler<dim, spacedim> &dof_handler,
     const IndexSet                  &locally_relevant_dofs)
{
  AffineConstraints<double> constraints(dof_handler.locally_owned_dofs(),
                                        locally_relevant_dofs);

  std::vector<types::global_dof_index> dof_indices(
    dof_handler.get_fe().n_dofs_per_face());

  for (const auto &cell : dof_handler.active_cell_iterators())
    for (const auto face : cell->face_indices())
      if (cell->is_locally_owned() && !cell->at_boundary(face) &&
          cell->material_id() != cell->neighbor(face)->material_id())
        {
          cell->face(face)->get_dof_indices(dof_indices);
          for (const auto i : dof_indices)
            if (constraints.is_constrained(i) == false)
              constraints.constrain_dof_to_zero(i);
        }

  const auto a = collect_lines(constraints, dof_handler.n_dofs());
  a.print(deallog.get_file_stream());

  constraints.make_consistent_in_parallel(dof_handler.locally_owned_dofs(),
                                          locally_relevant_dofs,
                                          dof_handler.get_mpi_communicator());

  const auto b = collect_lines(constraints, dof_handler.n_dofs());
  b.print(deallog.get_file_stream());
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  const int dim = 2;

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(tria, 2);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_Q<dim>(1));

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->center()[0] > 0.5 && cell->center()[1] > 0.5)
      cell->set_material_id(1);

  const IndexSet locally_active_dofs =
    DoFTools::extract_locally_active_dofs(dof_handler);
  test(dof_handler, locally_active_dofs);

  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);
  test(dof_handler, locally_relevant_dofs);
}
