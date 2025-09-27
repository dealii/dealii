// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// distribute dofs on a shared and distributed triangulation and compare
// Note: this doesn't pass for all meshes since, for some complicated refinement
//       schemes, cells may be traversed in different order in a distributed
//       triangulation as opposed to a shared triangulation.

#include <deal.II/base/index_set.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

template <int dim>
void
compare_meshes(DoFHandler<dim> &shared_dof_handler,
               DoFHandler<dim> &distributed_dof_handler)
{
  FE_Q<dim> fe(2);

  shared_dof_handler.distribute_dofs(fe);
  distributed_dof_handler.distribute_dofs(fe);

  IndexSet shared_dofs      = shared_dof_handler.locally_owned_dofs();
  IndexSet distributed_dofs = distributed_dof_handler.locally_owned_dofs();
  Assert(shared_dofs == distributed_dofs, ExcInternalError());
  shared_dofs.print(deallog.get_file_stream());

  std::vector<IndexSet> shared_dofs_per_proc =
    Utilities::MPI::all_gather(MPI_COMM_WORLD,
                               shared_dof_handler.locally_owned_dofs());
  std::vector<IndexSet> distributed_dofs_per_proc =
    Utilities::MPI::all_gather(MPI_COMM_WORLD,
                               distributed_dof_handler.locally_owned_dofs());
  for (unsigned int i = 0; i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
       ++i)
    Assert(shared_dofs_per_proc[i] == distributed_dofs_per_proc[i],
           ExcInternalError());

  typename DoFHandler<dim>::active_cell_iterator
    cell = distributed_dof_handler.begin_active(),
    endc = distributed_dof_handler.end();
  for (; cell != endc; ++cell)
    {
      if (cell->subdomain_id() == numbers::artificial_subdomain_id)
        continue;

      typename DoFHandler<dim>::active_cell_iterator dof_shared_cell(
        &shared_dof_handler.get_triangulation(),
        cell->level(),
        cell->index(),
        &shared_dof_handler);

      std::vector<types::global_dof_index> distributed_cell_dofs(
        fe.dofs_per_cell);
      std::vector<types::global_dof_index> shared_cell_dofs(fe.dofs_per_cell);
      cell->get_dof_indices(distributed_cell_dofs);
      dof_shared_cell->get_dof_indices(shared_cell_dofs);
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        Assert(distributed_cell_dofs[i] == shared_cell_dofs[i],
               ExcInternalError());
    }
}



template <int dim>
void
test()
{
  parallel::shared::Triangulation<dim> shared_tria(
    MPI_COMM_WORLD,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::limit_level_difference_at_vertices),
    true,
    typename parallel::shared::Triangulation<dim>::Settings(
      parallel::shared::Triangulation<dim>::partition_zorder));
  DoFHandler<dim> shared_dof_handler(shared_tria);

  parallel::distributed::Triangulation<dim> distributed_tria(
    MPI_COMM_WORLD, Triangulation<dim>::limit_level_difference_at_vertices);
  DoFHandler<dim> distributed_dof_handler(distributed_tria);

  GridGenerator::subdivided_hyper_cube(shared_tria, 3, -1, 1);
  GridGenerator::subdivided_hyper_cube(distributed_tria, 3, -1, 1);

  compare_meshes(shared_dof_handler, distributed_dof_handler);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();

  deallog << "OK" << std::endl;
}
