// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// distribute dofs on a shared triangulation. Tests the change
// from coin_flip to smallest proc index method of distribution
// along partition interface

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
test()
{
  parallel::shared::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::limit_level_difference_at_vertices),
    true,
    typename parallel::shared::Triangulation<dim>::Settings(
      parallel::shared::Triangulation<dim>::partition_zorder));
  DoFHandler<dim> dof_handler(tria);
  GridGenerator::subdivided_hyper_cube(tria, 3, -1, 1);
  {
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      {
        if (cell->index() == 0 || cell->index() == 5 || cell->index() == 6)
          cell->set_refine_flag();
      }
    tria.execute_coarsening_and_refinement();
  }

  FE_Q<dim> fe(2);
  dof_handler.distribute_dofs(fe);

  deallog << "Number of locally owned dofs: "
          << dof_handler.n_locally_owned_dofs() << std::endl;

  std::vector<IndexSet> shared_dofs_per_proc =
    Utilities::MPI::all_gather(MPI_COMM_WORLD,
                               dof_handler.locally_owned_dofs());
  for (unsigned int i = 0; i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
       ++i)
    shared_dofs_per_proc[i].print(deallog.get_file_stream());

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      if (cell->subdomain_id() == numbers::artificial_subdomain_id)
        continue;

      std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
      deallog << "cell" << cell->index() << " has dofs: ";
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        deallog << local_dof_indices[i] << ' ';
      deallog << std::endl;
    }
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
