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


// distribute mg dofs for a shared triangulation and
// check level_ghost_owners is set properly

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
write_dof_data(DoFHandler<dim> &dof_handler)
{
  FE_Q<dim> fe(2);

  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  unsigned int n_levels = dof_handler.get_triangulation().n_global_levels();
  for (unsigned int lvl = 0; lvl < n_levels; ++lvl)
    {
      std::vector<IndexSet> dof_index_per_proc =
        Utilities::MPI::all_gather(MPI_COMM_WORLD,
                                   dof_handler.locally_owned_mg_dofs(lvl));
      for (unsigned int i = 0; i < dof_index_per_proc.size(); ++i)
        dof_index_per_proc[i].print(deallog);


      typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin(lvl),
                                              endc = dof_handler.end(lvl);
      for (; cell != endc; ++cell)
        {
          if (cell->level_subdomain_id() == numbers::artificial_subdomain_id)
            continue;

          std::vector<types::global_dof_index> local_mg_dof_indices(
            fe.dofs_per_cell);
          cell->get_mg_dof_indices(local_mg_dof_indices);
          deallog << "proc " << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                  << ", "
                  << "cell " << cell->id() << ", mg_dof_indices: ";
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            deallog << local_mg_dof_indices[i] << ' ';
          deallog << std::endl;
        }
    }
}



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
      parallel::shared::Triangulation<dim>::partition_zorder |
      parallel::shared::Triangulation<dim>::construct_multigrid_hierarchy));

  DoFHandler<dim> dof_handler(tria);

  GridGenerator::hyper_cube(tria, -1, 1);

  tria.refine_global(1);

  for (typename Triangulation<dim>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell)
    {
      if (dim == 2)
        if (cell->center()[0] < 0)
          cell->set_refine_flag();
      if (dim == 3)
        if (cell->center()[0] < 0 && cell->center()[1] < 0)
          cell->set_refine_flag();
    }
  tria.execute_coarsening_and_refinement();

  deallog << "proc " << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
          << ", level_ghost_owners: ";
  for (std::set<types::subdomain_id>::iterator it =
         tria.level_ghost_owners().begin();
       it != tria.level_ghost_owners().end();
       ++it)
    deallog << *it << ' ';
  deallog << std::endl;

  write_dof_data(dof_handler);
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
}
