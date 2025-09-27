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


// distribute mg dofs on a shared and distributed mesh and compare
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
  shared_dof_handler.distribute_mg_dofs();

  distributed_dof_handler.distribute_dofs(fe);
  distributed_dof_handler.distribute_mg_dofs();

  unsigned int n_levels =
    distributed_dof_handler.get_triangulation().n_global_levels();
  for (unsigned int lvl = 0; lvl < n_levels; ++lvl)
    {
      IndexSet shared_dofs = shared_dof_handler.locally_owned_mg_dofs(lvl);
      IndexSet distributed_dofs =
        distributed_dof_handler.locally_owned_mg_dofs(lvl);
      Assert(shared_dofs == distributed_dofs, ExcInternalError());

      typename DoFHandler<dim>::cell_iterator cell =
                                                distributed_dof_handler.begin(
                                                  lvl),
                                              endc =
                                                distributed_dof_handler.end(
                                                  lvl);
      for (; cell != endc; ++cell)
        {
          if (cell->level_subdomain_id() == numbers::artificial_subdomain_id)
            continue;

          typename Triangulation<dim>::cell_iterator tria_shared_cell =
            shared_dof_handler.get_triangulation().create_cell_iterator(
              cell->id());
          typename DoFHandler<dim>::cell_iterator dof_shared_cell(
            &shared_dof_handler.get_triangulation(),
            tria_shared_cell->level(),
            tria_shared_cell->index(),
            &shared_dof_handler);

          std::vector<types::global_dof_index> distributed_cell_dofs(
            fe.dofs_per_cell);
          std::vector<types::global_dof_index> shared_cell_dofs(
            fe.dofs_per_cell);
          cell->get_mg_dof_indices(distributed_cell_dofs);
          dof_shared_cell->get_mg_dof_indices(shared_cell_dofs);
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            Assert(distributed_cell_dofs[i] == shared_cell_dofs[i],
                   ExcInternalError());
        }
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
      parallel::shared::Triangulation<dim>::partition_zorder |
      parallel::shared::Triangulation<dim>::construct_multigrid_hierarchy));
  DoFHandler<dim> shared_dof_handler(shared_tria);

  parallel::distributed::Triangulation<dim> distributed_tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  DoFHandler<dim> distributed_dof_handler(distributed_tria);

  GridGenerator::hyper_cube(shared_tria);
  shared_tria.refine_global(2);
  for (typename Triangulation<dim>::active_cell_iterator cell =
         shared_tria.begin_active();
       cell != shared_tria.end();
       ++cell)
    {
      if (dim == 1)
        if (cell->center()[0] < 0.5)
          cell->set_refine_flag();
      if (dim == 2)
        if (cell->center()[0] < 0.5 && cell->center()[1] < 0.5)
          cell->set_refine_flag();
      if (dim == 3)
        if (cell->center()[0] < 0.5 && cell->center()[1] < 0.5 &&
            cell->center()[2] < 0.5)
          cell->set_refine_flag();
    }
  shared_tria.execute_coarsening_and_refinement();

  GridGenerator::hyper_cube(distributed_tria);
  distributed_tria.refine_global(2);
  for (typename Triangulation<dim>::active_cell_iterator cell =
         distributed_tria.begin_active();
       cell != distributed_tria.end();
       ++cell)
    {
      if (dim == 1)
        if (cell->is_locally_owned() && cell->center()[0] < 0.5)
          cell->set_refine_flag();
      if (dim == 2)
        if (cell->is_locally_owned() && cell->center()[0] < 0.5 &&
            cell->center()[1] < 0.5)
          cell->set_refine_flag();
      if (dim == 3)
        if (cell->is_locally_owned() && cell->center()[0] < 0.5 &&
            cell->center()[1] < 0.5 && cell->center()[2] < 0.5)
          cell->set_refine_flag();
    }
  distributed_tria.execute_coarsening_and_refinement();

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
