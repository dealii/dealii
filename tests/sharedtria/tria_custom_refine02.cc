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


// create a shared tria mesh and distribute with a custom function
// tests that the correct cells are set to artificial in multigrid

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim>
void
mypartition(parallel::shared::Triangulation<dim> &tria)
{
  std::vector<unsigned int> assignment = {
    0, 0, 1, 2, 0, 0, 2, 1, 0, 2, 2, 1, 2, 2, 0, 0};
  {
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    unsigned int index                                     = 0;
    for (; cell != endc; ++cell, ++index)
      cell->set_subdomain_id(assignment[index % 16]);
  }

  for (int lvl = tria.n_levels() - 1; lvl >= 0; --lvl)
    {
      typename parallel::shared::Triangulation<dim>::cell_iterator cell =
                                                                     tria.begin(
                                                                       lvl),
                                                                   endc =
                                                                     tria.end(
                                                                       lvl);
      for (; cell != endc; ++cell)
        {
          if (cell->is_active())
            cell->set_level_subdomain_id(cell->subdomain_id());
          else
            cell->set_level_subdomain_id(cell->child(0)->level_subdomain_id());
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
      parallel::shared::Triangulation<dim>::partition_custom_signal |
      parallel::shared::Triangulation<dim>::construct_multigrid_hierarchy));
  shared_tria.signals.post_refinement.connect(
    std::bind(&mypartition<dim>, std::ref(shared_tria)));


  GridGenerator::hyper_cube(shared_tria);
  shared_tria.refine_global(2);

  {
    deallog << "(CellId,subdomain_id) for each active cell:" << std::endl;
    typename Triangulation<dim>::active_cell_iterator cell = shared_tria
                                                               .begin_active(),
                                                      endc = shared_tria.end();
    for (; cell != endc; ++cell)
      if (cell->subdomain_id() != numbers::artificial_subdomain_id)
        deallog << '(' << cell->id().to_string() << ',' << cell->subdomain_id()
                << ')' << std::endl;
  }

  {
    deallog << "(CellId,level_subdomain_id) for each cell:" << std::endl;
    typename Triangulation<dim>::cell_iterator cell = shared_tria.begin(),
                                               endc = shared_tria.end();
    for (; cell != endc; ++cell)
      if (cell->level_subdomain_id() != numbers::artificial_subdomain_id)
        deallog << '(' << cell->id().to_string() << ','
                << cell->level_subdomain_id() << ')' << std::endl;
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
