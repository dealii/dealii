// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Check global level cell ids when construct_multigrid_hierarchy is enabled.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test(int n_refinements, MPI_Comm comm)
{
  parallel::distributed::Triangulation<dim> tria(
    comm,
    Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_refinements);

  for (unsigned int l = 0; l < tria.n_global_levels(); ++l)
    {
      deallog.push("level=" + std::to_string(l));
      for (const auto cell : tria.cell_iterators_on_level(l))
        if (cell->level_subdomain_id() !=
            dealii::numbers::artificial_subdomain_id)
          deallog << cell->id() << " -> " << cell->level_subdomain_id() << ' '
                  << cell->global_level_cell_index() << std::endl;

      const Utilities::MPI::Partitioner &part =
        *tria.global_level_cell_index_partitioner(l).lock();

      part.locally_owned_range().print(deallog);
      part.ghost_indices().print(deallog);

      deallog << std::endl;
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const int      n_refinements = 2;
  const MPI_Comm comm          = MPI_COMM_WORLD;

  {
    deallog.push("2d");
    test<2>(n_refinements, comm);
    deallog.pop();
  }
  if (false)
    {
      deallog.push("3d");
      test<3>(n_refinements, comm);
      deallog.pop();
    }
}
