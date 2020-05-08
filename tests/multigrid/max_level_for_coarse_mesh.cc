// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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


// check function MGTools::max_level_for_coarse_mesh()

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_cube(tria, 0, 1);
  tria.refine_global(1);
  for (unsigned int cycle = 0; cycle < 2; ++cycle)
    {
      for (typename parallel::distributed::Triangulation<
             dim>::active_cell_iterator cell = tria.begin_active();
           cell != tria.end();
           ++cell)
        for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
          {
            if (dim == 2)
              if (cell->vertex(v)[0] < 0.25 && cell->vertex(v)[1] < 0.25)
                {
                  cell->set_refine_flag();
                  break;
                }
            if (dim == 3)
              if (cell->vertex(v)[0] < 0.25 && cell->vertex(v)[1] < 0.25 &&
                  cell->vertex(v)[2] < 0.25)
                {
                  cell->set_refine_flag();
                  break;
                }
          }
      tria.execute_coarsening_and_refinement();
    }

  const unsigned int max_possible_level =
    MGTools::max_level_for_coarse_mesh(tria);
  deallog << "Max possible level for coarse mesh: " << max_possible_level
          << std::endl;
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
