// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

// this tests distribute_mg_dofs on a mesh with periodic boundary
// conditions. we used to forget to enforce the 2:1 level balance over
// periodic boundaries in artificial cells, which lead to the case where the
// coarser MG level did not provide the correct ghost layer. This test
// explicitly lists all cells on all processors (processor 1 used to miss a
// few cells in the implementation before January 2016).

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    dealii::Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube(tria, 3);

  // set periodic boundary conditions in all directions
  for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
       cell != tria.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        cell->face(f)->set_all_boundary_ids(f + 10);

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;
  GridTools::collect_periodic_faces(tria, 0 + 10, 1 + 10, 0, periodic_faces);
  GridTools::collect_periodic_faces(tria, 2 + 10, 3 + 10, 1, periodic_faces);
  tria.add_periodicity(periodic_faces);

  // adaptively refine into the lower left corner
  for (unsigned int i = 0; i < 2; ++i)
    {
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tria.begin_active();
           cell != tria.end();
           ++cell)
        if (cell->at_boundary(0) && cell->at_boundary(2) &&
            cell->is_locally_owned())
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
       cell != tria.end();
       ++cell)
    deallog << cell->id().to_string() << " " << cell->level_subdomain_id()
            << std::endl;

  if (0)
    {
      std::ofstream grid_output(
        ("out" + Utilities::to_string(myid) + ".svg").c_str());
      GridOut           grid_out;
      GridOutFlags::Svg flags;
      flags.label_level_subdomain_id = true;
      flags.coloring                 = GridOutFlags::Svg::level_subdomain_id;
      flags.convert_level_number_to_height = true;
      grid_out.set_flags(flags);

      grid_out.write_svg(tria, grid_output);
    }
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();

  return 0;
}
