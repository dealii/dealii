// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

// test level subdomain ids for periodic boundary conditions

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::subdivided_hyper_cube(tria, 3, 0., 1.);

  for (const auto &cell : tria.cell_iterators())
    for (unsigned int face_index : GeometryInfo<dim>::face_indices())
      {
        if (std::abs(cell->face(face_index)->center()(face_index / 2)) < 1e-12)
          cell->face(face_index)->set_all_boundary_ids(face_index);
        if (std::abs(cell->face(face_index)->center()(face_index / 2) - 1.) <
            1e-12)
          cell->face(face_index)->set_all_boundary_ids(face_index);
      }

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;
  for (unsigned int d = 0; d < dim; ++d)
    GridTools::collect_periodic_faces(static_cast<Triangulation<dim> &>(tria),
                                      2 * d,
                                      2 * d + 1,
                                      d,
                                      periodic_faces);

  tria.add_periodicity(periodic_faces);

  for (unsigned int i = 0; i < 2; ++i)
    {
      for (const auto &cell : tria.cell_iterators())
        if (cell->level_subdomain_id() == tria.locally_owned_subdomain())
          {
            deallog << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) << " "
                    << cell->id() << " neighbor subdomain ids: ";
            for (unsigned int f : GeometryInfo<dim>::face_indices())
              {
                deallog << cell->neighbor_or_periodic_neighbor(f)->id() << " ";
                if (cell->is_active())
                  deallog
                    << cell->neighbor_or_periodic_neighbor(f)->subdomain_id()
                    << " ";
                deallog << cell->neighbor_or_periodic_neighbor(f)
                             ->level_subdomain_id()
                        << "  ";
              }
            deallog << std::endl;
          }
      tria.refine_global();
    }
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  try
    {
      test<2>();
      test<3>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
