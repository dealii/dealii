// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

// test Triangulation::communicate_locally_moved_vertices() on a
// Triangulation with periodic faces.

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

template <int dim>
void
test()
{
  deallog << "Testing " << Utilities::int_to_string(dim, 1) << "D" << std::endl;

  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation, -1., 1., true);

  std::vector<GridTools::PeriodicFacePair<
    typename parallel::distributed::Triangulation<dim>::cell_iterator>>
    periodicity_vector;
  for (int i = 0; i < dim; ++i)
    GridTools::collect_periodic_faces(
      triangulation, 2 * i, 2 * i + 1, i, periodicity_vector);
  triangulation.add_periodicity(periodicity_vector);

  triangulation.refine_global(3);

  auto       cell = triangulation.begin_active();
  const auto endc = triangulation.end();

  std::map<unsigned int, Point<dim>> non_artificial_vertices_old;
  for (cell = triangulation.begin_active(); cell != endc; ++cell)
    if (!cell->is_artificial())
      for (const unsigned int vertex_no : GeometryInfo<dim>::vertex_indices())
        non_artificial_vertices_old[cell->vertex_index(vertex_no)] =
          cell->vertex(vertex_no);

  std::vector<bool>       vertex_moved(triangulation.n_vertices(), false);
  const std::vector<bool> locally_owned_vertices =
    GridTools::get_locally_owned_vertices(triangulation);
  for (cell = triangulation.begin_active(); cell != endc; ++cell)
    if (cell->is_locally_owned())
      for (const unsigned int vertex_no : GeometryInfo<dim>::vertex_indices())
        {
          const unsigned global_vertex_no = cell->vertex_index(vertex_no);
          if (!vertex_moved[global_vertex_no] &&
              locally_owned_vertices[global_vertex_no])
            {
              cell->vertex(vertex_no)(0) += 1.e-1;
              vertex_moved[global_vertex_no] = true;
            }
        }

  triangulation.communicate_locally_moved_vertices(vertex_moved);

  std::map<unsigned int, Point<dim>> non_artificial_vertices_new;
  for (cell = triangulation.begin_active(); cell != endc; ++cell)
    if (!cell->is_artificial())
      for (const unsigned int vertex_no : GeometryInfo<dim>::vertex_indices())
        {
          Point<dim> point = cell->vertex(vertex_no);
          point(0) -= 1.e-1;
          non_artificial_vertices_new[cell->vertex_index(vertex_no)] = point;
        }

  for (const auto &pair : non_artificial_vertices_new)
    if ((non_artificial_vertices_old[pair.first] - pair.second).norm_square() >
        1.e-6)
      {
        deallog << pair.first << ": " << non_artificial_vertices_old[pair.first]
                << " vs. " << pair.second << std::endl;
        AssertThrow(
          false,
          ExcMessage(
            "Some of the vertices on ghost cell were not moved correctly!"));
      }

  /*std::string
  filename("grid"+Utilities::int_to_string(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD),2)+".vtu");
  std::ofstream out(filename.c_str());
  GridOut grid_out;
  grid_out.write_vtu(triangulation, out);*/
  deallog << Utilities::int_to_string(dim, 1) << "D OK" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  MPILogInitAll                    log;
  test<2>();
  test<3>();
  return 0;
}
