// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2023 by the deal.II authors
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

// Test distributed_compute_intersection_locations() for face intersections
// (1D-2D).

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
std::vector<std::vector<BoundingBox<dim>>>
get_global_bboxes(const Triangulation<dim> &tria,
                  const Mapping<dim>       &mapping,
                  const unsigned int        rtree_level = 0)
{
  std::vector<dealii::BoundingBox<dim>> local_boxes;
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      local_boxes.emplace_back(mapping.get_bounding_box(cell));

  // create r-tree of bounding boxes
  const auto local_tree = pack_rtree(local_boxes);

  // compress r-tree to a minimal set of bounding boxes
  std::vector<std::vector<BoundingBox<dim>>> global_bboxes(1);
  global_bboxes[0] = extract_rtree_level(local_tree, rtree_level);

  return global_bboxes;
}

template <typename Container>
void
print_vertices(const Container &vertices)
{
  deallog << "{";
  unsigned int i = 0;
  for (; i < vertices.size() - 1; ++i)
    {
      deallog << "(" << vertices[i] << "), ";
    }
  deallog << "(" << vertices[i] << ")}";
}

template <int structdim, int dim>
void
do_test(const bool use_marked_vertices = false)
{
  // create triangulation
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  // create cache
  MappingQ1<dim>              mapping;
  const GridTools::Cache<dim> cache(tria, mapping);

  // create intersection requests
  std::vector<std::vector<Point<dim>>> intersection_requests;
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      Triangulation<dim> tria_other;
      GridGenerator::hyper_cube(tria_other, 0.25, 0.75);
      const auto cell = tria_other.begin_active();

      Triangulation<dim> tria_temp;
      GridGenerator::hyper_cube(tria_temp, 0.75, 0.9);
      Triangulation<dim> simplex_tria;
      GridGenerator::convert_hypercube_to_simplex_mesh(tria_temp, simplex_tria);
      const auto simplex_cell = tria.begin_active();

      for (const auto f : cell->face_indices())
        {
          const unsigned int      n_vertices = cell->face(f)->n_vertices();
          std::vector<Point<dim>> vertices(n_vertices);
          std::copy_n(mapping.get_vertices(cell, f).begin(),
                      n_vertices,
                      vertices.begin());

          intersection_requests.emplace_back(vertices);
        }
    }

  std::vector<bool> marked_vertices;
  if (use_marked_vertices) // mark vertices for x < 0.5
    {
      marked_vertices.resize(tria.n_vertices(), false);

      for (const auto &cell : tria.active_cell_iterators())
        for (const auto v : cell->vertex_indices())
          if (cell->vertex(v)[0] < 0.49)
            marked_vertices[cell->vertex_index(v)] = true;
    }

  auto intersection_location =
    GridTools::internal::distributed_compute_intersection_locations<structdim>(
      cache,
      intersection_requests,
      get_global_bboxes<dim>(tria, mapping),
      marked_vertices,
      1.0e-9);

  deallog << "Recv Components " << std::endl;
  for (const auto &rc : intersection_location.recv_components)
    {
      deallog << std::get<0>(rc) << "; " << std::get<1>(rc) << "; ";
      print_vertices(std::get<2>(rc));
      deallog << std::endl;
    }
  deallog << std::endl;

  deallog << "Send Components " << std::endl;
  for (const auto &sc : intersection_location.send_components)
    {
      deallog << "<" << std::get<0>(sc).first << ", " << std::get<0>(sc).second
              << ">; " << std::get<1>(sc) << "; " << std::get<2>(sc) << "; ";
      print_vertices(std::get<3>(sc));
      deallog << std::endl;
    }
  deallog << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;
  deallog.precision(8);

  deallog << "Face intersections 1D-2D with marked vertices" << std::endl;
  do_test<1, 2>(true);
  deallog << "Face intersections 1D-2D without marked vertices" << std::endl;
  do_test<1, 2>(false);
}
