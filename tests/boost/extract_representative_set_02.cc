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

// Extract a representative vector of BoundingBox objects from an RTree

#include <deal.II/base/patterns.h>

#include <deal.II/boost_adaptors/bounding_box.h>
#include <deal.II/boost_adaptors/point.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/numerics/rtree.h>

#include <algorithm>

#include "../tests.h"


std::string
print(const BoundingBox<2> &box)
{
  std::stringstream str;

  const auto p  = box.get_boundary_points();
  const auto p1 = p.first;
  const auto p2 = p.second;

  str << p1 << std::endl
      << p2[0] << " " << p1[1] << std::endl
      << p2 << std::endl
      << p1[0] << " " << p2[1] << std::endl
      << p1 << std::endl;
  return str.str();
}

void
test(const unsigned int n_refinements)
{
  parallel::distributed::Triangulation<2> tria(MPI_COMM_WORLD);

  GridGenerator::hyper_ball(tria);
  tria.refine_global(n_refinements);

  std::vector<BoundingBox<2>> all_boxes(tria.n_locally_owned_active_cells());
  unsigned int                i = 0;
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      all_boxes[i++] = cell->bounding_box();

  const auto tree  = pack_rtree(all_boxes);
  const auto boxes = extract_rtree_level(tree, 0);

  deallog << "Refinements: " << n_refinements << std::endl;
  deallog << "LEVEL 0:  N boxes: " << boxes.size() << std::endl;
  for (const auto &b : boxes)
    deallog << print(b) << std::endl;

  // Uncomment the following lines to generate the images of the documentation
  //
  // std::ofstream ofile(
  //   "output_" +
  //   std::to_string(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)) +
  //   ".gpl");
  //
  // std::ofstream all(
  //   "all_" + std::to_string(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
  //   +
  //   ".gpl");
  //
  // for (const auto &b : all_boxes)
  //   all << print(b);
}


int
main(int argc, char **argv)
{
  auto          mpi_init = Utilities::MPI::MPI_InitFinalize(argc, argv);
  MPILogInitAll log;

  for (unsigned int i = 0; i <= 3; ++i)
    test(i);
}
