// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019-2019 by the deal.II authors
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

// Test GridTools::build_global_description_tree


#include <deal.II/base/bounding_box.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>

#include <deal.II/boost_adaptors/bounding_box.h>
#include <deal.II/boost_adaptors/point.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/rtree.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/geometry/index/detail/serialization.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/signals2.hpp>

#include <algorithm>
#include <utility>

#include "../tests.h"


namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;


template <int dim, int spacedim>
void
test()
{
  unsigned int n_procs      = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  unsigned int current_proc = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int next_p       = (current_proc + n_procs + 1) % n_procs;

  // For this test we actually set up a simple system of bounding boxes
  // The triangulation is really used only to get the right mpi_communicator
  Point<spacedim> p1;
  Point<spacedim> p2;
  for (unsigned int d = 1; d < spacedim; ++d)
    {
      p1[d] = current_proc;
      p2[d] = current_proc + 1;
    }

  BoundingBox<spacedim> box({p1, p2});

  std::vector<BoundingBox<spacedim>> local_description(1, box);

  auto built_tree =
    GridTools::build_global_description_tree(local_description, MPI_COMM_WORLD);

  Point<spacedim> my_point;
  Point<spacedim> point_inside_d_1;
  Point<spacedim> outside_point;

  for (unsigned int d = 1; d < spacedim; ++d)
    {
      my_point[d]         = current_proc + 0.1 + d / 4.2;
      point_inside_d_1[d] = next_p + 0.25 + d / 5.1;
      outside_point[d]    = -current_proc - 2.6;
    }

  std::vector<std::pair<BoundingBox<spacedim>, unsigned int>> test_results;
  built_tree.query(bgi::intersects(outside_point),
                   std::back_inserter(test_results));
  if (test_results.size() != 0)
    deallog
      << "Point found inside a bounding box! It should be outside all of them!"
      << std::endl;

  built_tree.query(bgi::intersects(my_point), std::back_inserter(test_results));
  if (std::get<1>(test_results[0]) != current_proc)
    deallog << "Error: Point found inside wrong process: "
            << std::get<1>(test_results[0]) << std::endl;

  built_tree.query(bgi::intersects(point_inside_d_1),
                   std::back_inserter(test_results));
  if (std::get<1>(test_results[1]) != next_p)
    deallog << "Error: Point found inside wrong process: "
            << std::get<1>(test_results[1]) << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
#ifdef DEAL_II_WITH_MPI
  MPILogInitAll log;
#else
  initlog();
  deallog.push("0");
#endif

  deallog << "Deal.II build_global_description_tree:" << std::endl;

  test<2, 2>();
  test<3, 3>();

  deallog << "Test finished" << std::endl;
  return 0;
}
