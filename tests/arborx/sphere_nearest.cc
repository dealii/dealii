// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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


// Check ArborX wrapper: nearest spheres from  others points


#include <deal.II/arborx/access_traits.h>
#include <deal.II/arborx/bvh.h>

#include <deal.II/base/point.h>

#include "../tests.h"



void
test_2d()
{
  std::vector<Point<2>> points;

  unsigned int n_points_1d = 5;
  for (unsigned int i = 0; i < n_points_1d; ++i)
    {
      for (unsigned int j = 0; j < n_points_1d; ++j)
        {
          points.emplace_back(j, i);
        }
    }


  std::vector<std::pair<Point<2>, double>> query_spheres;
  query_spheres.push_back(std::make_pair(Point<2>(0.5, 0.5), 0.1));
  query_spheres.push_back(std::make_pair(Point<2>(1.5, 1.5), 0.1));
  query_spheres.push_back(std::make_pair(Point<2>(2.2, 2.2), 0.1));
  query_spheres.push_back(std::make_pair(Point<2>(2.6, 2.6), 0.1));


  ArborXWrappers::BVH                    bvh(points);
  ArborXWrappers::SphereNearestPredicate sph_nearest(query_spheres, 1);
  auto             indices_offsets = bvh.query(sph_nearest);
  std::vector<int> indices         = indices_offsets.first;
  std::vector<int> offsets         = indices_offsets.second;

  std::vector<int> indices_ref = {0, 6, 12, 18};
  std::vector<int> offsets_ref = {0, 1, 2, 3, 4};

  AssertThrow(indices.size() == indices_ref.size(), ExcInternalError());
  AssertThrow(offsets.size() == offsets_ref.size(), ExcInternalError());
  for (unsigned int i = 0; i < offsets.size() - 1; ++i)
    {
      for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
          // The indices associated to each query are not ordered.
          bool found = false;
          for (int k = offsets[i]; k < offsets[i + 1]; ++k)
            {
              if (indices[j] == indices_ref[k])
                {
                  found = true;
                  break;
                }
            }
          AssertThrow(found, ExcInternalError());
        }
    }
  for (unsigned int i = 0; i < offsets.size(); ++i)
    AssertThrow(offsets[i] == offsets_ref[i], ExcInternalError());

  deallog << "OK" << std::endl;
}


void
test_3d()
{
  std::vector<Point<3>> points;

  unsigned int n_points_1d = 5;
  for (unsigned int i = 0; i < n_points_1d; ++i)
    {
      for (unsigned int j = 0; j < n_points_1d; ++j)
        {
          for (unsigned int k = 0; k < n_points_1d; ++k)
            {
              points.emplace_back(k, j, i);
            }
        }
    }

  std::vector<std::pair<Point<3>, double>> query_spheres;
  query_spheres.push_back(std::make_pair(Point<3>(0.5, 0.5, 0.5), 0.1));
  query_spheres.push_back(std::make_pair(Point<3>(1.5, 1.5, 1.5), 0.1));
  query_spheres.push_back(std::make_pair(Point<3>(2.2, 2.2, 2.2), 0.1));
  query_spheres.push_back(std::make_pair(Point<3>(2.6, 2.6, 2.6), 0.1));


  ArborXWrappers::BVH                    bvh(points);
  ArborXWrappers::SphereNearestPredicate sph_nearest(query_spheres, 1);
  auto             indices_offsets = bvh.query(sph_nearest);
  std::vector<int> indices         = indices_offsets.first;
  std::vector<int> offsets         = indices_offsets.second;

  std::vector<int> indices_ref = {0, 31, 62, 93};
  std::vector<int> offsets_ref = {0, 1, 2, 3, 4};

  AssertThrow(indices.size() == indices_ref.size(), ExcInternalError());
  AssertThrow(offsets.size() == offsets_ref.size(), ExcInternalError());
  for (unsigned int i = 0; i < offsets.size() - 1; ++i)
    {
      for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
          // The indices associated to each query are not ordered.
          bool found = false;
          for (int k = offsets[i]; k < offsets[i + 1]; ++k)
            {
              if (indices[j] == indices_ref[k])
                {
                  found = true;
                  break;
                }
            }
          AssertThrow(found, ExcInternalError());
        }
    }
  for (unsigned int i = 0; i < offsets.size(); ++i)
    AssertThrow(offsets[i] == offsets_ref[i], ExcInternalError());
  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  // Initialize Kokkos
  Kokkos::initialize(argc, argv);

  initlog();

  // Initialize ArborX
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  // tests
  test_2d();
  test_3d();

  Kokkos::finalize();
}
