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


// Check ArborX wrapper: intersection of spheres with points


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
  query_spheres.push_back(std::make_pair(Point<2>(0.5, 0.5), 0.8));
  query_spheres.push_back(std::make_pair(Point<2>(1.5, 1.5), 0.8));
  query_spheres.push_back(std::make_pair(Point<2>(2.2, 2.2), 1.2));
  query_spheres.push_back(std::make_pair(Point<2>(2.6, 2.6), 0.9));


  ArborXWrappers::BVH                      bvh(points);
  ArborXWrappers::SphereIntersectPredicate sph_intersect(query_spheres);
  auto             indices_offsets = bvh.query(sph_intersect);
  std::vector<int> indices         = indices_offsets.first;
  std::vector<int> offsets         = indices_offsets.second;

  std::vector<int> indices_ref = {
    0, 5, 1, 6, 6, 11, 7, 12, 12, 17, 13, 18, 12, 17, 13, 18};
  std::vector<int> offsets_ref = {0, 4, 8, 12, 16};

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
  query_spheres.push_back(std::make_pair(Point<3>(0.5, 0.5, 0.5), 0.9));
  query_spheres.push_back(std::make_pair(Point<3>(1.5, 1.5, 1.5), 0.9));
  query_spheres.push_back(std::make_pair(Point<3>(2.2, 2.2, 2.2), 1.4));
  query_spheres.push_back(std::make_pair(Point<3>(2.6, 2.6, 2.6), 1.1));


  ArborXWrappers::BVH                      bvh(points);
  ArborXWrappers::SphereIntersectPredicate sph_intersect(query_spheres);
  auto             indices_offsets = bvh.query(sph_intersect);
  std::vector<int> indices         = indices_offsets.first;
  std::vector<int> offsets         = indices_offsets.second;

  std::vector<int> indices_ref = {0,  25, 5,  30, 1,  26, 6,  31, 31,
                                  56, 36, 61, 32, 57, 37, 62, 61, 57,
                                  37, 62, 87, 67, 92, 63, 88, 68, 93,
                                  62, 87, 67, 92, 63, 88, 68, 93};
  std::vector<int> offsets_ref = {0, 8, 16, 27, 35};

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
