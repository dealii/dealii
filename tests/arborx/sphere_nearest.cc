// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check ArborX wrapper: nearest spheres from  others points


#include <deal.II/arborx/access_traits.h>
#include <deal.II/arborx/bvh.h>

#include <deal.II/base/point.h>

#include "../tests.h"



void
test_1d()
{
  std::vector<Point<1, float>> points;

  unsigned int n_points_1d = 5;
  for (unsigned int i = 0; i < n_points_1d; ++i)
    points.emplace_back(i);

  std::vector<std::pair<Point<1, float>, float>> query_spheres;
  query_spheres.emplace_back(Point<1, float>(0.5), 0.1);
  query_spheres.emplace_back(Point<1, float>(1.5), 0.1);
  query_spheres.emplace_back(Point<1, float>(2.2), 0.1);
  query_spheres.emplace_back(Point<1, float>(2.6), 0.1);


#if ARBORX_VERSION_MAJOR < 2
  ArborXWrappers::BVH bvh(points);
#else
  std::vector<BoundingBox<1, float>> bounding_boxes;
  for (const auto &p : points)
    {
      bounding_boxes.emplace_back(std::pair(p, p));
    }
  ArborXWrappers::BVH<BoundingBox<1, float>> bvh(bounding_boxes);
#endif
  ArborXWrappers::SphereNearestPredicate sph_nearest(query_spheres, 1);
  auto             indices_offsets = bvh.query(sph_nearest);
  std::vector<int> indices         = indices_offsets.first;
  std::vector<int> offsets         = indices_offsets.second;

  std::vector<int> indices_ref = {0, 1, 2, 3};
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
test_2d()
{
  std::vector<Point<2, float>> points;

  unsigned int n_points_1d = 5;
  for (unsigned int i = 0; i < n_points_1d; ++i)
    {
      for (unsigned int j = 0; j < n_points_1d; ++j)
        {
          points.emplace_back(j, i);
        }
    }


  std::vector<std::pair<Point<2, float>, float>> query_spheres;
  query_spheres.push_back(std::make_pair(Point<2, float>(0.5, 0.5), 0.1));
  query_spheres.push_back(std::make_pair(Point<2, float>(1.5, 1.5), 0.1));
  query_spheres.push_back(std::make_pair(Point<2, float>(2.2, 2.2), 0.1));
  query_spheres.push_back(std::make_pair(Point<2, float>(2.6, 2.6), 0.1));


#if ARBORX_VERSION_MAJOR < 2
  ArborXWrappers::BVH bvh(points);
#else
  std::vector<BoundingBox<2, float>>         bounding_boxes;
  for (const auto &p : points)
    {
      bounding_boxes.emplace_back(std::pair(p, p));
    }
  ArborXWrappers::BVH<BoundingBox<2, float>> bvh(bounding_boxes);
#endif
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
  std::vector<Point<3, float>> points;

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

  std::vector<std::pair<Point<3, float>, float>> query_spheres;
  query_spheres.push_back(std::make_pair(Point<3, float>(0.5, 0.5, 0.5), 0.1));
  query_spheres.push_back(std::make_pair(Point<3, float>(1.5, 1.5, 1.5), 0.1));
  query_spheres.push_back(std::make_pair(Point<3, float>(2.2, 2.2, 2.2), 0.1));
  query_spheres.push_back(std::make_pair(Point<3, float>(2.6, 2.6, 2.6), 0.1));


#if ARBORX_VERSION_MAJOR < 2
  ArborXWrappers::BVH bvh(points);
#else
  std::vector<BoundingBox<3, float>>         bounding_boxes;
  for (const auto &p : points)
    {
      bounding_boxes.emplace_back(std::pair(p, p));
    }
  ArborXWrappers::BVH<BoundingBox<3, float>> bvh(bounding_boxes);
#endif
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
  initlog();

  // Initialize ArborX
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  // The 1D test hits a bug in clang:
  // https://github.com/llvm/llvm-project/issues/18060
#if defined(DEAL_II_HAVE_FP_EXCEPTIONS)
  {
    const int current_fe_except = fegetexcept();
    fedisableexcept(current_fe_except);
  }
#endif
  test_1d();
  test_2d();
  test_3d();
}
