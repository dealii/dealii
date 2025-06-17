// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check ArborX wrapper: intersection of bounding boxes


#include <deal.II/arborx/access_traits.h>
#include <deal.II/arborx/bvh.h>

#include <deal.II/base/bounding_box.h>

#include "../tests.h"

void
test_1d()
{
  std::vector<BoundingBox<1>> bounding_boxes;
  std::vector<Point<1>>       points;

  unsigned int n_points_1d = 5;
  for (unsigned int i = 0; i < n_points_1d; ++i)
    points.emplace_back(i);

  for (unsigned int i = 0; i < n_points_1d - 1; ++i)
    bounding_boxes.emplace_back(std::make_pair(points[i], points[i + 1]));

  Point<1> point_a(0.5);
  Point<1> point_b(1.5);
  Point<1> point_c(2.2);
  Point<1> point_d(2.6);

  std::vector<BoundingBox<1>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

#if ARBORX_VERSION_MAJOR < 2
  ArborXWrappers::BVH bvh(bounding_boxes);
#else
  ArborXWrappers::BVH<BoundingBox<1>> bvh(bounding_boxes);
#endif
  ArborXWrappers::BoundingBoxIntersectPredicate bb_intersect(
    query_bounding_boxes);
  auto             indices_offset = bvh.query(bb_intersect);
  std::vector<int> indices        = indices_offset.first;
  std::vector<int> offset         = indices_offset.second;

  std::vector<int> indices_ref = {0, 1, 2};
  std::vector<int> offset_ref  = {0, 2, 3};

  AssertDimension(indices.size(), indices_ref.size());
  AssertDimension(offset.size(), offset_ref.size());
  for (unsigned int i = 0; i < offset.size() - 1; ++i)
    {
      for (unsigned int j = offset[i]; j < offset[i + 1]; ++j)
        {
          // The indices associated to each query are not ordered.
          bool found = false;
          for (unsigned int k = offset[i]; k < offset[i + 1]; ++k)
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
  for (unsigned int i = 0; i < offset.size(); ++i)
    AssertThrow(offset[i] == offset_ref[i], ExcInternalError());

  deallog << "OK" << std::endl;
}

void
test_2d()
{
  std::vector<BoundingBox<2>> bounding_boxes;
  std::vector<Point<2>>       points;

  unsigned int n_points_1d = 5;
  for (unsigned int i = 0; i < n_points_1d; ++i)
    {
      for (unsigned int j = 0; j < n_points_1d; ++j)
        {
          points.emplace_back(j, i);
        }
    }

  for (unsigned int i = 0; i < n_points_1d - 1; ++i)
    {
      for (unsigned int j = 0; j < n_points_1d - 1; ++j)
        {
          unsigned int point_index = j + i * n_points_1d;
          bounding_boxes.push_back(
            std::make_pair(points[point_index],
                           points[point_index + n_points_1d + 1]));
        }
    }

  Point<2> point_a(0.5, 0.5);
  Point<2> point_b(1.5, 1.5);
  Point<2> point_c(2.2, 2.2);
  Point<2> point_d(2.6, 2.6);

  std::vector<BoundingBox<2>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

#if ARBORX_VERSION_MAJOR < 2
  ArborXWrappers::BVH bvh(bounding_boxes);
#else
  ArborXWrappers::BVH<BoundingBox<2>> bvh(bounding_boxes);
#endif
  ArborXWrappers::BoundingBoxIntersectPredicate bb_intersect(
    query_bounding_boxes);
  auto             indices_offset = bvh.query(bb_intersect);
  std::vector<int> indices        = indices_offset.first;
  std::vector<int> offset         = indices_offset.second;

  std::vector<int> indices_ref = {0, 1, 4, 5, 10};
  std::vector<int> offset_ref  = {0, 4, 5};

  AssertThrow(indices.size() == indices_ref.size(), ExcInternalError());
  AssertThrow(offset.size() == offset_ref.size(), ExcInternalError());
  for (unsigned int i = 0; i < offset.size() - 1; ++i)
    {
      for (unsigned int j = offset[i]; j < offset[i + 1]; ++j)
        {
          // The indices associated to each query are not ordered.
          bool found = false;
          for (unsigned int k = offset[i]; k < offset[i + 1]; ++k)
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
  for (unsigned int i = 0; i < offset.size(); ++i)
    AssertThrow(offset[i] == offset_ref[i], ExcInternalError());

  deallog << "OK" << std::endl;
}

void
test_3d()
{
  std::vector<BoundingBox<3>> bounding_boxes;
  std::vector<Point<3>>       points;

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

  for (unsigned int i = 0; i < n_points_1d - 1; ++i)
    {
      for (unsigned int j = 0; j < n_points_1d - 1; ++j)
        {
          for (unsigned int k = 0; k < n_points_1d - 1; ++k)
            {
              unsigned int point_index =
                k + j * n_points_1d + i * n_points_1d * n_points_1d;
              bounding_boxes.push_back(
                std::make_pair(points[point_index],
                               points[point_index + n_points_1d * n_points_1d +
                                      n_points_1d + 1]));
            }
        }
    }

  Point<3> point_a(0.5, 0.5, 0.5);
  Point<3> point_b(1.5, 1.5, 1.5);
  Point<3> point_c(2.2, 2.2, 2.2);
  Point<3> point_d(2.6, 2.6, 2.6);

  std::vector<BoundingBox<3>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

#if ARBORX_VERSION_MAJOR < 2
  ArborXWrappers::BVH bvh(bounding_boxes);
#else
  ArborXWrappers::BVH<BoundingBox<3>> bvh(bounding_boxes);
#endif
  ArborXWrappers::BoundingBoxIntersectPredicate bb_intersect(
    query_bounding_boxes);
  auto             indices_offset = bvh.query(bb_intersect);
  std::vector<int> indices        = indices_offset.first;
  std::vector<int> offset         = indices_offset.second;

  std::vector<int> indices_ref = {0, 1, 4, 5, 16, 17, 20, 21, 42};
  std::vector<int> offset_ref  = {0, 8, 9};

  AssertThrow(indices.size() == indices_ref.size(), ExcInternalError());
  AssertThrow(offset.size() == offset_ref.size(), ExcInternalError());
  for (unsigned int i = 0; i < offset.size() - 1; ++i)
    {
      for (unsigned int j = offset[i]; j < offset[i + 1]; ++j)
        {
          // The indices associated to each query are not ordered.
          bool found = false;
          for (unsigned int k = offset[i]; k < offset[i + 1]; ++k)
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
  for (unsigned int i = 0; i < offset.size(); ++i)
    AssertThrow(offset[i] == offset_ref[i], ExcInternalError());
  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  // Initialize Kokkos
  Kokkos::initialize(argc, argv);

  initlog();

  // tests
  test_1d();
  test_2d();
  test_3d();

  Kokkos::finalize();
}
