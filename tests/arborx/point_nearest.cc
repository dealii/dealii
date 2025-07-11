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


// Check ArborX wrapper: nearest points from bounding boxes and others points


#include <deal.II/arborx/access_traits.h>
#include <deal.II/arborx/bvh.h>

#include <deal.II/base/point.h>

#include "../tests.h"


void
test_bounding_box_2d()
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

  std::vector<Point<2>> query_points;
  query_points.emplace_back(0.5, 0.5);
  query_points.emplace_back(1.5, 1.5);
  query_points.emplace_back(2.2, 2.2);
  query_points.emplace_back(2.6, 2.6);


#if ARBORX_VERSION_MAJOR < 2
  ArborXWrappers::BVH bvh(bounding_boxes);
#else
  ArborXWrappers::BVH<BoundingBox<2>> bvh(bounding_boxes);
#endif
  ArborXWrappers::PointNearestPredicate pt_nearest(query_points, 1);
  auto                                  indices_offset = bvh.query(pt_nearest);
  std::vector<int>                      indices        = indices_offset.first;
  std::vector<int>                      offset         = indices_offset.second;

  std::vector<int> indices_ref = {0, 5, 10, 10};
  std::vector<int> offset_ref  = {0, 1, 2, 3, 4};

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
test_points_2d()
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


  std::vector<Point<2>> query_points;
  query_points.emplace_back(0.5, 0.5);
  query_points.emplace_back(1.5, 1.5);
  query_points.emplace_back(2.2, 2.2);
  query_points.emplace_back(2.6, 2.6);


#if ARBORX_VERSION_MAJOR < 2
  ArborXWrappers::BVH bvh(points);
#else
  ArborXWrappers::BVH<Point<2>>       bvh(points);
#endif
  ArborXWrappers::PointNearestPredicate pt_nearest(query_points, 1);
  auto                                  indices_offset = bvh.query(pt_nearest);
  std::vector<int>                      indices        = indices_offset.first;
  std::vector<int>                      offset         = indices_offset.second;

  std::vector<int> indices_ref = {0, 6, 12, 18};
  std::vector<int> offset_ref  = {0, 1, 2, 3, 4};

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
test_bounding_box_3d()
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

  std::vector<Point<3>> query_points;
  query_points.emplace_back(0.5, 0.5, 0.5);
  query_points.emplace_back(1.5, 1.5, 1.5);
  query_points.emplace_back(2.2, 2.2, 2.2);
  query_points.emplace_back(2.6, 2.6, 2.6);


#if ARBORX_VERSION_MAJOR < 2
  ArborXWrappers::BVH bvh(bounding_boxes);
#else
  ArborXWrappers::BVH<BoundingBox<3>> bvh(bounding_boxes);
#endif
  ArborXWrappers::PointNearestPredicate pt_nearest(query_points, 1);
  auto                                  indices_offset = bvh.query(pt_nearest);
  std::vector<int>                      indices        = indices_offset.first;
  std::vector<int>                      offset         = indices_offset.second;


  std::vector<int> indices_ref = {0, 21, 42, 42};
  std::vector<int> offset_ref  = {0, 1, 2, 3, 4};

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
test_points_3d()
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

  std::vector<Point<3>> query_points;
  query_points.emplace_back(0.5, 0.5, 0.5);
  query_points.emplace_back(1.5, 1.5, 1.5);
  query_points.emplace_back(2.2, 2.2, 2.2);
  query_points.emplace_back(2.6, 2.6, 2.6);


  ArborXWrappers::BVH                   bvh(points);
  ArborXWrappers::PointNearestPredicate pt_nearest(query_points, 1);
  auto                                  indices_offset = bvh.query(pt_nearest);
  std::vector<int>                      indices        = indices_offset.first;
  std::vector<int>                      offset         = indices_offset.second;

  std::vector<int> indices_ref = {0, 31, 62, 93};
  std::vector<int> offset_ref  = {0, 1, 2, 3, 4};

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
  initlog();

  // Initialize ArborX
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  // tests
  test_bounding_box_2d();
  test_bounding_box_3d();
  test_points_2d();
  test_points_3d();
}
