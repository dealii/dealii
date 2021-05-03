// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


// Check ArborX wrapper: nearest bounding boxes


#include <deal.II/arborx/access_traits.h>
#include <deal.II/arborx/bvh.h>

#include <deal.II/base/bounding_box.h>

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

  Point<2> point_a(-1.0, -1.0);
  Point<2> point_b(-0.5, -0.5);
  Point<2> point_c(5.2, 5.2);
  Point<2> point_d(5.6, 5.6);

  std::vector<BoundingBox<2>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

  ArborXWrappers::BVH                         bvh(bounding_boxes);
  ArborXWrappers::BoundingBoxNearestPredicate bb_nearest(query_bounding_boxes,
                                                         1);
  auto             indices_offset = bvh.query(bb_nearest);
  std::vector<int> indices        = indices_offset.first;
  std::vector<int> offset         = indices_offset.second;

  std::vector<int> indices_ref = {0, 15};
  std::vector<int> offset_ref  = {0, 1, 2};

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
test_point_2d()
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

  Point<3> point_a(0.5, 0.5, 0.5);
  Point<3> point_b(1.5, 1.5, 1.5);
  Point<3> point_c(2.2, 2.2, 2.2);
  Point<3> point_d(2.6, 2.6, 2.6);

  std::vector<BoundingBox<3>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

  ArborXWrappers::BVH                         bvh(points);
  ArborXWrappers::BoundingBoxNearestPredicate bb_nearest(query_bounding_boxes,
                                                         1);
  auto             indices_offset = bvh.query(bb_nearest);
  std::vector<int> indices        = indices_offset.first;
  std::vector<int> offset         = indices_offset.second;

  std::vector<int> indices_ref = {6, 12};
  std::vector<int> offset_ref  = {0, 1, 2};

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

  Point<3> point_a(-1.0, -1.0, -1.0);
  Point<3> point_b(-0.5, -0.5, -0.5);
  Point<3> point_c(10.2, 10.2, 10.2);
  Point<3> point_d(10.6, 10.6, 10.6);

  std::vector<BoundingBox<3>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

  ArborXWrappers::BVH                         bvh(bounding_boxes);
  ArborXWrappers::BoundingBoxNearestPredicate bb_nearest(query_bounding_boxes,
                                                         1);
  auto             indices_offset = bvh.query(bb_nearest);
  std::vector<int> indices        = indices_offset.first;
  std::vector<int> offset         = indices_offset.second;

  std::vector<int> indices_ref = {0, 63};
  std::vector<int> offset_ref  = {0, 1, 2};

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
test_point_3d()
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

  Point<3> point_a(0.5, 0.5, 0.5);
  Point<3> point_b(1.5, 1.5, 1.5);
  Point<3> point_c(2.2, 2.2, 2.2);
  Point<3> point_d(2.6, 2.6, 2.6);

  std::vector<BoundingBox<3>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

  ArborXWrappers::BVH                         bvh(points);
  ArborXWrappers::BoundingBoxNearestPredicate bb_nearest(query_bounding_boxes,
                                                         1);
  auto             indices_offset = bvh.query(bb_nearest);
  std::vector<int> indices        = indices_offset.first;
  std::vector<int> offset         = indices_offset.second;

  std::vector<int> indices_ref = {31, 62};
  std::vector<int> offset_ref  = {0, 1, 2};

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
  test_bounding_box_2d();
  test_bounding_box_3d();
  test_point_2d();
  test_point_3d();

  Kokkos::finalize();
}
