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


// Check ArborX wrapper: nearest bounding boxes


#include <deal.II/arborx/access_traits.h>
#include <deal.II/arborx/distributed_tree.h>

#include <deal.II/base/bounding_box.h>

#include "../tests.h"


void
test_bounding_box_2d()
{
  std::vector<BoundingBox<2>> bounding_boxes;
  std::vector<Point<2>>       points;

  unsigned int n_points_1d  = 5;
  const int    rank         = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const int    n_procs      = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int    offset_bb    = 2 * rank * n_points_1d;
  const int    offset_query = 2 * ((rank + 1) % n_procs) * n_points_1d;
  for (unsigned int i = 0; i < n_points_1d; ++i)
    {
      for (unsigned int j = 0; j < n_points_1d; ++j)
        {
          points.emplace_back(j + offset_bb, i + offset_bb);
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

  Point<2> point_a(-1.0 + offset_query, -1.0 + offset_query);
  Point<2> point_b(-0.5 + offset_query, -0.5 + offset_query);
  Point<2> point_c(5.2 + offset_query, 5.2 + offset_query);
  Point<2> point_d(5.6 + offset_query, 5.6 + offset_query);

  std::vector<BoundingBox<2>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

  ArborXWrappers::DistributedTree             distributed_tree(MPI_COMM_WORLD,
                                                   bounding_boxes);
  ArborXWrappers::BoundingBoxNearestPredicate bb_nearest(query_bounding_boxes,
                                                         1);
  auto indices_ranks_offsets = distributed_tree.query(bb_nearest);
  auto indices_ranks         = indices_ranks_offsets.first;
  auto offsets               = indices_ranks_offsets.second;

  std::vector<int> indices_ref = {0, 15};
  std::vector<int> offsets_ref = {0, 1, 2};

  AssertThrow(indices_ranks.size() == indices_ref.size(), ExcInternalError());
  AssertThrow(offsets.size() == offsets_ref.size(), ExcInternalError());
  for (unsigned int i = 0; i < offsets.size() - 1; ++i)
    {
      for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
          // The indices associated to each query are not ordered.
          bool found = false;
          for (int k = offsets[i]; k < offsets[i + 1]; ++k)
            {
              if ((indices_ranks[j].first == indices_ref[k]) &&
                  (indices_ranks[j].second == (rank + 1) % n_procs))
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
test_point_2d()
{
  std::vector<Point<2>> points;

  unsigned int n_points_1d  = 5;
  const int    rank         = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const int    n_procs      = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int    offset_pt    = 2 * rank * n_points_1d;
  const int    offset_query = 2 * ((rank + 1) % n_procs) * n_points_1d;
  for (unsigned int i = 0; i < n_points_1d; ++i)
    {
      for (unsigned int j = 0; j < n_points_1d; ++j)
        {
          points.emplace_back(j + offset_pt, i + offset_pt);
        }
    }

  Point<2> point_a(0.5 + offset_query, 0.5 + offset_query);
  Point<2> point_b(1.5 + offset_query, 1.5 + offset_query);
  Point<2> point_c(2.2 + offset_query, 2.2 + offset_query);
  Point<2> point_d(2.6 + offset_query, 2.6 + offset_query);

  std::vector<BoundingBox<2>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

  ArborXWrappers::DistributedTree distributed_tree(MPI_COMM_WORLD, points);
  ArborXWrappers::BoundingBoxNearestPredicate bb_nearest(query_bounding_boxes,
                                                         1);
  auto indices_ranks_offsets = distributed_tree.query(bb_nearest);
  auto indices_ranks         = indices_ranks_offsets.first;
  auto offsets               = indices_ranks_offsets.second;

  std::vector<int> indices_ref = {6, 12};
  std::vector<int> offsets_ref = {0, 1, 2};

  AssertThrow(indices_ranks.size() == indices_ref.size(), ExcInternalError());
  AssertThrow(offsets.size() == offsets_ref.size(), ExcInternalError());
  for (unsigned int i = 0; i < offsets.size() - 1; ++i)
    {
      for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
          // The indices associated to each query are not ordered.
          bool found = false;
          for (int k = offsets[i]; k < offsets[i + 1]; ++k)
            {
              if ((indices_ranks[j].first == indices_ref[k]) &&
                  (indices_ranks[j].second == (rank + 1) % n_procs))
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
test_bounding_box_3d()
{
  std::vector<BoundingBox<3>> bounding_boxes;
  std::vector<Point<3>>       points;

  unsigned int n_points_1d  = 5;
  const int    rank         = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const int    n_procs      = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int    offset_bb    = 20 * rank * n_points_1d;
  const int    offset_query = 20 * ((rank + 1) % n_procs) * n_points_1d;
  for (unsigned int i = 0; i < n_points_1d; ++i)
    {
      for (unsigned int j = 0; j < n_points_1d; ++j)
        {
          for (unsigned int k = 0; k < n_points_1d; ++k)
            {
              points.emplace_back(k + offset_bb, j + offset_bb, i + offset_bb);
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

  Point<3> point_a(-1.0 + offset_query,
                   -1.0 + offset_query,
                   -1.0 + offset_query);
  Point<3> point_b(-0.5 + offset_query,
                   -0.5 + offset_query,
                   -0.5 + offset_query);
  Point<3> point_c(10.2 + offset_query,
                   10.2 + offset_query,
                   10.2 + offset_query);
  Point<3> point_d(10.6 + offset_query,
                   10.6 + offset_query,
                   10.6 + offset_query);

  std::vector<BoundingBox<3>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

  ArborXWrappers::DistributedTree             distributed_tree(MPI_COMM_WORLD,
                                                   bounding_boxes);
  ArborXWrappers::BoundingBoxNearestPredicate bb_nearest(query_bounding_boxes,
                                                         1);
  auto indices_ranks_offsets = distributed_tree.query(bb_nearest);
  auto indices_ranks         = indices_ranks_offsets.first;
  auto offsets               = indices_ranks_offsets.second;

  std::vector<int> indices_ref = {0, 63};
  std::vector<int> offsets_ref = {0, 1, 2};

  AssertThrow(indices_ranks.size() == indices_ref.size(), ExcInternalError());
  AssertThrow(offsets.size() == offsets_ref.size(), ExcInternalError());
  for (unsigned int i = 0; i < offsets.size() - 1; ++i)
    {
      for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
          // The indices associated to each query are not ordered.
          bool found = false;
          for (int k = offsets[i]; k < offsets[i + 1]; ++k)
            {
              if ((indices_ranks[j].first == indices_ref[k]) &&
                  (indices_ranks[j].second == (rank + 1) % n_procs))
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
test_point_3d()
{
  std::vector<Point<3>> points;

  unsigned int n_points_1d  = 5;
  const int    rank         = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const int    n_procs      = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int    offset_pt    = 2 * rank * n_points_1d;
  const int    offset_query = 2 * ((rank + 1) % n_procs) * n_points_1d;
  for (unsigned int i = 0; i < n_points_1d; ++i)
    {
      for (unsigned int j = 0; j < n_points_1d; ++j)
        {
          for (unsigned int k = 0; k < n_points_1d; ++k)
            {
              points.emplace_back(k + offset_pt, j + offset_pt, i + offset_pt);
            }
        }
    }

  Point<3> point_a(0.5 + offset_query, 0.5 + offset_query, 0.5 + offset_query);
  Point<3> point_b(1.5 + offset_query, 1.5 + offset_query, 1.5 + offset_query);
  Point<3> point_c(2.2 + offset_query, 2.2 + offset_query, 2.2 + offset_query);
  Point<3> point_d(2.6 + offset_query, 2.6 + offset_query, 2.6 + offset_query);

  std::vector<BoundingBox<3>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

  ArborXWrappers::DistributedTree distributed_tree(MPI_COMM_WORLD, points);
  ArborXWrappers::BoundingBoxNearestPredicate bb_nearest(query_bounding_boxes,
                                                         1);
  auto indices_ranks_offsets = distributed_tree.query(bb_nearest);
  auto indices_ranks         = indices_ranks_offsets.first;
  auto offsets               = indices_ranks_offsets.second;

  std::vector<int> indices_ref = {31, 62};
  std::vector<int> offsets_ref = {0, 1, 2};

  AssertThrow(indices_ranks.size() == indices_ref.size(), ExcInternalError());
  AssertThrow(offsets.size() == offsets_ref.size(), ExcInternalError());
  for (unsigned int i = 0; i < offsets.size() - 1; ++i)
    {
      for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
          // The indices associated to each query are not ordered.
          bool found = false;
          for (int k = offsets[i]; k < offsets[i + 1]; ++k)
            {
              if ((indices_ranks[j].first == indices_ref[k]) &&
                  (indices_ranks[j].second == (rank + 1) % n_procs))
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
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  initlog();

  // tests
  test_bounding_box_2d();
  test_bounding_box_3d();
  test_point_2d();
  test_point_3d();
}
