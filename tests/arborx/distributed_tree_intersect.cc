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


// Check ArborX wrapper: intersection of bounding boxes


#include <deal.II/arborx/access_traits.h>
#include <deal.II/arborx/distributed_tree.h>

#include <deal.II/base/bounding_box.h>

#include "../tests.h"


void
test_1d()
{
  std::vector<BoundingBox<1>> bounding_boxes;
  std::vector<Point<1>>       points;

  const unsigned int n_points_1d = 5;
  const int          rank    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const int          n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int          offset_bb    = 2 * rank * n_points_1d;
  const int          offset_query = 2 * ((rank + 1) % n_procs) * n_points_1d;
  for (unsigned int i = 0; i < n_points_1d; ++i)
    points.emplace_back(i + offset_bb);

  for (unsigned int i = 0; i < n_points_1d - 1; ++i)
    bounding_boxes.emplace_back(std::make_pair(points[i], points[i + 1]));

  Point<1> point_a(0.5 + offset_query);
  Point<1> point_b(1.5 + offset_query);
  Point<1> point_c(2.2 + offset_query);
  Point<1> point_d(2.6 + offset_query);

  std::vector<BoundingBox<1>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

  ArborXWrappers::DistributedTree               distributed_tree(MPI_COMM_WORLD,
                                                   bounding_boxes);
  ArborXWrappers::BoundingBoxIntersectPredicate bb_intersect(
    query_bounding_boxes);
  auto indices_ranks_offset = distributed_tree.query(bb_intersect);
  auto indices_ranks        = indices_ranks_offset.first;
  auto offsets              = indices_ranks_offset.second;

  std::vector<int> indices_ref = {1, 0, 2};
  std::vector<int> offsets_ref = {0, 2, 3};

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
test_2d()
{
  std::vector<BoundingBox<2>> bounding_boxes;
  std::vector<Point<2>>       points;

  const unsigned int n_points_1d = 5;
  const int          rank    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const int          n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int          offset_bb    = 2 * rank * n_points_1d;
  const int          offset_query = 2 * ((rank + 1) % n_procs) * n_points_1d;
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

  Point<2> point_a(0.5 + offset_query, 0.5 + offset_query);
  Point<2> point_b(1.5 + offset_query, 1.5 + offset_query);
  Point<2> point_c(2.2 + offset_query, 2.2 + offset_query);
  Point<2> point_d(2.6 + offset_query, 2.6 + offset_query);

  std::vector<BoundingBox<2>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

  ArborXWrappers::DistributedTree               distributed_tree(MPI_COMM_WORLD,
                                                   bounding_boxes);
  ArborXWrappers::BoundingBoxIntersectPredicate bb_intersect(
    query_bounding_boxes);
  auto indices_ranks_offset = distributed_tree.query(bb_intersect);
  auto indices_ranks        = indices_ranks_offset.first;
  auto offsets              = indices_ranks_offset.second;

  std::vector<int> indices_ref = {0, 1, 4, 5, 10};
  std::vector<int> offsets_ref = {0, 4, 5};

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
test_3d()
{
  std::vector<BoundingBox<3>> bounding_boxes;
  std::vector<Point<3>>       points;

  unsigned int n_points_1d  = 5;
  const int    rank         = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const int    n_procs      = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int    offset_bb    = 2 * rank * n_points_1d;
  const int    offset_query = 2 * ((rank + 1) % n_procs) * n_points_1d;
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

  Point<3> point_a(0.5 + offset_query, 0.5 + offset_query, 0.5 + offset_query);
  Point<3> point_b(1.5 + offset_query, 1.5 + offset_query, 1.5 + offset_query);
  Point<3> point_c(2.2 + offset_query, 2.2 + offset_query, 2.2 + offset_query);
  Point<3> point_d(2.6 + offset_query, 2.6 + offset_query, 2.6 + offset_query);

  std::vector<BoundingBox<3>> query_bounding_boxes;
  query_bounding_boxes.emplace_back(std::make_pair(point_a, point_b));
  query_bounding_boxes.emplace_back(std::make_pair(point_c, point_d));

  ArborXWrappers::DistributedTree               distributed_tree(MPI_COMM_WORLD,
                                                   bounding_boxes);
  ArborXWrappers::BoundingBoxIntersectPredicate bb_intersect(
    query_bounding_boxes);
  auto indices_ranks_offset = distributed_tree.query(bb_intersect);
  auto indices_ranks        = indices_ranks_offset.first;
  auto offsets              = indices_ranks_offset.second;

  std::vector<int> indices_ref = {0, 1, 4, 5, 16, 17, 20, 21, 42};
  std::vector<int> offsets_ref = {0, 8, 9};

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
  test_1d();
  test_2d();
  test_3d();
}
