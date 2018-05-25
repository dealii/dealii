// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------
#include <deal.II/base/bounding_box.h>

DEAL_II_NAMESPACE_OPEN

template <int spacedim, typename Number>
bool
BoundingBox<spacedim, Number>::point_inside(
  const Point<spacedim, Number> &p) const
{
  for (unsigned int i = 0; i < spacedim; ++i)
    {
      // Bottom left-top right convention: the point is outside if it's smaller
      // than the first or bigger than the second boundary point The bounding
      // box is defined as a closed set
      if (std::numeric_limits<Number>::epsilon() *
              (std::abs(this->boundary_points.first[i] + p[i])) <
            this->boundary_points.first[i] - p[i] ||
          std::numeric_limits<Number>::epsilon() *
              (std::abs(this->boundary_points.second[i] + p[i])) <
            p[i] - this->boundary_points.second[i])
        return false;
    }
  return true;
}

template <int spacedim, typename Number>
void
BoundingBox<spacedim, Number>::merge_with(
  const BoundingBox<spacedim, Number> &other_bbox)
{
  for (unsigned int i = 0; i < spacedim; ++i)
    {
      this->boundary_points.first[i] = std::min(
        this->boundary_points.first[i], other_bbox.boundary_points.first[i]);
      this->boundary_points.second[i] = std::max(
        this->boundary_points.second[i], other_bbox.boundary_points.second[i]);
    }
}

template <int spacedim, typename Number>
NeighborType
BoundingBox<spacedim, Number>::get_neighbor_type(
  const BoundingBox<spacedim, Number> &other_bbox) const
{
  if (spacedim == 1)
    {
      // In dimension 1 if the two bounding box are neighbors
      // we can merge them
      if (this->point_inside(other_bbox.boundary_points.first) ||
          this->point_inside(other_bbox.boundary_points.second))
        return NeighborType::mergeable_neighbors;
      return NeighborType::not_neighbors;
    }
  else
    {
      std::vector<Point<spacedim, Number>> bbox1;
      bbox1.push_back(this->get_boundary_points().first);
      bbox1.push_back(this->get_boundary_points().second);
      std::vector<Point<spacedim, Number>> bbox2;
      bbox2.push_back(other_bbox.get_boundary_points().first);
      bbox2.push_back(other_bbox.get_boundary_points().second);

      // Step 1: testing if the boxes are close enough to intersect
      for (unsigned int d = 0; d < spacedim; ++d)
        if (bbox1[0][d] * (1 - std::numeric_limits<Number>::epsilon()) >
              bbox2[1][d] ||
            bbox2[0][d] * (1 - std::numeric_limits<Number>::epsilon()) >
              bbox1[1][d])
          return NeighborType::not_neighbors;

      // The boxes intersect: we need to understand now how they intersect.
      // We begin by computing the intersection:
      std::vector<double> intersect_bbox_min;
      std::vector<double> intersect_bbox_max;
      for (unsigned int d = 0; d < spacedim; ++d)
        {
          intersect_bbox_min.push_back(std::max(bbox1[0][d], bbox2[0][d]));
          intersect_bbox_max.push_back(std::min(bbox1[1][d], bbox2[1][d]));
        }

      // Finding the intersection's dimension

      unsigned int intersect_dim = spacedim;
      for (unsigned int d = 0; d < spacedim; ++d)
        if (std::abs(intersect_bbox_min[d] - intersect_bbox_max[d]) <=
            std::numeric_limits<Number>::epsilon() *
              (std::abs(intersect_bbox_min[d]) +
               std::abs(intersect_bbox_max[d])))
          --intersect_dim;

      if (intersect_dim == 0 || intersect_dim + 2 == spacedim)
        return NeighborType::simple_neighbors;

      // Checking the two mergeable cases: first if the boxes are aligned so
      // that they can be merged
      unsigned int not_align_1 = 0, not_align_2 = 0;
      bool         same_direction = true;
      for (unsigned int d = 0; d < spacedim; ++d)
        {
          if (std::abs(bbox2[0][d] - bbox1[0][d]) >
              std::numeric_limits<double>::epsilon() *
                (std::abs(bbox2[0][d]) + std::abs(bbox1[0][d])))
            ++not_align_1;
          if (std::abs(bbox1[1][d] - bbox2[1][d]) >
              std::numeric_limits<double>::epsilon() *
                (std::abs(bbox1[1][d]) + std::abs(bbox2[1][d])))
            ++not_align_2;
          if (not_align_1 != not_align_2)
            {
              same_direction = false;
              break;
            }
        }

      if (not_align_1 <= 1 && not_align_2 <= 1 && same_direction)
        return NeighborType::mergeable_neighbors;

      // Second: one box is contained/equal to the other
      if ((this->point_inside(bbox2[0]) && this->point_inside(bbox2[1])) ||
          (other_bbox.point_inside(bbox1[0]) &&
           other_bbox.point_inside(bbox1[1])))
        return NeighborType::mergeable_neighbors;

      // Degenerate and mergeable cases have been found, it remains:
      return NeighborType::attached_neighbors;
    }
}

template <int spacedim, typename Number>
const std::pair<Point<spacedim, Number>, Point<spacedim, Number>> &
BoundingBox<spacedim, Number>::get_boundary_points() const
{
  return this->boundary_points;
}

template <int spacedim, typename Number>
double
BoundingBox<spacedim, Number>::volume() const
{
  double vol = 1.0;
  for (unsigned int i = 0; i < spacedim; ++i)
    vol *= (this->boundary_points.second[i] - this->boundary_points.first[i]);
  return vol;
}

#include "bounding_box.inst"
DEAL_II_NAMESPACE_CLOSE
