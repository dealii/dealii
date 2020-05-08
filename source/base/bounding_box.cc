// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/geometry_info.h>

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
      this->boundary_points.first[i] =
        std::min(this->boundary_points.first[i],
                 other_bbox.boundary_points.first[i]);
      this->boundary_points.second[i] =
        std::max(this->boundary_points.second[i],
                 other_bbox.boundary_points.second[i]);
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
std::pair<Point<spacedim, Number>, Point<spacedim, Number>> &
BoundingBox<spacedim, Number>::get_boundary_points()
{
  return this->boundary_points;
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



template <int spacedim, typename Number>
Number
BoundingBox<spacedim, Number>::lower_bound(const unsigned int direction) const
{
  AssertIndexRange(direction, spacedim);

  return boundary_points.first[direction];
}



template <int spacedim, typename Number>
Number
BoundingBox<spacedim, Number>::upper_bound(const unsigned int direction) const
{
  AssertIndexRange(direction, spacedim);

  return boundary_points.second[direction];
}



template <int spacedim, typename Number>
Point<spacedim, Number>
BoundingBox<spacedim, Number>::center() const
{
  Point<spacedim, Number> point;
  for (unsigned int i = 0; i < spacedim; ++i)
    point[i] = .5 * (boundary_points.first[i] + boundary_points.second[i]);

  return point;
}



template <int spacedim, typename Number>
BoundingBox<1, Number>
BoundingBox<spacedim, Number>::bounds(const unsigned int direction) const
{
  AssertIndexRange(direction, spacedim);

  std::pair<Point<1, Number>, Point<1, Number>> lower_upper_bounds;
  lower_upper_bounds.first[0]  = lower_bound(direction);
  lower_upper_bounds.second[0] = upper_bound(direction);

  return BoundingBox<1, Number>(lower_upper_bounds);
}



template <int spacedim, typename Number>
Number
BoundingBox<spacedim, Number>::side_length(const unsigned int direction) const
{
  AssertIndexRange(direction, spacedim);

  return boundary_points.second[direction] - boundary_points.first[direction];
}



template <int spacedim, typename Number>
Point<spacedim, Number>
BoundingBox<spacedim, Number>::vertex(const unsigned int index) const
{
  AssertIndexRange(index, GeometryInfo<spacedim>::vertices_per_cell);

  const Point<spacedim> unit_cell_vertex =
    GeometryInfo<spacedim>::unit_cell_vertex(index);

  Point<spacedim, Number> point;
  for (unsigned int i = 0; i < spacedim; ++i)
    point[i] = boundary_points.first[i] + side_length(i) * unit_cell_vertex[i];

  return point;
}



template <int spacedim, typename Number>
BoundingBox<spacedim, Number>
BoundingBox<spacedim, Number>::child(const unsigned int index) const
{
  AssertIndexRange(index, GeometryInfo<spacedim>::max_children_per_cell);

  // Vertex closest to child.
  const Point<spacedim, Number> parent_vertex = vertex(index);
  const Point<spacedim, Number> parent_center = center();

  const Point<spacedim> upper_corner_unit_cell =
    GeometryInfo<spacedim>::unit_cell_vertex(
      GeometryInfo<spacedim>::vertices_per_cell - 1);

  const Point<spacedim> lower_corner_unit_cell =
    GeometryInfo<spacedim>::unit_cell_vertex(0);

  std::pair<Point<spacedim, Number>, Point<spacedim, Number>>
    child_lower_upper_corner;
  for (unsigned int i = 0; i < spacedim; ++i)
    {
      const double child_side_length = side_length(i) / 2;

      const double child_center = (parent_center[i] + parent_vertex[i]) / 2;

      child_lower_upper_corner.first[i] =
        child_center + child_side_length * (lower_corner_unit_cell[i] - .5);
      child_lower_upper_corner.second[i] =
        child_center + child_side_length * (upper_corner_unit_cell[i] - .5);
    }

  return BoundingBox<spacedim, Number>(child_lower_upper_corner);
}



template <int spacedim, typename Number>
BoundingBox<spacedim - 1, Number>
BoundingBox<spacedim, Number>::cross_section(const unsigned int direction) const
{
  AssertIndexRange(direction, spacedim);

  std::pair<Point<spacedim - 1, Number>, Point<spacedim - 1, Number>>
    cross_section_lower_upper_corner;
  for (unsigned int d = 0; d < spacedim - 1; ++d)
    {
      const int index_to_write_from =
        internal::coordinate_to_one_dim_higher<spacedim - 1>(direction, d);

      cross_section_lower_upper_corner.first[d] =
        boundary_points.first[index_to_write_from];

      cross_section_lower_upper_corner.second[d] =
        boundary_points.second[index_to_write_from];
    }

  return BoundingBox<spacedim - 1, Number>(cross_section_lower_upper_corner);
}



template <int dim, typename Number>
BoundingBox<dim, Number>
create_unit_bounding_box()
{
  std::pair<Point<dim, Number>, Point<dim, Number>> lower_upper_corner;
  for (unsigned int i = 0; i < dim; ++i)
    {
      lower_upper_corner.second[i] = 1;
    }
  return BoundingBox<dim, Number>(lower_upper_corner);
}


#include "bounding_box.inst"
DEAL_II_NAMESPACE_CLOSE
