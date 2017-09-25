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
bool BoundingBox<spacedim,Number>::point_inside (const Point<spacedim, Number> &p) const
{
  for (unsigned int i=0; i < spacedim; ++i)
    {
      //Bottom left-top right convention: the point is outside if it's smaller than the
      //first or bigger than the second boundary point
      //The bounding box is defined as a closed set
      if ( p[i] < this->boundary_points.first[i] || this->boundary_points.second[i] < p[i])
        return false;
    }
  return true;
}

template <int spacedim, typename Number>
void BoundingBox<spacedim,Number>::merge_with(const BoundingBox<spacedim,Number> &other_bbox)
{
  for (unsigned int i=0; i < spacedim; ++i)
    {
      this->boundary_points.first[i] = std::min( this->boundary_points.first[i], other_bbox.boundary_points.first[i] );
      this->boundary_points.second[i] = std::max( this->boundary_points.second[i], other_bbox.boundary_points.second[i] );
    }
}

template <int spacedim, typename Number>
const std::pair<Point<spacedim,Number>,Point<spacedim,Number>> &BoundingBox<spacedim,Number>::get_boundary_points () const
{
  return this->boundary_points;
}

template <int spacedim, typename Number>
double BoundingBox<spacedim,Number>::volume() const
{
  double vol = 1.0;
  for (unsigned int i=0; i < spacedim; ++i)
    vol *= ( this->boundary_points.second[i] - this->boundary_points.first[i] );
  return vol;
}

#include "bounding_box.inst"
DEAL_II_NAMESPACE_CLOSE
