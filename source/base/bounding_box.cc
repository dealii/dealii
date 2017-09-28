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
bool BoundingBox<spacedim,Number>::point_inside (const Point<spacedim, Number> &p, const double &err) const
{
  for (unsigned int i=0; i < spacedim; ++i)
    {
      //Bottom left-top right convention: the point is outside if it's smaller than the
      //first or bigger than the second boundary point
      //The bounding box is defined as a closed set
      if ( p[i] + err < this->boundary_points.first[i] || this->boundary_points.second[i] + err < p[i])
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
unsigned int BoundingBox < spacedim,Number >::is_neighbour(const BoundingBox<spacedim,Number> &other_bbox, const double &err) const
{
  if (spacedim == 1)
    {
      // In dimension 1 if the two bounding box are neighbours
      // we can merge them
      if ( this->point_inside(other_bbox.boundary_points.first,err)
           || this->point_inside(other_bbox.boundary_points.second,err) )
        return 2;
      return 0;
    }
  else
    {
      std::vector<bool > boundary_relations;

      // Check if the boundary points of one bounding box are owned by another
      boundary_relations.push_back(other_bbox.point_inside(this->boundary_points.first,err));
      boundary_relations.push_back(other_bbox.point_inside(this->boundary_points.second,err));

      boundary_relations.push_back(this->point_inside(other_bbox.boundary_points.first,err));
      boundary_relations.push_back(this->point_inside(other_bbox.boundary_points.second,err));

      std::vector<bool> all_false(4,false);

      if ( (boundary_relations[0] && boundary_relations[1]) || (boundary_relations[2] && boundary_relations[3]))
        // This happens if a bounding box is inside the other
        return 2;
      else if ( boundary_relations == all_false)
        return 0;

      //Translating the first point in 0
      Tensor< 1 ,spacedim , Number > boundary1 = this->boundary_points.second - this->boundary_points.first;
      std::vector< Tensor< 1 ,spacedim , Number > > boundary2 = {other_bbox.boundary_points.first-this->boundary_points.first,
                                                                 other_bbox.boundary_points.second-this->boundary_points.first
                                                                };

      // These variables keep how many of the projection of the various axis
      // have the following properties:
      // on_vertex_1 = combination of the vertices of bounding box1
      // inside_1 = inside point 1
      // this is done for the two boundary points boundary2[0,1]
      unsigned int on_vertex_1 = 0;
      unsigned int inside_1 = 0;
      for (unsigned int d=0; d<spacedim; ++d)
        {
          if ( std::abs( boundary2[0][d] ) < err ||
               std::abs( boundary2[0][d] - boundary1[d]) < err )
            ++on_vertex_1;
          else if ( -err < boundary2[0][d] && boundary2[0][d] < boundary1[d] + err )
            ++inside_1;
        }

      unsigned int on_vertex_2 = 0;
      unsigned int inside_2 = 0;
      for (unsigned int d=0; d<spacedim; ++d)
        {
          if ( std::abs( boundary2[1][d] ) < err || std::abs(boundary2[1][d] - boundary1[d] ) < err )
            ++on_vertex_2;
          else if ( -err < boundary2[1][d] && boundary2[1][d] < boundary1[d] + err )
            ++inside_2;
        }

      if ( (on_vertex_1 == spacedim && on_vertex_2 == spacedim -1 )
           || (on_vertex_2 == spacedim && on_vertex_1 == spacedim -1 ) )
        return 2;
      else if ( (on_vertex_1 == spacedim && (on_vertex_2 + inside_2 == 0))
                || (on_vertex_2 == spacedim && (on_vertex_1 + inside_1 == 0)) )
        //This happens if the two boxes share only a vertex
        return 0;
      else if ( spacedim == 3 )
        if ( (on_vertex_1 == spacedim - 1 && inside_1 == 1 && (on_vertex_2 + inside_2 == 0) )
             || (on_vertex_2 == spacedim - 1 && inside_2 == 1 && (on_vertex_1 + inside_1 == 0) ) )
          // This happens, in dimension 3, if the two boxes share only an edge
          return 0;
      return 1;
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
