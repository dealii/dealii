/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file DomainUtil.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "DomainUtil.hpp"
#include "MeshInterface.hpp"
#include "MsqVertex.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS { 
namespace DomainUtil { 

void bounding_box( const MsqVertex* coords,
                   size_t num_coords,
                   Vector3D& min,
                   Vector3D& max )
{
  min = max = coords[0];
  for (size_t i = 1; i < num_coords; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (coords[i][j] < min[j])
        min[j] = coords[i][j];
      if (coords[i][j] > max[j])
        max[j] = coords[i][j];
    }
  }
}

double max_box_extent( const MsqVertex* vertex_array, size_t num_vertices )
{
  Vector3D min, max;
  bounding_box( vertex_array, num_vertices, min, max );
  max -= min;
  return (max[0] >= max[1] && max[0] >= max[2]) ? max[0] : 
         (max[1] >= max[2]) ? max[1] : max[2];
}

void get_fixed_vertices( Mesh* mesh,
                         const Mesh::VertexHandle* verts,
                         size_t num_verts,
                         std::vector<Mesh::VertexHandle>& fixed_verts,
                         MsqError& err )
{
  std::vector<bool> fixed( num_verts );
  mesh->vertices_get_fixed_flag( verts, fixed, num_verts, err );
  MSQ_ERRRTN(err);
  for (size_t i = 0; i < num_verts; ++i) 
    if (fixed[i])
      fixed_verts.push_back(verts[i]);
}

bool non_colinear_vertices( const MsqVertex* verts,
                            size_t num_verts,
                            Vector3D coords_out[3],
                            double epsilon )
{
    // This function will attempt to find trhee non-colinear 
    // vertices from the input list.  Further, it will attempt 
    // to select three such vertices that are relatively far
    // apart so as to minimize rounding error in any calculation
    // using the results of this function.

    // Non-colinear, by definition, must be at least trhee unique points.
  if (num_verts < 3)
    return false;

    // Begin with the first input vertex
  size_t first_idx = 0;

    // Choose the vertex furthest from the initial one
  size_t second_idx = 1;
  double dist_sqr = (verts[first_idx] - verts[second_idx]).length_squared();
  for (size_t i = 2; i < num_verts; ++i) {
    double ds = (verts[second_idx] - verts[i]).length_squared();
    if (ds > dist_sqr) {
      dist_sqr = ds;
      second_idx = i;
    }
  }
    // fail if all vertices are coincident
  if (dist_sqr <= epsilon*epsilon) 
    return false;
  
    // re-select the first vertex as the one furthest from the second
  for (size_t i = 1; i < num_verts; ++i) {
    double ds = (verts[second_idx] - verts[i]).length_squared();
    if (ds > dist_sqr) {
      dist_sqr = ds;
      first_idx = i;
    }
  }
  
    // select the third vertex as the one furthest from the line formed
    // by the first two vertices
  Vector3D b = verts[first_idx];
  Vector3D m = verts[second_idx] - b;
  Vector3D mx = m * (1.0 / m.length_squared());
  dist_sqr = -1.0;
  size_t third_idx = 0;
  for (size_t i = 0; i < num_verts; ++i) {
    double t = mx % (verts[i] - b);
    double ds = ((b + t*m) - verts[i]).length_squared();
    if (ds > dist_sqr) {
      third_idx = i;
      dist_sqr = ds;
    }
  }
    // fail if all vertices are colinear
  if (dist_sqr <= epsilon*epsilon) 
    return false;
  
  coords_out[0] = verts[first_idx];
  coords_out[1] = verts[second_idx];
  coords_out[2] = verts[third_idx];
  return true;
}
    

bool non_coplanar_vertices( const MsqVertex* verts,
                            size_t num_verts,
                            Vector3D coords_out[4],
                            double epsilon )
{
    // This function will attempt to find four non-coplanar 
    // vertices from the input list.  Further, it will attempt 
    // to select four such vertices that are relatively far
    // apart so as to minimize rounding error in any calculation
    // using the results of this function.

    // Non-coplanar, by definition, must be at least four unique points.
  if (num_verts < 4)
    return false;

    // Get three non-colinear vertices
  if (!non_colinear_vertices( verts, num_verts, coords_out, epsilon ))
    return false;

    // The plane of the first three vertices:
  Vector3D norm = (coords_out[1] - coords_out[0]) * (coords_out[2] - coords_out[0]);
  norm /= norm.length();
  double d = -(norm % coords_out[0]);
  
    // Search for the fourth vertex that is furthest from the plane
    // of the first three
  double dist = -1.0;
  for (size_t i = 0; i < num_verts; ++i) {
    double disti = fabs( norm % verts[i] + d );
    if (disti > dist) {
      dist = disti;
      coords_out[3] = verts[i];
    }
  }

    // fail if all vertices are colinear
  return (dist > epsilon);
}

} // namespace DomainUtil  
} // namespace MESQUITE_NS
