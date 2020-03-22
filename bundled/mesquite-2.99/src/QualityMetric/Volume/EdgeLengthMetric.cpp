/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file EdgeLengthMetric.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "EdgeLengthMetric.hpp"
#include "PatchData.hpp"
#include "MsqError.hpp"
#include "Matrix3D.hpp"

namespace MESQUITE_NS {


EdgeLengthMetric::EdgeLengthMetric( double gamma )
  : mGamma( gamma ) {}

std::string EdgeLengthMetric::get_name() const
  { return "EdgeLength"; }

   //! 1 if metric should be minimized, -1 if metric should be maximized.
int EdgeLengthMetric::get_negate_flag() const
  { return 1; }

bool EdgeLengthMetric::evaluate( PatchData& pd, 
                 size_t handle, 
                 double& value, 
                 MsqError& err )
{
  MsqMeshEntity& e = pd.element_by_index( elem(handle) );
  const unsigned* vert_nums;
  vert_nums = TopologyInfo::edge_vertices( e.get_element_type(), edge(handle), err );
  MSQ_ERRZERO(err);
  size_t svi = e.get_vertex_index( vert_nums[0] );
  size_t evi = e.get_vertex_index( vert_nums[1] );
  Vector3D diff = pd.vertex_by_index( svi ) - pd.vertex_by_index( evi );
  double len_sqr = diff % diff - mGamma;
  if (len_sqr <= 0.0) {
    value = 0.0;
    return false; 
  }
  value = sqrt(len_sqr);
  return true;
}

bool EdgeLengthMetric::evaluate_with_gradient( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 MsqError& err )
{
  MsqMeshEntity& e = pd.element_by_index( elem(handle) );
  const unsigned* vert_nums;
  vert_nums = TopologyInfo::edge_vertices( e.get_element_type(), edge(handle), err );
  MSQ_ERRZERO(err);
  size_t svi = e.get_vertex_index( vert_nums[0] );
  size_t evi = e.get_vertex_index( vert_nums[1] );
  Vector3D diff = pd.vertex_by_index( svi ) - pd.vertex_by_index( evi );
  double val_sqr = diff % diff - mGamma;
  if (val_sqr <= 0.0) {
    value = 0.0;
    return false; 
  }
  value = sqrt(val_sqr);
  
  diff *= 1.0/value;
  indices.clear();
  gradient.clear();
  if (svi < pd.num_free_vertices())  {
    indices.push_back( svi );
    gradient.push_back( diff );
  }
  if (evi < pd.num_free_vertices()) {
    indices.push_back( evi ); 
    gradient.push_back( -diff );
  }
  
  return true;
}

} // namespace MESQUITE_NS
