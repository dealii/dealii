/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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
 
    (2007) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file AffineMapMetric.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "AffineMapMetric.hpp"
#include "MsqMatrix.hpp"
#include "ElementQM.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include "MappingFunction.hpp"
#include "WeightCalculator.hpp"
#include "TargetCalculator.hpp"
#include "TMetric.hpp"
#include "TargetMetricUtil.hpp"

#include <functional>
#include <algorithm>

namespace MESQUITE_NS {

const double TRI_XFORM_VALS[] = { 1.0, -1.0/sqrt(3.0), 0.0, 2.0/sqrt(3.0) };
MsqMatrix<2,2> TRI_XFORM( TRI_XFORM_VALS );

const double TET_XFORM_VALS[] = { 1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0),
                                  0.0,  2.0/sqrt(3.0), -1.0/sqrt(6.0),
                                  0.0,  0.0,            sqrt(3.0/2.0) };
MsqMatrix<3,3> TET_XFORM( TET_XFORM_VALS );
 
AffineMapMetric::AffineMapMetric( TargetCalculator* tc,
                                  WeightCalculator* wc,
                                  TMetric* target_metric ) 
  : targetCalc(tc),
    weightCalc(wc),
    targetMetric( target_metric )
{ }
 
AffineMapMetric::AffineMapMetric( TargetCalculator* tc,
                                  TMetric* target_metric ) 
  : targetCalc(tc),
    weightCalc(0),
    targetMetric( target_metric )
{ }
     

int AffineMapMetric::get_negate_flag( ) const { return 1; }

std::string AffineMapMetric::get_name() const
  { return std::string("AffineMap(") + targetMetric->get_name() + ')'; }

void AffineMapMetric::get_evaluations( PatchData& pd,
                                       std::vector<size_t>& handles,
                                       bool free,
                                       MsqError& err )
{
  get_sample_pt_evaluations( pd, handles, free, err );
}

void AffineMapMetric::get_element_evaluations( PatchData& pd,
                                               size_t elem,
                                               std::vector<size_t>& handles,
                                               MsqError& err )
{
  get_elem_sample_points( pd, elem, handles, err );
}

bool AffineMapMetric::evaluate( PatchData& pd, size_t handle, double& value, MsqError& err )
{
  Sample s = ElemSampleQM::sample( handle );
  size_t e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  EntityTopology type = elem.get_element_type();
  unsigned edim = TopologyInfo::dimension( type );
  const size_t* conn = elem.get_vertex_index_array();
  
    // This metric only supports sampling at corners, except for simplices.
    // If element is a simpex, then the Jacobian is constant over a linear 
    // element.  In this case, always evaluate at any vertex.
  //unsigned corner = s.number;
  if (s.dimension != 0) {
    if (type == TRIANGLE || type == TETRAHEDRON)
      /*corner = 0*/;
    else {
      MSQ_SETERR(err)("Invalid sample point for AffineMapMetric", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
  }
  
  bool rval;
  if (edim == 3) { // 3x3 or 3x2 targets ?
    Vector3D c[3] = { Vector3D(0,0,0), Vector3D(0,0,0), Vector3D(0,0,0) };
    unsigned n;
    const unsigned* adj = TopologyInfo::adjacent_vertices( type, s.number, n );
    c[0] = pd.vertex_by_index( conn[adj[0]] ) - pd.vertex_by_index( conn[s.number] );
    c[1] = pd.vertex_by_index( conn[adj[1]] ) - pd.vertex_by_index( conn[s.number] );
    c[2] = pd.vertex_by_index( conn[adj[2]] ) - pd.vertex_by_index( conn[s.number] );
    MsqMatrix<3,3> A;
    A.set_column( 0, MsqMatrix<3,1>(c[0].to_array()) );
    A.set_column( 1, MsqMatrix<3,1>(c[1].to_array()) );
    A.set_column( 2, MsqMatrix<3,1>(c[2].to_array()) );
    if (type == TETRAHEDRON)
      A = A * TET_XFORM;

    MsqMatrix<3,3> W;
    targetCalc->get_3D_target( pd, e, s, W, err ); MSQ_ERRZERO(err);
    rval = targetMetric->evaluate( A * inverse(W), value, err ); MSQ_ERRZERO(err);
  }
  else {
    Vector3D c[2] = { Vector3D(0,0,0), Vector3D(0,0,0) };
    unsigned n;
    const unsigned* adj = TopologyInfo::adjacent_vertices( type, s.number, n );
    c[0] = pd.vertex_by_index( conn[adj[0]] ) - pd.vertex_by_index( conn[s.number] );
    c[1] = pd.vertex_by_index( conn[adj[1]] ) - pd.vertex_by_index( conn[s.number] );
    MsqMatrix<3,2> App;
    App.set_column( 0, MsqMatrix<3,1>(c[0].to_array()) );
    App.set_column( 1, MsqMatrix<3,1>(c[1].to_array()) );
    
    MsqMatrix<3,2> Wp;
    targetCalc->get_surface_target( pd, e, s, Wp, err ); MSQ_ERRZERO(err);

    MsqMatrix<2,2> A, W;
    MsqMatrix<3,2> RZ;
    surface_to_2d( App, Wp, W, RZ );
    A = transpose(RZ) * App;
    if (type == TRIANGLE)
      A = A * TRI_XFORM;
    
    rval = targetMetric->evaluate( A*inverse(W), value, err ); MSQ_ERRZERO(err);
  }
  
    // apply target weight to value
  if (weightCalc) {
    double ck = weightCalc->get_weight( pd, e, s, err ); MSQ_ERRZERO(err);
    value *= ck;
  }
  return rval;
}

bool AffineMapMetric::evaluate_with_indices( PatchData& pd,
                                             size_t handle,
                                             double& value,
                                             std::vector<size_t>& indices,
                                             MsqError& err )
{
  Sample   s = ElemSampleQM::sample( handle );
  size_t   e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  EntityTopology type = elem.get_element_type();
  const size_t* conn = elem.get_vertex_index_array();
  
    // this metric only supports sampling at corners
  if (s.dimension != 0) {
    if (type != TRIANGLE && type != TETRAHEDRON) {
      MSQ_SETERR(err)("Invalid sample point for AffineMapMetric", MsqError::UNSUPPORTED_ELEMENT );
      return false;
    }
    s.dimension = 0;
    s.number = 0;
  }

  unsigned n;
  const unsigned* adj = TopologyInfo::adjacent_vertices( type, s.number, n );
  indices.clear();
  if (conn[s.number] < pd.num_free_vertices())
    indices.push_back(conn[s.number]);
  for (unsigned i = 0; i < n; ++i)
    if (conn[adj[i]] < pd.num_free_vertices())
      indices.push_back(conn[adj[i]]);
  
  return evaluate( pd, handle, value, err );
}

} // namespace Mesquite
