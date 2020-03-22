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


/** \file RefSizeTargetCalculator.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "RefSizeTargetCalculator.hpp"
#include "ReferenceMesh.hpp"
#include "PatchData.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

/* We now scale W_ideal such that det(W_ideal) == 1, so don't scale here 
static void init_scale_factors( double factors[MIXED] )
{
  for (int i = 0; i < MIXED; ++i)
    factors[i] = 0.0;
  
  factors[     TRIANGLE] = 0.21934566882541542;
  factors[QUADRILATERAL] = 0.25;
  factors[  TETRAHEDRON] = 0.08171340981758000;
  factors[      PYRAMID] = 0.06068647146341543;
  factors[        PRISM] = 0.07311522294180514;
  factors[   HEXAHEDRON] = 0.08333333333333333;
}
*/

RefSizeTargetCalculator::RefSizeTargetCalculator( 
                           ReferenceMesh* reference_mesh,
                           TargetCalculator* tc )
   : refMesh( reference_mesh ),
     scaledTargets( tc )
 {} //   { init_scale_factors( scaleFactor ); }
   
RefSizeTargetCalculator::RefSizeTargetCalculator( 
                           ReferenceMesh* reference_mesh  )
   : refMesh( reference_mesh ),
     scaledTargets( &defaultTargets )
 {} //   { init_scale_factors( scaleFactor );  }

double RefSizeTargetCalculator::average_edge_length( PatchData& pd, 
                                                     size_t element,
                                                     MsqError& err )
{  
  Vector3D coords[8];
  MsqMeshEntity& elem = pd.element_by_index( element );
  const size_t* conn = elem.get_vertex_index_array();
  size_t nvtx = elem.vertex_count();
  if (nvtx > (sizeof(coords)/sizeof(coords[0]))) {
    MSQ_SETERR(err)("Invalid element type", MsqError::UNSUPPORTED_ELEMENT );
    return false;
  }
  
  Mesh::VertexHandle handles[8];
  for (unsigned i = 0; i < nvtx; ++i)
    handles[i] = pd.get_vertex_handles_array()[conn[i]];
  
  refMesh->get_reference_vertex_coordinates( handles, nvtx, coords, err );
  MSQ_ERRZERO(err);

  EntityTopology type = elem.get_element_type();
  unsigned num_edges = TopologyInfo::edges( type );
  double len_sum = 0.0;
  for (unsigned i = 0; i < num_edges; ++i) {
    const unsigned* edge = TopologyInfo::edge_vertices( type, i );
    len_sum += (coords[edge[0]] - coords[edge[1]]).length();
  }
  return len_sum * (1.0/num_edges); // scaleFactor[type];
}
                                             
 
bool RefSizeTargetCalculator::get_3D_target( PatchData& pd, 
                                             size_t element,
                                             Sample sample,
                                             MsqMatrix<3,3>& W,
                                             MsqError& err )
{
  scaledTargets->get_3D_target( pd, element, sample, W, err );
  MSQ_ERRZERO(err);
 
  double f = average_edge_length( pd, element, err );
  MSQ_ERRZERO(err);
  W *= f;
  
  return true;
}

bool RefSizeTargetCalculator::get_surface_target( PatchData& pd, 
                                            size_t element,
                                            Sample sample,
                                            MsqMatrix<3,2>& W,
                                            MsqError& err )
{
  scaledTargets->get_surface_target( pd, element, sample, W, err );
  MSQ_ERRZERO(err);
 
  double f = average_edge_length( pd, element, err );
  MSQ_ERRZERO(err);
  W *= f;
  
  return true;
}

bool RefSizeTargetCalculator::get_2D_target( PatchData& pd, 
                                            size_t element,
                                            Sample sample,
                                            MsqMatrix<2,2>& W,
                                            MsqError& err )
{
  scaledTargets->get_2D_target( pd, element, sample, W, err );
  MSQ_ERRZERO(err);
 
  double f = average_edge_length( pd, element, err );
  MSQ_ERRZERO(err);
  W *= f;
  
  return true;
}



} // namespace MESQUITE_NS
