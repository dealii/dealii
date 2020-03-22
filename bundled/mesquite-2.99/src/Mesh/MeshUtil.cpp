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


/** \file MeshUtil.cpp
 *  \brief Implement MeshUtil class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MeshUtil.hpp"
#include "MeshInterface.hpp"
#include "EdgeIterator.hpp"
#include "SimpleStats.hpp"
#include "MsqError.hpp"
#include "TMPQualityMetric.hpp"
#include "MappingFunction.hpp"
#include "TargetCalculator.hpp"

#include <vector>
#include <algorithm>
#include <limits>
#include <sstream>

namespace MESQUITE_NS {


MeshUtil::~MeshUtil() { delete globalPatch; }

PatchData* MeshUtil::get_global_patch( MsqError& err ) 
{
  if (!globalPatch) {
    globalPatch = new PatchData;
    globalPatch->set_mesh( myMesh );
    if (mySettings)
      globalPatch->attach_settings( mySettings );
    globalPatch->fill_global_patch( err );
    if (MSQ_CHKERR(err)) {
      delete globalPatch;
      globalPatch = 0;
    }
  }
  return globalPatch;
}

void MeshUtil::edge_length_distribution( SimpleStats& results, MsqError& err )
{
  PatchData* pd = get_global_patch(err); MSQ_ERRRTN(err);
  EdgeIterator iter(pd, err); MSQ_ERRRTN(err);
  while (!iter.is_at_end()) {
    Vector3D diff = iter.start() - iter.end();
    results.add_squared( diff % diff );
    iter.step( err ); MSQ_ERRRTN(err);
  }
 
  if (results.empty()) { // no mesh
    MSQ_SETERR(err)("Mesh contains no elements", MsqError::INVALID_MESH);
  }
}

void MeshUtil::lambda_distribution( SimpleStats& results, MsqError& err )
{
  PatchData& pd = *get_global_patch(err); MSQ_ERRRTN(err);
  
  std::vector<size_t> handles;
  TMPQualityMetric::get_patch_evaluations( pd, handles, false, err );
  
  const size_t N = 64;
  size_t count;
  size_t indices[N];
  MsqVector<2> derivs2D[N];
  MsqVector<3> derivs3D[N];
  
  for (size_t i = 0; i < handles.size(); ++i) {
    double lambda;
    const Sample s = ElemSampleQM::sample( handles[i] );
    const size_t e = ElemSampleQM::  elem( handles[i] );
    EntityTopology type = pd.element_by_index(e).get_element_type();
    const NodeSet bits = pd.non_slave_node_set( e );

    if (TopologyInfo::dimension( type ) == 3)
    {
      const MappingFunction3D* mf = pd.get_mapping_function_3D( type );
      if (!mf) {
        MSQ_SETERR(err)( "No mapping function for element type", MsqError::UNSUPPORTED_ELEMENT );
        return;
      }
      
      MsqMatrix<3,3> W;
      mf->jacobian( pd, e, bits, s, indices, derivs3D, count, W, err ); MSQ_ERRRTN(err);
      assert(N >= count);
      lambda = TargetCalculator::size( W );
    }
    else 
    {
      assert( TopologyInfo::dimension( type ) == 2 );
      const MappingFunction2D* mf = pd.get_mapping_function_2D( type );
      if (!mf) {
        MSQ_SETERR(err)( "No mapping function for element type", MsqError::UNSUPPORTED_ELEMENT );
        return;
      }
      
      MsqMatrix<3,2> W;
      mf->jacobian( pd, e, bits, s, indices, derivs2D, count, W, err ); MSQ_ERRRTN(err);
      assert(N >= count);
      lambda = TargetCalculator::size( W );
    }

    results.add_value( lambda );
  }
  
  if (results.empty()) { // no mesh
    MSQ_SETERR(err)("Mesh contains no elements", MsqError::INVALID_MESH);
  }
}

bool MeshUtil::meshes_are_different(Mesh& mesh1, Mesh& mesh2, MsqError& err, double tol, bool do_print )
{
  //MsqError& err = (do_print ? MsqPrintError(cout) : MsqError());
    
  // mesh 1
  std::vector<Mesh::ElementHandle> elements1;
  std::vector<Mesh::VertexHandle> vertices1;
  std::vector<Mesh::VertexHandle> connectivity1;
  std::vector<size_t> offsets1;                                                                                      

  mesh1.get_all_elements(elements1, err);  MSQ_ERRZERO(err);
  mesh1.get_all_vertices(vertices1, err);  MSQ_ERRZERO(err);
  mesh1.elements_get_attached_vertices( arrptr(elements1), elements1.size(),
                                        connectivity1, offsets1, err ); MSQ_ERRZERO(err);
  std::vector<EntityTopology> types1(elements1.size());
  mesh1.elements_get_topologies( arrptr(elements1), arrptr(types1), elements1.size(), err ); MSQ_ERRZERO(err);

  // mesh 2
  std::vector<Mesh::ElementHandle> elements2;
  std::vector<Mesh::VertexHandle> vertices2;
  std::vector<Mesh::VertexHandle> connectivity2;
  std::vector<size_t> offsets2;                                                                                      

  mesh2.get_all_elements(elements2, err);  MSQ_ERRZERO(err);
  mesh2.get_all_vertices(vertices2, err);  MSQ_ERRZERO(err);
  mesh2.elements_get_attached_vertices( arrptr(elements2), elements2.size(),
                                        connectivity2, offsets2, err ); MSQ_ERRZERO(err);
  std::vector<EntityTopology> types2(elements2.size());
  mesh2.elements_get_topologies( arrptr(elements2), arrptr(types2), elements2.size(), err ); MSQ_ERRZERO(err);

#define MDRET(msg) do { if (do_print) std::cout << "MeshUtil::mesh_diff meshes are different: " << msg << std::endl; return true; } while(0)
  if (elements1.size() != elements2.size()) MDRET("elements sizes differ");
  if (vertices1.size() != vertices2.size()) MDRET("vertices sizes differ");
  if (offsets1.size() != offsets2.size()) MDRET("offsets sizes differ");
  if (connectivity1.size() != connectivity2.size()) MDRET("connectivity sizes differ");

  for (unsigned i=0; i < offsets1.size(); i++)
    {
      if (offsets1[i] != offsets2[i]) MDRET("offets differ");
    }

  for (unsigned i=0; i < vertices1.size(); i++)
    {
      MsqVertex vert1, vert2;
      mesh1.vertices_get_coordinates( &vertices1[i], &vert1, 1, err ); MSQ_ERRZERO(err);
      mesh2.vertices_get_coordinates( &vertices2[i], &vert2, 1, err ); MSQ_ERRZERO(err);
      double dist = Vector3D::distance_between(vert1, vert2);
      double v1 = vert1.length();
      double v2 = vert2.length();
      if ( dist > 0.5*(v1+v2)*tol ) {
        std::ostringstream ost;
        ost << "vertices coordinates differ more than tolerance [" << tol << "], vert1= " << vert1 << " vert2= " << vert2;
        MDRET(ost.str());
      }
    }

  for (unsigned i=0; i < connectivity1.size(); i++)
    {
      MsqVertex vert1, vert2;
      mesh1.vertices_get_coordinates( &connectivity1[i], &vert1, 1, err ); MSQ_ERRZERO(err);
      mesh2.vertices_get_coordinates( &connectivity2[i], &vert2, 1, err ); MSQ_ERRZERO(err);
      double dist = Vector3D::distance_between(vert1, vert2);
      double v1 = vert1.length();
      double v2 = vert2.length();
      if ( dist > 0.5*(v1+v2)*tol ) {
        std::ostringstream ost;
        ost << "connectivity coordinates differ more than tolerance [" << tol << "], vert1= " << vert1 << " vert2= " << vert2;
        MDRET(ost.str());
      }
    }

  return false;
}

} // namespace MESQUITE_NS
