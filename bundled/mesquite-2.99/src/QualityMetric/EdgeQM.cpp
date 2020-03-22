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


/** \file EdgeQM.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "EdgeQM.hpp"
#include "ElemSampleQM.hpp"
#include "PatchData.hpp"
#include "TopologyInfo.hpp"

namespace MESQUITE_NS {

EdgeQM::~EdgeQM()
{}


void EdgeQM::get_evaluations( PatchData& pd, 
                              std::vector<size_t>& handles,
                              bool free_vertices_only, 
                              MsqError& err )
{
  get_edge_evaluations( pd, handles, free_vertices_only, false, err );
}

void EdgeQM::get_single_pass( PatchData& pd, 
                              std::vector<size_t>& handles,
                              bool free_vertices_only, 
                              MsqError& err )
{
  get_edge_evaluations( pd, handles, free_vertices_only, true, err );
}


struct EdgeData {
  size_t endVtx;
  size_t adjElem;
  unsigned elemEdge;
};

inline bool operator<( const EdgeData& e1, const EdgeData& e2 )
  { return e1.endVtx < e2.endVtx; }
inline bool operator==( const EdgeData& e1, const EdgeData& e2 )
  { return e1.endVtx == e2.endVtx; }


void EdgeQM::get_edge_evaluations( PatchData& pd, 
                                   std::vector<size_t>& handles,
                                   bool free_vertices_only, 
                                   bool single_pass_evaluate,
                                   MsqError& err )
{
  std::vector<EdgeData> vtx_edges;
  size_t n_verts = free_vertices_only ? pd.num_free_vertices() : pd.num_nodes();
  size_t n_cutoff = single_pass_evaluate ? pd.num_nodes() : n_verts;
  handles.clear();

  for (size_t i = 0; i < n_verts; ++i) {
    if (pd.vertex_by_index(i).is_flag_set( MsqVertex::MSQ_PATCH_FIXED ))
      continue;

    vtx_edges.clear();

    size_t n_elems;
    const size_t* elems;
    elems = pd.get_vertex_element_adjacencies( i, n_elems, err );
    MSQ_ERRRTN(err);

    for (size_t j = 0; j < n_elems; ++j) {
      MsqMeshEntity& elem = pd.element_by_index(elems[j]);
      unsigned n_edges = TopologyInfo::edges( elem.get_element_type() );
      for (unsigned k = 0; k < n_edges; ++k) {
        const unsigned* edge = TopologyInfo::edge_vertices( elem.get_element_type(), k, err );
        MSQ_ERRRTN(err);

        size_t vtx1 = elem.get_vertex_index_array()[edge[0]];
        size_t vtx2 = elem.get_vertex_index_array()[edge[1]];
        size_t other;
        if (vtx1 == i)
          other = vtx2;
        else if (vtx2 == i)
          other = vtx1;
        else
          continue;
        
          // If !free_vertices_only, we'll visit every edge twice.  
          // The first check below ensures that we only add each edge
          // once.  The second check is never true unless free_vertices_only
          // is true and single_pass_evaluate is false.  In that case, it 
          // serves as an exception to the first rule for those cases in which 
          // we visit an edge only once.  For single_pass_evaluate (e.g.
          // BCD initialization or QualityAssessor) we want to avoid visiting
          // and edge more than once for every patch rather than just within
          // this patch.
        if (other > i || other > n_cutoff) {
          EdgeData d = { other, elems[j], k };
          vtx_edges.push_back(d);
        }
      } // end for each edge in element
    } // end for each element adjacent to vertex
    
    std::sort( vtx_edges.begin(), vtx_edges.end() );
    std::vector<EdgeData>::iterator it, end;
    end = std::unique( vtx_edges.begin(), vtx_edges.end() );
    for (it = vtx_edges.begin(); it != end; ++it) 
      handles.push_back( handle( it->elemEdge, it->adjElem ) );
  } // end for each (free) vertex
}

bool EdgeQM::evaluate_with_indices( PatchData& pd,
                                    size_t handle,
                                    double& value,
                                    std::vector<size_t>& indices,
                                    MsqError& err )
{
  const MsqMeshEntity& element = pd.element_by_index( elem(handle) );
  EntityTopology type = element.get_element_type(); 
  const unsigned* verts = TopologyInfo::edge_vertices( type, edge(handle) );
  const size_t* conn = element.get_vertex_index_array();
  indices.clear();
  if (conn[verts[0]] < pd.num_free_vertices())
    indices.push_back( conn[verts[0]] );
  if (conn[verts[1]] < pd.num_free_vertices())
    indices.push_back( conn[verts[1]] );
  return evaluate( pd, handle, value, err );
}

} // namespace MESQUITE_NS
