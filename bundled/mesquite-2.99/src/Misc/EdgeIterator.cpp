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


/** \file EdgeIterator.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "EdgeIterator.hpp"

#include <algorithm>

namespace MESQUITE_NS {


EdgeIterator::EdgeIterator( PatchData* p, MsqError& err )
  : patchPtr(p),
    vertIdx(0)
{
  p->generate_vertex_to_element_data();
  if (patchPtr->num_nodes()) {
    get_adjacent_vertices( err );
    if (adjIter == adjList.end()) {
      step(err);
      MSQ_ERRRTN(err);
    }
  }
  else
    adjIter = adjList.end();
}

void EdgeIterator::get_adjacent_vertices( MsqError& err )
{
  adjList.clear();
  
    // Get all adjacent elements
  size_t num_elem;
  const size_t* elems = patchPtr->get_vertex_element_adjacencies( vertIdx, num_elem, err );
  MSQ_ERRRTN(err);
  
    // Get all adjacent vertices from elements
  std::vector<size_t> elem_verts;
  for (size_t e = 0; e < num_elem; ++e)
  {
    MsqMeshEntity& elem = patchPtr->element_by_index(elems[e]);
    EntityTopology type = elem.get_element_type();
    size_t num_edges = TopologyInfo::edges( type );
    
    bool mid_edge, mid_face, mid_vol;
    TopologyInfo::higher_order( type, elem.node_count(), mid_edge, mid_face, mid_vol, err );
    MSQ_ERRRTN(err);
    
      // For each edge
    for (size_t d = 0; d < num_edges; ++d)
    {
      const unsigned* edge = TopologyInfo::edge_vertices( type, d, err );
      MSQ_ERRRTN(err);
      size_t vert1 = elem.get_vertex_index( edge[0] );
      size_t vert2 = elem.get_vertex_index( edge[1] );

      size_t mid = ~(size_t)0;
      if (mid_edge) {
        int p = TopologyInfo::higher_order_from_side( type, elem.node_count(), 1, d, err );
        MSQ_ERRRTN(err);
        mid = elem.get_vertex_index_array()[p];
      }

        // If this edge contains the input vertex (vert_idx)
        // AND the input vertex index is less than the 
        // other vertex (avoids iterating over this edge twice)
        // add it to the list.
      if (vert1 > vert2)
      {
        if (vert2 == vertIdx)
          adjList.push_back( Edge(vert1,mid) );
      } 
      else 
      {
        if (vert1 == vertIdx)
          adjList.push_back( Edge(vert2,mid) );
      }
    }
  }
  
    // Remove duplicates
  std::sort( adjList.begin(), adjList.end() );
  adjIter = std::unique( adjList.begin(), adjList.end() );
  adjList.resize( adjIter - adjList.begin() );
  adjIter = adjList.begin();
}



} // namespace MESQUITE_NS
