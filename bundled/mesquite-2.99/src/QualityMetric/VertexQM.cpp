/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file VertexQM.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "VertexQM.hpp"
#include "PatchData.hpp"
#include "ElemSampleQM.hpp"

namespace MESQUITE_NS {

VertexQM::~VertexQM() {}

void VertexQM::get_evaluations( PatchData& pd, 
                                std::vector<size_t>& handles,
                                bool free_vertices_only, 
                                MsqError& err )
{
  get_vertex_evaluations( pd, handles, free_vertices_only, err );
}

void VertexQM::get_vertex_evaluations( PatchData& pd, 
                                       std::vector<size_t>& handles,
                                       bool free_vertices_only, 
                                       MsqError& err )
{
  if (free_vertices_only) {
    handles.resize( pd.num_free_vertices() );
    for (size_t i = 0; i < pd.num_free_vertices(); ++i)
      handles[i] = i;
  }
  else {
    handles.clear();
    handles.reserve( pd.num_nodes() );
    for (size_t i = 0; i < pd.num_nodes(); ++i)
      if (!(pd.vertex_by_index(i).get_flags() & MsqVertex::MSQ_PATCH_FIXED))
        handles.push_back(i);
  }
}

void VertexQM::get_vertex_corner_handles( PatchData& pd, 
                                          size_t vtx_idx,
                                          std::vector<size_t>& handles,
                                          MsqError& err )
{
  size_t len;
  const size_t *elems = pd.get_vertex_element_adjacencies( vtx_idx, len, err );
  MSQ_ERRRTN(err);
  
  handles.resize(len);
  for (size_t i = 0; i < len; ++i) {
    const MsqMeshEntity& elem = pd.element_by_index( elems[i] );
    const size_t* verts = elem.get_vertex_index_array();
    const size_t* ptr = std::find( verts, verts+elem.node_count(), vtx_idx );
    unsigned idx = ptr - verts;
    handles[i] = ElemSampleQM::handle( Sample(0,idx), elems[i] );
  }
}

} // namespace Mesquite
