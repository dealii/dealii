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


/** \file ElementQM.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ElementQM.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {

ElementQM::~ElementQM() {}

void ElementQM::get_evaluations( PatchData& pd, 
                        std::vector<size_t>& handles,
                        bool free_vertices_only, 
                        MsqError& err )
{
  get_element_evaluations( pd, handles, free_vertices_only, err );
}

void ElementQM::get_element_evaluations( PatchData& pd, 
                        std::vector<size_t>& handles,
                        bool free_vertices_only, 
                        MsqError& err )
{
  size_t num_elem = pd.num_elements();
  if (!free_vertices_only) {
    handles.resize( num_elem );
    for (size_t i = 0; i < num_elem; ++i)
      handles[i] = i;
    return;
  }
  
  handles.clear();
  for (size_t i = 0; i < num_elem; ++i) {
      // check if element has any free vertices
    MsqMeshEntity& elem = pd.element_by_index(i);
    unsigned num_vtx = elem.node_count();
    size_t* vtx = elem.get_vertex_index_array();
    unsigned j;
    for (j = 0; j < num_vtx && vtx[j] >= pd.num_free_vertices(); ++j);
    if (j < num_vtx)
      handles.push_back(i);
  }
}

bool ElementQM::evaluate_with_indices( PatchData& pd,
                                       size_t handle,
                                       double& value,
                                       std::vector<size_t>& indices,
                                       MsqError& err )
{
  const MsqMeshEntity& elem = pd.element_by_index(handle);
  const size_t* vtx = elem.get_vertex_index_array();
  const unsigned n = elem.vertex_count();
  indices.clear();
  for (unsigned i = 0; i < n; ++i)
    if (vtx[i] < pd.num_free_vertices())
      indices.push_back( vtx[i] );
  
  bool rval = evaluate( pd, handle, value, err );
  return !MSQ_CHKERR(err) && rval;
}

} // namespace Mesquite
