/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
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
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */

#include "MeshImplData.hpp"
#include "TopologyInfo.hpp"
#include "MsqError.hpp"

// NOTE: If this is defined, the normal vertex query functions
//       will not return mid-nodes.
#undef SEPARATE_MID_NODES

namespace MESQUITE_NS {

const std::vector<size_t> dummy_list;
const Vector3D dummy_vtx;

size_t MeshImplData::num_vertices() const 
{
  size_t count = 0;
  for (std::vector<Vertex>::const_iterator iter = vertexList.begin();
       iter != vertexList.end(); ++iter)
#ifdef SEPARATE_MID_NODES
    if (iter->valid && iter->midcount < iter->adjacencies.size())
#else
    if (iter->valid)
#endif
      ++count;
  return count;
}

size_t MeshImplData::num_elements() const
  { return elementList.size() - deletedElementList.size(); }

size_t MeshImplData::num_vertex_uses() const
{
  size_t result = 0;
  for (std::vector<Element>::const_iterator iter = elementList.begin();
       iter != elementList.end(); ++iter)
  {
#ifdef SEPARATE_MID_NODES
    unsigned from_topo = TopologyInfo::corners( iter->topology );
    result += from_topo ? from_topo : iter->connectivity.size();
#else
    result += iter->connectivity.size();
#endif
  }
  return result;
}

const Vector3D& MeshImplData::get_vertex_coords( size_t index, MsqError& err ) const
{
  if (!is_vertex_valid( index ))
  {
    MSQ_SETERR(err)("Invalid vertex handle", MsqError::INVALID_ARG);
    return dummy_vtx;
  }
  
  return vertexList[index].coords;
}

bool MeshImplData::vertex_is_fixed( size_t index, MsqError& err ) const
{
  if (!is_vertex_valid( index ))
  {
    MSQ_SETERR(err)("Invalid vertex handle", MsqError::INVALID_ARG);
    return false;
  }
  
  return vertexList[index].fixed;
}

bool MeshImplData::vertex_is_slaved( size_t index, MsqError& err ) const
{
  if (!is_vertex_valid( index ))
  {
    MSQ_SETERR(err)("Invalid vertex handle", MsqError::INVALID_ARG);
    return false;
  }
  if (!have_slaved_flags())
  {
    MSQ_SETERR(err)("Slave flags not set", MsqError::INVALID_STATE);
    return false;
  }
  
  return vertexList[index].slaved;
}

void MeshImplData::fix_vertex( size_t index, bool flag, MsqError& err )
{
  if (!is_vertex_valid( index ))
  {
    MSQ_SETERR(err)("Invalid vertex handle", MsqError::INVALID_ARG);
    return;
  }
  
  vertexList[index].fixed = flag;
}

void MeshImplData::slave_vertex( size_t index, bool flag, MsqError& err )
{
  if (!is_vertex_valid( index ))
  {
    MSQ_SETERR(err)("Invalid vertex handle", MsqError::INVALID_ARG);
    return;
  }
  
  vertexList[index].slaved = flag;
  haveSlavedFlags = true;
}


unsigned char MeshImplData::get_vertex_byte( size_t index, MsqError& err ) const
{
  if (!is_vertex_valid( index ))
  {
    MSQ_SETERR(err)("Invalid vertex handle", MsqError::INVALID_ARG);
    return 0;
  }
  
  return vertexList[index].byte;
}

void MeshImplData::set_vertex_byte( size_t index, unsigned char value, MsqError& err )
{
  if (!is_vertex_valid( index ))
  {
    MSQ_SETERR(err)("Invalid vertex handle", MsqError::INVALID_ARG);
    return;
  }
  
  vertexList[index].byte = value;
}

EntityTopology MeshImplData::element_topology( size_t index, MsqError& err ) const
{
  if (!is_element_valid( index ))
  {
    MSQ_SETERR(err)("Invalid element handle", MsqError::INVALID_ARG);
    return MIXED;
  }

  return elementList[index].topology;
}

void MeshImplData::element_topology( size_t index, EntityTopology type, MsqError& err )
{
  if (!is_element_valid( index ))
  {
    MSQ_SETERR(err)("Invalid element handle", MsqError::INVALID_ARG);
    return ;
  }
  
  unsigned i, numvert;
  
  numvert = TopologyInfo::corners( elementList[index].topology );
  if (numvert) 
    for (i = numvert; i < elementList[index].connectivity.size(); ++i)
      --vertexList[elementList[index].connectivity[i]].midcount;

  elementList[index].topology = type;
  
  numvert = TopologyInfo::corners( elementList[index].topology );
  if (numvert) 
    for (i = numvert; i < elementList[index].connectivity.size(); ++i)
      ++vertexList[elementList[index].connectivity[i]].midcount;
}


const std::vector<size_t>& MeshImplData::element_connectivity( size_t index, MsqError& err ) const
{
  if (!is_element_valid( index ))
  {
    MSQ_SETERR(err)("Invalid element handle", MsqError::INVALID_ARG);
    return dummy_list;
  }

  return elementList[index].connectivity;
}

const std::vector<size_t>& MeshImplData::vertex_adjacencies( size_t index, MsqError& err ) const
{
  if (!is_vertex_valid( index ))
  {
    MSQ_SETERR(err)("Invalid vertex handle", MsqError::INVALID_ARG);
    return dummy_list;
  }
  
  return vertexList[index].adjacencies;
}

void MeshImplData::clear()
{
  vertexList.clear();
  elementList.clear();
  deletedVertexList.clear();
  deletedElementList.clear();
  haveSlavedFlags = false;
}

void MeshImplData::allocate_vertices( size_t count, MsqError& err )
{
  if (vertexList.size())
  {
    MSQ_SETERR(err)(MsqError::INVALID_STATE);
    return;
  }
  
  vertexList.resize(count);
}

void MeshImplData::allocate_elements( size_t count, MsqError& err )
{
  if (elementList.size())
  {
    MSQ_SETERR(err)(MsqError::INVALID_STATE);
    return;
  }
  
  elementList.resize( count );
}

void MeshImplData::set_vertex_coords( size_t index, 
                                      const Vector3D& coords,
                                      MsqError& err )
{
  if (!is_vertex_valid( index ))
  {
    MSQ_SETERR(err)("Invalid vertex handle", MsqError::INVALID_ARG);
    return;
  }
  
  vertexList[index].coords = coords;
}

void MeshImplData::reset_vertex( size_t index, 
                                 const Vector3D& coords,
                                 bool fixed, 
                                 MsqError& err )
{
  if (index >= vertexList.size())
  {
    MSQ_SETERR(err)("Invalid vertex handle", MsqError::INVALID_ARG);
    return;
  }
  
  Vertex& vert = vertexList[index];
  
  if (!vert.adjacencies.empty())
  {
    MSQ_SETERR(err)("Cannot overwrite referenced vertex", MsqError::INVALID_STATE);
    return;
  }
  
  vert.coords = coords;
  vert.fixed = fixed;
  vert.valid = true;
}

void MeshImplData::reset_element( size_t index,
                                  const std::vector<long>& vertices,
                                  EntityTopology topology,
                                  MsqError& err )
{
  clear_element( index, err );                   MSQ_ERRRTN(err);
  set_element( index, vertices, topology, err ); MSQ_ERRRTN(err);
}

void MeshImplData::reset_element( size_t index,
                                  const std::vector<size_t>& vertices,
                                  EntityTopology topology,
                                  MsqError& err )
{
  clear_element( index, err );                   MSQ_ERRRTN(err);
  set_element( index, vertices, topology, err ); MSQ_ERRRTN(err);
}


void MeshImplData::clear_element( size_t index, MsqError& err )
{  
  if (index >= elementList.size())
  {
    MSQ_SETERR(err)("Invalid element handle", MsqError::INVALID_ARG);
    return;
  }
  
  unsigned numvert = TopologyInfo::corners( elementList[index].topology );
  if (numvert) 
    for (unsigned i = numvert; i < elementList[index].connectivity.size(); ++i)
      --vertexList[elementList[index].connectivity[i]].midcount;

  std::vector<size_t>& conn = elementList[index].connectivity;
  for (std::vector<size_t>::iterator iter = conn.begin();
       iter != conn.end(); ++iter)
  {
    std::vector<size_t>& adj = vertexList[*iter].adjacencies;
    for (std::vector<size_t>::iterator iter2 = adj.begin();
         iter2 != adj.end(); ++iter2)
    {
      if (*iter2 == index)
      {
        adj.erase(iter2);
        break;
      }
    }
  }
  conn.clear();
}

void MeshImplData::set_element( size_t index,
                                const std::vector<long>& vertices,
                                EntityTopology topology,
                                MsqError& err )
{
  if (sizeof(long) == sizeof(size_t)) 
    set_element( index, *reinterpret_cast<const std::vector<size_t>*>(&vertices), topology, err );
  else {
    std::vector<size_t> conn(vertices.size());
    std::copy( vertices.begin(), vertices.end(), conn.begin() );
    set_element( index, conn, topology, err );
  }
}

void MeshImplData::set_element( size_t index,
                                const std::vector<size_t>& vertices,
                                EntityTopology topology,
                                MsqError& err )
{
  if (index >= elementList.size())
  {
    MSQ_SETERR(err)("Invalid element handle", MsqError::INVALID_ARG);
    return;
  }

  elementList[index].connectivity = vertices;
  elementList[index].topology = topology;
  
  for (std::vector<size_t>::const_iterator iter = vertices.begin();
       iter != vertices.end(); ++iter)
  {
    if (!is_vertex_valid( *iter ))
    {
      MSQ_SETERR(err)("Invalid vertex handle", MsqError::INVALID_ARG);
      return;
    }
    
    std::vector<size_t>& adj = vertexList[*iter].adjacencies;
    for (std::vector<size_t>::iterator iter2 = adj.begin();
         iter2 != adj.end(); ++iter2)
      if (*iter2 == index)
        return;
    
    adj.push_back( index );
  }
  
  unsigned numvert = TopologyInfo::corners( elementList[index].topology );
  if (numvert) 
    for (unsigned i = numvert; i < elementList[index].connectivity.size(); ++i)
      ++vertexList[elementList[index].connectivity[i]].midcount;
}

size_t MeshImplData::add_vertex( const Vector3D& coords, bool fixed, MsqError& err )
{
  size_t index;
  
  if (!deletedVertexList.empty())
  {
    index = deletedVertexList[deletedVertexList.size()-1];
    deletedVertexList.pop_back();
    reset_vertex( index, coords, fixed, err ); MSQ_ERRZERO(err);
  }
  else
  {
    index = vertexList.size();
    vertexList.push_back( Vertex(coords, fixed ) );
  }
  
  return index;
}

size_t MeshImplData::add_element( const std::vector<long>& vertices,
                                  EntityTopology topology,
                                  MsqError& err )
{
  size_t index;
  if (!deletedElementList.empty())
  {
    index = deletedElementList[deletedElementList.size()-1];
    deletedElementList.pop_back();
  }
  else
  {
    index = elementList.size();
    elementList.resize( elementList.size() + 1 );
  }

  set_element( index, vertices, topology, err ); MSQ_ERRZERO(err);
  return index;
}

size_t MeshImplData::add_element( const std::vector<size_t>& vertices,
                                  EntityTopology topology,
                                  MsqError& err )
{
  size_t index;
  if (!deletedElementList.empty())
  {
    index = deletedElementList[deletedElementList.size()-1];
    deletedElementList.pop_back();
  }
  else
  {
    index = elementList.size();
    elementList.resize( elementList.size() + 1 );
  }

  set_element( index, vertices, topology, err ); MSQ_ERRZERO(err);
  return index;
}

void MeshImplData::delete_vertex( size_t index, MsqError& err )
{
  if (!is_vertex_valid( index ))
  {
    MSQ_SETERR(err)("Invalid vertex handle", MsqError::INVALID_ARG);
    return;
  }
  
  vertexList[index].valid = false;
  deletedVertexList.push_back( index );
}

void MeshImplData::delete_element( size_t index, MsqError& err )
{
  clear_element( index, err ); MSQ_ERRRTN(err);
  deletedElementList.push_back( index );
}

void MeshImplData::copy_mesh( size_t* vertex_handle_array,
                              size_t* element_handle_array,
                              size_t* element_conn_offsets,
                              size_t* element_conn_indices )
{
  std::vector<size_t> vertex_map( vertexList.size() );
  size_t vh_index = 0;
  for (size_t v = 0; v < vertexList.size(); ++v)
  {
    if (vertexList[v].valid
#ifdef SEPARATE_MID_NODES
        && vertexList[v].midcount < vertexList[v].adjacencies.size()
#endif
       )
    {
      vertex_handle_array[vh_index] = v;
      vertex_map[v] = vh_index;
      ++vh_index;
    }
    else
    {
      vertex_map[v] = vertexList.size();
    }
  }
  
  size_t offset = 0;
  for (size_t e = 0; e < elementList.size(); ++e)
  {
    Element& elem = elementList[e];
    size_t cl;
#ifdef SEPARATE_MID_NODES
    cl = TopologyInfo::corners( elem.topology );
    if (!cl)
#endif
      cl = elem.connectivity.size();
    if (cl)
    {
      *element_handle_array = e;
      ++element_handle_array;
      
      *element_conn_offsets = offset;
      ++element_conn_offsets;
      offset += cl;
    
      std::vector<size_t>::iterator conn = elem.connectivity.begin();
      std::vector<size_t>::iterator end = conn + cl;
      while (conn != end)
      {
        *element_conn_indices = vertex_map[*conn];
        ++element_conn_indices;
        ++conn;
      }
    }
  }
  *element_conn_offsets = offset;
}

void MeshImplData::copy_higher_order( std::vector<size_t>& mid_nodes,
                                      std::vector<size_t>& vertices,
                                      std::vector<size_t>& vertex_indices,
                                      std::vector<size_t>& index_offsets,
                                      MsqError& err )
{
  mid_nodes.clear();
  vertices.clear();
  vertex_indices.clear();
  index_offsets.clear();
 
    // Create a map of from vertex handle to index in "vertices"
    // Use vertexList.size() to mean uninitialized.
  size_t v;
  std::vector<size_t> vert_map( vertexList.size() );
  for (v = 0; v < vertexList.size(); ++v)
    vert_map[v] = vertexList.size();

    // Loop over all mid-side vertices
  for (v = 0; v < vertexList.size(); ++v)
  {
    const Vertex& vert = vertexList[v];
    
      // Not a mid-side vertex, skip it
    if (!vert.valid || !vert.midcount)
      continue;
      
      // Populate "verts" with the handles of all adjacent corner vertices
    assert( vert.adjacencies.size() ); // shouldn't be able to fail if vert.midcount > 0
    int elem_indx = vert.adjacencies[0];
    Element& elem = elementList[elem_indx];
    
      // Find index of node in elem's connectivity list
    unsigned index;
    for (index = 0; index < elem.connectivity.size(); ++index)
      if (elem.connectivity[index] == v)
        break;
    if (index == elem.connectivity.size())
    {
      MSQ_SETERR(err)("Inconsistent data.", MsqError::INTERNAL_ERROR);
      return;
    }

      // Given the index in the element's connectivity list,
      // get the side of the element containing the mid-node.
    unsigned side_dim, side_num;
    TopologyInfo::side_number( elem.topology, elem.connectivity.size(),
                          index, side_dim, side_num, err ); MSQ_ERRRTN(err);

    if (!side_dim) // Not a mid-side node
    {
      MSQ_SETERR(err)(MsqError::INVALID_STATE,"Improperly connected mesh.");
      return;
    }

      // Get the adjacent corner vertices from the element side.
    unsigned num_corners;
    const unsigned* corner_indices = TopologyInfo::side_vertices( 
      elem.topology, side_dim, side_num, num_corners, err ); MSQ_ERRRTN(err);

      // Add the mid-side node to the output list
    mid_nodes.push_back( v );
      // Store offset at which the indices of the corner
      // vertices adjacent to this mid-side node will be 
      // stored in "vertex_indices".
    index_offsets.push_back( vertex_indices.size() );
      // For each adjacent corner vertex, if the vertex is not
      // already in "vertices" add it, and add the index to
      // the adjacency list for this mid-side node.
    for (unsigned i = 0; i < num_corners; ++i)
    {
      size_t vert_idx = elem.connectivity[corner_indices[i]];
      assert( is_vertex_valid(vert_idx) );
      
      if (vert_map[vert_idx] == vertexList.size())
      {
        vert_map[vert_idx] = vertices.size();
        vertices.push_back( vert_idx );
      }
      vertex_indices.push_back( vert_map[vert_idx] );
    }
  }
  index_offsets.push_back( vertex_indices.size() );
}

bool MeshImplData::is_mid_node( size_t index ) const
{
  return is_vertex_valid(index) && vertexList[index].midcount > 0;
}

bool MeshImplData::is_corner_node( size_t index ) const
{
  return is_vertex_valid(index) && vertexList[index].midcount < vertexList[index].adjacencies.size();
}



void MeshImplData::all_vertices( std::vector<size_t>& list, MsqError& ) const
{
  list.clear();
  for (size_t idx = 0; idx < vertexList.size(); ++idx)
    if (vertexList[idx].valid)
      list.push_back( idx );
}

void MeshImplData::all_elements( std::vector<size_t>& list, MsqError& err ) const
{
  list.clear();
  for (size_t idx = 0; idx < elementList.size(); ++idx)
    if (!elementList[idx].connectivity.empty())
      list.push_back( idx );
}

void MeshImplData::get_adjacent_elements( 
                        std::vector<size_t>::const_iterator node_iter,
                        std::vector<size_t>::const_iterator node_end,
                        std::vector<size_t>& elems, MsqError& err )
{
  if (node_iter == node_end || !is_vertex_valid( *node_iter ))
  {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    return;
  }
  
    // Get list of elements adjacent to first node
  elems = vertexList[*node_iter].adjacencies;
  
    // For each aditional node, intersect elems with elements adjacent to node
  for (++node_iter; node_iter != node_end; ++node_iter)
  {
    std::vector<size_t>::iterator elem_iter = elems.begin();
    while (elem_iter != elems.end())
    {
      std::vector<size_t>::const_iterator adj_iter = vertexList[*node_iter].adjacencies.begin();
      const std::vector<size_t>::const_iterator adj_end = vertexList[*node_iter].adjacencies.end();
      for (; adj_iter != adj_end; ++adj_iter)
        if (*elem_iter == *adj_iter)
          break;
      
      if (adj_iter == adj_end)
      {
        *elem_iter = elems[elems.size()-1];
        elems.pop_back();
      }
      else
      {
        ++elem_iter;
      }
    }
  }
}

bool MeshImplData::has_adjacent_elements( 
                        size_t elem,
                        const std::vector<size_t>& nodes,
                        MsqError& err )
{
  std::vector<size_t> adj_elems;
  const unsigned dim = TopologyInfo::dimension( elementList[elem].topology );
  get_adjacent_elements( nodes.begin(), nodes.end(), adj_elems, err );
  
  std::vector<size_t>::iterator iter;
  for (iter = adj_elems.begin(); iter != adj_elems.end(); ++iter)
    if (*iter != elem && 
        TopologyInfo::dimension( elementList[*iter].topology ) == dim )
      break;

  return iter != adj_elems.end();
}                 

void MeshImplData::skin( std::vector<size_t>& sides, MsqError& err ) 
{
  std::vector<size_t> side_nodes;
  
    // For each element in mesh
  for (size_t elem = 0; elem < elementList.size(); ++elem)
  {
    if (!is_element_valid(elem))
      continue;
    
      // For each side of the element, check if there
      // are any adjacent elements.
    const EntityTopology topo = elementList[elem].topology;
    std::vector<size_t>& conn = elementList[elem].connectivity;
    switch (topo)
    {
        // For normal elements (not poly****)
      default:
      {
        unsigned num = TopologyInfo::sides( topo );
        unsigned dim = TopologyInfo::dimension( topo ) - 1;
          // For each side
        for (unsigned side = 0; side < num; ++side)
        {
            // Get list of vertices defining the side
          unsigned count;
          const unsigned* indices = TopologyInfo::side_vertices( topo, dim, side, count, err );
          MSQ_ERRRTN(err);
          side_nodes.clear();
          for (unsigned k = 0; k < count; ++k)
            side_nodes.push_back( conn[indices[k]] );
          
            // If no adjacent element, add side to output list
          bool adj = has_adjacent_elements( elem, side_nodes, err );
          MSQ_ERRRTN(err);
          if ( !adj )
          {
            sides.push_back( elem );
            sides.push_back( side );
          }
        }
      }
      break;
      
      case POLYGON:
      {
        for (unsigned side = 0, next = 1; next < conn.size(); ++side, ++next)
        {
          side_nodes.clear();
          side_nodes.push_back( conn[side] );
          side_nodes.push_back( conn[next] );
          
            // If no adjacent element, add side to output list
          bool adj = has_adjacent_elements( elem, side_nodes, err );
          MSQ_ERRRTN(err);
          if ( !adj )
          {
            sides.push_back( elem );
            sides.push_back( side );
          }
        }
      }
      break;
      
      case POLYHEDRON:
      {
        for (unsigned side = 0; side < conn.size(); ++side)
        {
          side_nodes = elementList[conn[side]].connectivity;
          
            // If no adjacent element, add side to output list
          bool adj = has_adjacent_elements( elem, side_nodes, err );
          MSQ_ERRRTN(err);
          if ( !adj )
          {
            sides.push_back( elem );
            sides.push_back( side );
          }
        }
      }
      break;
    } // switch(topo)
  } // for (elementList)
}          


MeshImplVertIter::~MeshImplVertIter()
  {}

void MeshImplVertIter::restart()
{
  index = 0;
  if (!mesh->is_vertex_valid( index ))
    operator++();
}

void MeshImplVertIter::operator++()
{
  ++index;
  while (index < mesh->max_vertex_index() && 
         (!mesh->is_vertex_valid(index) ||
          !mesh->is_corner_node(index)))
    ++index;
}

Mesh::VertexHandle MeshImplVertIter::operator*() const
{
  return reinterpret_cast<Mesh::VertexHandle>(index);
}

bool MeshImplVertIter::is_at_end() const
{ 
  return index >= mesh->max_vertex_index();
}



MeshImplElemIter::~MeshImplElemIter()
  {}

void MeshImplElemIter::restart()
{
  index = 0;
  if (!mesh->is_element_valid( index ))
    operator++();
}

void MeshImplElemIter::operator++()
{
  ++index;
  while (index < mesh->max_element_index() && !mesh->is_element_valid(index))
    ++index;
}

Mesh::ElementHandle MeshImplElemIter::operator*() const
{
  return reinterpret_cast<Mesh::ElementHandle>(index);
}

bool MeshImplElemIter::is_at_end() const
{ 
  return index >= mesh->max_element_index();
}



} // namespace Mesquite

      
