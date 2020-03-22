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


/** \file MeshDecorator.cpp
 *  \brief Implementation of Mesquite::MeshDecorator class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MeshDecorator.hpp"

namespace MESQUITE_NS {

MeshDecorator::MeshDecorator() : myMesh(0) {}
MeshDecorator::MeshDecorator( Mesh* mesh ) : myMesh(mesh) {}

MeshDecorator::~MeshDecorator()
{ }

void MeshDecorator::set_mesh( Mesh* mesh )
  { myMesh = mesh; }

void MeshDecorator::vertices_get_coordinates( const VertexHandle vert_array[],
                                              MsqVertex* coordinates,
                                              size_t num_vtx,
                                              MsqError &err )
  { return get_mesh()->vertices_get_coordinates( vert_array, coordinates, num_vtx, err ); }
    
void MeshDecorator::vertex_set_coordinates( VertexHandle vertex,
                                            const Vector3D &coordinates,
                                            MsqError &err )
  { return get_mesh()->vertex_set_coordinates( vertex, coordinates, err ); }


//************ Operations on entire mesh ****************

int MeshDecorator::get_geometric_dimension( MsqError& err )
  { return get_mesh()->get_geometric_dimension( err ); }

void MeshDecorator::get_all_elements( std::vector<ElementHandle>& elements, MsqError& err )
  { get_mesh()->get_all_elements( elements, err ); }

void MeshDecorator::get_all_vertices( std::vector<VertexHandle>& vertices, MsqError& err )
  { get_mesh()->get_all_vertices( vertices, err ); }


//************ Vertex Properties ********************

void MeshDecorator::vertices_get_fixed_flag( const VertexHandle vert_array[], 
                                             std::vector<bool>& fixed_flag_array,
                                             size_t num_vtx, 
                                             MsqError &err )
  { get_mesh()->vertices_get_fixed_flag( vert_array, fixed_flag_array, num_vtx, err ); }

void MeshDecorator::vertices_get_slaved_flag( const VertexHandle vert_array[], 
                                              std::vector<bool>& flag_array,
                                              size_t num_vtx, 
                                              MsqError &err )
  { get_mesh()->vertices_get_slaved_flag( vert_array, flag_array, num_vtx, err ); }

void MeshDecorator::vertex_set_byte( VertexHandle vertex,
                                     unsigned char byte, 
                                     MsqError &err)
  { get_mesh()->vertex_set_byte( vertex, byte, err ); }

void MeshDecorator::vertices_set_byte( const VertexHandle *vert_array,
                                       const unsigned char *byte_array,
                                       size_t array_size, 
                                       MsqError &err )
  { get_mesh()->vertices_set_byte( vert_array, byte_array, array_size, err ); }

void MeshDecorator::vertex_get_byte( const VertexHandle vertex,
                                     unsigned char *byte, 
                                     MsqError &err )
  { get_mesh()->vertex_get_byte( vertex, byte, err ); }

void MeshDecorator::vertices_get_byte( const VertexHandle *vertex,
                                       unsigned char *byte_array,
                                       size_t array_size, 
                                       MsqError &err )
  { get_mesh()->vertices_get_byte( vertex, byte_array, array_size, err ); }
    
//**************** Vertex Topology *****************    

void MeshDecorator::vertices_get_attached_elements( 
                                const VertexHandle* vertex_array,
                                size_t num_vertex,
                                std::vector<ElementHandle>& elements,
                                std::vector<size_t>& offsets,
                                MsqError& err )
  { get_mesh()->vertices_get_attached_elements( vertex_array, num_vertex, elements, offsets, err ); }
    
//*************** Element Topology *************

void MeshDecorator::elements_get_attached_vertices(
                                   const ElementHandle *elem_handles,
                                   size_t num_elems,
                                   std::vector<VertexHandle>& vert_handles,
                                   std::vector<size_t>& offsets, 
                                   MsqError &err)
  { get_mesh()->elements_get_attached_vertices( elem_handles, num_elems, vert_handles, offsets, err ); }

void MeshDecorator::elements_get_topologies(
                                    const ElementHandle *element_handle_array,
                                    EntityTopology *element_topologies,
                                    size_t num_elements, MsqError &err)
  { get_mesh()->elements_get_topologies( element_handle_array, element_topologies, num_elements, err ); }

//***************  Tags  ***********

TagHandle MeshDecorator::tag_create( const std::string& tag_name,
                                     TagType type, unsigned length,
                                     const void* default_value,
                                     MsqError &err)
{ 
  return get_mesh()->tag_create( tag_name, type, length, default_value, err );
}

void MeshDecorator::tag_delete( TagHandle handle, MsqError& err )
  { get_mesh()->tag_delete( handle, err ); }

TagHandle MeshDecorator::tag_get( const std::string& name, MsqError& err )
  { return get_mesh()->tag_get( name, err ); }


void MeshDecorator::tag_properties( TagHandle handle,
                                    std::string& name_out,
                                    TagType& type_out,
                                    unsigned& length_out,
                                    MsqError& err )
  { get_mesh()->tag_properties( handle, name_out, type_out, length_out, err ); }

void MeshDecorator::tag_set_element_data( TagHandle handle,
                                          size_t num_elems,
                                          const ElementHandle* elem_array,
                                          const void* tag_data,
                                          MsqError& err )
  { get_mesh()->tag_set_element_data( handle, num_elems, elem_array, tag_data, err ); }

void MeshDecorator::tag_set_vertex_data ( TagHandle handle,
                                          size_t num_elems,
                                          const VertexHandle* node_array,
                                          const void* tag_data,
                                          MsqError& err )
  { get_mesh()->tag_set_vertex_data( handle, num_elems, node_array, tag_data, err ); }

void MeshDecorator::tag_get_element_data( TagHandle handle,
                                          size_t num_elems,
                                          const ElementHandle* elem_array,
                                          void* tag_data,
                                          MsqError& err )
  { get_mesh()->tag_get_element_data( handle, num_elems, elem_array, tag_data, err ); }

void MeshDecorator::tag_get_vertex_data ( TagHandle handle,
                                          size_t num_elems,
                                          const VertexHandle* node_array,
                                          void* tag_data,
                                          MsqError& err )
  { get_mesh()->tag_get_vertex_data( handle, num_elems, node_array, tag_data, err ); }
    
//**************** Memory Management ****************

void MeshDecorator::release_entity_handles( const EntityHandle *handle_array,
                                            size_t num_handles, 
                                            MsqError &err)
  { get_mesh()->release_entity_handles( handle_array, num_handles, err ); }

void MeshDecorator::release()
  { get_mesh()->release(); }


} // namespace Mesquite
