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


/** \file MeshDecorator.hpp
 *  \brief Definition of Mesquite::MeshDecorator class
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_MESH_DECORATOR_HPP
#define MSQ_MESH_DECORATOR_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"

namespace MESQUITE_NS {

/**\brief Utility class for implementing decorators for the Mesquite::Mesh interface
 *
 * This is a utility class that to assist with implementing decorators
 * for the Mesquite::Mesh interface.  This class implements a "dumb"
 * decorator that passes all operations to its underlying (decorated)
 * Mesh instance unchanged.  The intention is that useful decorators 
 * subclass this class, overriding only those functions for which they
 * modify the behavioir.
 */
class MESQUITE_EXPORT MeshDecorator : public Mesh
{
  private:
  
    Mesh* myMesh;          //< The actual Mesh this instance decorates
  
  protected:
  
    void set_mesh( Mesh* mesh );
  
  public:
  
    MeshDecorator();
    MeshDecorator( Mesh* real_mesh );
    
    virtual ~MeshDecorator();
    
    /**\brief Get the real Mesh instance */
    Mesh* get_mesh( ) const { return myMesh; }
    
    
//************ Operations on entire mesh ****************

    virtual int get_geometric_dimension(MsqError &err);

    virtual void get_all_elements( std::vector<ElementHandle>& elements,
                                   MsqError& err );

    virtual void get_all_vertices( std::vector<VertexHandle>& vertices,
                                   MsqError& err );

//************ Vertex Properties ********************

    virtual void vertices_get_fixed_flag( const VertexHandle vert_array[], 
                                          std::vector<bool>& fixed_flag_array,
                                          size_t num_vtx, 
                                          MsqError &err );

    virtual void vertices_get_slaved_flag( const VertexHandle vert_array[], 
                                           std::vector<bool>& slaved_flag_array,
                                           size_t num_vtx, 
                                           MsqError &err );

    virtual void vertices_get_coordinates( const VertexHandle vert_array[],
                                           MsqVertex* coordinates,
                                           size_t num_vtx,
                                           MsqError &err );

    virtual void vertex_set_coordinates( VertexHandle vertex,
                                         const Vector3D &coordinates,
                                         MsqError &err );

    virtual void vertex_set_byte( VertexHandle vertex,
                                  unsigned char byte, 
                                  MsqError &err);

    virtual void vertices_set_byte( const VertexHandle *vert_array,
                                    const unsigned char *byte_array,
                                    size_t array_size, 
                                    MsqError &err );

    virtual void vertex_get_byte( const VertexHandle vertex,
                                  unsigned char *byte, 
                                  MsqError &err );

    virtual void vertices_get_byte( const VertexHandle *vertex,
                                    unsigned char *byte_array,
                                    size_t array_size, 
                                    MsqError &err );
    
//**************** Vertex Topology *****************    

    virtual void vertices_get_attached_elements( 
                         const VertexHandle* vertex_array,
                         size_t num_vertex,
                         std::vector<ElementHandle>& elements,
                         std::vector<size_t>& offsets,
                         MsqError& err );
    
//*************** Element Topology *************

    virtual void elements_get_attached_vertices(
                                   const ElementHandle *elem_handles,
                                   size_t num_elems,
                                   std::vector<VertexHandle>& vert_handles,
                                   std::vector<size_t>& offsets, 
                                   MsqError &err);
    

    virtual void elements_get_topologies(const ElementHandle *element_handle_array,
                                         EntityTopology *element_topologies,
                                         size_t num_elements, MsqError &err);

//***************  Tags  ***********

    virtual TagHandle tag_create( const std::string& tag_name,
                                  TagType type, unsigned length,
                                  const void* default_value,
                                  MsqError &err);

    virtual void tag_delete( TagHandle handle, MsqError& err );

    virtual TagHandle tag_get( const std::string& name, 
                               MsqError& err );

    virtual void tag_properties( TagHandle handle,
                                 std::string& name_out,
                                 TagType& type_out,
                                 unsigned& length_out,
                                 MsqError& err );

    virtual void tag_set_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       const void* tag_data,
                                       MsqError& err );

    virtual void tag_set_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       const void* tag_data,
                                       MsqError& err );

    virtual void tag_get_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       void* tag_data,
                                       MsqError& err );

    virtual void tag_get_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       void* tag_data,
                                       MsqError& err );

    
//**************** Memory Management ****************

    virtual void release_entity_handles( const EntityHandle *handle_array,
                                         size_t num_handles, 
                                         MsqError &err);

    virtual void release();

};


} // namespace Mesquite

#endif
