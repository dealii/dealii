/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov  ,
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */
/*!
  \file   ParallelMeshImpl.hpp
  \brief  

  \author Darryl Melander
  \author Thomas Leurent
  \author Jason Kraftcheck
  \author Martin Isenburg
  \date   2003-04-17
*/

#ifndef PARALLEL_MESQUITE_MESH_IMPL_HPP
#define PARALLEL_MESQUITE_MESH_IMPL_HPP

#include "ParallelMeshInterface.hpp"

namespace MESQUITE_NS
{
  /*!  \class ParallelMeshImpl
  \brief ParallelMeshImpl is a Mesquite implementation of the ParallelMesh
  interface. It inherits all of the implementation from MeshImpl and only
  implements any additional functionality.
  */

  class MESQUITE_EXPORT ParallelMeshImpl : public ParallelMesh
  {
  public:

    ParallelMeshImpl(Mesh* myMesh, const char * gid_name=0, const char * pid_name=0);

    void set_global_id_tag(const char * name, MsqError& err);
    void set_processor_id_tag(const char * name, MsqError& err);

    void vertices_set_global_id(const VertexHandle vert_array[],
				size_t gid[],
				size_t num_vtx,
				MsqError& err);
    
    void vertices_set_processor_id(const VertexHandle vert_array[],
				   int pid[],
				   size_t num_vtx,
				   MsqError& err);
    
//**************** Inherited Methods from ParallelMesh ******************************

    /*! Get global ids for given vertices.
     */
    virtual void vertices_get_global_id(const VertexHandle vert_array[],
                                        size_t gid[],
                                        size_t num_vtx,
                                        MsqError& err);
         
    /*! Get processor ids for given vertices.
     */
    virtual void vertices_get_processor_id(const VertexHandle vert_array[],
                                           int pid[],
                                           size_t num_vtx,
                                           MsqError& err);
     
//**************** Inherited Methods from Mesh ******************************
    
      // Returns whether this mesh lies in a 2D or 3D coordinate system.
    virtual int get_geometric_dimension(MsqError &err) ;
    
    
    /** \brief Get all elements in mesh
     *
     * Get the handles of every element in the active mesh.
     */ 
    virtual void get_all_elements( std::vector<ElementHandle>& handles,
                                   MsqError& err );
    
    /** \brief Get all vertices in mesh
     *
     * Get the handles of every vertex in the active mesh
     */
    virtual void get_all_vertices( std::vector<VertexHandle>& vertices,
                                   MsqError& err );

//************ Vertex Properties ********************
      //! Returns true or false, indicating whether the vertex
      //! is allowed to be repositioned.  True indicates that the vertex
      //! is fixed and cannot be moved.  Note that this is a read-only
      //! property; this flag can't be modified by users of the
      //! Mesquite::Mesh interface.
    virtual void vertices_get_fixed_flag( const VertexHandle vert_array[], 
                                          std::vector<bool>& fixed_flag_array,
                                          size_t num_vtx, 
                                          MsqError &err);

    void vertices_set_fixed_flag( const VertexHandle vert_array[], 
                                  const bool fixed_flag_array[],
                                  size_t num_vtx, 
                                  MsqError &err);

    virtual void vertices_get_slaved_flag( const VertexHandle vert_array[], 
                                           std::vector<bool>& slaved_flag_array,
                                           size_t num_vtx, 
                                           MsqError &err );
    
      // Get/set location of a vertex
    virtual void vertices_get_coordinates(const Mesh::VertexHandle vert_array[],
                                          Mesquite::MsqVertex* coordinates,
                                          size_t num_vtx,
                                          MsqError &err);
    virtual void vertex_set_coordinates(VertexHandle vertex,
                                        const Vector3D &coordinates,
                                        MsqError &err);
    
      // Each vertex has a byte-sized flag that can be used to store
      // flags.  This byte's value is neither set nor used by the mesh
      // implementation.  It is intended to be used by Mesquite algorithms.
      // Until a vertex's byte has been explicitly set, its value is 0.
    virtual void vertex_set_byte (VertexHandle vertex,
                                  unsigned char byte,
                                  MsqError &err);
    virtual void vertices_set_byte (const VertexHandle *vert_array,
                                    const unsigned char *byte_array,
                                    size_t array_size,
                                    MsqError &err);
    
      // Retrieve the byte value for the specified vertex or vertices.
      // The byte value is 0 if it has not yet been set via one of the
      // *_set_byte() functions.
    virtual void vertex_get_byte(const VertexHandle vertex,
                                 unsigned char *byte,
                                 MsqError &err);
    virtual void vertices_get_byte(const VertexHandle *vertex,
                                   unsigned char *byte_array,
                                   size_t array_size,
                                   MsqError &err);
    
//**************** Vertex Topology *****************    

      /** \brief get elements adjacent to vertices
       *
       * Get adjacency data for vertices
       *
       *\param vertex_array    Array of vertex handles specifying the
       *                       list of vertices to retrieve adjacency
       *                       data for.
       *\param num_vertex      Number of vertex handles in #vertex_array
       *\param elements     The array in which to place the handles of
       *                       elements adjacent to the input vertices.
       *\param offsets    For each vertex in #vertex_array, the
       *                       value in the corresponding position in this
       *                       array is the index into #elem_array at
       *                       which the adjacency list begins for that
       *                       vertex.
       */
    virtual void vertices_get_attached_elements( 
                         const VertexHandle* vertex_array,
                         size_t num_vertex,
                         std::vector<ElementHandle>& elements,
                         std::vector<size_t>& offsets,
                         MsqError& err );
    
//*************** Element Topology *************
    
      /** \brief Get element connectivity
       *
       * Get the connectivity (ordered list of vertex handles) for
       * each element in the input array.
       *
       *\param elem_handles  The array of element handles for which to
       *                     retrieve the connectivity list.
       *\param num_elems     The length of #elem_handles
       *\param vert_handles  Array in which to place the vertex handles
       *                     in each elements connectivity.
       *\param offsets       For each element in #elem_handles, the
       *                     value in the same position in this array
       *                     is the index into #vert_handles at which
       *                     the connectivity list for that element begins.
       */
    virtual void elements_get_attached_vertices(
                                   const ElementHandle *elem_handles,
                                   size_t num_elems,
                                   std::vector<VertexHandle>& vert_handles,
                                   std::vector<size_t>& offsets, 
                                   MsqError &err);

    
      // Returns the topologies of the given entities.  The "entity_topologies"
      // array must be at least "num_elements" in size.
    virtual void elements_get_topologies(const ElementHandle *element_handle_array,
                                         EntityTopology *element_topologies,
                                         size_t num_elements,
                                         MsqError &err);

    
//*************** Tags  ***********

      /** \brief Create a tag
       *
       * Create a user-defined data type that can be attached
       * to any element or vertex in the mesh.  For an opaque or
       * undefined type, use type=BYTE and length=sizeof(..).
       *
       * \param tag_name  A unique name for the data object
       * \param type      The type of the data
       * \param length    Number of values per entity (1->scalar, >1 ->vector)
       * \param default_value Default value to assign to all entities - may be NULL
       * \return - Handle for tag definition 
       */
    virtual TagHandle tag_create( const std::string& tag_name,
                                  TagType type, unsigned length,
                                  const void* default_value,
                                  MsqError &err);
     
      /** \brief Remove a tag and all corresponding data
       *
       * Delete a tag.
       */
    virtual void tag_delete( TagHandle handle, MsqError& err );
    
    
      /** \brief Get handle for existing tag, by name. */
    virtual TagHandle tag_get( const std::string& name, 
                               MsqError& err );
     
      /** \brief Get properites of tag
       *
       * Get data type and number of values per entity for tag.
       * \param handle     Tag to get properties of.
       * \param name_out   Passed back tag name.
       * \param type_out   Passed back tag type.
       * \param length_out Passed back number of values per entity.
       */
    virtual void tag_properties( TagHandle handle,
                                 std::string& name_out,
                                 TagType& type_out,
                                 unsigned& length_out,
                                 MsqError& err );
    
      /** \brief Set tag values on elements
       * 
       * Set the value of a tag for a list of mesh elements.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of elements for which to set the tag value.
       * \param tag_data   Tag data for each element, contiguous in memory.
       *                   This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_set_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       const void* tag_data,
                                       MsqError& err );

      /** \brief Set tag values on vertices
       * 
       * Set the value of a tag for a list of mesh vertices.
       * \param handle     The tag 
       * \param num_elems  Length of node_array
       * \param node_array Array of vertices for which to set the tag value.
       * \param tag_data   Tag data for each element, contiguous in memory.
       *                   This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_set_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       const void* tag_data,
                                       MsqError& err );
    
    
      /** \brief Get tag values on elements
       * 
       * Get the value of a tag for a list of mesh elements.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of elements for which to get the tag value.
       * \param tag_data   Return buffer in which to copy tag data, contiguous 
       *                   in memory.  This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_get_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       void* tag_data,
                                       MsqError& err );
    
      /** \brief Get tag values on vertices.
       * 
       * Get the value of a tag for a list of mesh vertices.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of vertices for which to get the tag value.
       * \param tag_data   Return buffer in which to copy tag data, contiguous 
       *                   in memory.  This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_get_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       void* tag_data,
                                       MsqError& err );

//**************** Memory Management ****************
      // Tells the mesh that the client is finished with a given
      // entity handle.  
    virtual void release_entity_handles( const EntityHandle *handle_array,
                                         size_t num_handles,
                                         MsqError &err );
    
      // Instead of deleting a Mesh when you think you are done,
      // call release().  In simple cases, the implementation could
      // just call the destructor.  More sophisticated implementations
      // may want to keep the Mesh object to live longer than Mesquite
      // is using it.
    virtual void release();

  private:
    Mesh* myMesh;
    void* gid_tag;
    void* pid_tag;
  };
}

#endif
