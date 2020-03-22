/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file ArrayMesh.hpp
 *  \brief Access mesh stored in arrays
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_ARRAY_MESH_HPP
#define MSQ_ARRAY_MESH_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"

namespace MESQUITE_NS {

class ArrayMesh : public Mesh
{
  public:
  
      /** Create a Mesquite::Mesh instance that wraps application-provided
       *  arrays.  
       *
       * Note:  An instance of this class will reference the passed in 
       *        arrays.  It will not copy them.  
       *
       *\param coords_per_vertex Dimension of the mesh (2 or 3)
       *\param num_vertices      Number of vertices in the mesh
       *\param interleaved_vertex_coords Vertex coordinates.  Ordered as
       *                         [x1, y1, z1, x2, y2, z2, ...]
       *\param vertex_fixed_flags One value per vertex.  Zero if vertex is
       *                         free, one if the position is fixed.
       *\param num_elements      Number of elements in the mesh
       *\param element_type      The type of the elements
       *\param element_connectivity_array Element connectivity, specified
       *                         as vertex indices such that the location
       *                         of the vertex coordinates in vertex_coords
       *                         is at 3 times the value in this array.
       *\param one_based_conn_indices Use one-based (Fortran) array indexing.
       *\param nodes_per_element Number of nodes in each element.  If not
       *                         specified, number of nodes in a linear
       *                         element with the type 'element_type' is
       *                         assumed.
       *\param vertex_slaved_flags One value per vertex.  Zero if vertex is
       *                         free, one if the vertex is slaved to the
       *                         logical position of the element mapping/shape
       *                         function..
       */
    ArrayMesh( int coords_per_vertex,
               unsigned long num_vertices,
               double* interleaved_vertex_coords,
               const int* vertex_fixed_flags,
               unsigned long num_elements,
               EntityTopology element_type,
               const unsigned long* element_connectivity_array,
               bool one_based_conn_indices = false,
               unsigned nodes_per_element = 0,
               const int* vertex_slaved_flags = 0 );
               
      /** Create a Mesquite::Mesh instance that wraps application-provided
       *  arrays.  
       *
       * Note:  An instance of this class will reference the passed in 
       *        arrays.  It will not copy them.  
       *
       *\param coords_per_vertex Dimension of the mesh (2 or 3)
       *\param num_vertices      Number of vertices in the mesh
       *\param interleaved_vertex_coords Vertex coordinates.  Ordered as
       *                         [x1, y1, z1, x2, y2, z2, ...]
       *\param vertex_fixed_flags One value per vertex.  Zero if vertex is
       *                         free, one if the poistion is fixed.
       *\param num_elements      Number of elements in the mesh
       *\param element_types     The topological type for each element.
       *\param element_connectivity_array Element connectivity, specified
       *                         as vertex indices such that the location
       *                         of the vertex coordinates in vertex_coords
       *                         is at 3 times the value in this array.
       *\param element_connectivity_offsets An optional array of length one greater
       *                         than num_elements.  Each entry other than
       *                         the last should contain the value of the
       *                         index into element_connectivity_array at
       *                         which the connectivity data for the 
       *                         corresponding element begins.  The last
       *                         entry in the array must be the next-to-last
       *                         entry plus the length of the connectivity
       *                         data for the last element, such that the
       *                         length of the connectivity data for any
       *                         element 'e' can be calculated with:
       *                           n = element_connectivity_offsets[e+1] - element_connectivity_offsets[e]
       *
       *                         If this array is not specified, then it
       *                         will be assumed that the length of the
       *                         connectivity data for each element is the
       *                         number of corners necessary to represent
       *                         its topological type (that it has no higher-order
       *                         nodes.)
       *\param one_based_conn_indices Use one-based (Fortran) array indexing.
       *\param vertex_slaved_flags One value per vertex.  Zero if vertex is
       *                         free, one if the vertex is slaved to the
       *                         logical position of the element mapping/shape
       *                         function..
       */
    ArrayMesh( int coords_per_vertex,
               unsigned long num_vertices,
               double* interleaved_vertex_coords,
               const int* vertex_fixed_flags,
               unsigned long num_elements,
               const EntityTopology* element_types,
               const unsigned long* element_connectivity_array,
               const unsigned long* element_connectivity_offsets = NULL,
               bool one_based_conn_indices = false,
               const int* vertex_slaved_flags = 0 );


    ArrayMesh();
    
    ~ArrayMesh();
    
    void clear_mesh();
    void set_mesh( int coords_per_vertex,
               unsigned long num_vertices,
               double* interleaved_vertex_coords,
               const int* vertex_fixed_flags,
               unsigned long num_elements,
               EntityTopology element_type,
               const unsigned long* element_connectivity_array,
               bool one_based_conn_indices = false,
               unsigned nodes_per_element = 0,
               const int* vertex_slaved_flags = 0 );
    
    /**\brief Give mesquite access to per-entity application data via a tag
     *
     * Allow mesquite to access application-owned per-element and/or per-vertex
     * data by assigning a tag name to said data.  Mesquite will be allowed
     * to read application-owned data.  Any attempt by Mesquite to modify 
     * the data will result in an internal error condition.
     *
     *\param tag_name  Tag name through which Mesquite can access the data
     *\param data_type The type of the data
     *\param vals_per_entity Number of values of type \c data_type that each entity has
     *                       Data is assumed to be interleaved such that all values
     *                       associated with a a single entity are adjacent in memory
     *\param vertex_data     A pointer to the region of memory in which per-vertex
     *                       values are stored.  May be NULL if per-vertex data is
     *                       not available.  Values are assumed to be in the same 
     *                       order as the vertex coordinates passed to the constructor
     *                       or \c set_mesh.
     *\param vertex_data     A pointer to the region of memory in which per-element
     *                       values are stored.  May be NULL if per-element data is
     *                       not available.  Values are assumed to be in the same 
     *                       order as the element connetivity passed to the constructor
     *                       or \c set_mesh.
     *\param default_value   Value to return for all vertices and/or all elements
     *                       if the cooresponding data array is null.  May be
     *                       NULL if no default value.
     */
    TagHandle add_read_only_tag_data( const char* tag_name,
                                      TagType data_type,
                                      int vals_per_entity,
                                      const void* vertex_data,
                                      const void* element_data,
                                      const void* default_value,
                                      MsqError& err );
    
    /**\brief Give mesquite access to per-entity application data via a tag
     *
     * Allow mesquite to access and/or set application-owned per-element 
     * and/or per-vertex data by assigning a tag name to said data.  
     *
     *\param tag_name  Tag name through which Mesquite can access the data
     *\param data_type The type of the data
     *\param vals_per_entity Number of values of type \c data_type that each entity has
     *                       Data is assumed to be interleaved such that all values
     *                       associated with a a single entity are adjacent in memory
     *\param vertex_data     A pointer to the region of memory in which per-vertex
     *                       values are stored.  May be NULL if per-vertex data is
     *                       not available.  Values are assumed to be in the same 
     *                       order as the vertex coordinates passed to the constructor
     *                       or \c set_mesh.
     *\param vertex_data     A pointer to the region of memory in which per-element
     *                       values are stored.  May be NULL if per-element data is
     *                       not available.  Values are assumed to be in the same 
     *                       order as the element connetivity passed to the constructor
     *                       or \c set_mesh.
     *\param default_value   Value to return for all vertices and/or all elements
     *                       if the cooresponding data array is null.  May be
     *                       NULL if no default value.
     */
    TagHandle add_writable_tag_data( const char* tag_name,
                                     TagType tag_data_type,
                                     int vals_per_entity,
                                     void* vertex_data,
                                     void* element_data,
                                     const void* default_value,
                                     MsqError& err );
    
    
    virtual int get_geometric_dimension( MsqError& err );

    virtual void get_all_elements( std::vector<ElementHandle>& elements,
                                   MsqError& err );
    virtual void get_all_vertices( std::vector<VertexHandle>& vertices,
                                   MsqError& err );
    
    virtual VertexIterator* vertex_iterator(MsqError &err);
    virtual ElementIterator* element_iterator(MsqError &err);

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
    
    virtual void vertices_get_attached_elements( 
                         const VertexHandle* vertex_array,
                         size_t num_vertex,
                         std::vector<ElementHandle>& elements,
                         std::vector<size_t>& offsets,
                         MsqError& err );
    
    virtual void elements_get_attached_vertices(
                                   const ElementHandle *elem_handles,
                                   size_t num_elems,
                                   std::vector<VertexHandle>& vert_handles,
                                   std::vector<size_t>& offsets, 
                                   MsqError &err);
    
    virtual void elements_get_topologies(const ElementHandle *element_handle_array,
                                         EntityTopology *element_topologies,
                                         size_t num_elements, MsqError &err);

    
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

    
    
    virtual void release_entity_handles( const EntityHandle *handle_array,
                                         size_t num_handles, 
                                         MsqError &err);
    
    virtual void release();

  private:
    
    inline const unsigned long* elem_verts( size_t elem_index, int& num_vertex ) const;
    
    void build_vertex_adjacency_list();
    
    int mDimension;                  //!< Coordinates per vertex
    unsigned long vertexCount;       //!< Number of vertices
    double* coordArray;              //!< Interleaved vertex coordinates
    const int* fixedFlags;           //!< Vertex fixed flags
    const int* slavedFlags;          //!< Vertex slaved flags
    unsigned char* vertexByteArray;  //!< Vertex bytes
    
    unsigned long elementCount;      //!< Number of elements
    const unsigned long* connArray;  //!< Element connectivity
    const unsigned long* connOffsets;//!< Offsets into connectivity array
                                     //!< for each element.  If NULL, then
                                     //!< all elements are of the same type
                                     //!< and have the same number of vertices.
    unsigned long* allocConnOffsets; //!< Same as connOffsets if allocated
                                     //!< by constructor.  NULL if connOffsets
                                     //!< is either NULL or application-provided
                                     //!< data.
    EntityTopology elementType;      //!< Type for all elements if connOffsets is NULL
    const EntityTopology* elementTypes; //!< Type for each element type if connOffsets is not NULL
    unsigned nodesPerElement;        //!< Nodes per element if connOffsets is NULL
    bool oneBasedArrays;             //!< FORTRAN-style array indexing
    
    unsigned long* vertexAdjacencyList;
    unsigned long* vertexAdjacencyOffsets;
    
    static unsigned bytes( TagType type );
    
    static void fill( unsigned char* buffer, const unsigned char* value, size_t size, size_t count );
    
    struct Tag {
      char* name;             //!< Tag name (null-terminated string)
      TagType type;          //!< tag data type
      unsigned size;         //!< number of *bytes* per entity
      bool owned;            //!< true if memory for tag storage is owned by this class
      unsigned char* vtxWritePtr;     //!< pointer to writable tag data (NULL if read-only)
      const unsigned char* vtxReadPtr;//!< pointer to tag data
      unsigned char* eleWritePtr;     //!< pointer to writable tag data (NULL if read-only)
      const unsigned char* eleReadPtr;//!< pointer to tag data
      unsigned char* defaultValue;    //!< Default value
      Tag* next;             //!< Linked-list next pointer
    };
    
    Tag* allocate_tag( const char* name, 
                       bool owned,
                       TagType type, 
                       unsigned size, 
                       const void* vertex_ro_data,
                       void* vertex_rw_data,
                       const void* element_ro_data,
                       void* element_rw_data,
                       const void* default_value,
                       MsqError& err );
    
    Tag* tagList;
    
    bool valid() const;
};



} // namespace Mesquite

#endif
