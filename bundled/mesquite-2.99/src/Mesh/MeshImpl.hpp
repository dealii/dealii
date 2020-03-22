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
  \file   MeshImpl.hpp
  \brief  

  \author Darryl Melander
  \author Thomas Leurent
  \author Jason Kraftcheck
  \date   2003-04-17
*/

#ifndef MESQUITE_MESH_IMPL_HPP
#define MESQUITE_MESH_IMPL_HPP

#include "MeshInterface.hpp"

#include <map>
#include <iosfwd>

namespace MESQUITE_NS
{
  class FileTokenizer;
  class MeshImplData;
  class MeshImplTags;
  struct TagDescription;
  
  /*!  \class MeshImpl

  \brief MeshImpl is a Mesquite implementation of the Mesh interface. 
  Applications can also provide their own implementation of the interface.
    
  MeshImpl can read in mesh files in VTK format and ExodusII format. 
  */
  class MESQUITE_EXPORT MeshImpl : public Mesquite::Mesh
  {
  public:
//********* Functions that are NOT inherited ************

    MeshImpl();
    
    virtual ~MeshImpl();

    /**\brief Initialize mesh by copying data from arrays.
     *
     *\param num_vertex  Number of vertices in the mesh
     *\param num_elem    Number of elements in the mesh
     *\param entity_topology Type of all elements in the mesh
     *\param fixed        Value of fixed flag for each vertex
     *\param coords      Interleaved vertex coordinates.
     *\param conn        Element corner vertices specified as indices into
     *                   vertex list
     */
    MeshImpl(int num_vertex, int num_elem, 
	     EntityTopology entity_topology, 
	     const bool *fixed, const double *coords, const int *conn);
    
    /**\brief Initialize mesh by copying data from arrays.
     *
     *\param num_vertex  Number of vertices in the mesh
     *\param num_elem    Number of elements in the mesh
     *\param entity_topologies The types of each element in the mesh
     *\param fixed       Value of fixed flag for each vertex
     *\param coords      Interleaved vertex coordinates.
     *\param conn        Element corner vertices specified as indices into
     *                   vertex list
     */
    MeshImpl(int num_vertex, int num_elem, 
	     const EntityTopology *element_topologies, 
	     const bool *fixed, const double *coords, const int *conn);
    
    /**\brief Read mesh from VTK file format version 3.0 or earlier */
    void read_vtk(const char* in_filename, Mesquite::MsqError &err);
    
    /**\brief Write mesh to VTK file format version 3.0 */
    void write_vtk(const char* out_filename, Mesquite::MsqError &err);
    
    /**\brief Read mesh from ExodusII file */
    void read_exodus(const char* in_filename, Mesquite::MsqError &err);
    
    /**\brief Write mesh to ExodusII file */
    void write_exodus(const char* out_filename, Mesquite::MsqError &err);
    
    /**\brief Set the value returned by vertices_get_fixed_flag for all vertices */
    void set_all_fixed_flags(bool value, MsqError& err);
    
    /**\brief Set the value returned by vertices_get_slaved_flag for all vertices
     *
     * Set value for vertices_get_slaved_flag for all corners nodes to false and
     * for all other nodes to the passed value.
     *
     *\param value Value to set on mid-edge, mid-face, or mid-region nodes.
     *
     *\NOTE You must change slave vertex mode to Settings::SLAVED_FLAG for the values
     *      returned by vertices_get_slaved_flag to be used during optimization.
     */
    void set_all_slaved_flags(bool value, MsqError& err);
    
    /**\brief Set values for vertices_get_fixed_flag and vertices_get_slaved_flag
     * on skin vertices
     *
     * Set flag values for vertices on the skin (i.e. boundary) of the mesh.
     * Does not modify flags on iterior vertices.  Call \c set_all_fixed _flags
     * and \c set_all_slaved_flags *before* calling this function to set values
     * on interior vertices.
     *\param corner_fixed_flag  Value for vertices_get_fixed_flag for vertices at element corners
     *\param midnode_fixed_flag  Value for vertices_get_fixed_flag for non-corner 
     *                           vertices (i.e. mid-face or mid-element nodes)
     *\param midnode_slaved_flag  Value for vertices_get_slaved_flag for non-corner 
     *                            vertices (i.e. mid-face or mid-element nodes)
     *
     *\NOTE You must change slave vertex mode to Settings::SLAVED_FLAG for the values
     *      returned by vertices_get_slaved_flag to be used during optimization.
     */
    void set_skin_flags( bool corner_fixed_flag, 
                         bool midnode_fixed_flag,
                         bool midnode_slaved_flag,
                         MsqError& err );
    
    
    /**\brief Find the vertices in the skin (i.e bounary) of the mesh and
     * mark them as 'fixed'.
     *\param clear_existing  If true only skin vertices are marked as fixed.
     *                       If false, skin vertices will be marked as fixed
     *                       in addition to any other vertices already marked
     *                       as fixed.
     */
    void mark_skin_fixed( MsqError& err, bool clear_existing = true );
                                     
//********* Functions that ARE inherited ************
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
    
      // Returns a pointer to an iterator that iterates over the
      // set of all vertices in this mesh.  The calling code should
      // delete the returned iterator when it is finished with it.
      // If vertices are added or removed from the Mesh after obtaining
      // an iterator, the behavior of that iterator is undefined.
    virtual VertexIterator* vertex_iterator(MsqError &err);
    
      // Returns a pointer to an iterator that iterates over the
      // set of all top-level elements in this mesh.  The calling code should
      // delete the returned iterator when it is finished with it.
      // If elements are added or removed from the Mesh after obtaining
      // an iterator, the behavior of that iterator is undefined.
    virtual ElementIterator* element_iterator(MsqError &err);

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
                                  const std::vector<bool>& fixed_flag_array,
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
    
      // Remove all data
    void clear();

  protected:    
    
    /** Coordinate values per vertex */
    int numCoords;
    
    MeshImplData* myMesh;
    MeshImplTags* myTags;

  private:

//**************** VTK Parsing ****************

      /** Read a data block from the file */
    void vtk_read_dataset( FileTokenizer& file, MsqError& err );
    
      /** Read structured point mesh */
    void vtk_read_structured_points( FileTokenizer& file, MsqError& err );
      /** Read structured grid mesh */
    void vtk_read_structured_grid  ( FileTokenizer& file, MsqError& err );
      /** Read rectilinear grid structured mesh */
    void vtk_read_rectilinear_grid ( FileTokenizer& file, MsqError& err );
      /** Read polydata mesh */
    void vtk_read_polydata         ( FileTokenizer& file, MsqError& err );
      /** Read unstructured mesh */
    void vtk_read_unstructured_grid( FileTokenizer& file, MsqError& err );
      /** Read file-level field data */
    void vtk_read_field            ( FileTokenizer& file, MsqError& err );
    
      /** Helper function for vtk_read_polydata() - reads polygon subsection */
    void vtk_read_polygons( FileTokenizer& file, MsqError& err );
      /** Helper function for readers of structured mesh - create elements */
    void vtk_create_structured_elems( const long* dims, MsqError& err );
    
      /** Read attribute data for vertices */
    void vtk_read_point_data( FileTokenizer& file, MsqError& err );
      /** Read attribute data for elements */
    void vtk_read_cell_data ( FileTokenizer& file, MsqError& err );
      /** Store data read in vtk_read_point_data into mesh tags */
    void vtk_store_point_data( const void* data, TagDescription& desc, MsqError& );
      /** Store data read in vtk_read_cell_data into mesh tags */
    void vtk_store_cell_data( const void* data, TagDescription& desc, MsqError& );
      /** Read actual data for both vtk_read_point_data() and vtk_read_cell_data() 
       *  Initializes all fields of passed TagDescription
       *\return NULL if field data, otherwise pointer to malloc'd data.
       */
    void* vtk_read_attrib_data( FileTokenizer& file, 
                                long num_data_to_read, 
                                TagDescription& tag_out,
                                MsqError& err );
      /** Read a 2-D array of data of the specified type from the file 
       *  Initializes size and type fields of passed TagDescroption */
    void* vtk_read_typed_data( FileTokenizer& file, int type,
                               size_t per_elem, size_t num_elem,
                               TagDescription& tag_out,
                               MsqError& err );
      /** Read field data 
       *\param count expected number of tuples, or zero if not known
       */
    void* vtk_read_field_data( FileTokenizer& file, size_t count,
                               size_t field_count,
                               const std::string& field_name, 
                               TagDescription& tag, 
                               MsqError& err );
    
      /** Read scalar attribute data  
       *  Initializes size and type fields of passed TagDescroption */
    void* vtk_read_scalar_attrib ( FileTokenizer& file, long count, 
                                   TagDescription& tag_out, MsqError& err );
      /** Read color attribute data  
       *  Initializes size and type fields of passed TagDescroption */
    void* vtk_read_color_attrib  ( FileTokenizer& file, long count, 
                                   TagDescription& tag_out, MsqError& err );
      /** Read vector or normal attribute data  
       *  Initializes size and type fields of passed TagDescroption */
    void* vtk_read_vector_attrib ( FileTokenizer& file, long count, 
                                   TagDescription& tag_out, MsqError& err );
      /** Read texture attribute data  
       *  Initializes size and type fields of passed TagDescroption */
    void* vtk_read_texture_attrib( FileTokenizer& file, long count, 
                                   TagDescription& tag_out, MsqError& err );
      /** Read tensor (3x3 matrix) data  
       *  Initializes size and type fields of passed TagDescroption */
    void* vtk_read_tensor_attrib ( FileTokenizer& file, long count, 
                                   TagDescription& tag_out, MsqError& err );

      /** Write tag data to VTK attributes */
    void vtk_write_attrib_data( std::ostream& file,
                                const TagDescription& desc,
                                const void* data, size_t count,
                                MsqError& err ) const;

      /** Convert tag data stored on vertices to boolean values.
       *  Deletes tag data.
       *  Passes back empty result vector if no tag.
       */
    void tag_to_bool( const char* tag_name, 
                      std::vector<bool>& vertex_vals,
                      MsqError& err  );

//**************** End VTK Parsing ****************
  };
}

#endif
