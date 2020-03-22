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

#ifndef MESQUITE_MESH_IMPL_TAGS_HPP
#define MESQUITE_MESH_IMPL_TAGS_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"

#include <vector>

namespace MESQUITE_NS {

  struct TagDescription {
      // The VTK attribute type the data was read from, or NONE if not VTK data.
      // This property is kept only so that when a vtk file is read and subseqquently
      // written, it can be preserved for other potential readers that actually care
      // about it.
    enum VtkType { NONE = 0, SCALAR, COLOR, VECTOR, NORMAL, TEXTURE, TENSOR, FIELD };

    std::string name;   //!< Tag name 
    Mesh::TagType type;     //!< Tag data type
    VtkType vtkType;        //!< Attribute type from VTK file
    size_t size;            //!< Size of tag data (sizeof(type)*array_length)
    std::string member; //!< Field member name for 1-member fields.

    inline TagDescription( std::string n, 
                           Mesh::TagType t, 
                           VtkType v, 
                           size_t s,
                           std::string m )
      : name(n), type(t), vtkType(v), size(s), member(m) {}

    inline TagDescription( )
      : type(Mesh::BYTE), vtkType(NONE), size(0) {}

    inline bool operator==(const TagDescription& o) const
      { return name == o.name && type == o.type && vtkType == o.vtkType && size == o.size; }
    inline bool operator!=(const TagDescription& o) const
      { return name != o.name || type != o.type || vtkType != o.vtkType || size != o.size; }
  };

/**\class MeshImplTags
 *
 * Store tags and tag data for Mesquite's native mesh representation.
 * Stores for each tag: properties, element data, and vertex data.
 * The tag element and vertex data sets are maps between some element
 * or vertex index and a tag value.
 */
class MeshImplTags {
  public:

  ~MeshImplTags() { clear(); }

  /** \class TagData
   * Store data for a single tag
   */
  struct TagData  {

      //! tag meta data
    const TagDescription desc;

      //! per-element data, or NULL if none has been set.
    void* elementData;

      //! number of entries in elementData
    size_t elementCount;

      //! per-vertex data, or NULL if none has been set.
    void* vertexData;
    
      //! number of entries in vertexData
    size_t vertexCount;

      //! Default value for tag
    void* defaultValue;

      /** \brief Construct tag
       *\param name Tag name
       *\param type Tag data type
       *\param length Tag array length (1 for scalar/non-array)
       *\param default_val Default value for tag
       *\param vtk_type Attribute type in VTK file
       */
    inline TagData( const std::string& name, 
             Mesh::TagType type, unsigned length,
             void* default_val = 0,
             TagDescription::VtkType vtk_type = TagDescription::NONE,
             const std::string& field_member = "")
      : desc(name, type, vtk_type, length*size_from_tag_type(type), field_member),
        elementData(0), elementCount(0),
        vertexData(0), vertexCount(0), 
        defaultValue(default_val) {}

      /** \brief Construct tag
       *\param desc Tag description object
       */
    inline TagData( const TagDescription& descr )
      : desc(descr), elementData(0), elementCount(0),
        vertexData(0), vertexCount(0),
        defaultValue(0) {}

    ~TagData();
  };
  
    /** \brief Get the size of the passed data type */
  static size_t size_from_tag_type( Mesh::TagType type );

    /** \brief Clear all data */
  void clear();

    /** \brief Get tag index from name */
  size_t handle( const std::string& name, MsqError& err ) const;
  
    /** \brief Get tag properties */
  const TagDescription& properties( size_t tag_handle, MsqError& err ) const;
  
    /** \brief Create a new tag 
     * 
     * Create a new tag with the passed properties
     *\param name Tag name (must be unique)
     *\param type Tag data type
     *\param length Number of values in tag (array length, 1 for scalar)
     *\param defval Optional default value for tag
     */
  size_t create( const std::string& name,
                 Mesh::TagType type,
                 unsigned length,
                 const void* defval,
                 MsqError& err );
                 
    /** \brief Create a new tag 
     * 
     * Create a new tag with the passed properties
     */
  size_t create( const TagDescription& desc,
                 const void* defval,
                 MsqError& err );

    /**\brief Remove a tag */
  void destroy( size_t tag_index, MsqError& err );
  
    /**\brief Set tag data on elements */
  void set_element_data( size_t tag_handle,
                         size_t num_indices,
                         const size_t* elem_indices,
                         const void* tag_data,
                         MsqError& err );
  
    /**\brief Set tag data on vertices */
  void set_vertex_data( size_t tag_handle,
                        size_t num_indices,
                        const size_t* elem_indices,
                        const void* tag_data,
                        MsqError& err );

    /**\brief Get tag data on elements */
  void get_element_data( size_t tag_handle,
                         size_t num_indices,
                         const size_t* elem_indices,
                         void* tag_data,
                         MsqError& err ) const;
  
    /**\brief Get tag data on vertices */
  void get_vertex_data( size_t tag_handle,
                        size_t num_indices,
                        const size_t* elem_indices,
                        void* tag_data,
                        MsqError& err ) const;
  
  /**\class TagIterator
   *
   * Iterate over list of valid tag handles 
   */
  class TagIterator 
  {
    public:
      TagIterator() : tags(0), index(0) {}
      TagIterator( MeshImplTags* d, size_t i ) : tags(d), index(i) {}
      size_t operator*() const  { return index+1; }
      TagIterator operator++();
      TagIterator operator--();
      TagIterator operator++(int);
      TagIterator operator--(int);
      bool operator==(TagIterator other) const { return index == other.index; }
      bool operator!=(TagIterator other) const { return index != other.index; }
    private:
      MeshImplTags* tags;
      size_t index;
  };
  TagIterator tag_begin();
  TagIterator tag_end()   { return TagIterator(this,tagList.size()); }

    /**\brief Check if any vertices have tag */
  bool tag_has_vertex_data( size_t index, MsqError& err ) ;
    /**\brief Check if any elements have tag */
  bool tag_has_element_data( size_t index, MsqError& err ) ;
  
  private:

  friend class MeshImplTags::TagIterator;  
  
    std::vector<TagData*> tagList;
}; // class MeshImplTags

} // namespace Mesquite

#endif

    
