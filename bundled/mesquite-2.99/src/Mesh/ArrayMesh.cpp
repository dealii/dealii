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


/** \file ArrayMesh.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ArrayMesh.hpp"
#include "TopologyInfo.hpp"
#include "MsqError.hpp"
#include "MsqVertex.hpp"

#include <iostream>

namespace MESQUITE_NS {

class IndexIterator : public EntityIterator
{
public:
  IndexIterator( size_t mStart, size_t mEnd )
    : mStart(mStart), mEnd(mEnd), mCurrent(mStart) {}
  virtual ~IndexIterator() {}
  virtual void restart() { mCurrent = mStart; }
  virtual Mesh::EntityHandle operator*() const
    { return (Mesh::EntityHandle)mCurrent; }
  virtual void operator++() { ++mCurrent; }
  virtual bool is_at_end() const { return mEnd - mCurrent <= 1; }
private:
  size_t mStart, mEnd, mCurrent;
};

ArrayMesh::ArrayMesh()
  : mDimension( 0 ),
    vertexCount( 0 ),
    coordArray( 0 ),
    fixedFlags( 0 ),
    slavedFlags( 0 ),
    vertexByteArray( 0 ),
    elementCount( 0 ),
    connArray( 0 ),
    connOffsets( 0 ),
    allocConnOffsets( 0 ),
    elementType( MIXED ),
    elementTypes( 0 ),
    nodesPerElement( 0 ),
    oneBasedArrays( false ),
    vertexAdjacencyList(0),
    vertexAdjacencyOffsets(0),
    tagList(0)
{}

ArrayMesh::ArrayMesh( int coords_per_vertex,
                      unsigned long num_vertices,
                      double* interleaved_vertex_coords,
                      const int* vertex_fixed_flags,
                      unsigned long num_elements,
                      EntityTopology element_type,
                      const unsigned long* element_connectivity_array,
                      bool one_based_conn_indices,
                      unsigned nodes_per_element,
                      const int* vertex_slaved_flags ) 
  : mDimension( coords_per_vertex ),
    vertexCount( num_vertices ),
    coordArray( interleaved_vertex_coords ),
    fixedFlags( vertex_fixed_flags ),
    slavedFlags( vertex_slaved_flags ),
    vertexByteArray( new unsigned char[num_vertices + one_based_conn_indices] ),
    elementCount( num_elements ),
    connArray( element_connectivity_array ),
    connOffsets( 0 ),
    allocConnOffsets( 0 ),
    elementType( element_type ),
    elementTypes( 0 ),
    nodesPerElement( nodes_per_element ),
    oneBasedArrays( one_based_conn_indices ),
    vertexAdjacencyList(0),
    vertexAdjacencyOffsets(0),
    tagList(0)
{
  if (oneBasedArrays) {
    coordArray -= mDimension;
    --fixedFlags;
  }
  
  if (nodesPerElement < 2)
    nodesPerElement = TopologyInfo::corners( element_type );
    
  assert(valid());
  memset( vertexByteArray, 0, num_vertices + one_based_conn_indices );
}

ArrayMesh::ArrayMesh( int coords_per_vertex,
                      unsigned long num_vertices,
                      double* interleaved_vertex_coords,
                      const int* vertex_fixed_flags,
                      unsigned long num_elements,
                      const EntityTopology* element_types,
                      const unsigned long* element_connectivity_array,
                      const unsigned long* element_connectivity_offsets,
                      bool one_based_conn_indices,
                      const int* vertex_slaved_flags ) 
  : mDimension( coords_per_vertex ),
    vertexCount( num_vertices ),
    coordArray( interleaved_vertex_coords ),
    fixedFlags( vertex_fixed_flags ),
    slavedFlags( vertex_slaved_flags ),
    vertexByteArray( new unsigned char[num_vertices + one_based_conn_indices] ),
    elementCount( num_elements ),
    connArray( element_connectivity_array ),
    connOffsets( element_connectivity_offsets ),
    allocConnOffsets( 0 ),
    elementType( MIXED ),
    elementTypes( element_types ),
    nodesPerElement( 0 ),
    oneBasedArrays( one_based_conn_indices ),
    vertexAdjacencyList(0),
    vertexAdjacencyOffsets(0),
    tagList(0)
{
  if (oneBasedArrays) {
    coordArray -= mDimension;
    --fixedFlags;
    if (element_connectivity_offsets)
      --connArray;
  }
  
  if (!element_connectivity_offsets) {
    connOffsets = allocConnOffsets = new unsigned long[num_elements+1];
    allocConnOffsets[0] = 0;
    for (unsigned long i = 1; i <= num_elements; ++i)
      allocConnOffsets[i] = allocConnOffsets[i-1] + TopologyInfo::corners( elementTypes[i-1] );
  }
  
  assert(valid());
  memset( vertexByteArray, 0, num_vertices + one_based_conn_indices );
}

bool ArrayMesh::valid() const
{
  unsigned long off = oneBasedArrays ? 1 : 0;
  for (unsigned long i = off; i < vertexCount+off; ++i) {
    if (fixedFlags[i] != 0 && fixedFlags[i] != 1) {
      std::cerr << "Invalid vertex fixed flag at index " << i << std::endl;
      return false;
    }
  }
  
  for (unsigned long i = 0; i < elementCount * nodesPerElement; ++i) {
    unsigned long j = connArray[i] - oneBasedArrays;
    if (j >= vertexCount) {
      std::cerr << "Invalid connectivity index at index " << j 
                << "(element " << j/elementCount << " node " << j%elementCount 
                << ')' << std::endl;
      return false;
    }
  }
  
  return true;
}

void ArrayMesh::clear_mesh()
{
  delete [] vertexByteArray;
  delete [] vertexAdjacencyList;
  delete [] vertexAdjacencyOffsets;
  delete [] allocConnOffsets;
  while (tagList) {
    Tag* dead = tagList;
    tagList = tagList->next;
    delete [] dead->name;
    delete [] dead->defaultValue;
    if (dead->owned) {
      delete [] dead->vtxWritePtr;
      delete [] dead->eleWritePtr;
    }
    delete dead;
  }
  mDimension = 0;
  vertexCount = 0;
  coordArray = 0;
  connOffsets = 0;
  allocConnOffsets = 0;
  fixedFlags = 0;
  vertexByteArray = 0;
  elementCount = 0;
  elementType = MIXED;
  elementTypes = 0;
  nodesPerElement = 0;
  oneBasedArrays = false;
  vertexAdjacencyList = 0;
  vertexAdjacencyOffsets = 0;
}

void ArrayMesh::set_mesh( int coords_per_vertex,
                          unsigned long num_vertices,
                          double* interleaved_vertex_coords,
                          const int* vertex_fixed_flags,
                          unsigned long num_elements,
                          EntityTopology element_type,
                          const unsigned long* element_connectivity_array,
                          bool one_based_conn_indices,
                          unsigned nodes_per_element,
                          const int* vertex_slaved_flags ) 
{
  clear_mesh();
  mDimension = coords_per_vertex;
  vertexCount = num_vertices;
  coordArray = interleaved_vertex_coords;
  fixedFlags = vertex_fixed_flags;
  slavedFlags = vertex_slaved_flags;
  elementCount = num_elements;
  connArray = element_connectivity_array;
  elementType = element_type;
  oneBasedArrays = one_based_conn_indices;
  
  if (oneBasedArrays) {
    coordArray -= mDimension;
    --fixedFlags;
  }
  
  if (nodes_per_element < 2)
    nodesPerElement = TopologyInfo::corners( element_type );
  else
    nodesPerElement = nodes_per_element;
    
  vertexByteArray = new unsigned char[num_vertices + one_based_conn_indices];
  assert(valid());
  memset( vertexByteArray, 0, num_vertices + one_based_conn_indices );
}
  
  

ArrayMesh::~ArrayMesh()
{
  clear_mesh();
}

inline const unsigned long* ArrayMesh::elem_verts( size_t e, int& n ) const
{
  assert( e < elementCount );
  if (connOffsets) {
    n = connOffsets[e+1] - connOffsets[e];
    return connArray + connOffsets[e];
  }
  else {
    n = nodesPerElement;
    return connArray + nodesPerElement*e;
  }
}

int ArrayMesh::get_geometric_dimension( MsqError& )
  { return mDimension; }

void ArrayMesh::get_all_elements( std::vector<ElementHandle>& elements, MsqError& )
{
  elements.resize( elementCount );
  for (unsigned long i = 0; i < elementCount; ++i)
    elements[i] = (Mesh::ElementHandle)i;
}

void ArrayMesh::get_all_vertices( std::vector<VertexHandle>& vertices, MsqError& )
{
  vertices.resize( vertexCount );
  for (unsigned long i = 0; i < vertexCount; ++i)
    vertices[i] = (Mesh::VertexHandle)(i+oneBasedArrays);
}

VertexIterator* ArrayMesh::vertex_iterator( MsqError& )
  { return new IndexIterator( oneBasedArrays, vertexCount + oneBasedArrays ); }

ElementIterator* ArrayMesh::element_iterator( MsqError& err )
  { return new IndexIterator( 0, elementCount ); }

void ArrayMesh::vertices_get_fixed_flag( const VertexHandle vert_array[], 
                                         std::vector<bool>& fixed_flag_array,
                                         size_t num_vtx, 
                                         MsqError & )
{
  fixed_flag_array.resize( num_vtx );
  const size_t* indices = (const size_t*)vert_array;
  for (size_t i = 0; i < num_vtx; ++i) {
    assert(indices[i] < vertexCount);
    fixed_flag_array[i] = !!fixedFlags[indices[i]];
  }
}

void ArrayMesh::vertices_get_slaved_flag( const VertexHandle* vert_array, 
                                          std::vector<bool>& slaved_flags,
                                          size_t num_vtx, 
                                          MsqError &err )
{
  if (!slavedFlags) {
    MSQ_SETERR(err)("No data provided to ArrayMesh for Settings::SLAVE_FLAG", 
                    MsqError::INVALID_STATE);
    return ;
  }
  
  slaved_flags.resize( num_vtx );
  const size_t* indices = (const size_t*)vert_array;
  for (size_t i = 0; i < num_vtx; ++i) {
    assert(indices[i] < vertexCount);
    slaved_flags[i] = !!slavedFlags[indices[i]];
  }
}

void ArrayMesh::vertices_get_coordinates( const VertexHandle vert_array[],
                                           MsqVertex* coordinates,
                                           size_t num_vtx,
                                           MsqError &err )
{
  const size_t* indices = (const size_t*)vert_array;
  if (mDimension == 3) 
    for (size_t i = 0; i < num_vtx; ++i) {
      assert( indices[i] < vertexCount+oneBasedArrays );
      coordinates[i].set( coordArray+3*indices[i] );
    }
  else if (mDimension == 2) 
    for (size_t i =0; i < num_vtx; ++i) {
      assert( indices[i] < vertexCount+oneBasedArrays );
      coordinates[i].set( coordArray[2*indices[i]], coordArray[2*indices[i]+1], 0.0 );
    }
  else
    MSQ_SETERR(err)(MsqError::INVALID_STATE);
}

void ArrayMesh::vertex_set_coordinates( VertexHandle vert,
                                          const Vector3D& coordinates,
                                           MsqError &err )
{
  size_t i = (size_t)vert;
  assert( i < vertexCount+oneBasedArrays );
  if (mDimension == 3) 
    coordinates.get_coordinates(coordArray+3*i);
  else if (mDimension == 2) {
    coordArray[2*i] = coordinates[0];
    coordArray[2*i+1] = coordinates[1];
  }
  else
    MSQ_SETERR(err)(MsqError::INVALID_STATE);
}

void ArrayMesh::vertex_set_byte( VertexHandle vertex, unsigned char byte, MsqError &)
{
  assert( (size_t)vertex < vertexCount+oneBasedArrays );
  vertexByteArray[(size_t)vertex] = byte;
}

void ArrayMesh::vertices_set_byte( const VertexHandle *vert_array,
                                   const unsigned char *byte_array,
                                   size_t array_size, 
                                   MsqError &err )
{
  const size_t* indices = (const size_t*)vert_array;
  for (size_t i = 0; i < array_size; ++i) {
    assert( indices[i] < vertexCount+oneBasedArrays );
    vertexByteArray[indices[i]] = byte_array[i];
  }
}

void ArrayMesh::vertex_get_byte( VertexHandle vertex, unsigned char* byte, MsqError &)
{
  assert( (size_t)vertex < vertexCount+oneBasedArrays );
  *byte = vertexByteArray[(size_t)vertex];
}

void ArrayMesh::vertices_get_byte( const VertexHandle *vert_array,
                                   unsigned char *byte_array,
                                   size_t array_size, 
                                   MsqError &err )
{
  const size_t* indices = (const size_t*)vert_array;
  for (size_t i = 0; i < array_size; ++i) {
    assert( indices[i] < vertexCount+oneBasedArrays );
    byte_array[i] = vertexByteArray[indices[i]];
  }
}

void ArrayMesh::vertices_get_attached_elements( 
                         const VertexHandle* vertex_array,
                         size_t num_vertex,
                         std::vector<ElementHandle>& elements,
                         std::vector<size_t>& offsets,
                         MsqError& )
{
  const size_t* indices = (const size_t*)vertex_array;
  if (!vertexAdjacencyList)
    build_vertex_adjacency_list();
  
  elements.clear();
  offsets.resize( num_vertex + 1 );
  for (size_t i = 0; i < num_vertex; ++i) {
    offsets[i] = elements.size();
    assert( indices[i] < vertexCount+oneBasedArrays );
    for (size_t j = vertexAdjacencyOffsets[indices[i]];
         j < vertexAdjacencyOffsets[indices[i]+1]; ++j)
      elements.push_back( (ElementHandle)vertexAdjacencyList[j] );
  }
  offsets[num_vertex] = elements.size();
}

void ArrayMesh::elements_get_attached_vertices(
                                   const ElementHandle *elem_handles,
                                   size_t num_elems,
                                   std::vector<VertexHandle>& vert_handles,
                                   std::vector<size_t>& offsets, 
                                   MsqError &)
{
  const size_t* indices = (const size_t*)elem_handles;
  offsets.resize( num_elems + 1);
  vert_handles.clear();
  for (size_t i = 0; i < num_elems; ++i) {
    assert( indices[i] < elementCount );
    int count;
    const unsigned long* conn = elem_verts( indices[i], count );
    size_t prev_size = vert_handles.size();
    offsets[i] = prev_size;
    vert_handles.resize( prev_size + count );
    std::copy( conn, conn+count, (size_t*)(&vert_handles[prev_size]) );
  }
  offsets[num_elems] = vert_handles.size();
}

void ArrayMesh::elements_get_topologies( const ElementHandle *handles,
                                         EntityTopology *element_topologies,
                                         size_t num_elements, MsqError& )
{
  const size_t* indices = (const size_t*)handles;
  if (elementType == MIXED) 
    for (size_t i = 0; i < num_elements; ++i) {
      assert( indices[i] < elementCount );
      element_topologies[i] = elementTypes[indices[i]];
    }
  else 
    for (size_t i = 0; i < num_elements; ++i) {
      assert( indices[i] < elementCount );
      element_topologies[i] = elementType;
    }
}

void ArrayMesh::release_entity_handles( const EntityHandle*, size_t, MsqError& err )
{
  MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED);
}

void ArrayMesh::release() {}

void ArrayMesh::build_vertex_adjacency_list()
{
  delete [] vertexAdjacencyList;
  delete [] vertexAdjacencyOffsets;
  vertexAdjacencyOffsets = new unsigned long[vertexCount+oneBasedArrays+1];
  
    // for each vertex, store the number of elements the previous
    // vertex occurs in.
  memset( vertexAdjacencyOffsets, 0, sizeof(unsigned long)*(vertexCount+oneBasedArrays+1) );
  for (size_t i = 0; i < elementCount; ++i) {
    int n;
    const unsigned long* conn = elem_verts( i, n );
    for (int j = 0; j < n; ++j)
      ++vertexAdjacencyOffsets[conn[j]+1];
  }
  
    // convert vertexAdjacencyOffsets from a shifted list of counts
    // to a list of offsts
  for (size_t i = 1; i <= vertexCount+oneBasedArrays; ++i)
    vertexAdjacencyOffsets[i] += vertexAdjacencyOffsets[i-1];
    
    // allocate space and populate with reverse connectivity
  vertexAdjacencyList = new unsigned long[vertexAdjacencyOffsets[vertexCount+oneBasedArrays]];
  for (size_t i = 0; i < elementCount; ++i) {
    int n;
    const unsigned long* conn = elem_verts( i, n );
    for (int j = 0; j < n; ++j)
      vertexAdjacencyList[vertexAdjacencyOffsets[conn[j]]++] = i;
  }
  
  for (size_t i = vertexCount+oneBasedArrays; i > 0; --i)
    vertexAdjacencyOffsets[i] = vertexAdjacencyOffsets[i-1];
  vertexAdjacencyOffsets[0] = 0; 
}

unsigned ArrayMesh::bytes( TagType type ) {
  switch (type) {
    case BYTE: case BOOL: return 1;
    case INT: return sizeof(int);
    case DOUBLE: return sizeof(double);
    case HANDLE: return sizeof(EntityHandle);
  }
  return 0;
}

void ArrayMesh::fill( unsigned char* buffer,
                      const unsigned char* value,
                      size_t size,
                      size_t count )
{
  if (!value) {
    memset( buffer, 0, size*count );
    return;
  }
  
  unsigned char* const end = buffer + size*count;
  for (unsigned char* iter = buffer; iter != end; iter += size)
    memcpy( iter, value, size );
}

ArrayMesh::Tag* ArrayMesh::allocate_tag( const char* name, 
                                         bool owned,
                                         TagType type, 
                                         unsigned size,
                                         const void* vertex_ro_data,
                                         void* vertex_rw_data,
                                         const void* element_ro_data,
                                         void* element_rw_data,
                                         const void* default_value,
                                         MsqError& err )
{
    // check if name is already in use
  for (Tag* iter = tagList; iter; iter = iter->next) {
    if (!strcmp( iter->name, name )) {
      MSQ_SETERR(err)(MsqError::TAG_ALREADY_EXISTS);
      return 0;
    }
  }
    
    // allocate object
  Tag* result = new Tag;
  
    // initialize members
  result->type = type;
  result->size = size*bytes(type);
  result->owned = owned;
  result->name = new char[strlen(name)+1];
  strcpy( result->name, name );

  result->vtxWritePtr = reinterpret_cast<unsigned char*>(vertex_rw_data);
  if (vertex_rw_data) 
    result->vtxReadPtr = reinterpret_cast<unsigned char*>(vertex_rw_data);
  else
    result->vtxReadPtr = reinterpret_cast<const unsigned char*>(vertex_ro_data);

  result->eleWritePtr = reinterpret_cast<unsigned char*>(element_rw_data);
  if (element_rw_data) 
    result->eleReadPtr = reinterpret_cast<unsigned char*>(element_rw_data);
  else
    result->eleReadPtr = reinterpret_cast<const unsigned char*>(element_ro_data);
  
  if (default_value) {
    result->defaultValue = new unsigned char[result->size];
    memcpy( result->defaultValue, default_value, result->size );
  }
  else {
    result->defaultValue = 0;
  }
  
    // prepend to tag list
  result->next = tagList;
  tagList = result;
  
  return result;
}

TagHandle ArrayMesh::add_read_only_tag_data( const char* tag_name,
                                             TagType data_type,
                                             int vals_per_entity,
                                             const void* vertex_data,
                                             const void* element_data,
                                             const void* default_value,
                                             MsqError& err )
{
  Tag* tag = allocate_tag( tag_name, false,
                           data_type, vals_per_entity, 
                           vertex_data, 0,
                           element_data, 0,
                           default_value, err );
  MSQ_ERRZERO(err);
  return reinterpret_cast<TagHandle>(tag);
}

TagHandle ArrayMesh::add_writable_tag_data( const char* tag_name,
                                            TagType data_type,
                                            int vals_per_entity,
                                            void* vertex_data,
                                            void* element_data,
                                            const void* default_value,
                                            MsqError& err )
{
  Tag* tag = allocate_tag( tag_name, false,
                           data_type, vals_per_entity, 
                           0, vertex_data,
                           0, element_data,
                           default_value, err );
  MSQ_ERRZERO(err);
  return reinterpret_cast<TagHandle>(tag);
}
  
TagHandle ArrayMesh::tag_create( const std::string& tag_name, 
                                 TagType data_type, 
                                 unsigned size, 
                                 const void* default_value, 
                                 MsqError &err)
{
  Tag* tag = allocate_tag( tag_name.c_str(), true,
                           data_type, size, 
                           0, 0,
                           0, 0,
                           default_value, err );
  MSQ_ERRZERO(err);
  return reinterpret_cast<TagHandle>(tag);
}

void ArrayMesh::tag_delete( TagHandle handle, MsqError& err )
{
  Tag* ptr = reinterpret_cast<Tag*>(handle);
  // find previous tag pointer in list
  if (!tagList) {
    MSQ_SETERR(err)("Invalid tag handle", MsqError::TAG_NOT_FOUND );
    return;
  }
  Tag* prev = 0;
  for (prev = tagList; prev && prev->next != ptr; prev = prev->next);
  if (!prev && tagList != ptr) {
    MSQ_SETERR(err)("Invalid tag handle", MsqError::TAG_NOT_FOUND );
    return;
  }
  delete [] ptr->name;
  delete [] ptr->defaultValue;
  if (ptr->owned) {
    delete [] ptr->vtxWritePtr;
    delete [] ptr->eleWritePtr;
  }
  if (prev)
    prev->next = ptr->next;
  else
    tagList = ptr->next;
  delete ptr;
}


TagHandle ArrayMesh::tag_get( const std::string& name, MsqError& err )
{
  for (Tag* iter = tagList; iter; iter = iter->next)
    if (name == iter->name)
      return reinterpret_cast<TagHandle>(iter);
  MSQ_SETERR(err)(MsqError::TAG_NOT_FOUND);
  return 0;
}

void ArrayMesh::tag_properties( TagHandle handle, 
                                std::string& name, 
                                TagType& type, 
                                unsigned& size, 
                                MsqError& err )
{
  const Tag* ptr = reinterpret_cast<const Tag*>(handle);
  name = ptr->name;
  type = ptr->type;
  size = ptr->size/bytes(ptr->type);
}

void ArrayMesh::tag_set_element_data( TagHandle handle, 
                                      size_t count, 
                                      const ElementHandle* entities, 
                                      const void* data, 
                                      MsqError& err )
{
  Tag* tag = reinterpret_cast<Tag*>(handle);
  if (!tag->eleWritePtr) {
    if (!tag->owned) {
      MSQ_SETERR(err)("Attempt to set non-writeable (application owned) tag data", 
                      MsqError::TAG_ALREADY_EXISTS);
      return;
    }
    else {
      assert(!tag->eleReadPtr);
      tag->eleReadPtr = tag->eleWritePtr = new unsigned char[elementCount*tag->size];
      fill( tag->eleWritePtr, tag->defaultValue, tag->size, elementCount );
    }
  }
  const unsigned char* ptr = reinterpret_cast<const unsigned char*>(data);
    // as we're working with user-supplied arrays, make sure we're not
    // memcpying overlapping regions
  assert(ptr+tag->size*elementCount <= tag->eleWritePtr ||
         tag->eleWritePtr+tag->size*elementCount <= ptr);
  for (size_t i = 0; i < count; ++i) {
    size_t idx = reinterpret_cast<size_t>(entities[i]);
    memcpy( tag->eleWritePtr + idx*tag->size,
            ptr + i*tag->size,
            tag->size );
  }
}
  
void ArrayMesh::tag_set_vertex_data ( TagHandle handle, 
                                      size_t count, 
                                      const VertexHandle* entities, 
                                      const void* data, 
                                      MsqError& err )
{
  Tag* tag = reinterpret_cast<Tag*>(handle);
  if (!tag->vtxWritePtr) {
    if (!tag->owned) {
      MSQ_SETERR(err)("Attempt to set non-writeable (application owned) tag data", 
                      MsqError::TAG_ALREADY_EXISTS);
      return;
    }
    else {
      assert(!tag->vtxReadPtr);
      tag->vtxReadPtr = tag->vtxWritePtr = new unsigned char[vertexCount*tag->size];
      fill( tag->vtxWritePtr, tag->defaultValue, tag->size, vertexCount );
    }
  }
  const unsigned char* ptr = reinterpret_cast<const unsigned char*>(data);
    // as we're working with user-supplied arrays, make sure we're not
    // memcpying overlapping regions
  assert(ptr+tag->size*vertexCount <= tag->vtxWritePtr ||
         tag->vtxWritePtr+tag->size*vertexCount <= ptr);
  for (size_t i = 0; i < count; ++i) {
    size_t idx = reinterpret_cast<size_t>(entities[i]) - oneBasedArrays;
    memcpy( tag->vtxWritePtr + idx*tag->size,
            ptr + i*tag->size,
            tag->size );
  }
}

void ArrayMesh::tag_get_element_data( TagHandle handle, 
                                      size_t count, 
                                      const ElementHandle* entities, 
                                      void* data, 
                                      MsqError& err )
{
  unsigned char* ptr = reinterpret_cast<unsigned char*>(data);
  const Tag* tag = reinterpret_cast<const Tag*>(handle);
  if (tag->eleReadPtr) {
    for (size_t i = 0; i < count; ++i) {
      size_t idx = reinterpret_cast<size_t>(entities[i]);
      memcpy( ptr + i*tag->size,
              tag->eleReadPtr + idx*tag->size,
              tag->size );
    }
  }
  else if (tag->defaultValue) {
    fill( ptr, tag->defaultValue, tag->size, count );
  }
  else {
    MSQ_SETERR(err)(MsqError::TAG_NOT_FOUND);
  }
}    

void ArrayMesh::tag_get_vertex_data ( TagHandle handle, 
                                      size_t count, 
                                      const VertexHandle* entities, 
                                      void* data, 
                                      MsqError& err )
{
  unsigned char* ptr = reinterpret_cast<unsigned char*>(data);
  const Tag* tag = reinterpret_cast<const Tag*>(handle);
  if (tag->vtxReadPtr) {
    for (size_t i = 0; i < count; ++i) {
      size_t idx = reinterpret_cast<size_t>(entities[i]) - oneBasedArrays;
      memcpy( ptr + i*tag->size,
              tag->vtxReadPtr + idx*tag->size,
              tag->size );
    }
  }
  else if (tag->defaultValue) {
    fill( ptr, tag->defaultValue, tag->size, count );
  }
  else {
    MSQ_SETERR(err)(MsqError::TAG_NOT_FOUND);
  }
}

  
} // namespace Mesquite
