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

#include "MeshImplTags.hpp"
#include "MsqError.hpp"
#include <assert.h>
#include <stdlib.h>

namespace MESQUITE_NS {


MeshImplTags::TagData::~TagData() 
{
  if (elementData) 
    free(elementData);
  if (vertexData)
    free(vertexData);
  if (defaultValue)
    free(defaultValue);
}


void MeshImplTags::clear()
{
  for (std::vector<TagData*>::iterator iter = tagList.begin();
       iter != tagList.end(); ++iter)
    if (*iter)
      delete *iter;
  
  tagList.clear();
}

size_t MeshImplTags::size_from_tag_type( Mesh::TagType type )
{
  switch( type ) {
    case Mesh::BYTE:   return 1;
    case Mesh::BOOL:   return sizeof(bool);
    case Mesh::DOUBLE: return sizeof(double);
    case Mesh::INT:    return sizeof(int);
    case Mesh::HANDLE: return sizeof(void*);
    default: assert(0); return 0;
  }
}

size_t MeshImplTags::create( const std::string& name,
                             Mesh::TagType type,
                             unsigned length,
                             const void* defval,
                             MsqError& err )
{
  size_t h = handle( name, err );
  if (h)
  {
    MSQ_SETERR(err)(name, MsqError::TAG_ALREADY_EXISTS);
    return 0;
  }
  
  if (length == 0 || size_from_tag_type(type) == 0)
  {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    return 0;
  }
  
  TagData* tag = new TagData( name, type, length );
  h = tagList.size();
  tagList.push_back(tag);
  
  if (defval)
  {
    tag->defaultValue = malloc( tag->desc.size );
    memcpy( tag->defaultValue, defval, tag->desc.size );
  }
  
  return h+1;
}

size_t MeshImplTags::create( const TagDescription& desc,
                             const void* defval,
                             MsqError& err )
{
  size_t h = handle( desc.name.c_str(), err );
  if (h)
  {
    MSQ_SETERR(err)(desc.name.c_str(), MsqError::TAG_ALREADY_EXISTS);
    return 0;
  }
  
  err.clear();
  if (desc.size == 0 || (desc.size % size_from_tag_type(desc.type)) != 0)
  {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    return 0;
  }
  
  TagData* tag = new TagData( desc );
  h = tagList.size();
  tagList.push_back(tag);
  
  if (defval)
  {
    tag->defaultValue = malloc( tag->desc.size );
    memcpy( tag->defaultValue, defval, tag->desc.size );
  }
  
  return h+1;
}

void MeshImplTags::destroy( size_t tag_index, MsqError& err )
{
  --tag_index;
  if (tag_index >= tagList.size() || 0 == tagList[tag_index])
  {
    MSQ_SETERR(err)(MsqError::TAG_NOT_FOUND);
    return ;
  }
  
  delete tagList[tag_index];
  tagList[tag_index] = 0;
}

size_t MeshImplTags::handle( const std::string& name, MsqError& err ) const
{
  for (size_t i = 0; i < tagList.size(); ++i)
    if (tagList[i] && tagList[i]->desc.name == name)
      return i+1;
      
  return 0;
}

const TagDescription& MeshImplTags::properties( size_t tag_index, MsqError& err ) const
{
  static TagDescription dummy_desc;
  --tag_index;
  
  if (tag_index >= tagList.size() || !tagList[tag_index])
  {
    MSQ_SETERR(err)("Invalid tag handle", MsqError::INVALID_ARG);
    return dummy_desc;
  }
  
  return tagList[tag_index]->desc;
}


void MeshImplTags::set_element_data( size_t tag_index,
                                     size_t num_indices,
                                     const size_t* index_array,
                                     const void* values,
                                     MsqError& err )
{
  size_t i;
  char* data;
  --tag_index;
  if (tag_index >= tagList.size() || !tagList[tag_index])
  {
    MSQ_SETERR(err)("Invalid tag handle", MsqError::INVALID_ARG);
    return;
  }
  
  TagData* tag = tagList[tag_index];
  
    // Get highest element index
  size_t total = tag->elementCount;
  for (i = 0; i < num_indices; ++i)
    if (index_array[i] >= total)
      total = index_array[i] + 1;
  
    // If need more space
  if (total > tag->elementCount)
  {
      // allocate more space
    tag->elementData = realloc( tag->elementData, tag->desc.size * total );
      // if a default value, initialize new space with it
    if (tag->defaultValue)
    {
      data = ((char*)tag->elementData) + tag->elementCount * tag->desc.size;
      for (i = tag->elementCount; i < total; ++i)
      {
        memcpy( data, tag->defaultValue, tag->desc.size );
        data += tag->desc.size;
      }
    }
    else
    {
      memset( (char*)tag->elementData + tag->elementCount * tag->desc.size, 0, 
              (total - tag->elementCount) * tag->desc.size );
    }
    tag->elementCount = total;
  }
  
    // Store passed tag values
  data = (char*)tag->elementData;
  const char* iter = (const char*)values;
  for (i = 0; i < num_indices; ++i)
  {
    memcpy( data + index_array[i]*tag->desc.size, iter, tag->desc.size );
    iter += tag->desc.size;
  }
}

void MeshImplTags::get_element_data( size_t tag_index,
                                     size_t num_indices,
                                     const size_t* index_array,
                                     void* values,
                                     MsqError& err ) const
{
  --tag_index;
  if (tag_index >= tagList.size() || !tagList[tag_index])
  {
    MSQ_SETERR(err)("Invalid tag handle", MsqError::INVALID_ARG);
    return;
  }
  
  TagData* tag = tagList[tag_index];
  
  char* iter = (char*)values;
  const char* data = (const char*)tag->elementData;
  
  for (size_t i = 0; i < num_indices; ++i)
  {
    const void* ptr;
    size_t index = index_array[i];
    if (index >= tag->elementCount)
    {
      ptr = tag->defaultValue;
      if (!ptr)
      {
        MSQ_SETERR(err)(MsqError::TAG_NOT_FOUND);
        return;
      }
    }
    else
    {
      ptr = data + index * tag->desc.size;
    }
    
    memcpy( iter, ptr, tag->desc.size );
    iter += tag->desc.size;
  }
}

void MeshImplTags::set_vertex_data( size_t tag_index,
                                    size_t num_indices,
                                    const size_t* index_array,
                                    const void* values,
                                    MsqError& err )
{
  size_t i;
  char* data;
  --tag_index;
  if (tag_index >= tagList.size() || !tagList[tag_index])
  {
    MSQ_SETERR(err)("Invalid tag handle", MsqError::INVALID_ARG);
    return;
  }
  
  TagData* tag = tagList[tag_index];
  
    // Get highest element index
  size_t total = tag->vertexCount;
  for (i = 0; i < num_indices; ++i)
    if (index_array[i] >= total)
      total = index_array[i] + 1;
  
    // If need more space
  if (total > tag->vertexCount)
  {
      // allocate more space
    tag->vertexData = realloc( tag->vertexData, tag->desc.size * total );
      // if a default value, initialize new space with it
    if (tag->defaultValue)
    {
      data = ((char*)tag->vertexData) + tag->vertexCount * tag->desc.size;
      for (i = tag->vertexCount; i < total; ++i)
      {
        memcpy( data, tag->defaultValue, tag->desc.size );
        data += tag->desc.size;
      }
    }
    else
    {
      memset( (char*)tag->vertexData + tag->vertexCount * tag->desc.size, 0, 
              (total - tag->vertexCount) * tag->desc.size );
    }
    tag->vertexCount = total;
  }
  
    // Store passed tag values
  data = (char*)tag->vertexData;
  const char* iter = (const char*)values;
  for (i = 0; i < num_indices; ++i)
  {
    memcpy( data + index_array[i]*tag->desc.size, iter, tag->desc.size );
    iter += tag->desc.size;
  }
}

void MeshImplTags::get_vertex_data( size_t tag_index,
                                    size_t num_indices,
                                    const size_t* index_array,
                                    void* values,
                                    MsqError& err ) const
{
  --tag_index;
  if (tag_index >= tagList.size() || !tagList[tag_index])
  {
    MSQ_SETERR(err)("Invalid tag handle", MsqError::INVALID_ARG);
    return;
  }
  
  TagData* tag = tagList[tag_index];
  
  char* iter = (char*)values;
  const char* data = (const char*)tag->vertexData;
  
  for (size_t i = 0; i < num_indices; ++i)
  {
    const void* ptr;
    size_t index = index_array[i];
    if (index >= tag->vertexCount)
    {
      ptr = tag->defaultValue;
      if (!ptr)
      {
        MSQ_SETERR(err)(MsqError::TAG_NOT_FOUND);
        return;
      }
    }
    else
    {
      ptr = data + index * tag->desc.size;
    }
    
    memcpy( iter, ptr, tag->desc.size );
    iter += tag->desc.size;
  }
}

bool MeshImplTags::tag_has_vertex_data( size_t tag_index, MsqError& err ) 
{
  --tag_index;
  if (tag_index >= tagList.size() || !tagList[tag_index])
  {
    MSQ_SETERR(err)("Invalid tag handle", MsqError::INVALID_ARG);
    return false;
  }
  
  TagData* tag = tagList[tag_index];
  return 0 != tag->vertexData || tag->defaultValue;
}  

bool MeshImplTags::tag_has_element_data( size_t tag_index, MsqError& err ) 
{
  --tag_index;
  if (tag_index >= tagList.size() || !tagList[tag_index])
  {
    MSQ_SETERR(err)("Invalid tag handle", MsqError::INVALID_ARG);
    return false;
  }
  
  TagData* tag = tagList[tag_index];
  return 0 != tag->elementData || tag->defaultValue;
}

MeshImplTags::TagIterator MeshImplTags::tag_begin()
{
  size_t index = 0;
  while (index < tagList.size() && tagList[index] == NULL)
    ++index;
  return TagIterator( this, index );
}

MeshImplTags::TagIterator MeshImplTags::TagIterator::operator++()
{
  ++index;
  while (index < tags->tagList.size() && NULL == tags->tagList[index])
    ++index;
  return TagIterator( tags, index );
}

MeshImplTags::TagIterator MeshImplTags::TagIterator::operator--()
{
  --index;
  while (index < tags->tagList.size() && NULL == tags->tagList[index])
    --index;
  return TagIterator( tags, index );
}

MeshImplTags::TagIterator MeshImplTags::TagIterator::operator++(int)
{
  size_t old = index;
  ++index;
  while (index < tags->tagList.size() && NULL == tags->tagList[index])
    ++index;
  return TagIterator( tags, old );
}

MeshImplTags::TagIterator MeshImplTags::TagIterator::operator--(int)
{
  size_t old = index;
  --index;
  while (index < tags->tagList.size() && NULL == tags->tagList[index])
    --index;
  return TagIterator( tags, old );
}


} //namespace Mesquite
