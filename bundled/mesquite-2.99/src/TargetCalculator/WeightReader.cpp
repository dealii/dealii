/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file WeightReader.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "WeightReader.hpp"
#include "PatchData.hpp"
#include "MsqError.hpp"
#include "ElemSampleQM.hpp"
#include <sstream>

namespace MESQUITE_NS {

static TagHandle get_tag( Mesh* mesh,
                          unsigned num_doubles,
                          const std::string& base_name,
                          MsqError& err )
{
  std::ostringstream str;
  str << base_name << num_doubles;
  
  TagHandle handle = mesh->tag_get( str.str().c_str(), err ); MSQ_ERRZERO(err);
  
    // double check tag type
  std::string temp_name;
  Mesh::TagType temp_type;
  unsigned temp_length;
  mesh->tag_properties( handle, temp_name, temp_type, temp_length, err );
  MSQ_ERRZERO(err);

  if (temp_type != Mesh::DOUBLE || temp_length != num_doubles)
  {
    MSQ_SETERR(err)( MsqError::TAG_ALREADY_EXISTS,
                    "Mismatched type or length for existing tag \"%s\"",
                     str.str().c_str() );
  }
  
  return handle;
}
   
  
  

WeightReader::WeightReader( std::string name )
  : tagBaseName(name) {}

WeightReader::~WeightReader()
{
}

double WeightReader::get_weight( PatchData &pd,
                                 size_t element,
                                 Sample sample,
                                 MsqError& err )
{
  WeightReaderData& data = get_data( pd );
  
    // calculate index of sample in array 
  NodeSet all_samples = pd.get_samples( element );
  unsigned offset = all_samples.num_before( sample );

  if (!data.weights.empty() && data.elementIndex == element) {
    assert(offset < data.weights.size());
    return data.weights[offset];
  }
  const unsigned num_samples = all_samples.num_nodes();
  const unsigned handle_idx = num_samples - 1;
  
    // get the tag handle
  const TagHandle INVALID_HANDLE = (TagHandle)-1;
  if (data.handles.size() <= handle_idx)
    data.handles.resize( handle_idx + 1, INVALID_HANDLE );
  TagHandle& tag_handle = data.handles[handle_idx];
  if (tag_handle == INVALID_HANDLE) {
    tag_handle = get_tag( pd.get_mesh(),
                          num_samples,
                          tagBaseName.c_str(),
                          err );
    MSQ_ERRZERO(err);
    assert(tag_handle != INVALID_HANDLE);
  }
  
    // get the tag data
  data.weights.resize( num_samples );
  pd.get_mesh()->tag_get_element_data( tag_handle, 1, 
                                       pd.get_element_handles_array() + element,
                                       &data.weights[0],
                                       err );
  if (MSQ_CHKERR(err)) {
    data.weights.clear();
    return false;
  }
  
  data.elementIndex = element;
  
  assert(offset < num_samples);
  return data.weights[offset];
}

  
  
  
void WeightReader::notify_patch_destroyed( WeightReaderData& data )
{
  data.handles.clear();
  data.weights.clear();
}

void WeightReader::notify_new_patch( PatchData&, WeightReaderData& data )
{
  data.weights.clear();
}

void WeightReader::notify_sub_patch( PatchData& pd, 
                                     WeightReaderData& data,
                                     PatchData& subpatch,
                                     const size_t* ,
                                     const size_t* ,
                                     MsqError& err )
{
  WeightReaderData& other = get_data(subpatch);
  if (other.handles.empty()) 
    other.handles = data.handles;
}

} // namespace Mesquite
