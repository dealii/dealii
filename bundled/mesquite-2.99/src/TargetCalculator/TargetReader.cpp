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


/** \file TargetReader.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TargetReader.hpp"
#include "PatchData.hpp"
#include "MsqError.hpp"
#include "MsqMatrix.hpp"
#include "ElemSampleQM.hpp"
#include <sstream>

namespace MESQUITE_NS {

static TagHandle get_tag( Mesh* mesh,
                          unsigned num_matrices,
                          unsigned dimension,
                          const char* base_name,
                          bool orient_surface,
                          MsqError& err )
{
  unsigned matrix_size;
  if (dimension == 2 && !orient_surface)
    matrix_size = 4;
  else 
    matrix_size = dimension*3;
  unsigned num_doubles = num_matrices * matrix_size;
  std::ostringstream str;
  str << base_name << num_doubles;
  
  TagHandle handle = mesh->tag_get( str.str().c_str(), err ); MSQ_ERRZERO(err);
  
    // check tag type
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
   
  
  

TargetReader::TargetReader( bool orient2d, std::string name )
  : tagBaseName(name),
    orient2D(orient2d) {}

TargetReader::~TargetReader()
{
}

bool TargetReader::get_3D_target( PatchData &pd,
                                  size_t element,
                                  Sample sample,
                                  MsqMatrix<3,3>& W_out,
                                  MsqError& err )
{
    // calculate index of sample in array 
  NodeSet all_samples = pd.get_samples( element );
  unsigned offset = all_samples.num_before( sample );

  int dim = TopologyInfo::dimension( pd.element_by_index(element).get_element_type() );
  if (dim == 2) {
    MSQ_SETERR(err)("Attempt to read 3D target for surface element", MsqError::INVALID_STATE );
    return false;
  }
  
  TargetReaderData& data = get_data( pd );
  if (!data.targets3D.empty() && data.elementIndex == element) {
    assert(offset < data.targets3D.size());
    W_out = data.targets3D[offset];
    return true;
  }
  const unsigned num_samples = all_samples.num_nodes();
  const unsigned handle_idx = num_samples - 1;
  
    // get the tag handle
  const size_t INVALID_HANDLE = (size_t)-1;
  if (data.handles3D.size() <= handle_idx)
    data.handles3D.resize( handle_idx + 1, (TagHandle)INVALID_HANDLE );
  TagHandle& tag_handle = data.handles3D[handle_idx];
  if (tag_handle == (TagHandle)INVALID_HANDLE) {
    tag_handle = get_tag( pd.get_mesh(),
                          num_samples, 3,
                          tagBaseName.c_str(),
                          orient2D,
                          err );
    MSQ_ERRZERO(err);
    assert(tag_handle != (TagHandle)INVALID_HANDLE);
  }
  
    // get the tag data
  data.targets3D.resize( num_samples );
  pd.get_mesh()->tag_get_element_data( tag_handle, 1, 
                                       pd.get_element_handles_array() + element,
                                       &data.targets3D[0],
                                       err );
  if (MSQ_CHKERR(err)) {
    data.targets3D.clear();
    return false;
  }

  data.elementIndex = element;
  W_out = data.targets3D[offset];
  return true;
}

bool TargetReader::get_2D_target( PatchData &pd,
                                  size_t element,
                                  Sample sample,
                                  MsqMatrix<2,2>& W_out,
                                  MsqError& err )
{
    // which type of 2D target do we actually have
  if (orient2D) {
    MSQ_SETERR(err)("Incorrect surface mesh target type", MsqError::INTERNAL_ERROR );
    return false;
  }
  
    // calculate index of sample in array 
  NodeSet all_samples = pd.get_samples( element );
  unsigned offset = all_samples.num_before( sample );

  int dim = TopologyInfo::dimension( pd.element_by_index(element).get_element_type() );
  if (dim == 3) {
    MSQ_SETERR(err)("Attempt to read surface target for region element", MsqError::INVALID_STATE );
    return false;
  }
  
  TargetReaderData& data = get_data( pd );
  if (!data.targets2D.empty() && data.elementIndex == element) {
    assert(offset < data.targets2D.size());
    W_out = data.targets2D[offset];
    return true;
  }
  const unsigned num_samples = all_samples.num_nodes();
  const unsigned handle_idx = num_samples - 1;
  
    // get the tag handle
  const size_t INVALID_HANDLE = (size_t)-1;
  if (data.handles2D.size() <= handle_idx)
    data.handles2D.resize( handle_idx + 1, (TagHandle)INVALID_HANDLE );
  TagHandle& tag_handle = data.handles2D[handle_idx];
  if (tag_handle == (TagHandle)INVALID_HANDLE) {
    tag_handle = get_tag( pd.get_mesh(),
                          num_samples, 2,
                          tagBaseName.c_str(),
                          orient2D,
                          err );
    MSQ_ERRZERO(err);
    assert(tag_handle != (TagHandle)INVALID_HANDLE);
  }
  
    // get the tag data
  data.targets2D.resize( num_samples );
  pd.get_mesh()->tag_get_element_data( tag_handle, 1, 
                                       pd.get_element_handles_array() + element,
                                       &data.targets2D[0],
                                       err );
  if (MSQ_CHKERR(err)) {
    data.targets2D.clear();
    return false;
  }
  
  data.elementIndex = element;
  W_out = data.targets2D[offset];
  return true;
}

bool TargetReader::get_surface_target( PatchData &pd,
                                  size_t element,
                                  Sample sample,
                                  MsqMatrix<3,2>& W_out,
                                  MsqError& err )
{
    // which type of 2D target do we actually have
  if (!orient2D) {
    MSQ_SETERR(err)("Incorrect surface mesh target type", MsqError::INTERNAL_ERROR );
    return false;
  }
  
    // calculate index of sample in array 
  NodeSet all_samples = pd.get_samples( element );
  unsigned offset = all_samples.num_before( sample );

  int dim = TopologyInfo::dimension( pd.element_by_index(element).get_element_type() );
  if (dim == 3) {
    MSQ_SETERR(err)("Attempt to read surface target for region element", MsqError::INVALID_STATE );
    return false;
  }
  
  TargetReaderData& data = get_data( pd );
  if (!data.targetsSurface.empty() && data.elementIndex == element) {
    assert(offset < data.targetsSurface.size());
    W_out = data.targetsSurface[offset];
    return true;
  }
  const unsigned num_samples = all_samples.num_nodes();
  const unsigned handle_idx = num_samples - 1;
  
    // get the tag handle
  const size_t INVALID_HANDLE = (size_t)-1;
  if (data.handles2D.size() <= handle_idx)
    data.handles2D.resize( handle_idx + 1, (TagHandle)INVALID_HANDLE );
  TagHandle& tag_handle = data.handles2D[handle_idx];
  if (tag_handle == (TagHandle)INVALID_HANDLE) {
    tag_handle = get_tag( pd.get_mesh(),
                          num_samples, 2,
                          tagBaseName.c_str(),
                          orient2D,
                          err );
    MSQ_ERRZERO(err);
    assert(tag_handle != (TagHandle)INVALID_HANDLE);
  }
  
    // get the tag data
  data.targetsSurface.resize( num_samples );
  pd.get_mesh()->tag_get_element_data( tag_handle, 1, 
                                       pd.get_element_handles_array() + element,
                                       &data.targetsSurface[0],
                                       err );
  if (MSQ_CHKERR(err)) {
    data.targetsSurface.clear();
    return false;
  }
  
  data.elementIndex = element;
  W_out = data.targetsSurface[offset];
  return true;
}
  
void TargetReader::notify_patch_destroyed( TargetReaderData& data )
{
  data.handles2D.clear();
  data.handles3D.clear();
  data.targets3D.clear();
  data.targets2D.clear();
  data.targetsSurface.clear();
}

void TargetReader::notify_new_patch( PatchData&, TargetReaderData& data )
{
  data.targets3D.clear();
  data.targets2D.clear();
  data.targetsSurface.clear();
}

void TargetReader::notify_sub_patch( PatchData& pd, 
                                     TargetReaderData& data,
                                     PatchData& subpatch,
                                     const size_t*,
                                     const size_t*,
                                     MsqError& err )
{
  TargetReaderData& other = get_data(subpatch);
  if (other.handles2D.empty() && other.handles3D.empty()) {
    other.handles2D = data.handles2D;
    other.handles3D = data.handles3D;
  }
}

} // namespace Mesquite
