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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov,
    kraftche@cae.wisc.edu
   
  ***************************************************************** */
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 26-March-08 at 10:26:21
//  LAST-MOD: 15-Nov-04 by kraftche@cae.wisc.edu
//
/*! \file ParallelMeshImpl.cpp

\brief This files contains a parallel mesh implementation that can be used
to run mesquite by default.

    \author Jason Kraftcheck
    \author Martin Isenburg
    \date 2008-3-26
 */

#include "ParallelMeshImpl.hpp"
#include "MeshImplData.hpp"
#include "MeshImplTags.hpp"
#include "MsqError.hpp"
#include "MsqDebug.hpp"

namespace MESQUITE_NS
{

ParallelMeshImpl::ParallelMeshImpl(Mesh* myMesh, const char * gid_name, const char * pid_name)
{
  MsqError err;

  this->myMesh = myMesh;
  this->helper = 0;

  if (gid_name)
    gid_tag = myMesh->tag_get(gid_name, err);
  else
    gid_tag = 0;

  if (pid_name)
    pid_tag = myMesh->tag_get(pid_name, err);
  else
    pid_tag = 0;
}

void ParallelMeshImpl::set_global_id_tag(const char * name, MsqError& err)
{
  gid_tag = myMesh->tag_get(name, err);
}

void ParallelMeshImpl::set_processor_id_tag(const char * name, MsqError& err)
{
  pid_tag = myMesh->tag_get(name, err);
}

//**************** Parallel Methods ******************************

void ParallelMeshImpl::vertices_get_global_id(const VertexHandle vert_array[],
					      size_t gid[],
					      size_t num_vtx,
					      MsqError& err)
{
  if (gid_tag)
  {
    myMesh->tag_get_vertex_data(gid_tag, num_vtx, vert_array, gid, err);
    MSQ_CHKERR(err);
  }
  else
  {
    MSQ_SETERR(err)( "Parallel mesh does not have Global IDs.", MsqError::INVALID_STATE);
  }
}

void ParallelMeshImpl::vertices_set_global_id(const VertexHandle vert_array[],
					      size_t gid[],
					      size_t num_vtx,
					      MsqError& err)
{
  if (gid_tag == 0)
  {
    const char GLOBAL_ID_NAME[] = "GLOBAL_ID";

    int default_gid = -1;
    gid_tag = tag_create( GLOBAL_ID_NAME, HANDLE, 1, &default_gid, err );
      // the 'HANDLE' is the type of data to store
      // the '1' is for one value per vertex
      // NULL for no default value, if you want them all
      // initialized to something, pass in a pointer to an int
      // with the value.
    MSQ_CHKERR(err);
  }

  tag_set_vertex_data( gid_tag, num_vtx, vert_array, gid, err );
  MSQ_CHKERR(err);
}
     
void ParallelMeshImpl::vertices_get_processor_id(const VertexHandle vert_array[],
						 int pid[],
						 size_t num_vtx,
						 MsqError& err)
{
  if (pid_tag)
  {
    tag_get_vertex_data( pid_tag, num_vtx, vert_array, pid, err );
    MSQ_CHKERR(err);
  }
  else
  {
    MSQ_SETERR(err)( "Parallel mesh does not have Processor IDs.", MsqError::INVALID_STATE);
  }
}

void ParallelMeshImpl::vertices_set_processor_id(const VertexHandle vert_array[],
						 int pid[],
						 size_t num_vtx,
						 MsqError& err)
{
  if (pid_tag == 0)
  {
    const char PROCESSOR_ID_NAME[] = "PROCESSOR_ID";

    int default_pid = -1;
    pid_tag = tag_create( PROCESSOR_ID_NAME, INT, 1, &default_pid, err );
      // the 'INT' is the type of data to store
      // the '1' is for one value per vertex
      // NULL for no default value, if you want them all
      // initialized to something, pass in a pointer to an int
      // with the value.
    MSQ_CHKERR(err);
  }

  tag_set_vertex_data( pid_tag, num_vtx, vert_array, pid, err );
  MSQ_CHKERR(err); 
}

int ParallelMeshImpl::get_geometric_dimension(MsqError& err)
{
  return myMesh->get_geometric_dimension(err);
}

void ParallelMeshImpl::get_all_elements(std::vector<ElementHandle>& elems,
					MsqError& err)
{
  myMesh->get_all_elements(elems, err);
}

void ParallelMeshImpl::get_all_vertices(std::vector<VertexHandle>& verts,
					MsqError& err)
{
  myMesh->get_all_vertices(verts, err);
}

void ParallelMeshImpl::vertices_get_fixed_flag(const VertexHandle vert_array[],
					       std::vector<bool>& flag_array,
					       size_t num_vtx,
					       MsqError& err)
{
  myMesh->vertices_get_fixed_flag(vert_array,
				  flag_array,
				  num_vtx,
				  err);
}

void ParallelMeshImpl::vertices_get_coordinates(const Mesh::VertexHandle vert_array[],
						MsqVertex* coordinates,
						size_t num_vtx,
						MsqError& err)
{
  myMesh->vertices_get_coordinates(vert_array,
				   coordinates,
				   num_vtx,
				   err);
}

void ParallelMeshImpl::vertices_get_slaved_flag(const VertexHandle vert_array[],
					       std::vector<bool>& flag_array,
					       size_t num_vtx,
					       MsqError& err)
{
  myMesh->vertices_get_slaved_flag(vert_array,
				   flag_array,
				   num_vtx,
				   err);
}

void ParallelMeshImpl::vertex_set_coordinates(VertexHandle vertex,
					      const Vector3D &coordinates,
					      MsqError& err)
{
  myMesh->vertex_set_coordinates(vertex,
				 coordinates,
				 err);
}

void ParallelMeshImpl::vertex_set_byte(VertexHandle vertex,
				       unsigned char byte,
				       MsqError& err)
{
  myMesh->vertex_set_byte(vertex,
			  byte,
			  err);
}

void ParallelMeshImpl::vertices_set_byte(const VertexHandle *vert_array,
					 const unsigned char *byte_array,
					 size_t array_size,
					 MsqError& err)
{
  myMesh->vertices_set_byte(vert_array,
			    byte_array,
			    array_size,
			    err);
}

void ParallelMeshImpl::vertex_get_byte(const VertexHandle vertex,
				       unsigned char *byte,
				       MsqError &err)
{
  myMesh->vertex_get_byte(vertex,
			  byte,
			  err);
}

void ParallelMeshImpl::vertices_get_byte(const VertexHandle *vert_array,
					 unsigned char *byte_array,
					 size_t array_size,
					 MsqError& err)
{
  myMesh->vertices_get_byte(vert_array,
			    byte_array,
			    array_size,
			    err);
}

void ParallelMeshImpl::vertices_get_attached_elements(const VertexHandle* vertices,
						      size_t num_vertices,
						      std::vector<ElementHandle>& elements,
						      std::vector<size_t>& offsets,
						      MsqError& err)
{
  myMesh->vertices_get_attached_elements(vertices,
					 num_vertices,
					 elements,
					 offsets,
					 err);
}

void ParallelMeshImpl::elements_get_attached_vertices(const ElementHandle *elements,
						      size_t num_elems,
						      std::vector<VertexHandle>& vertices,
						      std::vector<size_t>& offsets,
						      MsqError &err) 
{
  myMesh->elements_get_attached_vertices(elements,
					 num_elems,
					 vertices,
					 offsets,
					 err);
}

void ParallelMeshImpl::elements_get_topologies(const ElementHandle *element_handle_array,
					       EntityTopology *element_topologies,
					       size_t num_elements,
					       MsqError& err)
{
  myMesh->elements_get_topologies(element_handle_array,
				  element_topologies,
				  num_elements,
				  err);
}

TagHandle ParallelMeshImpl::tag_create(const std::string& name,
				       TagType type,
				       unsigned length,
				       const void* defval,
				       MsqError& err)
{
  return myMesh->tag_create(name,
			    type,
			    length,
			    defval,
			    err);
}

void ParallelMeshImpl::tag_delete(TagHandle handle, MsqError& err)
{
  myMesh->tag_delete(handle, err);
}

TagHandle ParallelMeshImpl::tag_get(const std::string& name, MsqError& err)
{
  return myMesh->tag_get(name, err);
}

void ParallelMeshImpl::tag_properties(TagHandle handle,
				      std::string& name,
				      TagType& type,
				      unsigned& length,
				      MsqError& err)
{
  myMesh->tag_properties(handle,
			 name,
			 type,
			 length,
			 err);
}

void ParallelMeshImpl::tag_set_element_data(TagHandle handle,
					    size_t num_elems,
					    const ElementHandle* elem_array,
					    const void* values,
					    MsqError& err)
{
  myMesh->tag_set_element_data(handle,
			       num_elems,
			       elem_array,
			       values,
			       err);
}


void ParallelMeshImpl::tag_get_element_data(TagHandle handle,
					    size_t num_elems,
					    const ElementHandle* elem_array,
					    void* values,
					    MsqError& err)
{
  myMesh->tag_get_element_data(handle,
			       num_elems,
			       elem_array,
			       values,
			       err);
}

void ParallelMeshImpl::tag_set_vertex_data(TagHandle handle,
					   size_t num_verts,
					   const VertexHandle* vert_array,
					   const void* values,
					   MsqError& err)
{
  myMesh->tag_set_vertex_data(handle,
			      num_verts,
			      vert_array,
			      values,
			      err);
}

void ParallelMeshImpl::tag_get_vertex_data(TagHandle handle,
					   size_t num_verts,
					   const VertexHandle* vert_array,
					   void* values,
					   MsqError& err)
{
  myMesh->tag_get_vertex_data(handle,
			      num_verts,
			      vert_array,
			      values,
			      err);
}

void ParallelMeshImpl::release_entity_handles(const EntityHandle* handle_array,
				      size_t num_handles,
				      MsqError &err)

{
  myMesh->release_entity_handles(handle_array,
				 num_handles,
				 err);
}

void ParallelMeshImpl::release()
{
  myMesh->release();
}

} // namespace Mesquite
