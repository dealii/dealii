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

#ifndef MSQ_GLOBAL_PATCH_CPP
#define MSQ_GLOBAL_PATCH_CPP

/** \file GlobalPatch.cpp
 *  \brief 
 *  \author Jason Kraftcheck
 */

#include "GlobalPatch.hpp"
#include "MeshInterface.hpp"
#include "MsqError.hpp"
#include <assert.h>

namespace MESQUITE_NS {

const PatchSet::PatchHandle GLOBAL_PATCH_HANDLE = 0;

GlobalPatch::~GlobalPatch() {}

void GlobalPatch::get_patch_handles( std::vector<PatchHandle>& patch_handles_out,
                                     MsqError& )
{
  patch_handles_out.resize(1);
  patch_handles_out[0] = GLOBAL_PATCH_HANDLE;
}

void GlobalPatch::get_patch( PatchHandle patch_handle,
                             std::vector<Mesh::ElementHandle>& elem_handles_out,
                             std::vector<Mesh::VertexHandle>& free_vertices_out,
                             MsqError& err )
{
  free_vertices_out.clear();
  assert(GLOBAL_PATCH_HANDLE == patch_handle);
  get_mesh()->get_all_elements( elem_handles_out, err ); MSQ_ERRRTN(err);
  //get_mesh()->get_all_vertices( free_vertices_out, err ); MSQ_ERRRTN(err);
}

} // namespace Mesquite

#endif
