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


/** \file ReferenceMesh.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ReferenceMesh.hpp"
#include "MsqError.hpp"
#include "MeshInterface.hpp"
#include <algorithm>

namespace MESQUITE_NS {

ReferenceMeshInterface::~ReferenceMeshInterface() {}

ReferenceMesh::~ReferenceMesh() {}

void ReferenceMesh::get_reference_vertex_coordinates( 
                           const Mesh::VertexHandle* vertices,
                           size_t num_vertices,
                           Vector3D* coordinates_out,
                           MsqError& err )
{
  tmpStorage.resize( num_vertices );
  mMesh->vertices_get_coordinates( vertices, arrptr(tmpStorage), num_vertices, err );
  MSQ_ERRRTN(err);
  std::copy( tmpStorage.begin(), tmpStorage.end(), coordinates_out );
}
    

} // namespace Mesquite
