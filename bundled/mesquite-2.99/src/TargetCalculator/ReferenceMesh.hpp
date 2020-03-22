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


/** \file ReferenceMesh.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_REFERENCE_MESH_HPP
#define MSQ_REFERENCE_MESH_HPP

#include "Mesquite.hpp"
#include "MsqVertex.hpp"
#include "MeshInterface.hpp"

namespace MESQUITE_NS {

class Mesh;

/**\brief Class for accessing reference element topology.
 *
 * This class is provided in order to access reference elements
 * for target calculators, where there is a static 1-to-1 mapping
 * between vertices in the active mesh and vertices in the 
 * reference mesh.
 */
class MESQUITE_EXPORT ReferenceMeshInterface
{
public:
  virtual ~ReferenceMeshInterface();

  virtual void get_reference_vertex_coordinates( 
                           const Mesh::VertexHandle* vertices,
                           size_t num_vertices,
                           Vector3D* coordinates_out,
                           MsqError& err ) = 0;
};

class ReferenceMesh : public ReferenceMeshInterface
{
public:
  MESQUITE_EXPORT ReferenceMesh( Mesh* mesh ) : mMesh(mesh) {}
  MESQUITE_EXPORT virtual ~ReferenceMesh();
  
  MESQUITE_EXPORT virtual 
  void get_reference_vertex_coordinates( 
                           const Mesh::VertexHandle* vertices,
                           size_t num_vertices,
                           Vector3D* coordinates_out,
                           MsqError& err );
private:
  Mesh* mMesh;
  std::vector<MsqVertex> tmpStorage;
};


} // namespace Mesquite

#endif
