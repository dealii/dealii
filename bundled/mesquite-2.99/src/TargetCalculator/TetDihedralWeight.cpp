/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TetDihedralWeight.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TetDihedralWeight.hpp"
#include "PatchData.hpp"
#include "ReferenceMesh.hpp"
#include <algorithm>

namespace MESQUITE_NS {

static inline double
da( double dot )
{ return 180 - (180/M_PI)*acos( dot ); }

double TetDihedralWeight::get_weight( PatchData& pd, 
                                      size_t element,
                                      Sample ,
                                      MsqError& err )
{
  const double eps = 1e-10;

  MsqMeshEntity &elem = pd.element_by_index( element );
  if (elem.get_element_type() != TETRAHEDRON) {
    MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT);
    return 0.0;
  }

  const size_t *indices = elem.get_vertex_index_array();
  Vector3D v01, v02, v31, v32;
  if (refMesh) {
    const Mesh::VertexHandle* vtx_hdl = pd.get_vertex_handles_array();
    Mesh::VertexHandle handles[] = { vtx_hdl[indices[0]],
                                     vtx_hdl[indices[1]],
                                     vtx_hdl[indices[2]],
                                     vtx_hdl[indices[3]] };
    Vector3D coords[4];
    refMesh->get_reference_vertex_coordinates( handles, 4, coords, err );
    MSQ_ERRZERO(err);

    v01 = coords[1] - coords[0];
    v02 = coords[2] - coords[0];
    v31 = coords[1] - coords[3];
    v32 = coords[2] - coords[3];
  }
  else {  
    const MsqVertex* coords = pd.get_vertex_array();
    v01 = coords[indices[1]] - coords[indices[0]];
    v02 = coords[indices[2]] - coords[indices[0]];
    v31 = coords[indices[1]] - coords[indices[3]];
    v32 = coords[indices[2]] - coords[indices[3]];
  }
  
  Vector3D n012 = v02 * v01;
  Vector3D n013 = v31 * v01;
  Vector3D n023 = v02 * v32;
  Vector3D n123 = v31 * v32;
  
    // normalize face vectors.
  double l012 = n012.length();
  double l013 = n013.length();
  double l023 = n023.length();
  double l123 = n123.length();
  n012 *= (l012 < eps) ? 0.0 : 1.0/l012;
  n013 *= (l013 < eps) ? 0.0 : 1.0/l013;
  n023 *= (l023 < eps) ? 0.0 : 1.0/l023;
  n123 *= (l123 < eps) ? 0.0 : 1.0/l123;
  
    // calculate dihedral handles for each edge
  double ds[] = { da(n012 % n013),
                  da(n012 % n123),
                  da(n012 % n023),
                  da(n013 % n023),
                  da(n013 % n123),
                  da(n023 % n123) };
  
    // calculate weight from max dihedral handle
  double d = *std::max_element( ds, ds+6 );
  return 1/(1 + exp(-mA*(d - mCutoff)));
}


} // namespace MESQUITE_NS
