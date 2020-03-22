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


/** \file JacobianCalculator.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "JacobianCalculator.hpp"
#include "MappingFunction.hpp"
#include "MsqError.hpp"
#include "TopologyInfo.hpp"

namespace MESQUITE_NS {

void JacobianCalculator::get_Jacobian_2D( const MappingFunction2D* mf,
                                          NodeSet ho_bits,
                                          Sample location,
                                          const Vector3D* verts,
                                          size_t num_type_vert,
                                          MsqMatrix<3,2>& J_out,
                                          MsqError& err )
{
  size_t num_vtx = 0;
  mf->derivatives( location, ho_bits, mIndices, mDerivs2D, num_vtx, err ); MSQ_ERRRTN(err);
  mf->convert_connectivity_indices( num_type_vert, mIndices, num_vtx, err ); MSQ_ERRRTN(err);
  const MsqVector<2>* d = mDerivs2D;
  const size_t* const e = mIndices + num_vtx;
  Vector3D c[2] = {Vector3D(0,0,0), Vector3D(0,0,0)};
  for (const size_t* i = mIndices; i != e; ++i, ++d) {
    c[0] += (*d)[0] * verts[*i];
    c[1] += (*d)[1] * verts[*i];
  }
  J_out.set_column( 0, MsqMatrix<3,1>(c[0].to_array()) );
  J_out.set_column( 1, MsqMatrix<3,1>(c[1].to_array()) );
}

void JacobianCalculator::get_Jacobian_3D( const MappingFunction3D* mf,
                                          NodeSet ho_bits,
                                          Sample location,
                                          const Vector3D* verts,
                                          size_t num_type_vert,
                                          MsqMatrix<3,3>& J_out,
                                          MsqError& err )
{
  size_t num_vtx = 0;
  mf->derivatives( location, ho_bits, mIndices, mDerivs3D, num_vtx, err ); MSQ_ERRRTN(err);
  mf->convert_connectivity_indices( num_type_vert, mIndices, num_vtx, err ); MSQ_ERRRTN(err);
  const MsqVector<3>* d = mDerivs3D;
  const size_t* const e = mIndices + num_vtx;
  Vector3D c[3] = {Vector3D(0,0,0), Vector3D(0,0,0), Vector3D(0,0,0)};
  for (const size_t* i = mIndices; i != e; ++i, ++d) {
    c[0] += (*d)[0] * verts[*i];;
    c[1] += (*d)[1] * verts[*i];;
    c[2] += (*d)[2] * verts[*i];;
  }
  J_out.set_column( 0, MsqMatrix<3,1>(c[0].to_array()) );
  J_out.set_column( 1, MsqMatrix<3,1>(c[1].to_array()) );
  J_out.set_column( 2, MsqMatrix<3,1>(c[2].to_array()) );
}


} // namespace Mesquite
