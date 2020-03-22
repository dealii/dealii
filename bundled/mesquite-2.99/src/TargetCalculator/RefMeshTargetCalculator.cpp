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


/** \file RefMeshTargetCalculator.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "PatchData.hpp"
#include "MsqError.hpp"
#include "ReferenceMesh.hpp"
#include "ElemSampleQM.hpp"

namespace MESQUITE_NS {

RefMeshTargetCalculator::~RefMeshTargetCalculator() {}

bool RefMeshTargetCalculator::get_3D_target( PatchData& pd, 
                                             size_t element,
                                             Sample sample,
                                             MsqMatrix<3,3>& W_out,
                                             MsqError& err )
{
  get_refmesh_Jacobian_3D( refMesh, pd, element, sample, W_out, err );
  return !MSQ_CHKERR(err);
}

bool RefMeshTargetCalculator::get_surface_target( PatchData& pd, 
                                                    size_t element,
                                                    Sample sample,
                                                    MsqMatrix<3,2>& W_out,
                                                    MsqError& err )
{
  get_refmesh_Jacobian_2D( refMesh, pd, element, sample, W_out, err );
  return !MSQ_CHKERR(err);
}

bool RefMeshTargetCalculator::get_2D_target( PatchData& pd, 
                                             size_t element,
                                             Sample sample,
                                             MsqMatrix<2,2>& W_out,
                                             MsqError& err )
{
  MsqMatrix<3,2> W_orient;
  bool valid = get_surface_target( pd, element, sample, W_orient, err );
  if (MSQ_CHKERR(err) || !valid) return false;
  
  MsqMatrix<3,2> V;
  MsqMatrix<2,2> Q, delta;
  double lambda;
  valid = factor_surface( W_orient, lambda, V, Q, delta, err );
  W_out = lambda * Q * delta;
  return !MSQ_CHKERR(err) && valid;
}

} // namespace Mesquite
