/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file IdealShapeTarget.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "IdealShapeTarget.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {

IdealShapeTarget::~IdealShapeTarget() {}

bool IdealShapeTarget::get_3D_target( PatchData& pd, 
                                      size_t element,
                                      Sample sample,
                                      MsqMatrix<3,3>& W_out,
                                      MsqError& err )
{
  MsqMeshEntity& elem = pd.element_by_index(element);
  ideal_shape_3D( elem.get_element_type(),
                  sample, pd, W_out, err );
  return !MSQ_CHKERR(err);
}

bool IdealShapeTarget::get_surface_target( PatchData& pd, 
                                           size_t element,
                                           Sample sample,
                                           MsqMatrix<3,2>& W_out,
                                           MsqError& err )
{
  MsqMatrix<2,2> W_2d;
  bool rval = get_2D_target(pd, element, sample, W_2d, err );
  if (MSQ_CHKERR(err) || !rval)
    return false;
  
  W_out(0,0) = W_2d(0,0); W_out(0,1) = W_2d(0,1);
  W_out(1,0) = W_2d(1,0); W_out(1,1) = W_2d(1,1);
  W_out(2,0) = 0.0;       W_out(2,1) = 0.0;
  return true;
}

bool IdealShapeTarget::get_2D_target( PatchData& pd, 
                                      size_t element,
                                      Sample sample,
                                      MsqMatrix<2,2>& W_out,
                                      MsqError& err )
{
  MsqMeshEntity& elem = pd.element_by_index(element);
  ideal_shape_2D( elem.get_element_type(),
                  sample, pd, W_out, err );
  return !MSQ_CHKERR(err);
}

bool IdealShapeTarget::have_surface_orient() const
  { return false; }


} // namespace MESQUITE_NS
