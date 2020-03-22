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


/** \file ScalarAddQualityMetric.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ScalarAddQualityMetric.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"
#include "MsqError.hpp"
#include <sstream>

namespace MESQUITE_NS {

std::string ScalarAddQualityMetric::get_name() const
{
  std::ostringstream str;
  str << mMetric->get_name() << "+" << mOffset;
  return str.str();
}

void ScalarAddQualityMetric::get_evaluations( PatchData& pd, 
                                              std::vector<size_t>& handles, 
                                              bool free_vertices_only,
                                              MsqError& err )
{ 
  mMetric->get_evaluations( pd, handles, free_vertices_only, err );
  MSQ_CHKERR(err);
}

bool ScalarAddQualityMetric::evaluate( PatchData& pd, size_t handle, double& value, MsqError& err )
{
  bool rval = mMetric->evaluate( pd, handle, value, err );
  value += mOffset;
  return rval;
}

bool ScalarAddQualityMetric::evaluate_with_indices( PatchData& pd,
                                                    size_t handle,
                                                    double& value,
                                                    std::vector<size_t>& indices,
                                                    MsqError& err )
{
  bool rval = mMetric->evaluate_with_indices( pd, handle, value, indices, err );
  value += mOffset;
  return !MSQ_CHKERR(err) && rval;
}

bool ScalarAddQualityMetric::evaluate_with_gradient( PatchData& pd,
                                                   size_t handle,
                                                   double& value,
                                                   std::vector<size_t>& indices,
                                                   std::vector<Vector3D>& gradient,
                                                   MsqError& err )
{
  bool rval = mMetric->evaluate_with_gradient( pd, handle, value, indices, gradient, err );
  value += mOffset;
  return !MSQ_CHKERR(err) && rval;
}
  
  
bool ScalarAddQualityMetric::evaluate_with_Hessian( PatchData& pd,
                                                  size_t handle,
                                                  double& value,
                                                  std::vector<size_t>& indices,
                                                  std::vector<Vector3D>& gradient,
                                                  std::vector<Matrix3D>& Hessian,
                                                  MsqError& err )
{
  bool rval = mMetric->evaluate_with_Hessian( pd, handle, value, indices, gradient, Hessian, err );
  value += mOffset;
  return !MSQ_CHKERR(err) && rval;
}


} // namespace Mesquite
