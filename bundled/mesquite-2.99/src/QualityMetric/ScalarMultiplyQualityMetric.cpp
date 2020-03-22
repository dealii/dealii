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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   ScalarMultiplyQualityMetric.cpp
  \brief 

  \author Todd Munson
  \date   2004-12-21
*/

#include "ScalarMultiplyQualityMetric.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"
#include "MsqError.hpp"

using namespace Mesquite;

std::string ScalarMultiplyQualityMetric::get_name() const
{
  return std::string("Scale(") + mMetric->get_name() + ")";
}

void ScalarMultiplyQualityMetric::get_evaluations( PatchData& pd, 
                                              std::vector<size_t>& handles, 
                                              bool free_vertices_only,
                                              MsqError& err )
{ 
  mMetric->get_evaluations( pd, handles, free_vertices_only, err );
  MSQ_CHKERR(err);
}

bool ScalarMultiplyQualityMetric::evaluate( PatchData& pd, size_t handle, double& value, MsqError& err )
{
  bool rval = mMetric->evaluate( pd, handle, value, err );
  value *= mScale;
  return rval;
}

bool ScalarMultiplyQualityMetric::evaluate_with_indices( PatchData& pd,
                                                    size_t handle,
                                                    double& value,
                                                    std::vector<size_t>& indices,
                                                    MsqError& err )
{
  bool rval = mMetric->evaluate_with_indices( pd, handle, value, indices, err );
  value *= mScale;
  return !MSQ_CHKERR(err) && rval;
}

bool ScalarMultiplyQualityMetric::evaluate_with_gradient( PatchData& pd,
                                                   size_t handle,
                                                   double& value,
                                                   std::vector<size_t>& indices,
                                                   std::vector<Vector3D>& gradient,
                                                   MsqError& err )
{
  bool rval = mMetric->evaluate_with_gradient( pd, handle, value, indices, gradient, err );
  value *= mScale;
  for (std::vector<Vector3D>::iterator i = gradient.begin();
       i != gradient.end(); ++i)
    *i *= mScale;
  return !MSQ_CHKERR(err) && rval;
}
  
  
bool ScalarMultiplyQualityMetric::evaluate_with_Hessian( PatchData& pd,
                                                  size_t handle,
                                                  double& value,
                                                  std::vector<size_t>& indices,
                                                  std::vector<Vector3D>& gradient,
                                                  std::vector<Matrix3D>& Hessian,
                                                  MsqError& err )
{
  bool rval = mMetric->evaluate_with_Hessian( pd, handle, value, indices, gradient, Hessian, err );
  value *= mScale;
  for (std::vector<Vector3D>::iterator i = gradient.begin();
       i != gradient.end(); ++i)
    *i *= mScale;
  for (std::vector<Matrix3D>::iterator j = Hessian.begin();
       j != Hessian.end(); ++j)
    *j *= mScale;
  return !MSQ_CHKERR(err) && rval;
}
