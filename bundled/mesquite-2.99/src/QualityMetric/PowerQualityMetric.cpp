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
  \file   PowerQualityMetric.cpp
  \brief 

  \author Michael Brewer
  \date   April 1, 2003
*/

#include "PowerQualityMetric.hpp"
#include "MsqError.hpp"

using namespace Mesquite;

PowerQualityMetric::PowerQualityMetric(QualityMetric* qm,
                                       double pow_factor )
  : mMetric(*qm), mPower( pow_factor )
{}

PowerQualityMetric::~PowerQualityMetric() {}

std::string PowerQualityMetric::get_name() const
{
  return std::string("pow(") + mMetric.get_name() + ")";
}

int PowerQualityMetric::get_negate_flag() const
{
  return mPower.value() < 0 ? mMetric.get_negate_flag() : -mMetric.get_negate_flag();
}

void PowerQualityMetric::get_evaluations( PatchData& pd,
                                        std::vector<size_t>& handles,
                                        bool free_only,
                                        MsqError& err )
{
  mMetric.get_evaluations( pd, handles, free_only, err ); MSQ_CHKERR(err);
}

bool PowerQualityMetric::evaluate( PatchData& pd,
                                   size_t handle,
                                   double& value,
                                   MsqError& err )
{
  bool rval = mMetric.evaluate( pd, handle, value, err );
  value = mPower.raise(value);
  return !MSQ_CHKERR(err) && rval;
}

bool PowerQualityMetric::evaluate_with_indices( PatchData& pd,
                                                size_t handle,
                                                double& value,
                                                std::vector<size_t>& indices,
                                                MsqError& err )
{
  bool rval = mMetric.evaluate_with_indices( pd, handle, value, indices, err );
  value = mPower.raise(value);
  return !MSQ_CHKERR(err) && rval;
}


bool PowerQualityMetric::evaluate_with_gradient( PatchData& pd,
                                                 size_t handle,
                                                 double& value,
                                                 std::vector<size_t>& indices,
                                                 std::vector<Vector3D>& gradient,
                                                 MsqError& err )
{
  bool rval = mMetric.evaluate_with_gradient( pd, handle, value, indices, gradient, err );
  const double v = mPower.raise(value);
  const double g = fabs(value) > DBL_EPSILON ? mPower.value() * v / value : 0;
  value = v;
  for (std::vector<Vector3D>::iterator i = gradient.begin(); i != gradient.end(); ++i)
    *i *= g;
  return !MSQ_CHKERR(err) && rval;
}
  
  
bool PowerQualityMetric::evaluate_with_Hessian( PatchData& pd,
                                                size_t handle,
                                                double& value,
                                                std::vector<size_t>& indices,
                                                std::vector<Vector3D>& gradient,
                                                std::vector<Matrix3D>& Hessian,
                                                MsqError& err )
{
  indices.clear();
  bool rval = mMetric.evaluate_with_Hessian( pd, handle, value, indices, gradient, Hessian, err );
  const double v = mPower.raise(value);
  const double g = fabs(value) > DBL_EPSILON ? mPower.value() * v / value : 0.0;
  const double h = fabs(value) > DBL_EPSILON ? g * (mPower.value() - 1) / value : 0.0;
  value = v;
  Matrix3D op;
  unsigned idx = 0;
  unsigned n = indices.size();
  for (unsigned r = 0; r < n; ++r) {
    for (unsigned c = r; c < n; ++c) {
      Hessian[idx] *= g;
      op.outer_product( gradient[r], gradient[c] );
      op *= h;
      Hessian[idx] += op;
      ++idx;
    }
    gradient[r] *= g;
  }
  return !MSQ_CHKERR(err) && rval;
}
