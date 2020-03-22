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


/** \file ElementMaxQM.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ElementMaxQM.hpp"
#include "ElemSampleQM.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include <limits>

namespace MESQUITE_NS {

ElementMaxQM::ElementMaxQM( ElemSampleQM* metric )
    : mMetric(metric)
    {}
    
ElementMaxQM::~ElementMaxQM() {}

std::string ElementMaxQM::get_name() const
{
  std::string result("ElementMaxQM(");
  result += mMetric->get_name();
  result += ")";
  return result;
}


int ElementMaxQM::get_negate_flag() const
{
  return get_quality_metric()->get_negate_flag();
}

bool ElementMaxQM::evaluate( PatchData& pd, 
                              size_t handle, 
                              double& value, 
                              MsqError& err )
{
  ElemSampleQM* qm = get_quality_metric();
  mHandles.clear();
  qm->get_element_evaluations( pd, handle, mHandles, err ); MSQ_ERRFALSE(err);

  bool valid = true;
  double tmpval;
  bool tmpvalid;

  value = -1.e+100; // initialize max computation
  for (std::vector<size_t>::iterator h = mHandles.begin(); h != mHandles.end(); ++h) { 
    tmpvalid = qm->evaluate( pd, *h, tmpval, err );  // MSQ_ERRZERO(err);
    if (!tmpvalid)
    {
      value = +1.e+100;
      return false;   // if any handle within the element makes tmpvalid false, then valid is false, no matter what the other handles say
    }
    else if (tmpval > value)
      value = tmpval;
  }

  return valid;
}
/*
bool ElementMaxQM::evaluate_with_gradient( PatchData& pd, 
                                            size_t handle, 
                                            double& value, 
                                            std::vector<size_t>& indices,
                                            std::vector<Vector3D>& gradient,
                                            MsqError& err )
{
  ElemSampleQM* qm = get_quality_metric();
  mHandles.clear();
  qm->get_element_evaluations( pd, handle, mHandles, err ); MSQ_ERRFALSE(err);

  bool valid = true;
  double tmpval;
  bool tmpvalid;
  unsigned count = 0;
  std::vector<size_t>::iterator h, i, j;

  value = -std::numeric_limits<double>::maximum();
  for (h = mHandles.begin(); h != mHandles.end(); ++h) {
    mIndices.clear();
    mGrad.clear();
    tmpvalid = qm->evaluate_with_gradient( pd, *h, tmpval, mIndices, mGrad, err );
    MSQ_ERRZERO(err);
    if (!tmpvalid) {
      valid = false;
    }
      // new value greater than previous max value
    else if (tmpval - value > 1e-6) {
      indices = mIndices;
      gradient = mGrad;
      count = 1;
      value = tmpval;
    }
      // new value equal to previous max value
    else if (tmpval - value >= -1e-6) {
      ++count;
      for (i = mIndices.begin(); i != mIndices.end(); ++i) {
        j = std::find( indices.begin(), indices.end(), *i );
        if (j == indices.end()) {
          indices.push_back( *i );
          gradient.push_back( mGrad[i - mIndices.begin()] );
        }
        else {
          gradient[j - indices.begin()] += mGrad[i - mIndices.begin()];
        }
      }
    }
  }
  if (count > 1) {
    const double inv_count = 1.0 / count;
    for (std::vector<Vector3D>::iterator g = gradient.begin(); g != gradient.end(); ++g)
      *g *= inv_count;
  }
    
  return valid;
}
*/

} // namespace Mesquite
