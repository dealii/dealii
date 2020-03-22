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


/** \file VertexMaxQM.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "VertexMaxQM.hpp"
#include "ElemSampleQM.hpp"
#include "VertexQM.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include <algorithm>
#include <limits>

namespace MESQUITE_NS {

VertexMaxQM::VertexMaxQM( ElemSampleQM* metric )
    : mMetric(metric)
    {}
  
VertexMaxQM::~VertexMaxQM() {}

std::string VertexMaxQM::get_name() const
{
  std::string result("VertexMaxQM(");
  result += mMetric->get_name();
  result += ")";
  return result;
}

int VertexMaxQM::get_negate_flag() const
  { return mMetric->get_negate_flag(); }

bool VertexMaxQM::evaluate( PatchData& pd, 
                            size_t handle, 
                            double& value, 
                            MsqError& err )
{
  ElemSampleQM* qm = get_quality_metric();
  get_vertex_corner_handles( pd, handle, mHandles, err ); MSQ_ERRFALSE(err);

  bool valid = true;
  double tmpval;
  bool tmpvalid;

  value = -std::numeric_limits<double>::infinity();
  for (std::vector<size_t>::iterator h = mHandles.begin(); h != mHandles.end(); ++h) {
    tmpvalid = qm->evaluate( pd, *h, tmpval, err ); MSQ_ERRZERO(err);
    if (!tmpvalid) 
      valid = false;
    else if (tmpval > value)
      value = tmpval;
  }
    
  return valid;
}


bool VertexMaxQM::evaluate_with_indices( PatchData& pd, 
                                         size_t handle, 
                                         double& value, 
                                         std::vector<size_t>& indices,
                                         MsqError& err )
{
  ElemSampleQM* qm = get_quality_metric();
  get_vertex_corner_handles( pd, handle, mHandles, err ); MSQ_ERRFALSE(err);

  bool valid = true;
  double tmpval;
  bool tmpvalid;
  std::vector<size_t>::iterator i, e, h;

  value = -std::numeric_limits<double>::infinity();
  for (h = mHandles.begin(); h != mHandles.end(); ++h) {
    mIndices.clear();
    tmpvalid = qm->evaluate_with_indices( pd, *h, tmpval, mIndices, err );
    MSQ_ERRZERO(err);
    if (!tmpvalid) 
      valid = false;
    else if (tmpval > value)
      value = tmpval;
      
    size_t size = indices.size();
    e = indices.begin() + size;
    for (i = mIndices.begin(); i != mIndices.end(); ++i) {
      if (std::find( indices.begin(), e, *i ) == e) {
        indices.push_back(*i);
        e = indices.begin() + size;
      }
    }
  }
    
  return valid;
}

} // namespace Mesquite
