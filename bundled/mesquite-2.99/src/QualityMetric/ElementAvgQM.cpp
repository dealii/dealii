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
 
  ***************************************************************** */


/** \file ElementAvgQM.cpp
 *  \brief 
 *  \author Boyd Tidwell 
 */

#include "Mesquite.hpp"
#include "ElementAvgQM.hpp"
#include "ElemSampleQM.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include <limits>

namespace MESQUITE_NS {

ElementAvgQM::ElementAvgQM( ElemSampleQM* metric )
    : mMetric(metric)
    {}
    
ElementAvgQM::~ElementAvgQM() {}

std::string ElementAvgQM::get_name() const
{
  std::string result("ElementAvgQM(");
  result += mMetric->get_name();
  result += ")";
  return result;
}


int ElementAvgQM::get_negate_flag() const
{
  return get_quality_metric()->get_negate_flag();
}

bool ElementAvgQM::evaluate( PatchData& pd, 
                              size_t handle, 
                              double& value, 
                              MsqError& err )
{
  ElemSampleQM* qm = get_quality_metric();
  mHandles.clear();
  qm->get_element_evaluations( pd, handle, mHandles, err ); MSQ_ERRFALSE(err);

  bool valid = true;
  double tmpval;
  double accumulate = 0.0;
  int num_values = 0;
  bool tmpvalid;

  value = -std::numeric_limits<double>::infinity();
  for (std::vector<size_t>::iterator h = mHandles.begin(); h != mHandles.end(); ++h) 
  {
    tmpvalid = qm->evaluate( pd, *h, tmpval, err ); MSQ_ERRZERO(err);
    if (!tmpvalid) 
    {
      valid = false;
      break;
    }
    else
    { 
      accumulate += tmpval;    
      num_values++;
    }
  }
  if (valid)
    value = accumulate/num_values;
    
  return valid;
}

} // namespace Mesquite
