/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file MetricWeight.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MetricWeight.hpp"
#include "ElemSampleQM.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

MetricWeight::~MetricWeight() {}


double MetricWeight::get_weight( PatchData& pd, 
                                 size_t element,
                                 Sample sample,
                                 MsqError& err )
{
  size_t h = ElemSampleQM::handle( sample, element );
  double value;
  bool flag = mMetric->evaluate( pd, h, value, err );
  MSQ_ERRZERO(err);
  if (!flag) {
    MSQ_SETERR(err)("Invalid metric value canot be used as target weight", MsqError::INVALID_STATE);
    return 0.0;
  }
  
  return value;
}

} // namespace Mesquite
