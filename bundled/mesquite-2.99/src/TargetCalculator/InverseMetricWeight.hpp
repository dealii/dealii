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


/** \file InverseMetricWeight.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_INVERSE_METRIC_WEIGHT_HPP
#define MSQ_INVERSE_METRIC_WEIGHT_HPP

#include "Mesquite.hpp"
#include "WeightCalculator.hpp"

namespace MESQUITE_NS {

class ElemSampleQM;

/**\brief Use inverse of metric value as target weight */
class MESQUITE_EXPORT InverseMetricWeight : public WeightCalculator
{
public:
  InverseMetricWeight( ElemSampleQM* metric ) : mMetric(metric) {}

  virtual ~InverseMetricWeight();
  
  virtual double get_weight( PatchData& pd, 
                             size_t element,
                             Sample sample,
                             MsqError& err );
private:

  ElemSampleQM* mMetric;
};

} // namespace Mesquite

#endif
