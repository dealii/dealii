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


/** \file ElementAvgQM.hpp
 *  \brief 
 *  \author Boyd Tidwell 
 */

#ifndef MSQ_ELEMENT_AVG_QM_HPP
#define MSQ_ELEMENT_AVG_QM_HPP

#include "Mesquite.hpp"
#include "ElementQM.hpp"

namespace MESQUITE_NS {

class ElemSampleQM;

class ElementAvgQM : public ElementQM
{
public:

  MESQUITE_EXPORT ElementAvgQM( ElemSampleQM* metric );
  
  MESQUITE_EXPORT virtual ~ElementAvgQM();
  
  MESQUITE_EXPORT ElemSampleQM* get_quality_metric() const 
    { return mMetric; }
  
  MESQUITE_EXPORT virtual std::string get_name() const;

  MESQUITE_EXPORT virtual int get_negate_flag() const;

  MESQUITE_EXPORT virtual
  bool evaluate( PatchData& pd, 
                 size_t handle, 
                 double& value, 
                 MsqError& err );

private:

  ElemSampleQM* mMetric;
  mutable std::vector<size_t> mHandles;
};  


} // namespace Mesquite

#endif
