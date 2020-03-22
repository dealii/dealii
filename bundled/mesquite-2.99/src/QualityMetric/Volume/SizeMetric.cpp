/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file SizeMetric.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "SizeMetric.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {

SizeMetric::~SizeMetric() {}

std::string SizeMetric::get_name() const { return "Size"; }

int SizeMetric::get_negate_flag() const { return 1; }

bool SizeMetric::evaluate( PatchData& pd, 
                           size_t handle, 
                           double& value, 
                           MsqError& err )
{
  value = pd.element_by_index(handle).compute_unsigned_area( pd, err );
  return true;
}

} // namespace MESQUITE_NS
