/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TMetricBarrier.cpp
 *  \brief 
 *  \author Boyd Tidwell 
 */

#include "Mesquite.hpp"
#include "TMetricBarrier.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include <limits>

namespace MESQUITE_NS {

TMetricBarrier::~TMetricBarrier() {}

TMetricBarrier2D::~TMetricBarrier2D() {}
TMetricBarrier3D::~TMetricBarrier3D() {}

bool TMetricBarrier2D::evaluate( const MsqMatrix<3,3>&, 
                                double&,
                                MsqError& err )
{
  MSQ_SETERR(err)("2D target metric cannot be evaluated for volume elements",
                  MsqError::UNSUPPORTED_ELEMENT);
  return false;
}

bool TMetricBarrier3D::evaluate( const MsqMatrix<2,2>&,
                                double&,
                                MsqError& err )
{
  MSQ_SETERR(err)("2D target metric cannot be evaluated for volume elements",
                  MsqError::UNSUPPORTED_ELEMENT);
  return false;
}


} // namespace Mesquite

