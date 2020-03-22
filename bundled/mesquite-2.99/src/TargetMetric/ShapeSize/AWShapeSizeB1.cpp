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
 
    (2010) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file AWShapeSizeB1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "AWShapeSizeB1.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string AWShapeSizeB1::get_name() const
  { return "AWShapeSizeB1"; }

AWShapeSizeB1::~AWShapeSizeB1() {}

bool AWShapeSizeB1::evaluate( const MsqMatrix<2,2>& A, 
                              const MsqMatrix<2,2>& W, 
                              double& result, 
                              MsqError& err )
{ 
  const double alpha = det(A);
  if (AWMetric::invalid_determinant( alpha ))
  {
    MSQ_SETERR(err)( barrier_violated_msg_aw, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  result = sqr_Frobenius( A - 1/alpha * transpose_adj(A) * transpose(W) * W );
  return true;
}

bool AWShapeSizeB1::evaluate( const MsqMatrix<3,3>& A, 
                              const MsqMatrix<3,3>& W, 
                              double& result, 
                              MsqError& err )
{
  const double alpha = det(A);
  if (AWMetric::invalid_determinant( alpha ))
  {
    MSQ_SETERR(err)( barrier_violated_msg_aw, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  result = sqr_Frobenius( A - 1/alpha * transpose_adj(A) * transpose(W) * W );
  return true;
}


} // namespace Mesquite
