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


/** \file AWShape2DB1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "AWShape2DB1.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string AWShape2DB1::get_name() const
  { return "AWShape2DB1"; }

AWShape2DB1::~AWShape2DB1() {}

bool AWShape2DB1::evaluate( const MsqMatrix<2,2>& A, 
                            const MsqMatrix<2,2>& W, 
                            double& result, 
                            MsqError& err )
{
  const double alpha = det(A);
  const double omega = det(W);
  const double prod = alpha * omega;
  if (AWMetric::invalid_determinant( prod ))
  {
    MSQ_SETERR(err)( barrier_violated_msg_aw, MsqError::BARRIER_VIOLATED );
    return false;
  }
  
  result =  sqr_Frobenius(A * adj(W));
  result += sqr_Frobenius(W * adj(A));
  result *= 0.25/prod;
  result -= 1.0;
  return true;
}


} // namespace Mesquite
