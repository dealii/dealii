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


/** \file AWSizeB1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "AWSizeB1.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string AWSizeB1::get_name() const
  { return "AWSizeB1"; }

AWSizeB1::~AWSizeB1() {}

bool AWSizeB1::evaluate( const MsqMatrix<2,2>& A, 
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
  
  result = (alpha - omega);
  result *= result;
  result /= prod;
  return true;
}

bool AWSizeB1::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                   const MsqMatrix<2,2>& W,
                                   double& result,
                                   MsqMatrix<2,2>& deriv_wrt_A,
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

  result = (alpha - omega);
  result *= result;
  result /= prod;
  deriv_wrt_A = transpose_adj(A);
  deriv_wrt_A *= (alpha*alpha - omega*omega)/(alpha*prod);
  return true;
}

bool AWSizeB1::evaluate( const MsqMatrix<3,3>& A, 
                         const MsqMatrix<3,3>& W, 
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
  };
  
  result = (alpha - omega);
  result *= result;
  result /= prod;
  return true;
}

bool AWSizeB1::evaluate_with_grad( const MsqMatrix<3,3>& A,
                                   const MsqMatrix<3,3>& W,
                                   double& result,
                                   MsqMatrix<3,3>& deriv_wrt_A,
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

  result = (alpha - omega);
  result *= result;
  result /= prod;
  deriv_wrt_A = transpose_adj(A);
  deriv_wrt_A *= (alpha*alpha - omega*omega)/(alpha*prod);
  return true;
}

} // namespace Mesquite
