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


/** \file AWShape2DNB1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "AWShape2DNB1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"

namespace MESQUITE_NS {

std::string AWShape2DNB1::get_name() const
  { return "AWShape2DNB1"; }

AWShape2DNB1::~AWShape2DNB1() {}

bool AWShape2DNB1::evaluate( const MsqMatrix<2,2>& A, 
                             const MsqMatrix<2,2>& W, 
                             double& result, 
                             MsqError& )
{
  MsqMatrix<2,2> tmp = det(W) * transpose(A) * A;
  tmp -= det(A) * W * transpose(W);
  result = sqr_Frobenius( tmp );
  return true;
}

bool AWShape2DNB1::evaluate_with_grad( const MsqMatrix<2,2>& A,
                                       const MsqMatrix<2,2>& W,
                                       double& result,
                                       MsqMatrix<2,2>& deriv_wrt_A,
                                       MsqError& )
{
  const double alpha = det(A);
  const double omega = det(W);
  const MsqMatrix<2,2> AtA = transpose(A) * A;
  const MsqMatrix<2,2> WtW = transpose(W) * W;
  const MsqMatrix<2,2> AWt = A * transpose(W);
  const double alphaWtW = alpha * sqr_Frobenius(WtW);
  const double omegaAtA = omega * sqr_Frobenius(AtA);
  const double omegaAWt = omega * sqr_Frobenius(AWt);
  result  =  omega*omegaAtA;
  result +=  alpha*alphaWtW;
  result -=  2*alpha*omegaAWt;
  
  deriv_wrt_A =  4*omega*omega * A * transpose(A) * A;
  deriv_wrt_A -= 4*alpha*omega * AWt * W;
  deriv_wrt_A += 2*(alphaWtW - omegaAWt) * transpose_adj(A);
  return true;
}

} // namespace Mesquite
