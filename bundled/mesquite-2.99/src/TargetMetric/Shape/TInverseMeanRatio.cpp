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

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TInverseMeanRatio.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TInverseMeanRatio.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

std::string TInverseMeanRatio::get_name() const
  { return "TInverseMeanRatio"; }

TInverseMeanRatio::~TInverseMeanRatio() {}

bool TInverseMeanRatio::evaluate( const MsqMatrix<2,2>& T, 
                                   double& result, 
                                   MsqError& err )
{
  const double d = det( T );
  if (invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  else {
    result = sqr_Frobenius(T) / (2 * d) - 1;
    return true;
  }
}


bool TInverseMeanRatio::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                             double& result,
                                             MsqMatrix<2,2>& deriv_wrt_T,
                                             MsqError& err )
{
  const double d = det( T );
  if (invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    deriv_wrt_T = MsqMatrix<2,2>(0.0);
    return false;
  }
  else {
    result = sqr_Frobenius(T) / (2 * d);
    deriv_wrt_T = transpose_adj(T);
    deriv_wrt_T *= -result;
    deriv_wrt_T += T;
    deriv_wrt_T *= 1.0/d;
    result -= 1.0;
    return true;
  }
}


bool TInverseMeanRatio::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                             double& result,
                                             MsqMatrix<2,2>& dA,
                                             MsqMatrix<2,2> d2A[3],
                                             MsqError& err )
{
  const double d = det( T );
  if (invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    dA = d2A[0] = d2A[1] = d2A[2] = MsqMatrix<2,2>(0.0);
    return false;
  }
  else {
    const double inv_det = 1.0/d;
    result = sqr_Frobenius(T) * 0.5 * inv_det;
    
    const MsqMatrix<2,2> AT = transpose_adj(T);
    dA = AT;
    dA *= -result;
    dA += T;
    dA *= inv_det;
    
    const double p3 = -result * inv_det;
    const double p1 = -2.0 * p3 * inv_det;
    const double p2 = -inv_det * inv_det;
    const MsqMatrix<2,2> AT_T_op_00 = outer( AT.row(0), T.row(0));
    const MsqMatrix<2,2> AT_T_op_11 = outer( AT.row(1), T.row(1));
    d2A[0] = p1 * outer( AT.row(0), AT.row(0))
           + p2 * (AT_T_op_00 + transpose(AT_T_op_00));
    d2A[1] = p1 * outer( AT.row(0), AT.row(1)) 
           + p2 * (outer( AT.row(0), T.row(1))
	     + outer( T.row(0), AT.row(1) ));
    d2A[2] = p1 * outer( AT.row(1), AT.row(1)) 
           + p2 * (AT_T_op_11 + transpose(AT_T_op_11));

    d2A[0](0,0) += inv_det;
    d2A[0](1,1) += inv_det;
    d2A[1](0,1) += p3;
    d2A[1](1,0) -= p3;
    d2A[2](0,0) += inv_det;
    d2A[2](1,1) += inv_det;
    
    result -= 1.0;
    return true;
  }
}


bool TInverseMeanRatio::evaluate( const MsqMatrix<3,3>& T, 
                                   double& result, 
                                   MsqError& err )
{
  const double d = det( T );
  if (invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    return false;
  }
  else {
    const double det_cbrt = Mesquite::cbrt(d);
    result = sqr_Frobenius(T) / (3 * det_cbrt * det_cbrt) - 1;
    return true;
  }
}


bool TInverseMeanRatio::evaluate_with_grad( const MsqMatrix<3,3>& T,
                                             double& result,
                                             MsqMatrix<3,3>& deriv_wrt_T,
                                             MsqError& err )
{  
  const double d = det( T );
  if (invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    deriv_wrt_T = MsqMatrix<3,3>(0.0);
    return false;
  }

  const double inv_det = 1.0/d;
  const double inv_det_cbrt = Mesquite::cbrt(inv_det);
  const double inv_3_det_twothirds = inv_det_cbrt * inv_det_cbrt / 3.0;
  const double fnorm = sqr_Frobenius(T);
  result = fnorm * inv_3_det_twothirds - 1;
  deriv_wrt_T = transpose_adj(T);
  deriv_wrt_T *= -fnorm * inv_det / 3.0;
  deriv_wrt_T += T;
  deriv_wrt_T *= 2.0 * inv_3_det_twothirds;
  return true;
}


bool TInverseMeanRatio::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                             double& result,
                                             MsqMatrix<3,3>& dA,
                                             MsqMatrix<3,3> d2A[6],
                                             MsqError& err )
{
  const double d = det( T );
  if (invalid_determinant(d)) {
    MSQ_SETERR(err)( barrier_violated_msg, MsqError::BARRIER_VIOLATED );
    dA = MsqMatrix<3,3>(0.0);
    return false;
  }

  const double f0 = 1.0/d;
  const double c = Mesquite::cbrt(f0);
  const double f1 = 1.0/3.0 * c * c;
  const double f2 = sqr_Frobenius(T);
  result = f1 * f2;

  const double f3 = 2 * f1;
  const double f4 = result * (10.0/9.0) * f0 * f0;
  const double f5 = (1.0/3.0) * f0 * f3;
  const double f6 = 2 * f5;
  const double f7 = f2 * f5;

  const MsqMatrix<3,3> AT = transpose_adj(T);
  dA = AT;
  dA *= (-1.0/3.0) * f0 * f2;
  dA += T;
  dA *= f3;

  MsqMatrix<3,3> op;
  int i = 0;
  double s = 1;
  for (int r = 0; r < 3; ++r) {
    d2A[i] = outer( AT.row(r), AT.row(r) );
    d2A[i] *= f4;
    op = outer( AT.row(r), T.row(r) );
    op += transpose(op);
    op *= f6;
    d2A[i] -= op;
    
    d2A[i](0,0) += f3;
    d2A[i](1,1) += f3;
    d2A[i](2,2) += f3;
    
    ++i;

    for (int c = r+1; c < 3; ++c) {
      d2A[i] = outer( AT.row(r), AT.row(c) );
      d2A[i] *= f4;
      op = outer( AT.row(r), T.row(c) );
      op += outer( T.row(r), AT.row(c) );
      op *= f6;
      d2A[i] -= op;

      MsqMatrix<1,3> rt = T.row(3 - r - c);
      rt *= s * f7;
      d2A[i](0,1) -= rt(0,2);
      d2A[i](0,2) += rt(0,1);
      d2A[i](1,0) += rt(0,2);
      d2A[i](1,2) -= rt(0,0);
      d2A[i](2,0) -= rt(0,1);
      d2A[i](2,1) += rt(0,0);

      ++i;
      s = -s;
    }
  }

  result -= 1.0;
  return true;
}


} // namespace Mesquite
