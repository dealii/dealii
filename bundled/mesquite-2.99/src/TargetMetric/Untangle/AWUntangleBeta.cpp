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


/** \file AWUntangleBeta.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "AWUntangleBeta.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string AWUntangleBeta::get_name() const
  { return "AWUntangleBeta"; }

AWUntangleBeta::~AWUntangleBeta() {}

const int P = 3;

template <unsigned DIM> inline
bool AWUntangleBeta::eval( const MsqMatrix<DIM,DIM>& A, 
                           const MsqMatrix<DIM,DIM>& W, 
                           double& result)
{
  const double alpha = det(A);
  const double omega = det(W);
  double tmp = alpha - mGamma * omega;
  tmp = fabs(tmp) - tmp;
  result = tmp;
  for (int i = 1; i < P; ++i)
    result *= tmp;
  
  return true;
}

template <unsigned DIM> inline
bool AWUntangleBeta::grad( const MsqMatrix<DIM,DIM>& A, 
                           const MsqMatrix<DIM,DIM>& W, 
                           double& result, 
                           MsqMatrix<DIM,DIM>& deriv )
{
  const double alpha = det(A);
  const double omega = det(W);
  double tmp = 2 * (mGamma * omega - alpha);
  if (tmp < 0.0) {
    result = 0.0;
    deriv = MsqMatrix<DIM,DIM>(0.0);
    return true;
  }
  
  double prod = 1.0;
  for (int i = 1; i < P; ++i)
    prod *= tmp;
  result = prod * tmp;
  
  deriv = -2 *P * prod * transpose_adj(A);
  return true;  
}
/*
template <unsigned DIM> inline
bool AWUntangleBeta::hess( const MsqMatrix<DIM,DIM>& A, 
                           const MsqMatrix<DIM,DIM>& W, 
                           double& result, 
                           MsqMatrix<DIM,DIM>& deriv, 
                           MsqMatrix<DIM,DIM>* second )
{
  result = det(A) - det(W);
  deriv = transpose_adj(A);
  set_scaled_outer_product( second, 2.0, deriv );
  pluseq_scaled_2nd_deriv_of_det( second, 2*result, A );
  deriv *= 2*result;
  result *= result;
  return true;
}
*/
TMP_AW_TEMPL_IMPL_COMMON_NO2ND(AWUntangleBeta)

} // namespace Mesquite
