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


/** \file AWShapeOrientNB1.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "AWShapeOrientNB1.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string AWShapeOrientNB1::get_name() const
  { return "AWShapeOrientNB1"; }

AWShapeOrientNB1::~AWShapeOrientNB1() {}

template <unsigned DIM> static inline
bool eval( const MsqMatrix<DIM,DIM>& A, 
           const MsqMatrix<DIM,DIM>& W, 
           double& result)
{
  result = std::sqrt(sqr_Frobenius(A) * sqr_Frobenius(W));
  result -= A % W;
  result *= result;
  return true;
}

template <unsigned DIM> static inline
bool grad( const MsqMatrix<DIM,DIM>& A, 
           const MsqMatrix<DIM,DIM>& W, 
           double& result, 
           MsqMatrix<DIM,DIM>& deriv )
{
  const double nsW = sqr_Frobenius(W);
  const double nsA = sqr_Frobenius(A);
  const double nW = std::sqrt(nsW);
  const double nA = std::sqrt(nsA);
  const double dot = A % W;
  result = nA*nW - dot;
  result *= result;
  
  deriv  = nA * W;
  double tmp;
  if (divide(dot,nA,tmp))
    deriv += tmp*A;
  deriv *= -nW;
  deriv += nsW * A;
  deriv += dot * W;
  deriv *= 2.0;
  return true;
}

/*
template <unsigned DIM> static inline
bool hess( const MsqMatrix<DIM,DIM>& A, 
           const MsqMatrix<DIM,DIM>& W, 
           double& result, 
           MsqMatrix<DIM,DIM>& deriv, 
           MsqMatrix<DIM,DIM>* second )
{
}
*/

TMP_AW_TEMPL_IMPL_COMMON_NO2ND(AWShapeOrientNB1)

} // namespace Mesquite
