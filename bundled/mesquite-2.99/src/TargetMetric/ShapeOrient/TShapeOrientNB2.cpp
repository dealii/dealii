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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file TShapeOrientNB2.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TShapeOrientNB2.hpp"
#include "MsqMatrix.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string TShapeOrientNB2::get_name() const
  { return "TShapeOrientNB2"; }

TShapeOrientNB2::~TShapeOrientNB2() {}

template <unsigned DIM> static inline
bool eval( const MsqMatrix<DIM,DIM>& T, double& result )
{
  const double tr = trace(T);
  result = sqr_Frobenius( T ) - DimConst<DIM>::inv() * tr * fabs(tr);
  return true;
}


template <unsigned DIM> static inline
bool grad( const MsqMatrix<DIM,DIM>& T, 
           double& result, 
           MsqMatrix<DIM,DIM>& deriv_wrt_T )
{
  const double tr = trace(T);
  const double f = DimConst<DIM>::inv() * fabs(tr);
  result = sqr_Frobenius( T ) - f * tr;
  deriv_wrt_T = T;
  pluseq_scaled_I( deriv_wrt_T, -f );
  deriv_wrt_T *= 2;
  return true;
}

template <unsigned DIM> static inline
bool hess( const MsqMatrix<DIM,DIM>& T, 
           double& result, 
           MsqMatrix<DIM,DIM>& deriv_wrt_T, 
           MsqMatrix<DIM,DIM>* second_wrt_T )
{
  const double tr = trace(T);
  const double f = DimConst<DIM>::inv() * fabs(tr);
  result = sqr_Frobenius( T ) - f * tr;
  deriv_wrt_T = T;
  pluseq_scaled_I( deriv_wrt_T, -f );
  deriv_wrt_T *= 2;
  set_scaled_I( second_wrt_T, 2.0 );
  pluseq_scaled_outer_product_I_I( second_wrt_T, DimConst<DIM>::inv() * (tr < 0 ? 2 : -2) );
  return true;
}

TMP_T_TEMPL_IMPL_COMMON(TShapeOrientNB2)

} // namespace Mesquite
