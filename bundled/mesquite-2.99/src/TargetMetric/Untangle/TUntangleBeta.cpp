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


/** \file TUntangleBeta.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TUntangleBeta.hpp"
#include "TMPDerivs.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {


TUntangleBeta::~TUntangleBeta()
{}

std::string TUntangleBeta::get_name() const
  { return "untangle beta"; }


template <unsigned DIM> inline
bool TUntangleBeta::eval( const MsqMatrix<DIM,DIM>& T, 
                          double& result )
{
  double tau = det(T);
  double d = tau - mGamma;
  double f = fabs(d) - d;
  result = 0.125*f*f*f;
  return true;
}

template <unsigned DIM> inline
bool TUntangleBeta::grad( const MsqMatrix<DIM,DIM>& T, 
                          double& result, 
                          MsqMatrix<DIM,DIM>& deriv_wrt_T )
{
  double tau = det(T);
  if (tau < mGamma) {
    double d = mGamma - tau;
    result = d*d*d;
    deriv_wrt_T = -3*d*d*transpose_adj(T);
  }
  else {
    result = 0.0;
    deriv_wrt_T = MsqMatrix<DIM,DIM>(0.0);
  }
  return true;
}

template <unsigned DIM> inline
bool TUntangleBeta::hess( const MsqMatrix<DIM,DIM>& T, 
                          double& result, 
                          MsqMatrix<DIM,DIM>& deriv_wrt_T, 
                          MsqMatrix<DIM,DIM>* second_wrt_T )
{
  double tau = det(T);
  if (tau < mGamma) {
    const MsqMatrix<DIM,DIM> adjt = transpose_adj(T);
    double d = mGamma - tau;
    result = d*d*d;
    deriv_wrt_T = -3*d*d*adjt;
    set_scaled_outer_product( second_wrt_T, 6*d, adjt );
    pluseq_scaled_2nd_deriv_of_det( second_wrt_T, -3*d*d, T );
  }
  else {
    result = 0.0;
    deriv_wrt_T = MsqMatrix<DIM,DIM>(0.0);
    set_scaled_I( second_wrt_T, 0.0 ); // zero everything
  }
  return true;
}

TMP_T_TEMPL_IMPL_COMMON(TUntangleBeta)

} // namespace MESQUITE_NS
