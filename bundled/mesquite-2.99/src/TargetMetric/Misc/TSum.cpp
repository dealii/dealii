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


/** \file TSum.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TMetricBarrier.hpp"
#include "TSum.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include "TMPCommon.hpp"

namespace MESQUITE_NS {

std::string TSum::get_name() const
  { return mu1->get_name() + '+' + mu2->get_name(); }

TSum::~TSum() {}

template <unsigned DIM> inline
bool TSum::eval( const MsqMatrix<DIM,DIM>& T, 
                 double& result, 
                 MsqError& err )
{
  double val2;
  bool rval = mu1->evaluate( T, result, err );
  MSQ_ERRZERO(err);
  bool rval2 = mu2->evaluate( T, val2, err );
  MSQ_ERRZERO(err);
  result += val2;
  return rval && rval2;
}

template <unsigned DIM> inline
bool TSum::grad( const MsqMatrix<DIM,DIM>& T, 
                 double& result, 
                 MsqMatrix<DIM,DIM>& deriv, 
                 MsqError& err )
{
  double val2;
  MsqMatrix<DIM,DIM> grad2;
  bool rval = mu1->evaluate_with_grad( T, result, deriv, err );
  MSQ_ERRZERO(err);
  bool rval2 = mu2->evaluate_with_grad( T, val2, grad2, err );
  MSQ_ERRZERO(err);
  result += val2;
  deriv += grad2;
  return rval && rval2;
}

template <unsigned DIM> inline
bool TSum::hess( const MsqMatrix<DIM,DIM>& T, 
                 double& result, 
                 MsqMatrix<DIM,DIM>& deriv_wrt_T, 
                 MsqMatrix<DIM,DIM>* second_wrt_T,
                 MsqError& err )
{
  const int HL = (DIM*(DIM+1))/2;
  double val2;
  MsqMatrix<DIM,DIM> grad2, hess2[HL];

  bool rval = mu1->evaluate_with_hess( T, result, deriv_wrt_T, second_wrt_T, err );
  MSQ_ERRZERO(err);
  bool rval2 = mu2->evaluate_with_hess( T, val2, grad2, hess2, err );
  MSQ_ERRZERO(err);
  result += val2;
  deriv_wrt_T += grad2;
  for (int i = 0; i < HL; ++i)
    second_wrt_T[i] += hess2[i];
  return rval && rval2;
}

TMP_T_TEMPL_IMPL_COMMON_ERR(TSum)

} // namespace MESQUITE_NS
