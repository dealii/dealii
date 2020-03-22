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


/** \file AWUntangleBeta.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_AW_UNTANGLE_BETA_HPP
#define MSQ_AW_UNTANGLE_BETA_HPP

#include "Mesquite.hpp"
#include "AWMetric.hpp"

namespace MESQUITE_NS {


/** \f$ \left\{ |\alpha - \gamma \, \omega| - ( \alpha - \gamma \, \omega ) \right\}^p \f$ */
class AWUntangleBeta : public AWMetric
{
  private:
  
  double mGamma;

  public:
  
  MESQUITE_EXPORT
  AWUntangleBeta( double gamma = 0.5 ) : mGamma(0.5) {}

  MESQUITE_EXPORT virtual
  ~AWUntangleBeta();

  MESQUITE_EXPORT virtual
  std::string get_name() const;

  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<2,2>& A, 
                 const MsqMatrix<2,2>& W, 
                 double& result, 
                 MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<3,3>& A, 
                 const MsqMatrix<3,3>& W, 
                 double& result, 
                 MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_grad( const MsqMatrix<2,2>& A,
                           const MsqMatrix<2,2>& W,
                           double& result,
                           MsqMatrix<2,2>& deriv_wrt_A,
                           MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_grad( const MsqMatrix<3,3>& A,
                           const MsqMatrix<3,3>& W,
                           double& result,
                           MsqMatrix<3,3>& deriv_wrt_A,
                           MsqError& err );
/*
  MESQUITE_EXPORT virtual
  bool evaluate_with_hess( const MsqMatrix<2,2>& A,
                           const MsqMatrix<2,2>& W,
                           double& result,
                           MsqMatrix<2,2>& deriv_wrt_A,
                           MsqMatrix<2,2> second_wrt_A[3],
                           MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_hess( const MsqMatrix<3,3>& A,
                           const MsqMatrix<3,3>& W,
                           double& result,
                           MsqMatrix<3,3>& deriv_wrt_A,
                           MsqMatrix<3,3> second_wrt_A[6],
                           MsqError& err );
*/

private:
  template <unsigned D> inline
  bool eval( const MsqMatrix<D,D>& A, const MsqMatrix<D,D>& W, double& result );
  template <unsigned D> inline
  bool grad( const MsqMatrix<D,D>& A, const MsqMatrix<D,D>& W, double& result, MsqMatrix<D,D>& first );
//  template <unsigned D> inline
//  bool hess( const MsqMatrix<D,D>& T, double& result, MsqMatrix<D,D>& first, MsqMatrix<D,D>* second );
};



} // namespace Mesquite

#endif
