/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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


/** \file TUntangleMu.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_T_UNTANGLE_MU_HPP
#define MSQ_T_UNTANGLE_MU_HPP

#include "Mesquite.hpp"
#include "TMetric.hpp"

namespace MESQUITE_NS {

/**\brief Composite untangle metric
 *
 * This metric should be combined with TRel2DSize or TRel2DShapeSize
 * to produce a concrete untangle metric.
 *
 * \f$ \mu^\prime(T) = \frac{1}{8}(|d| - d)^3 \f$
 * \f$ d(T) = \sigma - \epsilon - \mu(T()) \f$
 *
 */
class TUntangleMu : public TMetric
{
private:
  TMetric* mBaseMetric;
  double mConstant;

public:

  TUntangleMu( TMetric* base, 
               double sigma = 1.0 ) 
    : mBaseMetric(base),
      mConstant(0.99*sigma) /* default epsilon is 0.01*sigma */
    {}

  TUntangleMu( TMetric* base, 
               double sigma,
               double epsilon ) 
    : mBaseMetric(base),
      mConstant(sigma-epsilon) 
    {}

  MESQUITE_EXPORT virtual
  ~TUntangleMu();

  MESQUITE_EXPORT virtual
  std::string get_name() const;

  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<2,2>& T, 
                 double& result, 
                 MsqError& err );
   
  MESQUITE_EXPORT virtual
  bool evaluate_with_grad( const MsqMatrix<2,2>& T,
                           double& result,
                           MsqMatrix<2,2>& deriv_wrt_T,
                           MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_hess( const MsqMatrix<2,2>& T,
                           double& result,
                           MsqMatrix<2,2>& deriv_wrt_T,
                           MsqMatrix<2,2> second_wrt_T[3],
                           MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<3,3>& T, 
                 double& result, 
                 MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_grad( const MsqMatrix<3,3>& T,
                           double& result,
                           MsqMatrix<3,3>& deriv_wrt_T,
                           MsqError& err );
  
  MESQUITE_EXPORT virtual
  bool evaluate_with_hess( const MsqMatrix<3,3>& T,
                           double& result,
                           MsqMatrix<3,3>& deriv_wrt_T,
                           MsqMatrix<3,3> second_wrt_T[6],
                           MsqError& err );

private:
  template <unsigned D> inline
  bool eval( const MsqMatrix<D,D>& T, double& result, MsqError& err );
  template <unsigned D> inline
  bool grad( const MsqMatrix<D,D>& T, double& result, MsqMatrix<D,D>& first, MsqError& err );
  template <unsigned D> inline
  bool hess( const MsqMatrix<D,D>& T, double& result, MsqMatrix<D,D>& first, MsqMatrix<D,D>* second, MsqError& err );
};


} // namespace MESQUITE_NS

#endif
