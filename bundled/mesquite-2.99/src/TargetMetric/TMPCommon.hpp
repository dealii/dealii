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


/** \file TMPCommon.hpp
 *  \brief Common utility stuff for implementing target metrics
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TMP_COMMON_HPP
#define MSQ_TMP_COMMON_HPP

#include "Mesquite.hpp"

namespace MESQUITE_NS {

/**\def TMP_T_TEMPL_IMPL_COMMON(N)
 * \brief A macro that declares virtual evaluation functions for a \c TMetric subclass \c N
 * 
 * IF a \c TMetric class provides template functions named \c eval, \c grad, 
 * and \c hess, then this macro can be used to provide the trivial implementations
 * of the 2D and 3D \c evaluate, \c evaluate_with_grad, and \c evaluate_with_hess virtual
 * member functions in terms of those macros.
 */

#define TMP_T_TEMPL_IMPL_DIM(N,D) \
bool N::evaluate( const MsqMatrix<D,D>& T, double& r, MsqError& ) \
  { return eval( T, r ); } \
bool N::evaluate_with_grad( const MsqMatrix<D,D>& T, double& r, MsqMatrix<D,D>& d1, MsqError& ) \
  { return grad( T, r, d1 ); } \
bool N::evaluate_with_hess( const MsqMatrix<D,D>& T, double& r, MsqMatrix<D,D>& d1, MsqMatrix<D,D>* d2, MsqError& ) \
  { return hess( T, r, d1, d2 ); }

#define TMP_T_TEMPL_IMPL_COMMON(N) \
  TMP_T_TEMPL_IMPL_DIM(N,2) \
  TMP_T_TEMPL_IMPL_DIM(N,3) 


/**\def TMP_T_TEMPL_IMPL_COMMON_ERR(N)
 * \brief A macro that declares virtual evaluation functions for a \c TMetric subclass \c N
 * 
 * IF a \c TMetric class provides template functions named \c eval, \c grad, 
 * and \c hess, then this macro can be used to provide the trivial implementations
 * of the 2D and 3D \c evaluate, \c evaluate_with_grad, and \c evaluate_with_hess virtual
 * member functions in terms of those macros.
 * 
 * The difference between this macro and \c TMP_AW_TEMPL_IMPL_COMMON is that this
 * variation expects the \c eval, \c grad, and \c hess to accept an \c MsqError
 * as their last argument.  This variation is typically used for metrics that are
 * a function of some other metric.
 */

#define TMP_T_TEMPL_IMPL_ERR_DIM(N,D) \
bool N::evaluate( const MsqMatrix<D,D>& T, double& r, MsqError& err ) \
  { return eval( T, r, err ); } \
bool N::evaluate_with_grad( const MsqMatrix<D,D>& T, double& r, MsqMatrix<D,D>& d1, MsqError& err ) \
  { return grad( T, r, d1, err ); } \
bool N::evaluate_with_hess( const MsqMatrix<D,D>& T, double& r, MsqMatrix<D,D>& d1, MsqMatrix<D,D>* d2, MsqError& err ) \
  { return hess( T, r, d1, d2, err ); }

#define TMP_T_TEMPL_IMPL_COMMON_ERR(N) \
  TMP_T_TEMPL_IMPL_ERR_DIM(N,2) \
  TMP_T_TEMPL_IMPL_ERR_DIM(N,3) 

/**\def TMP_AW_TEMPL_IMPL_COMMON(N)
 * \brief A macro that declares virtual evaluation functions for a \c AWMetric subclass \c N
 * 
 * IF a \c AWMetric class provides template functions named \c eval, \c grad, 
 * and \c hess, then this macro can be used to provide the trivial implementations
 * of the 2D and 3D \c evaluate, \c evaluate_with_grad, and \c evaluate_with_hess virtual
 * member functions in terms of those macros.
 */

#define TMP_AW_TEMPL_IMPL_DIM(N,D) \
bool N::evaluate( const MsqMatrix<D,D>& A, const MsqMatrix<D,D>& W, double& r, MsqError& ) \
  { return eval( A, W, r ); } \
bool N::evaluate_with_grad( const MsqMatrix<D,D>& A, const MsqMatrix<D,D>& W, double& r, MsqMatrix<D,D>& d1, MsqError& ) \
  { return grad( A, W, r, d1 ); } \
bool N::evaluate_with_hess( const MsqMatrix<D,D>& A, const MsqMatrix<D,D>& W, double& r, MsqMatrix<D,D>& d1, MsqMatrix<D,D>* d2, MsqError& ) \
  { return hess( A, W, r, d1, d2 ); }

#define TMP_AW_TEMPL_IMPL_COMMON(N) \
  TMP_AW_TEMPL_IMPL_DIM(N,2) \
  TMP_AW_TEMPL_IMPL_DIM(N,3) 

/**\def TMP_AW_TEMPL_IMPL_COMMON(N)
 * \brief Like TMP_AW_TEMPL_IMPL_COMMON, except no implementation of 2nd derivs
 */

#define TMP_AW_TEMPL_IMPL_NO2ND_DIM(N,D) \
bool N::evaluate( const MsqMatrix<D,D>& A, const MsqMatrix<D,D>& W, double& r, MsqError& ) \
  { return eval( A, W, r ); } \
bool N::evaluate_with_grad( const MsqMatrix<D,D>& A, const MsqMatrix<D,D>& W, double& r, MsqMatrix<D,D>& d1, MsqError& ) \
  { return grad( A, W, r, d1 ); } 

#define TMP_AW_TEMPL_IMPL_COMMON_NO2ND(N) \
  TMP_AW_TEMPL_IMPL_NO2ND_DIM(N,2) \
  TMP_AW_TEMPL_IMPL_NO2ND_DIM(N,3) 

/**\brief Dimension-specific constants
 * 
 * Provide constants that depend on a template dimension parameter.
 * For use in implementing targe metrics for which the target metric
 * implementation is templatized on the dimension of of the matrix.
 * 
 * In optimized code these should reduce to literal constants.
 */
template <unsigned D> struct DimConst {};
template <> struct DimConst<2> {
  static inline double sqrt() { return MSQ_SQRT_TWO; }
  static inline double inv()  { return 0.5; }
};
template <> struct DimConst<3> {
  static inline double sqrt() { return MSQ_SQRT_THREE; }
  static inline double inv()  { return MSQ_ONE_THIRD; }
};



} // namespace MESQUITE_NS

#endif
