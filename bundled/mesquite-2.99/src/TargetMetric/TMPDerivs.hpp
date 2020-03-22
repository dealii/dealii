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


/** \file TMPDerivs.hpp
 *  \brief Utility methods for calculating derivatives of TMP metrics
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TMP_DERIVS_HPP
#define MSQ_TMP_DERIVS_HPP

#include "Mesquite.hpp"
#include "MsqMatrix.hpp"

/** If defined, then \f$ (A \otimes B)_{i,j} = row(A,i)^t \otimes row(B,j) \f$
    otherwise \f$ (A \otimes B)_{i,j} = column(A,i) \otimes column(B,j)^t \f$
  */
#define MSQ_ROW_BASED_OUTER_PRODUCT

namespace MESQUITE_NS {

/**\brief \f$ R *= s \f$ */
template <unsigned D> inline
void hess_scale_t( MsqMatrix<D,D> R[D*(D+1)/2], double alpha );

inline
void hess_scale( MsqMatrix<3,3> R[6], double alpha )
{ hess_scale_t<3>( R, alpha ); }

inline
void hess_scale( MsqMatrix<2,2> R[3], double alpha )
{ hess_scale_t<2>( R, alpha ); }

/**\brief \f$ R = \alpha I_9 \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
inline
void set_scaled_I( MsqMatrix<3,3> R[6], double alpha );

/**\brief \f$ R += \alpha I_9 \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
inline
void pluseq_scaled_I( MsqMatrix<3,3> R[6], double alpha );

/**\brief \f$ R += \alpha I_9 \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
inline
void pluseq_scaled_I( MsqMatrix<3,3>& R, double alpha );

/**\brief \f$ R = \alpha I_4 \f$
 *
 *\param R The 3 blocks of the upper triangular portion of a 4x4
 *         symmetric matrix.
 */
inline
void set_scaled_I( MsqMatrix<2,2> R[3], double alpha );

/**\brief \f$ R += \alpha I_4 \f$
 *
 *\param R The 3 blocks of the upper triangular portion of a 4x4
 *         symmetric matrix.
 */
inline
void pluseq_scaled_I( MsqMatrix<2,2> R[3], double alpha );

/**\brief \f$ R += \alpha I_4 \f$
 *
 *\param R The 3 blocks of the upper triangular portion of a 4x4
 *         symmetric matrix.
 */
inline
void pluseq_scaled_I( MsqMatrix<2,2>& R, double alpha );

/**\brief \f$ R += \alpha \frac{\partial}{\partial T}det(T) \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
inline
void pluseq_scaled_2nd_deriv_of_det( MsqMatrix<3,3> R[6], 
                                     double alpha,
                                     const MsqMatrix<3,3>& T );

/**\brief \f$ R += \alpha \frac{\partial}{\partial T}det(T) \f$
 *
 *\param R The 3 blocks of the upper triangular portion of a 4x4
 *         symmetric matrix.
 */
inline
void pluseq_scaled_2nd_deriv_of_det( MsqMatrix<2,2> R[3], double alpha );
inline
void pluseq_scaled_2nd_deriv_of_det( MsqMatrix<2,2> R[3], double alpha, const MsqMatrix<2,2>& )
  { pluseq_scaled_2nd_deriv_of_det( R, alpha ); }

/**\brief \f$ R = \alpha \frac{\partial}{\partial T}det(T) \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
inline
void set_scaled_2nd_deriv_of_det( MsqMatrix<3,3> R[6],
                                  double alpha,
                                  const MsqMatrix<3,3>& T );

/**\brief \f$ R = \alpha \frac{\partial}{\partial T}det(T) \f$
 *
 *\param R The 3 blocks of the upper triangular portion of a 4x4
 *         symmetric matrix.
 */
inline
void set_scaled_2nd_deriv_of_det( MsqMatrix<2,2> R[3], double alpha );
inline
void set_scaled_2nd_deriv_of_det( MsqMatrix<2,2> R[3], double alpha, const MsqMatrix<2,2>& )
  { set_scaled_2nd_deriv_of_det( R, alpha ); }

/**\brief \f$ R += \alpha \frac{\partial^2}{\partial T^2}tr(adj T) \f$
 *
 * Note: 2nd deriv is a constant, so T is not passed as an argument.
 */
inline 
void pluseq_scaled_2nd_deriv_tr_adj( MsqMatrix<2,2> R[3], double alpha );

/**\brief \f$ R += \alpha \frac{\partial^2}{\partial T^2}tr(adj T) \f$
 *
 * Note: 2nd deriv is a constant, so T is not passed as an argument.
 */
inline 
void pluseq_scaled_2nd_deriv_tr_adj( MsqMatrix<3,3> R[6], double alpha );

/**\brief \f$ R += \alpha \frac{\partial^2}{\partial T^2}|adj T|^2 \f$
  */
inline 
void set_scaled_2nd_deriv_norm_sqr_adj( MsqMatrix<2,2> R[3], 
                                        double alpha,
                                        const MsqMatrix<2,2>& T );

/**\brief \f$ R += \alpha \frac{\partial^2}{\partial T^2}|adj T|^2 \f$
  */
inline 
void set_scaled_2nd_deriv_norm_sqr_adj( MsqMatrix<3,3> R[6], 
                                        double alpha,
                                        const MsqMatrix<3,3>& T );

/**\brief \f$ R += \alpha \left( M \otimes M \right) \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
template <unsigned D> inline
void pluseq_scaled_outer_product_t( MsqMatrix<D,D> R[D*(D+1)/2],
                                    double alpha,
                                    const MsqMatrix<D,D>& M );
inline
void pluseq_scaled_outer_product( MsqMatrix<3,3> R[6], 
                                  double alpha,
                                  const MsqMatrix<3,3>& M )
{ pluseq_scaled_outer_product_t<3>(R,alpha,M); }

inline
void pluseq_scaled_outer_product( MsqMatrix<2,2> R[3], 
                                  double alpha,
                                  const MsqMatrix<2,2>& M )
{ pluseq_scaled_outer_product_t<2>(R,alpha,M); }

/**\brief \f$ R = \alpha \left( M \otimes M \right) \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
template <unsigned D> inline
void set_scaled_outer_product_t( MsqMatrix<D,D> R[D*(D+1)/2],
                                 double alpha,
                                 const MsqMatrix<D,D>& M );

inline
void set_scaled_outer_product( MsqMatrix<3,3> R[6],
                               double alpha,
                               const MsqMatrix<3,3>& M )
{ set_scaled_outer_product_t<3>( R, alpha, M ); }

inline
void set_scaled_outer_product( MsqMatrix<2,2> R[3],
                               double alpha,
                               const MsqMatrix<2,2>& M )
{ set_scaled_outer_product_t<2>( R, alpha, M ); }

/**\brief \f$ R += \alpha \left( A \otimes B + B \otimes A \right) \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
inline
void pluseq_scaled_sum_outer_product( MsqMatrix<3,3> R[6],
                                      double alpha,
                                      const MsqMatrix<3,3>& A,
                                      const MsqMatrix<3,3>& B );

/**\brief \f$ R += \alpha \left( A \otimes B + B \otimes A \right) \f$
 *
 *\param R The 3 blocks of the upper triangular portion of a 4x4
 *         symmetric matrix.
 */
inline
void pluseq_scaled_sum_outer_product( MsqMatrix<2,2> R[3],
                                      double alpha,
                                      const MsqMatrix<2,2>& A,
                                      const MsqMatrix<2,2>& B );

/**\brief \f$ R = \alpha \left( A \otimes B + B \otimes A \right) \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
inline
void set_scaled_sum_outer_product( MsqMatrix<3,3> R[6],
                                   double alpha,
                                   const MsqMatrix<3,3>& A,
                                   const MsqMatrix<3,3>& B );

/**\brief \f$ R = \alpha \left( A \otimes B + B \otimes A \right) \f$
 *
 *\param R The 3 blocks of the upper triangular portion of a 4x4
 *         symmetric matrix.
 */
inline
void set_scaled_sum_outer_product( MsqMatrix<2,2> R[3],
                                   double alpha,
                                   const MsqMatrix<2,2>& A,
                                   const MsqMatrix<2,2>& B );

/**\brief \f$ R += \alpha (I \otimes I) \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
inline 
void pluseq_scaled_outer_product_I_I( MsqMatrix<3,3> R[6], double alpha );

/**\brief \f$ R += \alpha (I \otimes I) \f$
 *
 *\param R The 3 blocks of the upper triangular portion of a 4x4
 *         symmetric matrix.
 */
inline 
void pluseq_scaled_outer_product_I_I( MsqMatrix<2,2> R[3], double alpha );

/**\brief \f$ R += I \otimes A \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
inline 
void pluseq_I_outer_product( MsqMatrix<3,3> R[6], const MsqMatrix<3,3>& A );

/**\brief \f$ R += A \otimes I \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
inline 
void pluseq_outer_product_I( MsqMatrix<3,3> R[6], const MsqMatrix<3,3>& A );

/**\brief \f$ R += \alpha \left( I \otimes A + A \otimes I \right) \f$
 *
 *\param R The 6 blocks of the upper triangular portion of a 9x9
 *         symmetric matrix.
 */
template <unsigned D> inline
void pluseq_scaled_sum_outer_product_I_t( MsqMatrix<D,D> R[D*(D+1)/2],
                                        double alpha,
                                        const MsqMatrix<D,D>& A );
inline
void pluseq_scaled_sum_outer_product_I( MsqMatrix<3,3> R[6],
                                        double alpha,
                                        const MsqMatrix<3,3>& A_in )
{ pluseq_scaled_sum_outer_product_I_t<3>( R, alpha, A_in ); }

inline
void pluseq_scaled_sum_outer_product_I( MsqMatrix<2,2> R[2],
                                        double alpha,
                                        const MsqMatrix<2,2>& A_in )
{ pluseq_scaled_sum_outer_product_I_t<2>( R, alpha, A_in ); }


/**\brief \f$ \frac{\partial^2 f}{\partial (AZ)^2} \Rightarrow \frac{\partial^2 f}{\partial A^2} \f$
 *
 * Given the second derivatives of a function with respect to
 * a matrix product, calculate the derivatives of the function
 * with respect to the first matrix in the product.
 */
template <unsigned D> inline
void second_deriv_wrt_product_factor_t( MsqMatrix<D,D> R[D*(D+1)/2],
                                        const MsqMatrix<D,D>& Z );

inline
void second_deriv_wrt_product_factor( MsqMatrix<3,3> R[6],
                                      const MsqMatrix<3,3>& Z )
{ second_deriv_wrt_product_factor_t<3>( R, Z ); }

inline
void second_deriv_wrt_product_factor( MsqMatrix<2,2> R[3],
                                      const MsqMatrix<2,2>& Z )
{ second_deriv_wrt_product_factor_t<2>( R, Z ); }

/**\brief \f$ R = \alpha * \frac{\partial^2 \psi(T)}{\partial T^2} \f$
 *
 * \f$ \psi(T) = \sqrt{ |T|^2 + 2 \tau } \f$
 */
inline
void set_scaled_2nd_deriv_wrt_psi( MsqMatrix<2,2> R[3],
                                   const double alpha,
                                   const double psi,
                                   const MsqMatrix<2,2>& T );

/**\brief \f$  R = R + \alpha * Z \f$
 */
template <unsigned D> inline
void pluseq_scaled_t( MsqMatrix<D,D> R[D*(D+1)/2],
                      double alpha,
                      const MsqMatrix<D,D> Z[D*(D+1)/2] );


void set_scaled_I( MsqMatrix<3,3> R[6], double alpha )
{
  R[0] = R[3] = R[5] = MsqMatrix<3,3>(alpha);
  R[1] = R[2] = R[4] = MsqMatrix<3,3>(0.0);
}

void pluseq_scaled_I( MsqMatrix<3,3>& R, double alpha )
{
  R(0,0) += alpha;
  R(1,1) += alpha;
  R(2,2) += alpha;
}

void pluseq_scaled_I( MsqMatrix<2,2>& R, double alpha )
{
  R(0,0) += alpha;
  R(1,1) += alpha;
}

void pluseq_scaled_I( MsqMatrix<3,3> R[6], double alpha )
{
  pluseq_scaled_I( R[0], alpha );
  pluseq_scaled_I( R[3], alpha );
  pluseq_scaled_I( R[5], alpha );
}

void set_scaled_I( MsqMatrix<2,2> R[3], double alpha )
{
  R[0] = R[2] = MsqMatrix<2,2>(alpha);
  R[1]        = MsqMatrix<2,2>(0.0);
}

void pluseq_scaled_I( MsqMatrix<2,2> R[3], double alpha )
{
  pluseq_scaled_I( R[0], alpha );
  pluseq_scaled_I( R[2], alpha );
}

void pluseq_scaled_2nd_deriv_of_det( MsqMatrix<3,3> R[6],
                                     double alpha,
                                     const MsqMatrix<3,3>& T )
{
  MsqMatrix<3,3> A(T);
  A *= alpha;

  R[1](0,1) += A(2,2);
  R[1](1,0) -= A(2,2);
  R[1](0,2) -= A(2,1);
  R[1](2,0) += A(2,1);
  R[1](1,2) += A(2,0);
  R[1](2,1) -= A(2,0);
  
  R[2](0,1) -= A(1,2);
  R[2](1,0) += A(1,2);
  R[2](0,2) += A(1,1);
  R[2](2,0) -= A(1,1);
  R[2](1,2) -= A(1,0);
  R[2](2,1) += A(1,0);
 
  R[4](0,1) += A(0,2);
  R[4](1,0) -= A(0,2);
  R[4](0,2) -= A(0,1);
  R[4](2,0) += A(0,1);
  R[4](1,2) += A(0,0);
  R[4](2,1) -= A(0,0);
}

void set_scaled_2nd_deriv_of_det( MsqMatrix<3,3> R[6],
                                  double alpha,
                                  const MsqMatrix<3,3>& T )
{
  MsqMatrix<3,3> A(T);
  A *= alpha;

  R[0] = R[3] = R[5] = MsqMatrix<3,3>(0.0);

  R[1](0,0) = R[1](1,1) = R[1](2,2) = 0;
  R[1](0,1) =  A(2,2);
  R[1](1,0) = -A(2,2);
  R[1](0,2) = -A(2,1);
  R[1](2,0) =  A(2,1);
  R[1](1,2) =  A(2,0);
  R[1](2,1) = -A(2,0);
  
  R[2](0,0) = R[2](1,1) = R[2](2,2) = 0;
  R[2](0,1) = -A(1,2);
  R[2](1,0) =  A(1,2);
  R[2](0,2) =  A(1,1);
  R[2](2,0) = -A(1,1);
  R[2](1,2) = -A(1,0);
  R[2](2,1) =  A(1,0);
 
  R[4](0,0) = R[4](1,1) = R[4](2,2) = 0;
  R[4](0,1) =  A(0,2);
  R[4](1,0) = -A(0,2);
  R[4](0,2) = -A(0,1);
  R[4](2,0) =  A(0,1);
  R[4](1,2) =  A(0,0);
  R[4](2,1) = -A(0,0);
}

void pluseq_scaled_2nd_deriv_of_det( MsqMatrix<2,2> R[3], double alpha )
{
  R[1](0,1) += alpha;
  R[1](1,0) -= alpha;
}

void set_scaled_2nd_deriv_of_det( MsqMatrix<2,2> R[3], double alpha )
{
  R[0] = R[2] = MsqMatrix<2,2>(0.0);
  R[1](0,0) =  0.0;
  R[1](0,1) =  alpha;
  R[1](1,0) = -alpha;
  R[1](1,1) =  0.0;
  
}

void pluseq_scaled_2nd_deriv_tr_adj( MsqMatrix<2,2>*, double  )
{
  // 2nd derivative is zero
}

void pluseq_scaled_2nd_deriv_tr_adj( MsqMatrix<3,3> R[6], double alpha )
{
  R[1](0,1) += alpha;
  R[1](1,0) -= alpha;
  R[2](0,2) += alpha;
  R[2](2,0) -= alpha;
  R[4](1,2) += alpha;
  R[4](2,1) -= alpha;
}

void set_scaled_2nd_deriv_norm_sqr_adj( MsqMatrix<2,2> R[3], 
                                        double alpha,
                                        const MsqMatrix<2,2>&  )
{
  set_scaled_I( R, 2*alpha );
}

void set_scaled_2nd_deriv_norm_sqr_adj( MsqMatrix<3,3> R[6], 
                                        double alpha,
                                        const MsqMatrix<3,3>& T )
{
  set_scaled_outer_product( R, 1, T );
  double tmp01, tmp02, tmp12;
  //R[1] = 2*R[1] - transpose(R[1]);
  tmp01 = R[1](0,1);
  tmp02 = R[1](0,2);
  tmp12 = R[1](1,2);
  R[1](0,1) = 2 * R[1](0,1) - R[1](1,0);
  R[1](0,2) = 2 * R[1](0,2) - R[1](2,0);
  R[1](1,2) = 2 * R[1](1,2) - R[1](2,1);
  R[1](1,0) = 2 * R[1](1,0) - tmp01;
  R[1](2,0) = 2 * R[1](2,0) - tmp02;
  R[1](2,1) = 2 * R[1](2,1) - tmp12;
  //R[2] = 2*R[2] - transpose(R[1]);
  tmp01 = R[2](0,1);
  tmp02 = R[2](0,2);
  tmp12 = R[2](1,2);
  R[2](0,1) = 2 * R[2](0,1) - R[2](1,0);
  R[2](0,2) = 2 * R[2](0,2) - R[2](2,0);
  R[2](1,2) = 2 * R[2](1,2) - R[2](2,1);
  R[2](1,0) = 2 * R[2](1,0) - tmp01;
  R[2](2,0) = 2 * R[2](2,0) - tmp02;
  R[2](2,1) = 2 * R[2](2,1) - tmp12;
  //R[4] = 2*R[4] - transpose(R[1]);
  tmp01 = R[4](0,1);
  tmp02 = R[4](0,2);
  tmp12 = R[4](1,2);
  R[4](0,1) = 2 * R[4](0,1) - R[4](1,0);
  R[4](0,2) = 2 * R[4](0,2) - R[4](2,0);
  R[4](1,2) = 2 * R[4](1,2) - R[4](2,1);
  R[4](1,0) = 2 * R[4](1,0) - tmp01;
  R[4](2,0) = 2 * R[4](2,0) - tmp02;
  R[4](2,1) = 2 * R[4](2,1) - tmp12;
  
  const MsqMatrix<3,3> TpT = transpose(T) * T;
  R[0] -= TpT;
  R[3] -= TpT;
  R[5] -= TpT;
  
  const MsqMatrix<3,3> TTp = T * transpose(T);
  const double ns = sqr_Frobenius(T);
  R[0](0,0) += ns - TTp(0,0);
  R[0](1,1) += ns - TTp(0,0);
  R[0](2,2) += ns - TTp(0,0);
  R[1](0,0) +=    - TTp(0,1);
  R[1](1,1) +=    - TTp(0,1);
  R[1](2,2) +=    - TTp(0,1);
  R[2](0,0) +=    - TTp(0,2);
  R[2](1,1) +=    - TTp(0,2);
  R[2](2,2) +=    - TTp(0,2);
  R[3](0,0) += ns - TTp(1,1);
  R[3](1,1) += ns - TTp(1,1);
  R[3](2,2) += ns - TTp(1,1);
  R[4](0,0) +=    - TTp(1,2);
  R[4](1,1) +=    - TTp(1,2);
  R[4](2,2) +=    - TTp(1,2);
  R[5](0,0) += ns - TTp(2,2);
  R[5](1,1) += ns - TTp(2,2);
  R[5](2,2) += ns - TTp(2,2);
  
  alpha *= 2;
  R[0] *= alpha;
  R[1] *= alpha;
  R[2] *= alpha;
  R[3] *= alpha;
  R[4] *= alpha;
  R[5] *= alpha;
}

#ifdef MSQ_ROW_BASED_OUTER_PRODUCT
template <unsigned D>
void pluseq_scaled_outer_product_t( MsqMatrix<D,D> R[D*(D+1)/2],
                                  double alpha,
                                  const MsqMatrix<D,D>& M )
{
  MsqMatrix<D,D> aM(M);
  aM *= alpha;
  unsigned h = 0;
  for (unsigned i = 0; i < D; ++i)
    for (unsigned j = i; j < D; ++j)
      R[h++] += transpose(aM.row(i)) * M.row(j);
}
#else
template <unsigned D>
void pluseq_scaled_outer_product_t( MsqMatrix<D,D> R[D*(D+1)/2],
                                  double alpha,
                                  const MsqMatrix<D,D>& M )
{
  MsqMatrix<D,D> aM(transpose(M));
  aM *= alpha;
  unsigned h = 0;
  for (unsigned i = 0; i < D; ++i)
    for (unsigned j = i; j < D; ++j)
      R[h++] += M.column(i) * aM.row(0);
}
#endif

#ifdef MSQ_ROW_BASED_OUTER_PRODUCT
template <unsigned D>
void set_scaled_outer_product_t( MsqMatrix<D,D> R[D*(D+1)/2],
                               double alpha,
                               const MsqMatrix<D,D>& M )
{
  MsqMatrix<D,D> aM(M);
  aM *= alpha;
  unsigned h = 0;
  for (unsigned i = 0; i < D; ++i)
    for (unsigned j = i; j < D; ++j)
      R[h++] = transpose(aM.row(i)) * M.row(j);
}
#else
template <unsigned D>
void set_scaled_outer_product_t( MsqMatrix<D,D> R[D*(D+1)/2],
                               double alpha,
                               const MsqMatrix<D,D>& M )
{
  MsqMatrix<D,D> aM(transpose(M));
  aM *= alpha;
  unsigned h = 0;
  for (unsigned i = 0; i < D; ++i)
    for (unsigned j = i; j < D; ++j)
      R[h++] = M.column(i) * aM.row(j);
}
#endif


#ifdef MSQ_ROW_BASED_OUTER_PRODUCT
void pluseq_scaled_sum_outer_product( MsqMatrix<3,3> R[6],
                                      double alpha,
                                      const MsqMatrix<3,3>& A_in,
                                      const MsqMatrix<3,3>& B )
{
    // apply scalar first
  MsqMatrix<3,3> A(A_in), tmp;
  A *= alpha;

    // block 0,0
  tmp = transpose(A.row(0)) * B.row(0);
  R[0] += tmp;
  R[0] += transpose(tmp);
  
    // block 1,1
  tmp = transpose(A.row(1)) * B.row(1);
  R[3] += tmp;
  R[3] += transpose(tmp);
  
    // block 2,2
  tmp = transpose(A.row(2)) * B.row(2);
  R[5] += tmp;
  R[5] += transpose(tmp);
   
    // block 0,1
  R[1] += transpose(A.row(0)) * B.row(1) +
          transpose(B.row(0)) * A.row(1);
  
    // block 0,2
  R[2] += transpose(A.row(0)) * B.row(2) +
          transpose(B.row(0)) * A.row(2);
  
    // block 1,2
  R[4] += transpose(A.row(1)) * B.row(2) +
          transpose(B.row(1)) * A.row(2);
}
void set_scaled_sum_outer_product( MsqMatrix<3,3> R[6],
                                   double alpha,
                                   const MsqMatrix<3,3>& A_in,
                                   const MsqMatrix<3,3>& B )
{
    // apply scalar first
  MsqMatrix<3,3> A(A_in);
  A *= alpha;

    // block 0,0
  R[0] =  transpose(A.row(0)) * B.row(0);
  R[0] += transpose(R[0]);
  
    // block 1,1
  R[3] =  transpose(A.row(1)) * B.row(1);
  R[3] += transpose(R[3]);
  
    // block 2,2
  R[5] =  transpose(A.row(2)) * B.row(2);
  R[5] += transpose(R[5]);
   
    // block 0,1
  R[1] =  transpose(A.row(0)) * B.row(1);
  R[1] += transpose(B.row(0)) * A.row(1);
  
    // block 0,2
  R[2] =  transpose(A.row(0)) * B.row(2);
  R[2] += transpose(B.row(0)) * A.row(2);
  
    // block 1,2
  R[4] =  transpose(A.row(1)) * B.row(2);
  R[4] += transpose(B.row(1)) * A.row(2);
}
#else
void pluseq_scaled_sum_outer_product( MsqMatrix<3,3> R[6],
                                      double alpha,
                                      const MsqMatrix<3,3>& A_in,
                                      const MsqMatrix<3,3>& B )
{
    // apply scalar first
  MsqMatrix<3,3> A(A_in), tmp;
  A *= alpha;

    // block 0,0
  tmp = A.column(0) * transpose(B.column(0));
  R[0] += tmp;
  R[0] += transpose(tmp);
  
    // block 1,1
  tmp = A.column(1) * transpose(B.column(1));
  R[3] += tmp;
  R[3] += transpose(tmp);
  
    // block 2,2
  tmp = A.column(2) * transpose(B.column(2));
  R[5] += tmp;
  R[5] += transpose(tmp);
   
    // block 0,1
  R[1] += A.column(0) * transpose(B.column(1)) +
          B.column(0) * transpose(A.column(1));
  
    // block 0,2
  R[2] += A.column(0) * transpose(B.column(2)) +
          B.column(0) * transpose(A.column(2));
  
    // block 1,2
  R[4] += A.column(1) * transpose(B.column(2)) +
          B.column(1) * transpose(A.column(2));
}
void set_scaled_sum_outer_product( MsqMatrix<3,3> R[6],
                                   double alpha,
                                   const MsqMatrix<3,3>& A_in,
                                   const MsqMatrix<3,3>& B )
{
    // apply scalar first
  MsqMatrix<3,3> A(A_in);
  A *= alpha;

    // block 0,0
  R[0] = A.column(0) * transpose(B.column(0));
  R[0] += transpose(R[0]);
  
    // block 1,1
  R[3] = A.column(1) * transpose(B.column(1));
  R[3] += transpose(R[3]);
  
    // block 2,2
  R[5] = A.column(2) * transpose(B.column(2));
  R[5] += transpose(R[5]);
   
    // block 0,1
  R[1] =  A.column(0) * transpose(B.column(1));
  R[1] += B.column(0) * transpose(A.column(1));
  
    // block 0,2
  R[2] =  A.column(0) * transpose(B.column(2));
  R[2] += B.column(0) * transpose(A.column(2));
  
    // block 1,2
  R[4] =  A.column(1) * transpose(B.column(2));
  R[4] += B.column(1) * transpose(A.column(2));
}
#endif

#ifdef MSQ_ROW_BASED_OUTER_PRODUCT
void pluseq_scaled_sum_outer_product( MsqMatrix<2,2> R[3],
                                      double alpha,
                                      const MsqMatrix<2,2>& A_in,
                                      const MsqMatrix<2,2>& B )
{
    // apply scalar first
  MsqMatrix<2,2> A(A_in), tmp;
  A *= alpha;

    // block 0,0
  tmp = transpose(A.row(0)) * B.row(0);
  R[0] += tmp;
  R[0] += transpose(tmp);
  
    // block 1,1
  tmp = transpose(A.row(1)) * B.row(1);
  R[2] += tmp;
  R[2] += transpose(tmp);
   
    // block 0,1
  R[1] += transpose(A.row(0)) * B.row(1) +
          transpose(B.row(0)) * A.row(1);
}
void set_scaled_sum_outer_product( MsqMatrix<2,2> R[3],
                                   double alpha,
                                   const MsqMatrix<2,2>& A_in,
                                   const MsqMatrix<2,2>& B )
{
    // apply scalar first
  MsqMatrix<2,2> A(A_in);
  A *= alpha;

    // block 0,0
  R[0] =  transpose(A.row(0)) * B.row(0);
  R[0] += transpose(R[0]);
  
    // block 1,1
  R[2] =  transpose(A.row(1)) * B.row(1);
  R[2] += transpose(R[2]);
   
    // block 0,1
  R[1] =  transpose(A.row(0)) * B.row(1);
  R[1] += transpose(B.row(0)) * A.row(1);
}
#else
void pluseq_scaled_sum_outer_product( MsqMatrix<2,2> R[3],
                                      double alpha,
                                      const MsqMatrix<2,2>& A_in,
                                      const MsqMatrix<2,2>& B )
{
    // apply scalar first
  MsqMatrix<2,2> A(A_in), tmp;
  A *= alpha;

    // block 0,0
  tmp = A.column(0) * transpose(B.column(0));
  R[0] += tmp;
  R[0] += transpose(tmp);
  
    // block 1,1
  tmp = A.column(1) * transpose(B.column(1));
  R[2] += tmp;
  R[2] += transpose(tmp);
   
    // block 0,1
  R[1] += A.column(0) * transpose(B.column(1)) +
          B.column(0) * transpose(A.column(1));
}
void set_scaled_sum_outer_product( MsqMatrix<2,2> R[3],
                                   double alpha,
                                   const MsqMatrix<2,2>& A_in,
                                   const MsqMatrix<2,2>& B )
{
    // apply scalar first
  MsqMatrix<2,2> A(A_in);
  A *= alpha;

    // block 0,0
  R[0] =  A.column(0) * transpose(B.column(0));
  R[0] += transpose(R[0]);
  
    // block 1,1
  R[2] =  A.column(1) * transpose(B.column(1));
  R[2] += transpose(R[2]);
   
    // block 0,1
  R[1] =  A.column(0) * transpose(B.column(1));
  R[1] += B.column(0) * transpose(A.column(1));
}
#endif

void pluseq_scaled_outer_product_I_I( MsqMatrix<3,3> R[6], double alpha )
{
  R[0](0,0) += alpha;
  R[1](0,1) += alpha;
  R[2](0,2) += alpha;
  R[3](1,1) += alpha;
  R[4](1,2) += alpha;
  R[5](2,2) += alpha;
}

void pluseq_scaled_outer_product_I_I( MsqMatrix<2,2> R[3], double alpha )
{
  R[0](0,0) += alpha;
  R[1](0,1) += alpha;
  R[2](1,1) += alpha;
}

#ifdef MSQ_ROW_BASED_OUTER_PRODUCT
inline void pluseq_I_outer_product( MsqMatrix<3,3> R[6], const MsqMatrix<3,3>& A )
{
  R[0].add_row( 0, A.row(0) );
  R[1].add_row( 0, A.row(1) );
  R[2].add_row( 0, A.row(2) );
  R[3].add_row( 1, A.row(1) );
  R[4].add_row( 1, A.row(2) );
  R[5].add_row( 2, A.row(2) );
}
inline void pluseq_I_outer_product( MsqMatrix<2,2> R[3], const MsqMatrix<2,2>& A )
{
  R[0].add_row( 0, A.row(0) );
  R[1].add_row( 0, A.row(1) );
  R[2].add_row( 1, A.row(1) );
}
#else
inline void pluseq_I_outer_product( MsqMatrix<3,3> R[6], const MsqMatrix<3,3>& A )
{
  R[0].add_row( 0, transpose(A.column(0)) );
  R[1].add_row( 0, transpose(A.column(1)) );
  R[2].add_row( 0, transpose(A.column(2)) );
  R[3].add_row( 1, transpose(A.column(1)) );
  R[4].add_row( 1, transpose(A.column(2)) );
  R[5].add_row( 2, transpose(A.column(2)) );
}
inline void pluseq_I_outer_product( MsqMatrix<2,2> R[3], const MsqMatrix<2,2>& A )
{
  R[0].add_row( 0, transpose(A.column(0)) );
  R[1].add_row( 0, transpose(A.column(1)) );
  R[2].add_row( 1, transpose(A.column(1)) );
}
#endif

#ifdef MSQ_ROW_BASED_OUTER_PRODUCT
inline void pluseq_outer_product_I( MsqMatrix<3,3> R[6], const MsqMatrix<3,3>& A )
{
  R[0].add_column( 0, transpose(A.row(0)) );
  R[1].add_column( 1, transpose(A.row(0)) );
  R[2].add_column( 2, transpose(A.row(0)) );
  R[3].add_column( 1, transpose(A.row(1)) );
  R[4].add_column( 2, transpose(A.row(1)) );
  R[5].add_column( 2, transpose(A.row(2)) );
}
inline void pluseq_outer_product_I( MsqMatrix<2,2> R[3], const MsqMatrix<2,2>& A )
{
  R[0].add_column( 0, transpose(A.row(0)) );
  R[1].add_column( 1, transpose(A.row(0)) );
  R[2].add_column( 1, transpose(A.row(1)) );
}
#else
inline void pluseq_outer_product_I( MsqMatrix<3,3> R[6], const MsqMatrix<3,3>& A )
{
  R[0].add_column( 0, A.column(0) );
  R[1].add_column( 1, A.column(0) );
  R[2].add_column( 2, A.column(0) );
  R[3].add_column( 1, A.column(1) );
  R[4].add_column( 2, A.column(1) );
  R[5].add_column( 2, A.column(2) );
}
inline void pluseq_outer_product_I( MsqMatrix<2,2> R[3], const MsqMatrix<2,2>& A )
{
  R[0].add_column( 0, A.column(0) );
  R[1].add_column( 1, A.column(0) );
  R[2].add_column( 1, A.column(1) );
}
#endif

template <unsigned D>
void pluseq_scaled_sum_outer_product_I_t( MsqMatrix<D,D> R[D*(D+1)/2],
                                        double alpha,
                                        const MsqMatrix<D,D>& A_in )
{
    // apply scalar first
  MsqMatrix<D,D> A(A_in);
  A *= alpha;
  pluseq_I_outer_product( R ,A );
  pluseq_outer_product_I( R, A );
}

template <unsigned D>
void second_deriv_wrt_product_factor_t( MsqMatrix<D,D> R[D*(D+1)/2],
                                      const MsqMatrix<D,D>& Z )
{
  const MsqMatrix<D,D> Zt = transpose(Z);
  for (unsigned i = 0; i < D*(D+1)/2; ++i)
    R[i] = (Z * R[i]) * Zt;
}

void set_scaled_2nd_deriv_wrt_psi( MsqMatrix<2,2> R[3],
                                   const double alpha,
                                   const double psi,
                                   const MsqMatrix<2,2>& T )
{
  const double t = trace(T);
  const double f = alpha/(psi*psi*psi);
  const double s = T(0,1) - T(1,0);
  R[0](0,0) = R[1](0,1) = R[2](1,1) =  f*s*s;
  R[0](0,1) = R[0](1,0) = R[1](1,1) = -f*s*t;
  R[1](0,0) = R[2](0,1) = R[2](1,0) =  f*s*t;
  R[0](1,1) = R[2](0,0) = -(R[1](1,0) = -f*t*t);
}

template <unsigned D> inline
void pluseq_scaled( MsqMatrix<D,D> R[D*(D+1)/2],
                    double alpha,
                    const MsqMatrix<D,D> Z[D*(D+1)/2] )
{
  for (unsigned i = 0; i < D*(D+1)/2; ++i)
    R[i] += alpha * Z[i];
}


/**\brief \f$ R *= s */
template <unsigned D> inline
void hess_scale_t( MsqMatrix<D,D> R[D*(D+1)/2], double alpha )
{
  for (unsigned i = 0; i < D*(D+1)/2; ++i)
    R[i] *= alpha;
}

} // namespace MESQUITE_NS

#endif
