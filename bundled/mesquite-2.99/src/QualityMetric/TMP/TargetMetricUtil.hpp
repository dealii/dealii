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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TargetMetricUtil.hpp
 *  \brief A collection of utility code used by QualtiyMetrics
 *         composed of TMP Target Metrics
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_METRIC_UTIL_HPP
#define MSQ_TARGET_METRIC_UTIL_HPP

#include "Mesquite.hpp"
#include "SymMatrix3D.hpp"
#include <vector>
#include <assert.h>

namespace MESQUITE_NS {

template <unsigned R, unsigned C> class MsqMatrix;
template <unsigned C> class MsqVector;
class PatchData;
class MsqError;
class Vector3D;
class Matrix3D;


/**\brief Calculate R and Z such that \f$W\prime = Z^{-1} W\f$ and 
 *        \f$A\prime = (RZ)^{-1} A\f$
 *
 * Calculate the matrices required to transform the active and target
 * matrices from the 3x2 surface domain to a 2x2 2D domain.
 *\param A    Input: Element Jacobian matrix.
 *\param W_32 Input: Target Jacobian matrix.
 *\param W_22 Output: 2D Target matrix.
 *\param RZ   Output: Product of R and Z needed to calculate the 2D 
 *            element matrix.
 */
void surface_to_2d( const MsqMatrix<3,2>& A, 
                    const MsqMatrix<3,2>& W_32,
                    MsqMatrix<2,2>& W_22,
                    MsqMatrix<3,2>& RZ );
/*
void surface_to_2d( const MsqMatrix<3,2>& A_in,
                    const MsqMatrix<3,2>& W_in,
                    MsqMatrix<2,2>& A_out,
                    MsqMatrix<2,2>& W_out );
*/
void get_sample_pt_evaluations( PatchData& pd,
                                std::vector<size_t>& handles,
                                bool free,
                                MsqError& err );
                    
void get_elem_sample_points( PatchData& pd,
                             size_t elem,
                             std::vector<size_t>& handles,
                             MsqError& err );


/**\brief Calculate gradient from derivatives of mapping function terms
 *        and derivatives of target metric. */
template <int DIM> inline
void gradient( size_t num_free_verts,
               const MsqVector<DIM>* dNdxi,
               const MsqMatrix<3,DIM>& dmdA,
               std::vector<Vector3D>& grad )
{
  grad.clear();
  grad.resize( num_free_verts, Vector3D(0,0,0) );
  for (size_t i = 0; i < num_free_verts; ++i)
    grad[i] = Vector3D( (dmdA * dNdxi[i]).data() );
}

/**\brief Calculate Hessian from derivatives of mapping function terms
 *        and derivatives of target metric. */
template <int DIM, typename MAT> inline
void hessian( size_t num_free_verts,
              const MsqVector<DIM>* dNdxi,
              const MsqMatrix<DIM,DIM>* d2mdA2,
              MAT* hess )
{
  MsqMatrix<1,DIM> tmp[DIM][DIM];
  size_t h = 0; // index of current Hessian block

  for (size_t i = 0; i < num_free_verts; ++i) {
  
      // Populate TMP with vector-matrix procucts common
      // to terms of this Hessian row.
    const MsqMatrix<1,DIM>& gi = transpose(dNdxi[i]);
    switch (DIM) {
      case 3:
        tmp[0][2] = gi * d2mdA2[2];
        tmp[1][2] = gi * d2mdA2[4];
        tmp[2][0] = gi * transpose(d2mdA2[2]);
        tmp[2][1] = gi * transpose(d2mdA2[4]);
        tmp[2][2] = gi * d2mdA2[5];
     case 2:
        tmp[0][1] = gi * d2mdA2[1];
        tmp[1][0] = gi * transpose(d2mdA2[1]);
        tmp[1][1] = gi * d2mdA2[DIM];
      case 1:
        tmp[0][0] = gi * d2mdA2[0];
      case 0: 
        break;
      default: assert(false);
    }

      // Calculate Hessian diagonal block
    MAT& H = hess[h++];
    switch (DIM) {
      case 3:
        H(0,2) = H(2,0) = tmp[0][2] * transpose(gi);
        H(1,2) = H(2,1) = tmp[1][2] * transpose(gi);
        H(2,2) =          tmp[2][2] * transpose(gi);
      case 2:
        H(0,1) = H(1,0) = tmp[0][1] * transpose(gi);
        H(1,1) =          tmp[1][1] * transpose(gi);
      case 1:
        H(0,0) =          tmp[0][0] * transpose(gi);
      case 0: 
        break;
      default: assert(false);
    }
    
      // Calculate remainder of Hessian row
    for (size_t j = i+1; j < num_free_verts; ++j) {
      MAT& H = hess[h++];
      const MsqMatrix<DIM,1>& gj = dNdxi[j];
      switch (DIM) {
        case 3:
          H(0,2) = tmp[0][2] * gj;
          H(1,2) = tmp[1][2] * gj;
          H(2,0) = tmp[2][0] * gj;
          H(2,1) = tmp[2][1] * gj;
          H(2,2) = tmp[2][2] * gj;
        case 2:
          H(0,1) = tmp[0][1] * gj;
          H(1,0) = tmp[1][0] * gj;
          H(1,1) = tmp[1][1] * gj;
        case 1:
          H(0,0) = tmp[0][0] * gj;
        case 0: 
          break;
        default: assert(false);
      }
    }
  }
}

/**\brief Calculate Hessian from derivatives of mapping function terms
 *        and derivatives of target metric. */
template <int DIM> inline
void hessian_diagonal( size_t num_free_verts,
              const MsqVector<DIM>* dNdxi,
              const MsqMatrix<DIM,DIM>* d2mdA2,
              SymMatrix3D* diagonal )
{
  for (size_t i = 0; i < num_free_verts; ++i) {
    SymMatrix3D& H = diagonal[i];
    for (unsigned j = 0; j < ((DIM)*(DIM+1)/2); ++j)
      H[j] = transpose(dNdxi[i]) * d2mdA2[j] * dNdxi[i];
  }
}


#ifdef PRINT_INFO
template <int R, int C> inline 
void write_vect( char n, const MsqMatrix<R,C>& M )
{
  std::cout << "  " << n << ':';
  for (int c = 0; c < C; ++c) {
    std::cout << '[';
    for (int r = 0; r < R; ++r)
      std::cout << M(r,c) << ' ';
    std::cout << ']';
  }
  std::cout << std::endl;
}

template <int D> inline
void print_info( size_t elem, Sample sample,
                 const MsqMatrix<3,D>& A,
                 const MsqMatrix<3,D>& W,
                 const MsqMatrix<D,D>& T )
{
  std::cout << "Elem " << elem << " Dim " << sample.dimension << " Num " << sample.number << " :" << std::endl;
  write_vect<3,D>( 'A', A );
  write_vect<3,D>( 'W', W );
  write_vect<D,D>( 'T', T );
}
#endif

                    
} // namespace Mesquite

#endif
