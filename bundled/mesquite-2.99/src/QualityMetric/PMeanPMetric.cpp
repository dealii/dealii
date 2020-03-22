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


/** \file PMeanPMetric.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "PMeanPMetric.hpp"
#include "MsqError.hpp"
#include "QualityMetric.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"
#include "SymMatrix3D.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {

bool PMeanPMetric::average( PatchData& pd,
                            QualityMetric* metric,
                            const std::vector<size_t>& qm_handles,
                            double& value, 
                            MsqError& err )
{
  bool rval = true;
  value = 0;
  for (std::vector<size_t>::const_iterator i= qm_handles.begin();
       i != qm_handles.end(); ++i) {
    double mval;
    if (!metric->evaluate( pd, *i, mval, err )) 
      rval = false;
    MSQ_ERRZERO(err);
    value += P.raise(mval);
  }
  value /= qm_handles.size();
  return rval;
}


bool PMeanPMetric::average_with_indices( PatchData& pd, 
                                         QualityMetric* metric,
                                         const std::vector<size_t>& qm_handles,
                                         double& value, 
                                         std::vector<size_t>& indices,
                                         MsqError& err )
{
  indices.clear();
  
  bool rval = true;
  value = 0;
  for (std::vector<size_t>::const_iterator i = qm_handles.begin();
       i != qm_handles.end(); ++i) {
    double mval;
    mIndices.clear();
    if (!metric->evaluate_with_indices( pd, *i, mval, mIndices, err )) 
      rval = false;
    MSQ_ERRZERO(err);
    value += P.raise(mval);
    
    std::copy( mIndices.begin(), mIndices.end(), std::back_inserter(indices) );
  }
  std::sort( indices.begin(), indices.end() );
  indices.erase( std::unique( indices.begin(), indices.end() ), indices.end() );
  
  value /= qm_handles.size();
  return rval;
}

bool PMeanPMetric::average_with_gradient( PatchData& pd, 
                                          QualityMetric* metric,
                                          const std::vector<size_t>& qm_handles,
                                          double& value, 
                                          std::vector<size_t>& indices,
                                          std::vector<Vector3D>& gradient,
                                          MsqError& err )
{
  indices.clear();
  gradient.clear();

  std::vector<Vector3D>::iterator  g;
  std::vector<size_t>::iterator j, k;
  std::vector<size_t>::const_iterator i;
  
  bool rval = true;
  value = 0;
  for (i = qm_handles.begin(); i != qm_handles.end(); ++i) {
    
    double mval;
    mIndices.clear();
    mGrad.clear();
    if (!metric->evaluate_with_gradient( pd, *i, mval, mIndices, mGrad, err )) 
      rval = false;
    MSQ_ERRZERO(err);
    value += P.raise(mval);
    
    double p1val = P.value() * P1.raise( mval );
    for (j = mIndices.begin(), g = mGrad.begin(); j != mIndices.end(); ++j, ++g) {
      
      *g *= p1val;
      k = std::lower_bound( indices.begin(), indices.end(), *j );
      if (k == indices.end() || *k != *j) {
        k = indices.insert( k, *j );
        size_t idx = k - indices.begin();
        gradient.insert( gradient.begin() + idx, *g );
      }
      else {
        size_t idx = k - indices.begin();
        gradient[idx] += *g;
      }
    }
  }
  
  double inv_n = 1.0 / qm_handles.size();
  value *= inv_n;
  for (g = gradient.begin(); g != gradient.end(); ++g)
    *g *= inv_n;
  
  return rval;
}  

   
bool PMeanPMetric::average_with_Hessian_diagonal( PatchData& pd,
                                        QualityMetric* metric,
                                        const std::vector<size_t>& qm_handles,
                                        double& value,
                                        std::vector<size_t>& indices,
                                        std::vector<Vector3D>& gradient,
                                        std::vector<SymMatrix3D>& diagonal,
                                        MsqError& err )
{
    // clear temporary storage
  mIndices.clear();
  mGrad.clear();
  mDiag.clear();
  mOffsets.clear();
  mValues.clear();
  
    // Evaluate metric for all sample points,
    // accumulating indices, gradients, and Hessians
  bool rval = true;
  std::vector<size_t>::const_iterator q;
  for (q = qm_handles.begin(); q != qm_handles.end(); ++q) {
    double mval;
    indices.clear();
    gradient.clear();
    diagonal.clear();
    if (!metric->evaluate_with_Hessian_diagonal( pd, *q, mval, indices, gradient, diagonal, err )) 
      rval = false;
    MSQ_ERRZERO(err);
    
    mValues.push_back( mval );
    mOffsets.push_back( mIndices.size() );
    std::copy( indices.begin(), indices.end(), std::back_inserter(mIndices) );
    std::copy(gradient.begin(),gradient.end(), std::back_inserter(mGrad) );
    std::copy(diagonal.begin(),diagonal.end(), std::back_inserter(mDiag) );
  }
  mOffsets.push_back( mIndices.size() );
  
    // Combine lists of free vertex indices, and update indices
    // in per-evaluation lists to point into the combined gradient
    // and Hessian lists.
  indices = mIndices;
  std::sort( indices.begin(), indices.end() );
  indices.erase( std::unique( indices.begin(), indices.end() ), indices.end() );
  std::vector<size_t>::iterator i, j;
  for (i = mIndices.begin(); i != mIndices.end(); ++i) {
    j = std::lower_bound( indices.begin(), indices.end(), *i );
    assert( *j == *i );
    *i = j - indices.begin();
  }
  
    // Allocate space and zero output gradient and Hessian lists
  const size_t n = indices.size();
  const size_t m = mValues.size();
  const double inv_n = 1.0 / m;
  double v;
  gradient.clear();
  gradient.resize( n, Vector3D(0,0,0) );
  diagonal.clear();
  diagonal.resize( n, SymMatrix3D(0.0) );
  
    // Average values, gradients, Hessians
  value = 0.0;
  for (size_t k = 0; k < m; ++k) { // for each metric evaluate
    if (P.value() == 1.0) {
      v = P.raise( mValues[k] );
        // for each gradient (or Hessian) for the local metric evaluation
      for (size_t r = mOffsets[k]; r < mOffsets[k+1]; ++r) {
        const size_t nr = mIndices[r];

        mDiag[r] *= inv_n;
        diagonal[nr] += mDiag[r];
        mGrad[r] *= inv_n;
        gradient[nr] += mGrad[r];
      }
    }
    else {
      const double r2 = P2.raise(mValues[k]);
      const double r1 = r2 * mValues[k];
      v = r1 * mValues[k];
      const double h = r2 * P.value() * (P.value() - 1) * inv_n;
      const double g = r1 * P.value() * inv_n; 
        // for each gradient (or Hessian) for the local metric evaluation
      for (size_t r = mOffsets[k]; r < mOffsets[k+1]; ++r) {
        const size_t nr = mIndices[r];

        mDiag[r] *= g;
        diagonal[nr] += mDiag[r];
        diagonal[nr] += h * outer( mGrad[r] );
        mGrad[r] *= g;
        gradient[nr] += mGrad[r];
      }
    }
    value += v;
  }
  
  value *= inv_n;
  
  return rval;
}

bool PMeanPMetric::average_with_Hessian( PatchData& pd, 
                                         QualityMetric* metric,
                                         const std::vector<size_t>& qm_handles,
                                         double& value, 
                                         std::vector<size_t>& indices,
                                         std::vector<Vector3D>& gradient,
                                         std::vector<Matrix3D>& Hessian,
                                         MsqError& err )
{
    // clear temporary storage
  mIndices.clear();
  mGrad.clear();
  mHess.clear();
  mOffsets.clear();
  mValues.clear();
  
    // Evaluate metric for all sample points,
    // accumulating indices, gradients, and Hessians
  bool rval = true;
  std::vector<size_t>::const_iterator q;
  for (q = qm_handles.begin(); q != qm_handles.end(); ++q) {
    double mval;
    indices.clear();
    gradient.clear();
    Hessian.clear();
    if (!metric->evaluate_with_Hessian( pd, *q, mval, indices, gradient, Hessian, err )) 
      rval = false;
    MSQ_ERRZERO(err);
    
    mValues.push_back( mval );
    mOffsets.push_back( mIndices.size() );
    std::copy( indices.begin(), indices.end(), std::back_inserter(mIndices) );
    std::copy(gradient.begin(),gradient.end(), std::back_inserter(mGrad) );
    std::copy( Hessian.begin(), Hessian.end(), std::back_inserter(mHess) );
  }
  mOffsets.push_back( mIndices.size() );
  
    // Combine lists of free vertex indices, and update indices
    // in per-evaluation lists to point into the combined gradient
    // and Hessian lists.
  indices = mIndices;
  std::sort( indices.begin(), indices.end() );
  indices.erase( std::unique( indices.begin(), indices.end() ), indices.end() );
  std::vector<size_t>::iterator i, j;
  for (i = mIndices.begin(); i != mIndices.end(); ++i) {
    j = std::lower_bound( indices.begin(), indices.end(), *i );
    assert( *j == *i );
    *i = j - indices.begin();
  }
  
    // Allocate space and zero output gradient and Hessian lists
  const size_t N = indices.size();
  const size_t m = mValues.size();
  const double inv_n = 1.0 / m;
  double v;
  gradient.clear();
  gradient.resize( N, Vector3D(0,0,0) );
  Hessian.clear();
  Hessian.resize( N*(N+1)/2, Matrix3D(0.0) );
  
    // Average values, gradients, Hessians
  Matrix3D outer;
  value = 0.0;
  std::vector<Matrix3D>::iterator met_hess_iter = mHess.begin();
  for (size_t k = 0; k < m; ++k) { // for each metric evaluate
    if (P.value() == 1.0) {
      v = P.raise( mValues[k] );
        // for each gradient (or Hessian row) for the local metric evaluation
      for (size_t r = mOffsets[k]; r < mOffsets[k+1]; ++r) {
        const size_t nr = mIndices[r];
          // for each column of the local metric Hessian
        for (size_t c = r; c < mOffsets[k+1]; ++c) {
          const size_t nc = mIndices[c];
          *met_hess_iter *= inv_n;
          if (nr <= nc)
            Hessian[N*nr - nr*(nr+1)/2 + nc] += *met_hess_iter;
          else
            Hessian[N*nc - nc*(nc+1)/2 + nr].plus_transpose_equal( *met_hess_iter );
          ++met_hess_iter;
        }
        mGrad[r] *= inv_n;
        gradient[nr] += mGrad[r];
      }
    }
    else {
      const double r2 = P2.raise(mValues[k]);
      const double r1 = r2 * mValues[k];
      v = r1 * mValues[k];
      const double h = r2 * P.value() * (P.value() - 1) * inv_n;
      const double g = r1 * P.value() * inv_n; 
        // for each gradient (or Hessian row) for the local metric evaluation
      for (size_t r = mOffsets[k]; r < mOffsets[k+1]; ++r) {
        const size_t nr = mIndices[r];
          // for each column of the local metric Hessian
        for (size_t c = r; c < mOffsets[k+1]; ++c) {
          const size_t nc = mIndices[c];
          *met_hess_iter *= g;
          outer.outer_product( mGrad[r], mGrad[c] );
          outer *= h;
          outer += *met_hess_iter;
          if (nr <= nc)
            Hessian[N*nr - nr*(nr+1)/2 + nc] += outer;
          else
            Hessian[N*nc - nc*(nc+1)/2 + nr].plus_transpose_equal( outer );
          ++met_hess_iter;
        }
        mGrad[r] *= g;
        gradient[nr] += mGrad[r];
      }
    }
    value += v;
  }
  
  value *= inv_n;
  
  return rval;
}  

} // namespace Mesquite
