/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   MultiplyQualityMetric.cpp
  \brief 

  \author Michael Brewer
  \date   2002-05-09
*/

#include "MultiplyQualityMetric.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"
#include "MsqError.hpp"

using namespace Mesquite;

MultiplyQualityMetric::MultiplyQualityMetric(QualityMetric* qm1, QualityMetric* qm2, MsqError &err)
  : metric1(*qm1), metric2(*qm2)
{
  if (qm1->get_metric_type() != qm2->get_metric_type() ||
      qm1->get_negate_flag() != qm2->get_negate_flag())
  {
    MSQ_SETERR(err)("Incompatible metrics", MsqError::INVALID_ARG);
  }
}

MultiplyQualityMetric::~MultiplyQualityMetric() {}

QualityMetric::MetricType MultiplyQualityMetric::get_metric_type() const
{
  return metric1.get_metric_type();
}

std::string MultiplyQualityMetric::get_name() const
{
  return metric1.get_name() + "*" + metric2.get_name();
}

int MultiplyQualityMetric::get_negate_flag() const
{
  return metric1.get_negate_flag();
}

void MultiplyQualityMetric::get_evaluations( PatchData& pd,
                                        std::vector<size_t>& handles,
                                        bool free_only,
                                        MsqError& err )
{
  metric1.get_evaluations( pd, handles, free_only, err );
  MSQ_ERRRTN(err);
  metric2.get_evaluations( pd, mHandles, free_only, err );
  MSQ_ERRRTN(err);
  if (handles != mHandles) {
    MSQ_SETERR(err)("Incompatible metrics", MsqError::INVALID_STATE);
  }
}

bool MultiplyQualityMetric::evaluate( PatchData& pd,
                                 size_t handle,
                                 double& value,
                                 MsqError& err )
{
  double val1, val2;
  bool rval1, rval2;
  rval1 = metric1.evaluate( pd, handle, val1, err ); MSQ_ERRZERO(err);
  rval2 = metric2.evaluate( pd, handle, val2, err ); MSQ_ERRZERO(err);
  value = val1 * val2;
  return rval1 && rval2;
}

bool MultiplyQualityMetric::evaluate_with_indices( PatchData& pd,
                                              size_t handle,
                                              double& value,
                                              std::vector<size_t>& indices,
                                              MsqError& err )
{
  double val1, val2;
  bool rval1, rval2;
  rval1 = metric1.evaluate_with_indices( pd, handle, val1, indices1, err ); MSQ_ERRZERO(err);
  rval2 = metric2.evaluate_with_indices( pd, handle, val2, indices2, err ); MSQ_ERRZERO(err);
  
  indices.clear();
  std::sort( indices1.begin(), indices1.end() );
  std::sort( indices2.begin(), indices2.end() );
  std::set_union( indices1.begin(), indices1.end(),
                  indices2.begin(), indices2.end(),
                  std::back_inserter( indices ) );
  
  value = val1 * val2;
  return rval1 && rval2;
}

bool MultiplyQualityMetric::evaluate_with_gradient( PatchData& pd,
                                               size_t handle,
                                               double& value,
                                               std::vector<size_t>& indices,
                                               std::vector<Vector3D>& gradient,
                                               MsqError& err )
{
  std::vector<size_t>::iterator i;
  size_t j;
  double val1, val2;
  bool rval1, rval2;
  rval1 = metric1.evaluate_with_gradient( pd, handle, val1, indices1, grad1, err ); MSQ_ERRZERO(err);
  rval2 = metric2.evaluate_with_gradient( pd, handle, val2, indices2, grad2, err ); MSQ_ERRZERO(err);
  
  indices.resize( indices1.size() + indices2.size() );
  i = std::copy( indices1.begin(), indices1.end(), indices.begin() );
  std::copy( indices2.begin(), indices2.end(), i );
  std::sort( indices.begin(), indices.end() );
  indices.erase( std::unique( indices.begin(), indices.end() ), indices.end() );
  
  gradient.clear();
  gradient.resize( indices.size(), Vector3D(0.0) );
  for (j = 0; j < indices1.size(); ++j)
  {
    i = std::lower_bound( indices.begin(), indices.end(), indices1[j] );
    size_t k = i - indices.begin();
    gradient[k] += val2 * grad1[j];
  }
  for (j = 0; j < indices2.size(); ++j)
  {
    i = std::lower_bound( indices.begin(), indices.end(), indices2[j] );
    size_t k = i - indices.begin();
    gradient[k] += val1 * grad2[j];
  }
  
  value = val1 * val2;
  return rval1 && rval2;
}

bool MultiplyQualityMetric::evaluate_with_Hessian( PatchData& pd,
                                              size_t handle,
                                              double& value,
                                              std::vector<size_t>& indices,
                                              std::vector<Vector3D>& gradient,
                                              std::vector<Matrix3D>& Hessian,
                                              MsqError& err )
{
  std::vector<size_t>::iterator i;
  size_t j, r, c, n, h;
  double val1, val2;
  bool rval1, rval2;
  rval1 = metric1.evaluate_with_Hessian( pd, handle, val1, indices1, grad1, Hess1, err ); MSQ_ERRZERO(err);
  rval2 = metric2.evaluate_with_Hessian( pd, handle, val2, indices2, grad2, Hess2, err ); MSQ_ERRZERO(err);
    // merge index lists
  indices.resize( indices1.size() + indices2.size() );
  i = std::copy( indices1.begin(), indices1.end(), indices.begin() );
  std::copy( indices2.begin(), indices2.end(), i );
  std::sort( indices.begin(), indices.end() );
  indices.erase( std::unique( indices.begin(), indices.end() ), indices.end() );
    // calculate grads and convert index lists to indices into output list
  gradient.clear();
  gradient.resize( indices.size(), Vector3D(0.0) );
  for (j = 0; j < indices1.size(); ++j)
  {
    i = std::lower_bound( indices.begin(), indices.end(), indices1[j] );
    indices1[j] = i - indices.begin();
    gradient[indices1[j]] += val2 * grad1[j];
  }
  for (j = 0; j < indices2.size(); ++j)
  {
    i = std::lower_bound( indices.begin(), indices.end(), indices2[j] );
    indices2[j] = i - indices.begin();
    gradient[indices2[j]] += val1 * grad2[j];
  }
    // allocate space for hessians, and zero it
  const size_t N = indices.size();
  Hessian.clear();
  Hessian.resize( N * (N+1) / 2, Matrix3D(0.0) );
    // add hessian terms from first metric
  n = indices1.size();
  h = 0; 
  for (r = 0; r < n; ++r) 
  {
    const size_t nr = indices1[r];
    for (c = r; c < n; ++c)
    {
      const size_t nc = indices1[c];
      Hess1[h] *= val2;
      if (nr <= nc)
        Hessian[N*nr - nr*(nr+1)/2 + nc] += Hess1[h];
      else
        Hessian[N*nc - nc*(nc+1)/2 + nr].plus_transpose_equal( Hess1[h] );
      ++h;
    }
  }
    // add hessian terms from second metric
  n = indices2.size();
  h = 0; 
  for (r = 0; r < n; ++r) 
  {
    const size_t nr = indices2[r];
    for (c = r; c < n; ++c)
    {
      const size_t nc = indices2[c];
      Hess2[h] *= val1;
      if (nr <= nc)
        Hessian[N*nr - nr*(nr+1)/2 + nc] += Hess2[h];
      else
        Hessian[N*nc - nc*(nc+1)/2 + nr].plus_transpose_equal( Hess2[h] );
      ++h;
    }
  }
    // add gradient outer products
  n = indices1.size();
  size_t m = indices2.size();
  Matrix3D outer;
  for (r = 0; r < n; ++r)
  {
    const size_t nr = indices1[r];
    for (c = 0; c < m; ++c)
    {
      const size_t nc = indices2[c];
      outer.outer_product( grad1[r], grad2[c] );
      if (nr == nc) 
        Hessian[N*nr - nr*(nr+1)/2 + nc] += outer.plus_transpose_equal(outer);
      else if (nr < nc)
        Hessian[N*nr - nr*(nr+1)/2 + nc] += outer;
      else
        Hessian[N*nc - nc*(nc+1)/2 + nr].plus_transpose_equal(outer);
    }
  }
  
  value = val1 * val2;
  return rval1 && rval2;
}

