/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file AddQualityMetric.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "AddQualityMetric.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

AddQualityMetric::AddQualityMetric( QualityMetric* qm1, 
                                    QualityMetric* qm2,
                                    MsqError& err )
  : metric1(*qm1), 
    metric2(*qm2)
{
  if (qm1->get_metric_type() != qm2->get_metric_type() ||
      qm1->get_negate_flag() != qm2->get_negate_flag())
  {
    MSQ_SETERR(err)("Incompatible metrics", MsqError::INVALID_ARG);
  }
}

AddQualityMetric::~AddQualityMetric() {}

QualityMetric::MetricType AddQualityMetric::get_metric_type() const
{
  return metric1.get_metric_type();
}

std::string AddQualityMetric::get_name() const
{
  return std::string("Sum(") + metric1.get_name() + ", " + metric2.get_name() + ")";
}

int AddQualityMetric::get_negate_flag() const
{
  return metric1.get_negate_flag();
}

void AddQualityMetric::get_evaluations( PatchData& pd,
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

bool AddQualityMetric::evaluate( PatchData& pd,
                                 size_t handle,
                                 double& value,
                                 MsqError& err )
{
  double val1, val2;
  bool rval1, rval2;
  rval1 = metric1.evaluate( pd, handle, val1, err ); MSQ_ERRZERO(err);
  rval2 = metric2.evaluate( pd, handle, val2, err ); MSQ_ERRZERO(err);
  value = val1 + val2;
  return rval1 && rval2;
}

bool AddQualityMetric::evaluate_with_indices( PatchData& pd,
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
  
  value = val1 + val2;
  return rval1 && rval2;
}

bool AddQualityMetric::evaluate_with_gradient( PatchData& pd,
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
    gradient[k] += grad1[j];
  }
  for (j = 0; j < indices2.size(); ++j)
  {
    i = std::lower_bound( indices.begin(), indices.end(), indices2[j] );
    size_t k = i - indices.begin();
    gradient[k] += grad2[j];
  }
  
  value = val1 + val2;
  return rval1 && rval2;
}

bool AddQualityMetric::evaluate_with_Hessian( PatchData& pd,
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
    indices1[j] = i - indices.begin();
    gradient[indices1[j]] += grad1[j];
  }
  for (j = 0; j < indices2.size(); ++j)
  {
    i = std::lower_bound( indices.begin(), indices.end(), indices2[j] );
    indices2[j] = i - indices.begin();
    gradient[indices2[j]] += grad2[j];
  }

  const size_t N = indices.size();
  Hessian.clear();
  Hessian.resize( N * (N+1) / 2, Matrix3D(0.0) );
  
  n = indices1.size();
  h = 0; 
  for (r = 0; r < n; ++r) 
  {
    const size_t nr = indices1[r];
    for (c = r; c < n; ++c)
    {
      const size_t nc = indices1[c];
      if (nr <= nc)
        Hessian[N*nr - nr*(nr+1)/2 + nc] += Hess1[h++];
      else
        Hessian[N*nc - nc*(nc+1)/2 + nr].plus_transpose_equal( Hess1[h++] );
    }
  }
  
  n = indices2.size();
  h = 0; 
  for (r = 0; r < n; ++r) 
  {
    const size_t nr = indices2[r];
    for (c = r; c < n; ++c)
    {
      const size_t nc = indices2[c];
      if (nr <= nc)
        Hessian[N*nr - nr*(nr+1)/2 + nc] += Hess2[h++];
      else
        Hessian[N*nc - nc*(nc+1)/2 + nr].plus_transpose_equal( Hess2[h++] );
    }
  }
  
  value = val1 + val2;
  return rval1 && rval2;
}


} // namespace Mesquite
