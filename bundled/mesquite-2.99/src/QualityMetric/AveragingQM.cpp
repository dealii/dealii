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
    
    (2006) kraftche@cae.wisc.edu  
   
  ***************************************************************** */
/*!
  \file   AveragingQM.cpp
  \brief  

  \author Michael Brewer
  \author Thomas Leurent
  \author Jason Kraftcheck
  \date   2002-05-14
*/

#include "AveragingQM.hpp"
#include "MsqVertex.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqDebug.hpp"
#include "MsqTimer.hpp"
#include "PatchData.hpp"


namespace MESQUITE_NS {
        
double AveragingQM::average_corner_gradients( EntityTopology type,
                                  uint32_t fixed_vertices,
                                  unsigned num_corner,
                                  double corner_values[],
                                  const Vector3D corner_grads[],
                                  Vector3D vertex_grads[],
                                  MsqError& err )
{
  const unsigned num_vertex = TopologyInfo::corners( type );
  const unsigned dim = TopologyInfo::dimension(type);
  const unsigned per_vertex = dim+1;
  
  unsigned i, j, num_adj;
  const unsigned *adj_idx, *rev_idx;
  
    // NOTE: This function changes the corner_values array such that
    //       it contains the gradient coefficients.
  double avg = average_metric_and_weights( corner_values, num_corner, err );
  MSQ_ERRZERO(err);

  for (i = 0; i < num_vertex; ++i)
  {
    if (fixed_vertices & (1<<i))  // skip fixed vertices
      continue;
    
    adj_idx = TopologyInfo::adjacent_vertices( type, i, num_adj );
    rev_idx = TopologyInfo::reverse_vertex_adjacency_offsets( type, i, num_adj );
    if (i < num_corner) // not all vertices are corners (e.g. pyramid)
      vertex_grads[i] = corner_values[i] * corner_grads[per_vertex*i];
    else
      vertex_grads[i] = 0;
    for (j = 0; j < num_adj; ++j)
    {
      const unsigned v = adj_idx[j], c = rev_idx[j]+1;
      if (v >= num_corner) // if less corners than vertices (e.g. pyramid apex)
        continue;
      vertex_grads[i] += corner_values[v] * corner_grads[per_vertex*v+c];
    }
  }

  return avg;
}

/**\brief Iterate over only diagonal blocks of element corner Hessian data 
 *
 * Given concatenation of corner Hessian data for an element, iterate
 * over only the diagonal terms for each corner.  This class allows 
 * common code to be used to generate Hessian diagonal blocks from either
 * the diagonal blocks for each corner or the full Hessian data for each 
 * corner, where this class is used for the latter.
 */
class CornerHessDiagIterator 
{
private:
  const Matrix3D* cornerHess;    //!< Current location in concatenated Hessian data.
  const EntityTopology elemType; //!< Element topology for Hessian data
  unsigned mCorner;              //!< The element corner for which cornerHess 
                                 //!< is pointing into the corresponding Hessian data.
  unsigned mStep;                //!< Amount to step to reach next diagonal block.
public:
  CornerHessDiagIterator( const Matrix3D* corner_hessians, 
                          EntityTopology elem_type )
    : cornerHess(corner_hessians), 
      elemType(elem_type), 
      mCorner(0)
    {
      TopologyInfo::adjacent_vertices( elemType, mCorner, mStep );
      ++mStep;
    }
  
  SymMatrix3D operator*() const 
    { return cornerHess->upper(); }
  
  CornerHessDiagIterator& operator++()
    { 
      cornerHess += mStep;
      if (!--mStep) {
        TopologyInfo::adjacent_vertices( elemType, ++mCorner, mStep );
        ++mStep;
      }
      return *this;
    }
  
  CornerHessDiagIterator operator++(int)
    { CornerHessDiagIterator copy(*this); operator++(); return copy; }
};

template <typename HessIter>
static inline 
double sum_corner_diagonals( EntityTopology type,
                             unsigned num_corner,
                             const double corner_values[],
                             const Vector3D corner_grads[],
                             HessIter corner_diag_blocks,
                             Vector3D vertex_grads[],
                             SymMatrix3D vertex_hessians[] )
{
  unsigned i, n, r, R, idx[4];
  const unsigned* adj_list;
  double avg = 0.0;
  
    // calculate mean
  for (i = 0; i < num_corner; ++i)
    avg += corner_values[i];

  const Vector3D* grad = corner_grads;
  HessIter hess = corner_diag_blocks;
  for (i = 0; i < num_corner; ++i)
  {
    adj_list = TopologyInfo::adjacent_vertices( type, i, n );
    idx[0] = i;
    idx[1] = adj_list[0];
    idx[2] = adj_list[1];
    idx[3] = adj_list[2%n]; // %n so don't read off end if 2D

    for (r = 0; r <= n; ++r) 
    {
      R = idx[r];
      vertex_grads[R] += *grad;
      vertex_hessians[R] += *hess;
      ++grad;
      ++hess;
    }
  }
  return avg;
}

template <typename HessIter>
static inline 
double sum_sqr_corner_diagonals( EntityTopology type,
                                 unsigned num_corner,
                                 const double corner_values[],
                                 const Vector3D corner_grads[],
                                 HessIter corner_diag_blocks,
                                 Vector3D vertex_grads[],
                                 SymMatrix3D vertex_hessians[] )
{
  unsigned i, n, r, R, idx[4];
  const unsigned* adj_list;
  double v, avg = 0.0;
  
    // calculate mean
  for (i = 0; i < num_corner; ++i)
    avg += corner_values[i]*corner_values[i];

  const Vector3D* grad = corner_grads;
  HessIter hess = corner_diag_blocks;
  for (i = 0; i < num_corner; ++i)
  {
    adj_list = TopologyInfo::adjacent_vertices( type, i, n );
    idx[0] = i;
    idx[1] = adj_list[0];
    idx[2] = adj_list[1];
    idx[3] = adj_list[2%n]; // %n so don't read off end if 2D
    ++n;

    v = 2.0*corner_values[i];
    for (r = 0; r < n; ++r) 
    {
      R = idx[r];
      vertex_grads[R] += v * *grad;
      vertex_hessians[R] += 2.0 * outer( *grad );
      vertex_hessians[R] += v * *hess;
      ++grad;
      ++hess;
    }
  }
  return avg;
}

template <typename HessIter>
static inline 
double pmean_corner_diagonals( EntityTopology type,
                               unsigned num_corner,
                               const double corner_values[],
                               const Vector3D corner_grads[],
                               HessIter corner_diag_blocks,
                               Vector3D vertex_grads[],
                               SymMatrix3D vertex_hessians[],
                               double p )
{
  const unsigned N = TopologyInfo::corners( type );
  unsigned i, n, r, R, idx[4];
  const unsigned* adj_list;
  double m = 0.0, nm;
  double gf[8], hf[8];
  double inv = 1.0/num_corner;
  assert(num_corner <= 8);
  
    // calculate mean
  for (i = 0; i < num_corner; ++i)
  {
    nm = pow(corner_values[i], p);
    m += nm;
    
    gf[i] = inv * p * nm / corner_values[i];
    hf[i] = (p-1) * gf[i] / corner_values[i];
  }
  nm = inv * m;
  
  const Vector3D* grad = corner_grads;
  HessIter hess = corner_diag_blocks;
  for (i = 0; i < num_corner; ++i)
  {
    adj_list = TopologyInfo::adjacent_vertices( type, i, n );
    idx[0] = i;
    idx[1] = adj_list[0];
    idx[2] = adj_list[1];
    idx[3] = adj_list[2%n]; // %n so don't read off end if 2D
    ++n;

    for (r = 0; r < n; ++r) 
    {
      R = idx[r];
      vertex_grads[R] += gf[i] * *grad;
      vertex_hessians[R] += hf[i] * outer( *grad );
      vertex_hessians[R] += gf[i] * *hess;
      ++grad;
      ++hess;
    }
  }
  
  m = pow( nm, 1.0/p );
  gf[0] = m / (p * nm );
  hf[0] = (1.0/p - 1) * gf[0] / nm;
  for (r = 0; r < N; ++r) 
  {
    vertex_hessians[r] *= gf[0];
    vertex_hessians[r] += hf[0] * outer( vertex_grads[r] );
    vertex_grads[r] *= gf[0];
  }
  
  return m;
}
                                          

template <typename HessIter>
static inline
double average_corner_diagonals( EntityTopology type,
                                 QualityMetric::AveragingMethod method,
                                 unsigned num_corner,
                                 const double corner_values[],
                                 const Vector3D corner_grads[],
                                 HessIter corner_diag_blocks,
                                 Vector3D vertex_grads[],
                                 SymMatrix3D vertex_hessians[],
                                 MsqError& err )
{
  unsigned i;
  double avg, inv;
  
    // Zero gradients and Hessians
  const unsigned num_vertex = TopologyInfo::corners( type );
  for (i = 0; i < num_vertex; ++i) {
    vertex_grads[i].set(0.0);
    vertex_hessians[i] = SymMatrix3D(0.0);
  }
  
  switch (method)
  {
  case QualityMetric::SUM:
    avg = sum_corner_diagonals( type, num_corner, corner_values, 
                               corner_grads, corner_diag_blocks,
                               vertex_grads, vertex_hessians );
    break;
  
  case QualityMetric::LINEAR:
    avg = sum_corner_diagonals( type, num_corner, corner_values, 
                               corner_grads, corner_diag_blocks,
                               vertex_grads, vertex_hessians );
    inv = 1.0/num_corner;
    avg *= inv;
    for (i = 0; i < num_vertex; ++i) {
      vertex_grads[i] *= inv;
      vertex_hessians[i] *= inv;
    }
    break;

  case QualityMetric::SUM_SQUARED:
    avg = sum_sqr_corner_diagonals( type, num_corner, corner_values, 
                                   corner_grads, corner_diag_blocks,
                                   vertex_grads, vertex_hessians );
    break;

  case QualityMetric::RMS:
    avg = pmean_corner_diagonals( type, num_corner, corner_values, 
                                 corner_grads, corner_diag_blocks,
                                 vertex_grads, vertex_hessians,
                                 2.0 );
    break;

  case QualityMetric::HARMONIC:
    avg = pmean_corner_diagonals( type, num_corner, corner_values, 
                                 corner_grads, corner_diag_blocks,
                                 vertex_grads, vertex_hessians,
                                 -1.0 );
    break;

  case QualityMetric::HMS:
    avg = pmean_corner_diagonals( type, num_corner, corner_values, 
                                 corner_grads, corner_diag_blocks,
                                 vertex_grads, vertex_hessians,
                                 -2.0 );
    break;

  default:
    MSQ_SETERR(err)("averaging method not available.",MsqError::INVALID_STATE);
    return 0.0;
  }
    
  return avg;
}

double AveragingQM::average_corner_hessian_diagonals( EntityTopology element_type,
                                              uint32_t ,
                                              unsigned num_corners,
                                              const double corner_values[],
                                              const Vector3D corner_grads[],
                                              const Matrix3D corner_hessians[],
                                              Vector3D vertex_grads[],
                                              SymMatrix3D vertex_hessians[],
                                              MsqError& err )
{
  return average_corner_diagonals( element_type, avgMethod, num_corners,
                                   corner_values, corner_grads,
                                   CornerHessDiagIterator( corner_hessians, element_type ),
                                   vertex_grads, vertex_hessians, err );
}

double AveragingQM::average_corner_hessian_diagonals( EntityTopology element_type,
                                              uint32_t ,
                                              unsigned num_corners,
                                              const double corner_values[],
                                              const Vector3D corner_grads[],
                                              const SymMatrix3D corner_hess_diag[],
                                              Vector3D vertex_grads[],
                                              SymMatrix3D vertex_hessians[],
                                              MsqError& err )
{
  return average_corner_diagonals( element_type, avgMethod, num_corners,
                                   corner_values, corner_grads,
                                   corner_hess_diag,
                                   vertex_grads, vertex_hessians, err );
}

static inline double sum_corner_hessians( EntityTopology type,
                                          unsigned num_corner,
                                          const double corner_values[],
                                          const Vector3D corner_grads[],
                                          const Matrix3D corner_hessians[],
                                          Vector3D vertex_grads[],
                                          Matrix3D vertex_hessians[] )
{
  const unsigned N = TopologyInfo::corners( type );
  unsigned i, n, r, c, R, C, idx[4];
  const unsigned* adj_list;
  double avg = 0.0;
  
    // calculate mean
  for (i = 0; i < num_corner; ++i)
    avg += corner_values[i];

  const Vector3D* grad = corner_grads;
  const Matrix3D* hess = corner_hessians;
  for (i = 0; i < num_corner; ++i)
  {
    adj_list = TopologyInfo::adjacent_vertices( type, i, n );
    idx[0] = i;
    idx[1] = adj_list[0];
    idx[2] = adj_list[1];
    idx[3] = adj_list[2%n]; // %n so don't read off end if 2D

    for (r = 0; r <= n; ++r) 
    {
      R = idx[r];
      vertex_grads[R] += *grad;
      ++grad;
      for (c = r; c <= n; ++c) 
      {
        C = idx[c];
        if (R <= C)
          vertex_hessians[N*R - R*(R+1)/2 + C] += *hess;
        else 
          vertex_hessians[N*C - C*(C+1)/2 + R].plus_transpose_equal(*hess);
        ++hess;
      }
    }
  }
  return avg;
}

static inline double sum_sqr_corner_hessians( EntityTopology type,
                                              unsigned num_corner,
                                              const double corner_values[],
                                              const Vector3D corner_grads[],
                                              const Matrix3D corner_hessians[],
                                              Vector3D vertex_grads[],
                                              Matrix3D vertex_hessians[] )
{
  const unsigned N = TopologyInfo::corners( type );
  unsigned i, n, r, c, R, C, idx[4];
  const unsigned* adj_list;
  double v, avg = 0.0;
  Matrix3D op;
  
    // calculate mean
  for (i = 0; i < num_corner; ++i)
    avg += corner_values[i]*corner_values[i];

  const Vector3D* grad = corner_grads;
  const Matrix3D* hess = corner_hessians;
  for (i = 0; i < num_corner; ++i)
  {
    adj_list = TopologyInfo::adjacent_vertices( type, i, n );
    idx[0] = i;
    idx[1] = adj_list[0];
    idx[2] = adj_list[1];
    idx[3] = adj_list[2%n]; // %n so don't read off end if 2D
    ++n;

    v = 2.0*corner_values[i];
    for (r = 0; r < n; ++r) 
    {
      R = idx[r];
      vertex_grads[R] += v*grad[r];
      for (c = r; c < n; ++c) 
      {
        C = idx[c];
        op.outer_product( 2.0*grad[r], grad[c] );
        op += v * *hess;
        if (R <= C)
          vertex_hessians[N*R - R*(R+1)/2 + C] += op;
        else 
          vertex_hessians[N*C - C*(C+1)/2 + R].plus_transpose_equal(op);
        ++hess;
      }
    }
    grad += n;
  }
  return avg;
}

static inline double pmean_corner_hessians( EntityTopology type,
                                            unsigned num_corner,
                                            const double corner_values[],
                                            const Vector3D corner_grads[],
                                            const Matrix3D corner_hessians[],
                                            Vector3D vertex_grads[],
                                            Matrix3D vertex_hessians[],
                                            double p )
{
  const unsigned N = TopologyInfo::corners( type );
  unsigned i, n, r, c, R, C, idx[4];
  const unsigned* adj_list;
  double m = 0.0, nm;
  Matrix3D op;
  double gf[8], hf[8];
  double inv = 1.0/num_corner;
  assert(num_corner <= 8);
  
    // calculate mean
  for (i = 0; i < num_corner; ++i)
  {
    nm = pow(corner_values[i], p);
    m += nm;
    
    gf[i] = inv * p * nm / corner_values[i];
    hf[i] = (p-1) * gf[i] / corner_values[i];
  }
  nm = inv * m;
  
  const Vector3D* grad = corner_grads;
  const Matrix3D* hess = corner_hessians;
  for (i = 0; i < num_corner; ++i)
  {
    adj_list = TopologyInfo::adjacent_vertices( type, i, n );
    idx[0] = i;
    idx[1] = adj_list[0];
    idx[2] = adj_list[1];
    idx[3] = adj_list[2%n]; // %n so don't read off end if 2D
    ++n;

    for (r = 0; r < n; ++r) 
    {
      R = idx[r];
      vertex_grads[R] += gf[i]*grad[r];
      for (c = r; c < n; ++c) 
      {
        C = idx[c];
        op.outer_product( grad[r], grad[c] );
        op *= hf[i];
        op += gf[i] * *hess;
        if (R <= C)
          vertex_hessians[N*R - R*(R+1)/2 + C] += op;
        else 
          vertex_hessians[N*C - C*(C+1)/2 + R].plus_transpose_equal(op);
        ++hess;
      }
    }
    grad += n;
  }
  
  m = pow( nm, 1.0/p );
  gf[0] = m / (p * nm );
  hf[0] = (1.0/p - 1) * gf[0] / nm;
  for (r = 0; r < N; ++r) 
  {
    for (c = r; c < N; ++c)
    {
      op.outer_product( vertex_grads[r], vertex_grads[c] );
      op *= hf[0];
      *vertex_hessians *= gf[0];
      *vertex_hessians += op;
      ++vertex_hessians;
    }
    vertex_grads[r] *= gf[0];
  }
  
  return m;
}


double AveragingQM::average_corner_hessians( EntityTopology type,
                                     uint32_t ,
                                     unsigned num_corner,
                                     const double corner_values[],
                                     const Vector3D corner_grads[],
                                     const Matrix3D corner_hessians[],
                                     Vector3D vertex_grads[],
                                     Matrix3D vertex_hessians[],
                                     MsqError& err )
{
  unsigned i;
  double avg, inv;
  
    // Zero gradients and Hessians
  const unsigned num_vertex = TopologyInfo::corners( type );
  for (i = 0; i < num_vertex; ++i)
    vertex_grads[i].set(0.0);
  const unsigned num_hess = num_vertex * (num_vertex+1) / 2;
  for (i = 0; i < num_hess; ++i)
    vertex_hessians[i].zero();
  
  switch (avgMethod)
  {
  case QualityMetric::SUM:
    avg = sum_corner_hessians( type, num_corner, corner_values, 
                               corner_grads, corner_hessians,
                               vertex_grads, vertex_hessians );
    break;
  
  case QualityMetric::LINEAR:
    avg = sum_corner_hessians( type, num_corner, corner_values, 
                               corner_grads, corner_hessians,
                               vertex_grads, vertex_hessians );
    inv = 1.0/num_corner;
    avg *= inv;
    for (i = 0; i < num_vertex; ++i)
      vertex_grads[i] *= inv;
    for (i = 0; i < num_hess; ++i)
      vertex_hessians[i] *= inv;
    break;

  case QualityMetric::SUM_SQUARED:
    avg = sum_sqr_corner_hessians( type, num_corner, corner_values, 
                                   corner_grads, corner_hessians,
                                   vertex_grads, vertex_hessians );
    break;

  case QualityMetric::RMS:
    avg = pmean_corner_hessians( type, num_corner, corner_values, 
                                 corner_grads, corner_hessians,
                                 vertex_grads, vertex_hessians,
                                 2.0 );
    break;

  case QualityMetric::HARMONIC:
    avg = pmean_corner_hessians( type, num_corner, corner_values, 
                                 corner_grads, corner_hessians,
                                 vertex_grads, vertex_hessians,
                                 -1.0 );
    break;

  case QualityMetric::HMS:
    avg = pmean_corner_hessians( type, num_corner, corner_values, 
                                 corner_grads, corner_hessians,
                                 vertex_grads, vertex_hessians,
                                 -2.0 );
    break;

  default:
    MSQ_SETERR(err)("averaging method not available.",MsqError::INVALID_STATE);
    return 0.0;
  }
    
  return avg;
}

double AveragingQM::average_metric_and_weights( double metrics[],
                                                  int count, 
                                                  MsqError& err )
{
  static bool min_max_warning = false;
  double avg = 0.0;
  int i, tmp_count;
  double f;
  
  switch (avgMethod)
  {
  
  case QualityMetric::MINIMUM:
    if (!min_max_warning) {
      MSQ_DBGOUT(1) <<
        "Minimum and maximum not continuously differentiable.\n"
        "Element of subdifferential will be returned.\n";
      min_max_warning = true;
    }

    avg = metrics[0];
    for (i = 1; i < count; ++i)
      if (metrics[i] < avg)
        avg = metrics[i];
    
    tmp_count = 0;
    for (i = 0; i < count; ++i)
    {
      if( metrics[i] - avg <= MSQ_MIN )
      {
        metrics[i] = 1.0;
        ++tmp_count;
      }
      else
      {
        metrics[i] = 0.0;
      }
    }
    
    f = 1.0 / tmp_count;
    for (i = 0; i < count; ++i)
      metrics[i] *= f;
      
    break;

  
  case QualityMetric::MAXIMUM:
    if (!min_max_warning) {
      MSQ_DBGOUT(1) <<
        "Minimum and maximum not continuously differentiable.\n"
        "Element of subdifferential will be returned.\n";
      min_max_warning = true;
    }

    avg = metrics[0];
    for (i = 1; i < count; ++i)
      if (metrics[i] > avg)
        avg = metrics[i];
    
    tmp_count = 0;
    for (i = 0; i < count; ++i)
    {
      if( avg - metrics[i] <= MSQ_MIN )
      {
        metrics[i] = 1.0;
        ++tmp_count;
      }
      else
      {
        metrics[i] = 0.0;
      }
    }
    
    f = 1.0 / tmp_count;
    for (i = 0; i < count; ++i)
      metrics[i] *= f;
      
    break;

  
  case QualityMetric::SUM:
    for (i = 0; i < count; ++i)
    {
      avg += metrics[i];
      metrics[i] = 1.0;
    }
      
    break;

  
  case QualityMetric::SUM_SQUARED:
    for (i = 0; i < count; ++i)
    {
      avg += (metrics[i]*metrics[i]);
      metrics[i] *= 2;
    }
      
    break;

  
  case QualityMetric::LINEAR:
    f = 1.0 / count;
    for (i = 0; i < count; ++i)
    {
      avg += metrics[i];
      metrics[i] = f;
    }
    avg *= f;
      
    break;

  
  case QualityMetric::GEOMETRIC:
    avg = 1.0;
    for (i = 0; i < count; ++i)
      avg *= metrics[i];
    avg = pow( avg, 1.0/count );

    f = avg / count;
    for (i = 0; i < count; ++i)
      metrics[i] = f / metrics[i];
      
    break;

  
  case QualityMetric::RMS:
    for (i = 0; i < count; ++i)
      avg += metrics[i] * metrics[i];
    avg = sqrt( avg / count );
    
    f = 1. / (avg*count);
    for (i = 0; i < count; ++i)
      metrics[i] *= f;
      
    break;

  
  case QualityMetric::HARMONIC:
    for (i = 0; i < count; ++i)
      avg += 1.0 / metrics[i];
    avg = count / avg;
  
    for (i = 0; i < count; ++i)
      metrics[i] = (avg * avg) / (count * metrics[i] * metrics[i]);
      
    break;

  
  case QualityMetric::HMS:
    for (i = 0; i < count; ++i)
      avg += 1. / (metrics[i] * metrics[i]);
    avg = sqrt( count / avg );
    
    f = avg*avg*avg / count;
    for (i = 0; i < count; ++i)
      metrics[i] = f / (metrics[i] * metrics[i] * metrics[i]);
      
    break;

  
  default:
    MSQ_SETERR(err)("averaging method not available.",MsqError::INVALID_STATE);
  }
  
  return avg;
}

   

   /*! 
     average_metrics takes an array of length num_value and averages the
     contents using averaging method 'method'.
   */
double AveragingQM::average_metrics( const double metric_values[],
                                     int num_values, MsqError &err)
{
    //MSQ_MAX needs to be made global?
  //double MSQ_MAX=1e10;
  double total_value=0.0;
  double temp_value=0.0;
  int i=0;
  int j=0;
    //if no values, return zero
  if (num_values<=0){
    return 0.0;
  }

  switch(avgMethod){
    case QualityMetric::GEOMETRIC:
       total_value=1.0;
       for (i=0;i<num_values;++i){
         total_value*=(metric_values[i]);
       }
       total_value=pow(total_value, 1.0/num_values);
       break;

    case QualityMetric::HARMONIC:
         //ensure no divide by zero, return zero
       for (i=0;i<num_values;++i){
         if(metric_values[i]<MSQ_MIN){
           if(metric_values[i]>MSQ_MIN){
             return 0.0;
           }
         }
         total_value+=(1/metric_values[i]);
       }
         //ensure no divide by zero, return MSQ_MAX_CAP
       if(total_value<MSQ_MIN){
         if(total_value>MSQ_MIN){
           return MSQ_MAX_CAP;
         }
       }
       total_value=num_values/total_value;
       break;

    case QualityMetric::LINEAR:
       for (i=0;i<num_values;++i){
         total_value+=metric_values[i];
       }
       total_value/= (double) num_values;
       break;

    case QualityMetric::MAXIMUM:
       total_value = metric_values[0];
       for (i = 1; i < num_values; ++i){
         if (metric_values[i] > total_value){
           total_value = metric_values[i];
         }
       }
       break;

    case QualityMetric::MINIMUM:
       total_value = metric_values[0];
       for (i = 1; i < num_values; ++i){
         if (metric_values[i] < total_value) {
           total_value = metric_values[i];
         }
       }
       break;

    case QualityMetric::RMS:
       for (i=0;i<num_values;++i){
         total_value+=(metric_values[i]*metric_values[i]);
       }
       total_value/= (double) num_values;
       total_value=sqrt(total_value);
       break;

    case QualityMetric::HMS:
       //ensure no divide by zero, return zero
       for (i=0;i<num_values;++i){
         if (metric_values[i]*metric_values[i] < MSQ_MIN) {
           return 0.0;
         }
         total_value += (1.0/(metric_values[i]*metric_values[i]));
       }

       //ensure no divide by zero, return MSQ_MAX_CAP
       if (total_value < MSQ_MIN) {
         return MSQ_MAX_CAP;
       }
       total_value = sqrt(num_values/total_value);
       break;

    case QualityMetric::STANDARD_DEVIATION:
       total_value=0;
       temp_value=0;
       for (i=0;i<num_values;++i){
         temp_value+=metric_values[i];
         total_value+=(metric_values[i]*metric_values[i]);
       }
       temp_value/= (double) num_values;
       temp_value*=temp_value;
       total_value/= (double) num_values;
       total_value=total_value-temp_value;
       break;

    case QualityMetric::SUM:
       for (i=0;i<num_values;++i){
         total_value+=metric_values[i];
       }
       break;

    case QualityMetric::SUM_SQUARED:
       for (i=0;i<num_values;++i){
         total_value+= (metric_values[i]*metric_values[i]);
       }
       break;

    case QualityMetric::MAX_MINUS_MIN:
       //total_value used to store the maximum
         //temp_value used to store the minimum
       temp_value=MSQ_MAX_CAP;
       for (i=0;i<num_values;++i){
         if(metric_values[i]<temp_value){
           temp_value=metric_values[i];
         }
         if(metric_values[i]>total_value){
           total_value=metric_values[i];
         }
       }

       total_value-=temp_value;
       break;

    case QualityMetric::MAX_OVER_MIN:
         //total_value used to store the maximum
         //temp_value used to store the minimum
       temp_value=MSQ_MAX_CAP;
       for (i=0;i<num_values;++i){
         if(metric_values[i]<temp_value){
           temp_value=metric_values[i];
         }
         if(metric_values[i]>total_value){
           total_value=metric_values[i];
         }
       }

         //ensure no divide by zero, return MSQ_MAX_CAP
       if (fabs(temp_value) < MSQ_MIN) {
         return MSQ_MAX_CAP;
       }
       total_value/=temp_value;
       break;

    case QualityMetric::SUM_OF_RATIOS_SQUARED:
       for (j=0;j<num_values;++j){
         //ensure no divide by zero, return MSQ_MAX_CAP
         if (fabs(metric_values[j]) < MSQ_MIN) {
           return MSQ_MAX_CAP;
         }
         for (i=0;i<num_values;++i){
           total_value+=((metric_values[i]/metric_values[j])*
                         (metric_values[i]/metric_values[j]));
         }
       }
       total_value/=(double)(num_values*num_values);
       break;

    default:
         //Return error saying Averaging Method mode not implemented
       MSQ_SETERR(err)("Requested Averaging Method Not Implemented", MsqError::NOT_IMPLEMENTED);
       return 0;
  }
  return total_value;
}



} // namespace Mesquite

