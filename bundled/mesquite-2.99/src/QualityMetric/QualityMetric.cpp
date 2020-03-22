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
    kraftche@cae.wisc.edu
   
  ***************************************************************** */
/*!
  \file   QualityMetric.cpp
  \brief  

  \author Michael Brewer
  \author Thomas Leurent
  \date   2002-05-14
  \author Jason Kraftcheck
  \date   2006-04-20
*/

#include "QualityMetric.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {

void QualityMetric::initialize_queue( MeshDomainAssoc* ,
                                      const Settings* ,
                                      MsqError&  )
{}
  
void QualityMetric::get_single_pass( PatchData& pd, 
                        std::vector<size_t>& handles,
                        bool free_vertices_only, 
                        MsqError& err )
{
  get_evaluations( pd, handles, free_vertices_only, err );
}


static inline 
double get_delta_C( const PatchData& pd,
                    const std::vector<size_t>& indices,
                    MsqError& err )
{
  if (indices.empty()) {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    return 0.0;
  }

  std::vector<size_t>::const_iterator beg, iter, iter2, end;
  
  std::vector<size_t> tmp_vect;
  if (indices.size() == 1u) {
    pd.get_adjacent_vertex_indices( indices.front(), tmp_vect, err );
    MSQ_ERRZERO(err);
    assert(!tmp_vect.empty());
    tmp_vect.push_back( indices.front() );
    beg = tmp_vect.begin();
    end = tmp_vect.end();
  }
  else {
    beg = indices.begin();
    end = indices.end();
  }
  
  double min_dist_sqr = HUGE_VAL, sum_dist_sqr = 0.0;
  for (iter = beg; iter != end; ++iter) {
    for (iter2 = iter+1; iter2 != end; ++iter2) {
      Vector3D diff = pd.vertex_by_index(*iter);
      diff -= pd.vertex_by_index(*iter2);
      double dist_sqr = diff % diff;
      if (dist_sqr < min_dist_sqr)
        min_dist_sqr = dist_sqr;
      sum_dist_sqr += dist_sqr;
    }
  }
  
  return 3e-6*sqrt(min_dist_sqr) + 5e-7*sqrt(sum_dist_sqr/(end-beg));
}
    

bool QualityMetric::evaluate_with_gradient( PatchData& pd,
                              size_t handle,
                              double& value,
                              std::vector<size_t>& indices,
                              std::vector<Vector3D>& gradient,
                              MsqError& err )
{
  indices.clear();
  bool valid = evaluate_with_indices( pd, handle, value, indices, err);
  if (MSQ_CHKERR(err) || !valid)
    return false;
  if (indices.empty())
    return true;

    // get initial pertubation amount
  double delta_C = finiteDiffEps;
  if (!haveFiniteDiffEps) {
    delta_C = get_delta_C( pd, indices, err ); MSQ_ERRZERO(err);
    if (keepFiniteDiffEps) {
      finiteDiffEps = delta_C;
      haveFiniteDiffEps = true;
    }
  }
  const double delta_inv_C = 1.0/delta_C;
  const int reduction_limit = 15;

  gradient.resize( indices.size() );
  for (size_t v=0; v<indices.size(); ++v) 
  {
    const Vector3D pos = pd.vertex_by_index(indices[v]);
    
    /* gradient in the x, y, z direction */
    for (int j=0;j<3;++j) 
    {
      double delta = delta_C;
      double delta_inv = delta_inv_C;
      double metric_value;
      Vector3D delta_v( 0, 0, 0 );
      
        //perturb the node and calculate gradient.  The while loop is a
        //safety net to make sure the epsilon perturbation does not take
        //the element out of the feasible region.
      int counter = 0;
      for (;;) {
          // perturb the coordinates of the free vertex in the j direction
          // by delta       
        delta_v[j] = delta;
        pd.set_vertex_coordinates( pos+delta_v, indices[v], err ); MSQ_ERRZERO(err);

          //compute the function at the perturbed point location
        valid = evaluate( pd, handle, metric_value, err);
        if (!MSQ_CHKERR(err) && valid) 
          break;
          
        if (++counter >= reduction_limit) {
          MSQ_SETERR(err)("Perturbing vertex by delta caused an inverted element.",
                          MsqError::INTERNAL_ERROR);
          return false;
        }
        
        delta*=0.1;
        delta_inv*=10.;
      }
        // put the coordinates back where they belong
      pd.set_vertex_coordinates( pos, indices[v], err );
        // compute the numerical gradient
      gradient[v][j] = (metric_value - value) * delta_inv;
    } // for(j)
  } // for(indices)
  return true;
}
   

bool QualityMetric::evaluate_with_Hessian( PatchData& pd,
                              size_t handle,
                              double& value,
                              std::vector<size_t>& indices,
                              std::vector<Vector3D>& gradient,
                              std::vector<Matrix3D>& Hessian,
                              MsqError& err )
{
  indices.clear();
  gradient.clear();
  keepFiniteDiffEps = true;
  bool valid = evaluate_with_gradient( pd, handle, value, indices, gradient, err );
  keepFiniteDiffEps = false;
  if (MSQ_CHKERR(err) || !valid) {
    haveFiniteDiffEps = false;
    return false;
  }
  if (indices.empty()){
    haveFiniteDiffEps = false;
    return true;
  }

    // get initial pertubation amount
  double delta_C;
  if (haveFiniteDiffEps) {
    delta_C = finiteDiffEps;
  }
  else {
    delta_C = get_delta_C( pd, indices, err ); MSQ_ERRZERO(err);
  }
  assert(delta_C < 1e30);
  const double delta_inv_C = 1.0/delta_C;
  const int reduction_limit = 15;

  std::vector<Vector3D> temp_gradient( indices.size() );
  const int num_hess = indices.size() * (indices.size() + 1) / 2;
  Hessian.resize( num_hess );
  
  for (unsigned v = 0; v < indices.size(); ++v) {
    const Vector3D pos = pd.vertex_by_index(indices[v]);
    for (int j = 0; j < 3; ++j ) { // x, y, and z
      double delta = delta_C;
      double delta_inv = delta_inv_C;
      double metric_value;
      Vector3D delta_v(0,0,0);
      
        // find finite difference for gradient
      int counter = 0;
      for (;;) {
        delta_v[j] = delta;
        pd.set_vertex_coordinates( pos+delta_v, indices[v], err ); MSQ_ERRZERO(err);
        valid = evaluate_with_gradient( pd, handle, metric_value, indices, temp_gradient, err );
        if (!MSQ_CHKERR(err) && valid) 
          break;
        
        if (++counter >= reduction_limit) {
          MSQ_SETERR(err)("Algorithm did not successfully compute element's "
                           "Hessian.\n",MsqError::INTERNAL_ERROR);
          haveFiniteDiffEps = false;
          return false;
        }
        
        delta *= 0.1;
        delta_inv *= 10.0;
      }
      pd.set_vertex_coordinates( pos, indices[v], err ); MSQ_ERRZERO(err);
      
        //compute the numerical Hessian
      for (unsigned w = 0; w <= v; ++w) {
          //finite difference to get some entries of the Hessian
        Vector3D fd( temp_gradient[w] );
        fd -= gradient[w];
        fd *= delta_inv;
          // For the block at position w,v in a matrix, we need the corresponding index
          // (mat_index) in a 1D array containing only upper triangular blocks.
        unsigned sum_w = w*(w+1)/2; // 1+2+3+...+w
        unsigned mat_index = w*indices.size() + v - sum_w;
        Hessian[mat_index][0][j] = fd[0];
        Hessian[mat_index][1][j] = fd[1];
        Hessian[mat_index][2][j] = fd[2];
      }
    }  // for(j)
  } // for(indices)
  haveFiniteDiffEps = false;
  return true;
}

bool QualityMetric::evaluate_with_Hessian_diagonal( PatchData& pd,
                                size_t handle,
                                double& value,
                                std::vector<size_t>& indices,
                                std::vector<Vector3D>& gradient,
                                std::vector<SymMatrix3D>& Hessian_diagonal,
                                MsqError& err )
{
  bool rval = evaluate_with_Hessian( pd, handle, value, indices, gradient, tmpHess, err );
  if (MSQ_CHKERR(err) || !rval)
    return rval;
  size_t s = indices.size();
  Hessian_diagonal.resize( s );
  std::vector<Matrix3D>::const_iterator h = tmpHess.begin();
  for (size_t i = 0; i < indices.size(); ++i) {
    Hessian_diagonal[i] = h->upper();
    h += s--;
  }
  return rval;
}


uint32_t QualityMetric::fixed_vertex_bitmap( PatchData& pd,
                                           const MsqMeshEntity* elem,
                                           std::vector<size_t>& indices )
{
  indices.clear();
  uint32_t result = ~(uint32_t)0;
  unsigned num_vtx = elem->vertex_count();
  const size_t* vertices = elem->get_vertex_index_array();
  indices.clear();
  for (unsigned i = 0; i < num_vtx; ++i) {
    if (vertices[i] < pd.num_free_vertices()) {
      indices.push_back( vertices[i] );
      result &= ~(uint32_t)(1<<i);
    }
  }
  return result;
}
  

void QualityMetric::remove_fixed_gradients( EntityTopology elem_type,
                                          uint32_t fixed,
                                          std::vector<Vector3D>& grads )
{
  const unsigned num_vertex = TopologyInfo::corners( elem_type );
  unsigned r, w;
  for (r = 0; r < num_vertex && !(fixed & (1<<r)); ++r);
  for (w = r++; r < num_vertex; ++r) {
    if (!(fixed & (1<<r))) {
      grads[w] = grads[r];
      ++w;
    }
  }
  grads.resize(w);
}

void QualityMetric::remove_fixed_diagonals( EntityTopology type, 
                                          uint32_t fixed, 
                                          std::vector<Vector3D>& grads,
                                          std::vector<SymMatrix3D>& diags )
{
  const unsigned num_vertex = TopologyInfo::corners( type );
  unsigned r, w;
  for (r = 0; r < num_vertex && !(fixed & (1<<r)); ++r);
  for (w = r++; r < num_vertex; ++r) {
    if (!(fixed & (1<<r))) {
      grads[w] = grads[r];
      diags[w] = diags[r];
      ++w;
    }
  }
  grads.resize(w);
  diags.resize(w);
}

void QualityMetric::remove_fixed_hessians( EntityTopology elem_type,
                                         uint32_t fixed,
                                         std::vector<Matrix3D>& hessians )
{
  const unsigned num_vertex = TopologyInfo::corners( elem_type );
  unsigned r, c, i = 0, w = 0;
  for (r = 0; r < num_vertex; ++r) {
    if (fixed & (1<<r)) {
      i += num_vertex - r;
      continue;
    }
    for (c = r; c < num_vertex; ++c) {
      if (!(fixed & (1<<c))) {
        hessians[w] = hessians[i];
        ++w;
      }
      ++i;
    }
  }
  hessians.resize(w);
}      

double QualityMetric::weighted_average_metrics(const double coef[],
                                             const double metric_values[],
                                             const int& num_values, 
                                             MsqError &err)
{
  //MSQ_MAX needs to be made global?
  //double MSQ_MAX=1e10;
  double total_value=0.0;
  int i=0;
  //if no values, return zero
  if (num_values<=0){
    return 0.0;
  }

  for (i=0;i<num_values;++i){
    total_value += coef[i]*metric_values[i];
  }
  total_value /= (double) num_values;

  return total_value;
}
   
} // namespace Mesquite
