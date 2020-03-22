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


/** \file TMPQualityMetric.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#undef PRINT_INFO

#include "Mesquite.hpp"
#include "TMPQualityMetric.hpp"
#include "MsqMatrix.hpp"
#include "ElementQM.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include "MappingFunction.hpp"
#include "WeightCalculator.hpp"
#include "TargetCalculator.hpp"
#include "TargetMetricUtil.hpp"

#ifdef PRINT_INFO
#  include <iostream>
#endif

#include <functional>
#include <algorithm>

namespace MESQUITE_NS {

int TMPQualityMetric::get_negate_flag( ) const { return 1; }


void TMPQualityMetric::get_evaluations( PatchData& pd,
                                      std::vector<size_t>& handles,
                                      bool free,
                                      MsqError& err )
{
  get_sample_pt_evaluations( pd, handles, free, err );
}

void TMPQualityMetric::get_patch_evaluations( PatchData& pd,
                                      std::vector<size_t>& handles,
                                      bool free,
                                      MsqError& err )
{
  get_sample_pt_evaluations( pd, handles, free, err );
}

void TMPQualityMetric::get_element_evaluations( PatchData& pd,
                                              size_t elem,
                                              std::vector<size_t>& handles,
                                              MsqError& err )
{
  get_elem_sample_points( pd, elem, handles, err );
}

bool TMPQualityMetric::evaluate( PatchData& pd, size_t handle, double& value, MsqError& err )
{
  size_t num_idx;
  bool valid = evaluate_internal( pd, handle, value, mIndices, num_idx, err );
  if (MSQ_CHKERR(err) || !valid)
    return false;
  
    // apply target weight to value
  if (weightCalc) {
    const Sample s = ElemSampleQM::sample( handle );
    const size_t e = ElemSampleQM::  elem( handle );
    double ck = weightCalc->get_weight( pd, e, s, err ); MSQ_ERRZERO(err);
    value *= ck;
  }
  return true;
}


bool TMPQualityMetric::evaluate_with_indices( PatchData& pd,
                                              size_t handle,
                                              double& value,
                                              std::vector<size_t>& indices,
                                              MsqError& err )
{
  indices.resize( MAX_ELEM_NODES );
  size_t num_idx = 0;
  bool result = evaluate_internal( pd, handle, value, arrptr(indices), num_idx, err );
  if (MSQ_CHKERR(err) || !result)
    return false;

  indices.resize( num_idx );
  
    // apply target weight to value
  if (weightCalc) {
    const Sample s = ElemSampleQM::sample( handle );
    const size_t e = ElemSampleQM::  elem( handle );
    double ck = weightCalc->get_weight( pd, e, s, err ); MSQ_ERRZERO(err);
    value *= ck;
  }

  return true;
}

static void get_u_perp( const MsqVector<3>& u,
                        MsqVector<3>& u_perp )
{
  double a = sqrt(u[0]*u[0] + u[1]*u[1]);
  if (a < 1e-10) {
    u_perp[0] = 1.0;
    u_perp[1] = u_perp[2] = 0.0;
  }
  else {
    double b = -u[2]/a;
    u_perp[0] = u[0]*b;
    u_perp[1] = u[1]*b;
    u_perp[2] = a;
  }
}



/* Do transform M_hat = S_a M_{3x2}, M_{2x2} Theta^-1 M_hat
 * where the plane into which we are projecting is orthogonal
 * to the passed u vector.
 */
static inline bool
project_to_perp_plane(  MsqMatrix<3,2> J,
                        const MsqVector<3>& u,
                        const MsqVector<3>& u_perp,
                        MsqMatrix<2,2>& A,
                        MsqMatrix<3,2>& S_a_transpose_Theta )
{
  MsqVector<3> n_a = J.column(0) * J.column(1);
  double sc, len = length(n_a);
  if (!divide(1.0, len, sc))
    return false;
  n_a *= sc;
  double ndot = n_a % u;
  double sigma = (ndot < 0.0) ? -1 : 1;
  double cosphi = sigma * ndot;
  MsqVector<3> cross = n_a * u;
  double sinphi = length(cross);

  MsqMatrix<3,2> Theta;
  Theta.set_column(0,   u_perp);
  Theta.set_column(1, u*u_perp);

    // If columns of J are not in plane orthogonal to u, then
    // rotate J such that they are.
  if (sinphi > 1e-12) {
    MsqVector<3> m = sigma * cross;
    MsqVector<3> n = (1/sinphi) * m;
    MsqVector<3> p = (1-cosphi) * n;
    double s_a[] = 
      { p[0]*n[0] + cosphi, p[0]*n[1] - m[2],   p[0]*n[2] + m[1],
        p[1]*n[0] + m[2],   p[1]*n[1] + cosphi, p[1]*n[2] - m[0],
        p[2]*n[0] - m[1],   p[2]*n[1] + m[0],   p[2]*n[2] + cosphi };
    MsqMatrix<3,3> S_a(s_a);
    J = S_a * J;
    S_a_transpose_Theta = transpose(S_a) * Theta;
  } 
  else {
    S_a_transpose_Theta = Theta;
//    J *= sigma;
  }

    // Project to get 2x2 A from A_hat (which might be equal to J)
  A = transpose(Theta) * J;
  return true;
}

/* Do transform M_hat = S_a M_{3x2}, M_{2x2} Theta^-1 M_hat
 * where the plane into which we are projecting is the cross
 * product of the columns of M, such that S_a is I.  Use the
 * first column of M as u_perp.  
 *
 * Also pass back the cross product of the columns of M as u,
 * and the first column of M as u_perp, both normalized.
 */
static inline void
project_to_matrix_plane( const MsqMatrix<3,2>& M_in,
                         MsqMatrix<2,2>& M_out,
                         MsqVector<3>& u,
                         MsqVector<3>& u_perp )
{
  u = M_in.column(0) * M_in.column(1);
  u_perp = M_in.column(0);
  double len0 = length(u_perp);
  double u_len = length(u);
  double d_perp, d_u;
  if (!divide(1.0, len0, d_perp)) {
    // try the other column
    u_perp = M_in.column(1);
    len0 = length(u_perp);
    if (!divide(1.0, len0, d_perp)) {
      // matrix is all zeros
      u[0] = 0; u[1] = 0; u[2] = 1;
      u_perp[0] = 1; u_perp[1] = 0; u_perp[2] = 0;
      M_out = MsqMatrix<2,2>(0.0);
    }
    else {
      MsqMatrix<3,2> junk;
      get_u_perp( u_perp, u );
      project_to_perp_plane( M_in, u, u_perp, M_out, junk );
    }
  }
  else if (!divide(1.0, u_len, d_u )) {
    MsqMatrix<3,2> junk;
    get_u_perp( u_perp, u );
    project_to_perp_plane( M_in, u, u_perp, M_out, junk );
  }
  else { // the normal case (neither column is zero)
    u *= d_u;
    u_perp *= d_perp;

     // M_out = transpose(theta)*M_in
    M_out(0,0) = len0;
    M_out(0,1) = u_perp % M_in.column(1);
    M_out(1,0) = 0.0;
    M_out(1,1) = u_len / len0;
  }
}

bool
TMPQualityMetric::evaluate_surface_common( PatchData& pd,
                                           Sample s,
                                           size_t e,
                                           const NodeSet& bits,
                                           size_t* indices,
                                           size_t& num_indices,
                                           MsqVector<2>* derivs,
                                           MsqMatrix<2,2>& W,
                                           MsqMatrix<2,2>& A,
                                           MsqMatrix<3,2>& S_a_transpose_Theta,
                                           MsqError& err )
{
  EntityTopology type = pd.element_by_index( e ).get_element_type();

  const MappingFunction2D* mf = pd.get_mapping_function_2D( type );
  if (!mf) {
    MSQ_SETERR(err)( "No mapping function for element type", MsqError::UNSUPPORTED_ELEMENT );
    return false;
  }

  MsqMatrix<3,2> J;
  mf->jacobian( pd, e, bits, s, indices, derivs, num_indices, J, err );

    // If we have a 3x2 target matrix 
  if (targetCalc->have_surface_orient()) {
    MsqVector<3> u, u_perp;
    MsqMatrix<3,2> W_hat;
    targetCalc->get_surface_target( pd, e, s, W_hat, err ); MSQ_ERRZERO(err);
      // Use the cross product of the columns of W as the normal of the 
      // plane to work in (i.e. u.).  W should have been constructed such
      // that said cross product is in the direction of (n_s)_init.  And if
      // for some reason it as not, then using something other than said
      // cross product is likely to produce very wrong results.
    project_to_matrix_plane( W_hat, W, u, u_perp );
      // Do the transforms on A to align it with W and project into the plane.
    if (!project_to_perp_plane( J, u, u_perp, A, S_a_transpose_Theta ))
      return false;
  }
    // Otherwise if we have a 2x2 target matrix (i.e. the target does
    // not contain orientation information), project into the plane
    // tangent to J.
  else {
    MsqVector<3> u, u_perp;
    targetCalc->get_2D_target( pd, e, s, W, err ); MSQ_ERRZERO(err);
    project_to_matrix_plane( J, A, u, u_perp );
    S_a_transpose_Theta.set_column(0, u_perp);
    S_a_transpose_Theta.set_column(1, u*u_perp);
      // If the domain is set, adjust the sign of things correctly
      // for the case where the element is inverted with respect
      // to the domain.
    if (pd.domain_set()) {
      Vector3D n;
      pd.get_domain_normal_at_sample( e, s, n, err );
      MSQ_ERRZERO(err);
        // if sigma == -1
      if (Vector3D(u.data()) % n < 0.0) {
          // flip u
        u = -u;
          // S_a_transpose_Theta == Theta, because S_a == I here.
          // u_perp is unaffected by flipping u, so only the second
          // column of S_a_transpose_Theta and the second row of A
          // are flipped because u x u_perp will be flipped.
        S_a_transpose_Theta.set_column(1, -S_a_transpose_Theta.column(1) );
        A.set_row( 1, -A.row(1) );
      }
    }
  }
  
  return true;
}                    

void TMPQualityMetric::weight( PatchData& pd,
                               Sample sample,
                               size_t elem,
                               int num_idx,
                               double& value,
                               Vector3D* grad,
                               SymMatrix3D* diag,
                               Matrix3D* hess,
                               MsqError& err )
{
  if (!weightCalc)
    return;
  
  double ck = weightCalc->get_weight( pd, elem, sample, err ); MSQ_ERRRTN(err);
  value *= ck;
  if (grad) {
    for (int i = 0; i < num_idx; ++i)
      grad[i] *= ck;
  }
  if (diag) {
    for (int i = 0; i < num_idx; ++i)
      diag[i] *= ck;
  }
  if (hess) {
    const int n = num_idx * (num_idx+1) / 2;
    for (int i = 0; i < n; ++i)
      hess[i] *= ck;
  }
}

void TMPQualityMetric::initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                         const Settings* settings,
                                         MsqError& err )
{
  targetCalc->initialize_queue( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);
  if (weightCalc) {
    weightCalc->initialize_queue( mesh_and_domain, settings, err ); 
    MSQ_ERRRTN(err);
  }
}


} // namespace Mesquite
