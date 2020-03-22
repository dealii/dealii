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

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TrustRegion.cpp
 *  \brief Port Todd Munson's trust region solver to Mesquite
 *  \author Jason Kraftcheck (Mesquite Port)
 */

#include "Mesquite.hpp"
#include "TrustRegion.hpp"
#include "MsqDebug.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"

#define USE_FN_PC1    // Use 1st preconditioner from Todd's code
                      // (alternate is whatever is in MsqHessian already)
#undef  DO_STEEP_DESC // Jason's apparently broken hack to fall back to
                      // steepest descent search direction

namespace MESQUITE_NS {

// Force std::vector to release allocated memory
template <typename T>
static inline void free_vector( std::vector<T>& v )
{
  std::vector<T> temp;
  temp.swap( v );
}

std::string TrustRegion::get_name() const { return "TrustRegion"; }

PatchSet* TrustRegion::get_patch_set()
  { return PatchSetUser::get_patch_set(); }

TrustRegion::TrustRegion( ObjectiveFunction* of )
  : VertexMover( of ), PatchSetUser(true), mMemento(0)
{}

TrustRegion::~TrustRegion()
{
  delete mMemento;
  mMemento = 0;
}


void TrustRegion::initialize( PatchData& pd, MsqError& err )
{
  mMemento = pd.create_vertices_memento( err );
  MSQ_CHKERR(err);
}

void TrustRegion::initialize_mesh_iteration( PatchData& /*pd*/, MsqError& /*err*/ )
{ }

void TrustRegion::terminate_mesh_iteration( PatchData& /*pd*/, MsqError& /*err*/ )
{ }

void TrustRegion::cleanup()
{ 
    // release Memento
  delete mMemento;
  mMemento = 0;
    // release temporary array memory
  mHess.clear();
  free_vector(mGrad);
  free_vector(wVect);
  free_vector(zVect);
  free_vector(dVect);
  free_vector(pVect);
  free_vector(rVect);
  free_vector(preCond);
}

static inline void negate( Vector3D* out, const Vector3D* in, size_t nn )
{
  for (size_t i = 0; i < nn; ++i)
    out[i] = -in[i];
}

// Do v += s * x, where v and x are arrays of length n
static inline void plus_eq_scaled( Vector3D* v, double s, const Vector3D* x, size_t n )
{
  Vector3D* end = v + n;
  for (; v != end; ++v, ++x)
    *v += s * *x;
}

// Do v = s*v - x, where v and x are arrays of length n
static inline void times_eq_minus( Vector3D* v, double s, const Vector3D* x, size_t n )
{
  Vector3D* end = v + n;
  for (; v != end; ++v, ++x) {
    *v *= s;
    *v -= *x;
  }
}

void TrustRegion::compute_preconditioner( MsqError& err )
{
#ifndef USE_FN_PC1
  mHessian.calculate_preconditioner(err);
#else
  double dia;
  preCond.resize( mHess.size() );
  for (size_t i = 0; i < mHess.size(); ++i) {
    const Matrix3D& m = *mHess.get_block(i,i);
    dia = m[0][0] + m[1][1] + m[2][2];
    preCond[i] = dia < DBL_EPSILON ? 1.0 : 1.0/dia;
  }
#endif
}

void TrustRegion::apply_preconditioner( Vector3D* z, Vector3D* r, MsqError& err )
{
#ifndef USE_FN_PC1
  mHessian.apply_preconditioner( z, r, err );
#else
  for (size_t i = 0; i < preCond.size(); ++i) 
    z[i] = preCond[i] * r[i];
#endif
}

void TrustRegion::optimize_vertex_positions( PatchData& pd, MsqError& err )
{
  TerminationCriterion& term = *get_inner_termination_criterion();
  OFEvaluator& func = get_objective_function_evaluator();
  
  const double cg_tol = 1e-2;
  const double eta_1  = 0.01;
  const double eta_2  = 0.90;
  const double tr_incr = 10;
  const double tr_decr_def = 0.25;
  const double tr_decr_undef = 0.25;
  const double tr_num_tol = 1e-6;
  const int max_cg_iter = 10000;
  
  double radius = 1000;		/* delta*delta */

 const int nn = pd.num_free_vertices();
  wVect.resize(nn); Vector3D* w = arrptr(wVect);
  zVect.resize(nn); Vector3D* z = arrptr(zVect);
  dVect.resize(nn); Vector3D* d = arrptr(dVect);
  pVect.resize(nn); Vector3D* p = arrptr(pVect);
  rVect.resize(nn); Vector3D* r = arrptr(rVect);

  double norm_r, norm_g;
  double alpha, beta, kappa;
  double rz, rzm1;
  double dMp, norm_d, norm_dp1, norm_p;
  double obj, objn;

  int cg_iter;
  bool valid;

  mHess.initialize( pd, err );  //hMesh(mesh);
  valid = func.update( pd, obj, mGrad, mHess, err ); MSQ_ERRRTN(err);
  if (!valid) {
    MSQ_SETERR(err)("Initial objective function is not valid", MsqError::INVALID_MESH);
    return;
  }
  compute_preconditioner( err ); MSQ_ERRRTN(err);
  pd.recreate_vertices_memento( mMemento, err ); MSQ_ERRRTN(err);

  while (!term.terminate() && (radius > 1e-20)) {

    norm_r = length_squared(arrptr(mGrad), nn);
    norm_g = sqrt(norm_r);

    memset(d, 0, 3*sizeof(double)*nn);
    memcpy(r, arrptr(mGrad), nn*sizeof(Vector3D)); //memcpy(r, mesh->g, 3*sizeof(double)*nn);
    norm_g *= cg_tol;

    apply_preconditioner( z, r, err); MSQ_ERRRTN(err); //prec->apply(z, r, prec, mesh);
    negate(p, z, nn);
    rz = inner(r, z, nn);

    dMp    = 0;
    norm_p = rz;
    norm_d = 0;

    cg_iter = 0;
    while ((sqrt(norm_r) > norm_g) && 
#ifdef DO_STEEP_DESC
         (norm_d > tr_num_tol) && 
#endif
         (cg_iter < max_cg_iter)) 
    {
      ++cg_iter;

      memset(w, 0, 3*sizeof(double)*nn);
      //matmul(w, mHess, p); //matmul(w, mesh, p);
      mHess.product( w, p );

      kappa = inner(p, w, nn);
      if (kappa <= 0.0) {
        alpha = (sqrt(dMp*dMp+norm_p*(radius-norm_d))-dMp)/norm_p;
        plus_eq_scaled( d, alpha, p, nn );
	break;
      }

      alpha = rz / kappa;

      norm_dp1 = norm_d + 2.0*alpha*dMp + alpha*alpha*norm_p;
      if (norm_dp1 >= radius) {
        alpha = (sqrt(dMp*dMp+norm_p*(radius-norm_d))-dMp)/norm_p;
        plus_eq_scaled( d, alpha, p, nn );
	break;
      }

      plus_eq_scaled( d, alpha, p, nn );
      plus_eq_scaled( r, alpha, w, nn );
      norm_r = length_squared(r, nn);

      apply_preconditioner( z, r, err); MSQ_ERRRTN(err); //prec->apply(z, r, prec, mesh);

      rzm1 = rz;
      rz = inner(r, z, nn);
      beta = rz / rzm1;
      times_eq_minus( p, beta, z, nn );

      dMp = beta*(dMp + alpha*norm_p);
      norm_p = rz + beta*beta*norm_p;
      norm_d = norm_dp1;
    }

#ifdef DO_STEEP_DESC    
    if (norm_d <= tr_num_tol) {
      norm_g = length(arrptr(mGrad), nn);
      double ll = 1.0;
      if (norm_g < tr_num_tol)
        break;
      if (norm_g > radius)
        ll = radius / nurm_g;
      for (int i = 0; i < nn; ++i)
        d[i] = ll * mGrad[i];
    }
#endif

    alpha = inner( arrptr(mGrad), d, nn ); // inner(mesh->g, d, nn);

    memset(p, 0, 3*sizeof(double)*nn);
    //matmul(p, mHess, d); //matmul(p, mesh, d);
    mHess.product( p, d );
    beta = 0.5*inner(p, d, nn);
    kappa = alpha + beta;

    /* Put the new point into the locations */
    pd.move_free_vertices_constrained( d, nn, 1.0, err ); MSQ_ERRRTN(err);

    valid = func.evaluate( pd, objn, err ); 
    if (err.error_code() == err.BARRIER_VIOLATED)
      err.clear();  // barrier violated does not represent an actual error here
    MSQ_ERRRTN(err);

    if (!valid) {
      /* Function not defined at trial point */
      radius *= tr_decr_undef;
      pd.set_to_vertices_memento( mMemento, err ); MSQ_ERRRTN(err); 
      continue;
    }
      

    if ((fabs(kappa) <= tr_num_tol) && (fabs(objn - obj) <= tr_num_tol)) {
      kappa = 1;
    }
    else {
      kappa = (objn - obj) / kappa;
    }
    
    if (kappa < eta_1) {
      /* Iterate is unacceptable */
      radius *= tr_decr_def;
      pd.set_to_vertices_memento( mMemento, err ); MSQ_ERRRTN(err); 
      continue;
    }

      /* Iterate is acceptable */
    if (kappa >= eta_2) {
      /* Iterate is a very good step, increase radius */
      radius *= tr_incr;
      if (radius > 1e20) {
	radius = 1e20;
      }
    }

    func.update( pd, obj, mGrad, mHess, err );
    compute_preconditioner( err ); MSQ_ERRRTN(err);
    pd.recreate_vertices_memento( mMemento, err ); MSQ_ERRRTN(err);

    // checks stopping criterion 
    term.accumulate_patch( pd, err ); MSQ_ERRRTN(err);
    term.accumulate_inner( pd, objn, arrptr(mGrad), err ); MSQ_ERRRTN(err);
  }
}    

} // namespace Mesquite
