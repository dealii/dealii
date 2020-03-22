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


/** \file QuasiNewton.cpp
 *  \brief Port Todd Munson's quasi-Newton solver to Mesquite
 *  \author Jason Kraftcheck (Mesquite Port)
 */

#include "Mesquite.hpp"
#include "QuasiNewton.hpp"
#include "MsqDebug.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"

namespace MESQUITE_NS {

// Force std::vector to release allocated memory
template <typename T>
static inline void free_vector( std::vector<T>& v )
{
  std::vector<T> temp;
  temp.swap( v );
}

std::string QuasiNewton::get_name() const { return "QuasiNewton"; }

PatchSet* QuasiNewton::get_patch_set()
  { return PatchSetUser::get_patch_set(); }

QuasiNewton::QuasiNewton( ObjectiveFunction* of )
  : VertexMover( of ), PatchSetUser(true), mMemento(0)
{ }

QuasiNewton::~QuasiNewton()
{
  delete mMemento;
  mMemento = 0;
}


void QuasiNewton::initialize( PatchData& pd, MsqError& err )
{
  if (!mMemento) {
    mMemento = pd.create_vertices_memento( err );
    MSQ_CHKERR(err);
  }
}

void QuasiNewton::initialize_mesh_iteration( PatchData& /*pd*/, MsqError& /*err*/ )
{ }

void QuasiNewton::terminate_mesh_iteration( PatchData& /*pd*/, MsqError& /*err*/ )
{ }

void QuasiNewton::cleanup()
{ 
    // release memento
  delete mMemento;
  mMemento = 0;

    // release coordinates
  for (size_t i = 0; i < (sizeof(w)/sizeof(w[0])); ++i)   
    free_vector(w[i]);
    // release gradients
  for (size_t i = 0; i < (sizeof(v)/sizeof(v[0])); ++i) 
    free_vector(v[i]);
  
    // release Hessian memory
  free_vector(mHess);
  
    // release temporary array memory
  free_vector(x);
  free_vector(d);
}

// Do v += s * x, where v and x are arrays of length n
static inline void plus_eq_scaled( Vector3D* v, double s, const Vector3D* x, size_t n )
{
  Vector3D* end = v + n;
  for (; v != end; ++v, ++x)
    *v += s * *x;
}

void QuasiNewton::solve( Vector3D* z_arr, const Vector3D* v_arr ) const
{
  SymMatrix3D pd;
  
  const double small = DBL_EPSILON;
  const size_t nn = mHess.size();
  for (size_t i = 0; i < nn; ++i) {

      // ensure positive definite: perturb a bit if
      // diagonal values are zero.
    SymMatrix3D d = mHess[i];
    while (fabs(d[0]) < small || fabs(d[3]) < small || fabs(d[5]) < small)
      d += small;

      // factor
    pd[0] = 1.0 / d[0];
    pd[1] = d[1] * pd[0];
    pd[2] = d[2] * pd[0];
    
    pd[3] = 1.0 / (d[3] - d[1]*pd[1]);
    pd[5] = d[4] - d[2]*pd[1];
    pd[4] = pd[3] * pd[5];
    pd[5] = 1.0 / (d[5] - d[2]*pd[2] - pd[4]*pd[5]);
    
    if (pd[0] <= 0.0 || pd[3] <= 0.0 || pd[5] <= 0.0) {
      if (d[0] + d[3] + d[5] <= 0) {
          // switch to diagonal
        pd[0] = 1.0 / fabs(d[0]);
        pd[1] = 0.0;
        pd[2] = 0.0;
        pd[3] = 1.0 / fabs(d[3]);
        pd[4] = 0.0;
        pd[5] = 1.0 / fabs(d[5]);
      }
      else {
          // diagonal preconditioner
        pd[0] = pd[3] = pd[5] = 1.0 / (d[0] + d[3] + d[5]);
        pd[1] = pd[2] = pd[4] = 0.0;
      }
    }

      // solve
    const Vector3D& v = v_arr[i];
    Vector3D& z = z_arr[i];
    z[0] = v[0];
    z[1] = v[1] - pd[1]*z[0];
    z[2] = v[2] - pd[2]*z[0] - pd[4]*z[1];

    z[0] *= pd[0];
    z[1] *= pd[3];
    z[2] *= pd[5];

    z[1] -= pd[4]*z[2];
    z[0] -= pd[1]*z[1] + pd[2]*z[2];
  }
}
    

void QuasiNewton::optimize_vertex_positions( PatchData& pd, MsqError& err )
{
  TerminationCriterion& term = *get_inner_termination_criterion();
  OFEvaluator& func = get_objective_function_evaluator();
  
  const double sigma = 1e-4;
  const double beta0 = 0.25;
  const double beta1 = 0.80;
  const double tol1 = 1e-8;
  const double epsilon = 1e-10;

  double norm_r; //, norm_g;
  double alpha, beta;
  double obj, objn;

  size_t i;
  
    // Initialize stuff
  const size_t nn = pd.num_free_vertices();
  double a[QNVEC], b[QNVEC], r[QNVEC];
  for (i = 0; i < QNVEC; ++i)
    r[i] = 0;
  for (i = 0; i <= QNVEC; ++i) {
    v[i].clear();
    v[i].resize( nn, Vector3D(0.0) );
    w[i].clear();
    w[i].resize( nn, Vector3D(0.0) );
  }
  d.resize( nn );
  mHess.resize( nn );  //hMesh(mesh);

  bool valid = func.update( pd, obj, v[QNVEC], mHess, err ); MSQ_ERRRTN(err);
  if (!valid) {
    MSQ_SETERR(err)("Initial objective function is not valid", MsqError::INVALID_MESH);
    return;
  }

  while (!term.terminate()) {
    pd.recreate_vertices_memento( mMemento, err ); MSQ_ERRRTN(err);
    pd.get_free_vertex_coordinates( w[QNVEC] );

    x = v[QNVEC];
    for (i = QNVEC; i--; ) {
      a[i] = r[i] * inner( &(w[i][0]), arrptr(x), nn );
      plus_eq_scaled( arrptr(x), -a[i], &v[i][0], nn );
    }
     
    solve( arrptr(d), arrptr(x) );
  
    for (i = QNVEC; i--; ) {
      b[i] = r[i] * inner( &(v[i][0]), arrptr(d), nn );
      plus_eq_scaled( arrptr(d), a[i]-b[i], &(w[i][0]), nn );
    }
    
    alpha = -inner( &(v[QNVEC][0]), arrptr(d), nn );  /* direction is negated */
    if (alpha > 0.0) {
      MSQ_SETERR(err)("No descent.", MsqError::INVALID_MESH);
      return;
    }
   
    alpha *= sigma;
    beta = 1.0;
    
    pd.move_free_vertices_constrained( arrptr(d), nn, -beta, err ); MSQ_ERRRTN(err);
    valid = func.evaluate( pd, objn, v[QNVEC], err ); 
    if (err.error_code() == err.BARRIER_VIOLATED)             
      err.clear();  // barrier violated does not represent an actual error here
    MSQ_ERRRTN(err);
    if (!valid ||
        (obj - objn < -alpha*beta - epsilon &&
         length( &(v[QNVEC][0]), nn ) >= tol1)) {
      
      if (!valid)  // function not defined at trial point
        beta *= beta0;
      else  // unacceptable iterate
        beta *= beta1;
      
      for (;;) {
        if (beta < tol1) {
          pd.set_to_vertices_memento( mMemento, err ); MSQ_ERRRTN(err);
          MSQ_SETERR(err)("Newton step not good", MsqError::INTERNAL_ERROR);
          return;
        }
      
        pd.set_free_vertices_constrained( mMemento, arrptr(d), nn, -beta, err ); MSQ_ERRRTN(err);
        valid = func.evaluate( pd, objn, err );
        if (err.error_code() == err.BARRIER_VIOLATED)             
          err.clear();  // barrier violated does not represent an actual error here
        MSQ_ERRRTN(err);
        if (!valid) // function undefined at trial point
          beta *= beta0;
        else if (obj - objn < -alpha*beta - epsilon) // unacceptlable iterate
          beta *= beta1;
        else
          break;
      }
    }
    
    for (i = 0; i < QNVEC-1; ++i) {
      r[i] = r[i+1];
      w[i].swap( w[i+1] );
      v[i].swap( v[i+1] );
    }
    w[QNVEC-1].swap( w[0] );
    v[QNVEC-1].swap( v[0] );
    
    func.update( pd, obj, v[QNVEC], mHess, err ); MSQ_ERRRTN(err);
    norm_r = length_squared( &(v[QNVEC][0]), nn );
    //norm_g = sqrt(norm_r);

    // checks stopping criterion 
    term.accumulate_patch( pd, err ); MSQ_ERRRTN(err);
    term.accumulate_inner( pd, objn, &v[QNVEC][0], err ); MSQ_ERRRTN(err);
  }
}

} // namespace Mesquite
