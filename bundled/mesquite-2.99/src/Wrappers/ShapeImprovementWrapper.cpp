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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ShapeImprovementWrapper.cpp

Member functions of the Mesquite::ShapeImprovementWrapper class

  \author Michael Brewer
  \date   June 6, 2003
 */

#include "InstructionQueue.hpp"
#include "ShapeImprovementWrapper.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "ConjugateGradient.hpp"
#include "TerminationCriterion.hpp"
#include "IdealWeightInverseMeanRatio.hpp"
#include "FeasibleNewton.hpp"
#include "InstructionQueue.hpp"
#include "QualityAssessor.hpp"

namespace MESQUITE_NS {

const double DEF_UNT_BETA = 1e-8;
const double DEF_SUC_EPS = 1e-4;

/*! The consturctor allows for two values.  The first is a 
  time bound (in seconds) used as a termination criterion.  If
  this value is non-positive, no time bound will be set.
  By default, the value is set to zero and no time bound
  is used.  The second value is the tolerance for the gradient
  norm termination criteria.  The default value is 1.e-6.*/
ShapeImprovementWrapper::ShapeImprovementWrapper(MsqError& ,
                                                 double cpu_time,
                                                 double grad_norm,
                                                 int parallel_iterations) 
 : maxTime(cpu_time), 
   gradNorm(grad_norm),
   untBeta(DEF_UNT_BETA),
   successiveEps(DEF_SUC_EPS),
   parallelIterations(parallel_iterations)
{}

ShapeImprovementWrapper::ShapeImprovementWrapper(double cpu_time,
                                                 double grad_norm,
                                                 int parallel_iterations) 
 : maxTime(cpu_time), 
   gradNorm(grad_norm),
   untBeta(DEF_UNT_BETA),
   successiveEps(DEF_SUC_EPS),
   parallelIterations(parallel_iterations)
{}
 
 
void ShapeImprovementWrapper::run_wrapper( MeshDomainAssoc* mesh_and_domain,
                                           ParallelMesh* pmesh,
                                           Settings* settings,
                                           QualityAssessor* qa,
                                           MsqError& err )
{
    // Define an untangler
  UntangleBetaQualityMetric untangle_metric( untBeta );
  LPtoPTemplate untangle_func( 2, &untangle_metric );
  ConjugateGradient untangle_global( &untangle_func );
  TerminationCriterion untangle_inner, untangle_outer;
  untangle_global.use_global_patch();
  untangle_inner.add_absolute_quality_improvement( 0.0 );
  untangle_inner.add_absolute_successive_improvement( successiveEps );
  untangle_outer.add_iteration_limit( 1 );
  untangle_global.set_inner_termination_criterion( &untangle_inner );
  untangle_global.set_outer_termination_criterion( &untangle_outer );

    // define shape improver
  IdealWeightInverseMeanRatio inverse_mean_ratio;
  inverse_mean_ratio.set_averaging_method( QualityMetric::LINEAR );
  LPtoPTemplate obj_func( 2, &inverse_mean_ratio );
  FeasibleNewton feas_newt( &obj_func );
  TerminationCriterion term_inner, term_outer;
  feas_newt.use_global_patch();
  qa->add_quality_assessment( &inverse_mean_ratio );
  term_inner.add_absolute_gradient_L2_norm( gradNorm );
  term_inner.add_relative_successive_improvement( successiveEps );
  term_outer.add_iteration_limit( pmesh ? parallelIterations : 1 );
  feas_newt.set_inner_termination_criterion( &term_inner );
  feas_newt.set_outer_termination_criterion( &term_outer );

    // Apply CPU time limit to untangler
  if (maxTime > 0.0)
    untangle_inner.add_cpu_time( maxTime );
  
    // Run untangler
  InstructionQueue q1;
  Timer totalTimer;
  q1.set_master_quality_improver( &untangle_global, err ); MSQ_ERRRTN(err);
  q1.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q1.run_common( mesh_and_domain, pmesh, settings, err ); MSQ_ERRRTN(err);
  
    // If limited by CPU time, limit next step to remaning time
  if (maxTime > 0.0) {
    double remaining = maxTime - totalTimer.since_birth();
    if (remaining <= 0.0 ){
      MSQ_DBGOUT(2) << "Optimization is terminating without perfoming shape improvement." << std::endl;
      remaining = 0.0;
    }
    term_inner.add_cpu_time( remaining );
  }
  
    // Run shape improver
  InstructionQueue q2;
  q2.set_master_quality_improver( &feas_newt, err ); MSQ_ERRRTN(err);
  q2.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q2.run_common( mesh_and_domain, pmesh, settings, err ); MSQ_ERRRTN(err);
}

} // namespace Mesquite

  
