/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file UntangleWrapper.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "UntangleWrapper.hpp"
#include "MeshUtil.hpp"
#include "SimpleStats.hpp"
#include "MsqError.hpp"

#include "LambdaConstant.hpp"
#include "IdealShapeTarget.hpp"
#include "TQualityMetric.hpp"
#include "PMeanPTemplate.hpp"
#include "TerminationCriterion.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"
#include "QualityAssessor.hpp"
#include "InstructionQueue.hpp"
#include "ElementPMeanP.hpp"
#include "Instruction.hpp"

#include "TUntangleBeta.hpp"
#include "TUntangleMu.hpp"
#include "TSizeNB1.hpp"
#include "TShapeSize2DNB1.hpp"
#include "TShapeSize3DNB1.hpp"
#include "TMixed.hpp"

#include <memory>

const int NUM_INNER_ITERATIONS = 1;
const double DEFAULT_MOVEMENT_FACTOR = 0.001;
const bool CULLING_DEFAULT = true;
const bool JACOBI_DEFAULT = false;

namespace MESQUITE_NS {

UntangleWrapper::UntangleWrapper() 
  : qualityMetric( SIZE ),
    maxTime(-1),
    movementFactor( DEFAULT_MOVEMENT_FACTOR ),
    metricConstant( -1 ),
    maxIterations(-1),
    doCulling(CULLING_DEFAULT),
    doJacobi(JACOBI_DEFAULT)
{}

UntangleWrapper::UntangleWrapper(UntangleMetric m) 
  : qualityMetric( m ),
    maxTime(-1),
    movementFactor( DEFAULT_MOVEMENT_FACTOR ),
    metricConstant( -1 ),
    maxIterations(-1),
    doCulling(CULLING_DEFAULT),
    doJacobi(JACOBI_DEFAULT)
{}

UntangleWrapper::~UntangleWrapper()
{}

void UntangleWrapper::set_untangle_metric( UntangleMetric metric )
  { qualityMetric = metric; }

void UntangleWrapper::set_metric_constant( double value )
  { metricConstant = value; }

void UntangleWrapper::set_cpu_time_limit( double seconds )
  { maxTime = seconds; }

void UntangleWrapper::set_outer_iteration_limit( int maxIt )
  { maxIterations = maxIt; }

void UntangleWrapper::set_vertex_movement_limit_factor( double f )
  { movementFactor = f; }


void UntangleWrapper::run_wrapper( MeshDomainAssoc* mesh_and_domain,
                                   ParallelMesh* pmesh,
                                   Settings* settings,
                                   QualityAssessor* qa,
                                   MsqError& err )
{
  Instruction::initialize_vertex_byte( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);

  Mesh* mesh = mesh_and_domain->get_mesh();

    // get some global mesh properties
  SimpleStats edge_len, lambda;
  std::auto_ptr<MeshUtil> tool(new MeshUtil( mesh, settings ));
  tool->edge_length_distribution( edge_len, err ); MSQ_ERRRTN(err);
  tool->lambda_distribution( lambda, err ); MSQ_ERRRTN(err);
  tool.reset(0);
  
    // get target metrics from user perferences
  TSizeNB1 mu_size;
  TShapeSize2DNB1 mu_shape_2d;
  TShapeSize3DNB1 mu_shape_3d;
  TMixed mu_shape( &mu_shape_2d, &mu_shape_3d );
  std::auto_ptr<TMetric> mu;
  if (qualityMetric == BETA) {
    double beta = metricConstant;
    if (beta < 0) 
      beta = (lambda.average()*lambda.average())/20;
    //std::cout << "beta= " << beta << std::endl;
    mu.reset(new TUntangleBeta( beta ));
  }
  else {
    TMetric* sub = 0;
    if (qualityMetric == SIZE)
      sub = &mu_size;
    else 
      sub = &mu_shape;
    if (metricConstant >= 0) 
      mu.reset(new TUntangleMu( sub, metricConstant ));
    else 
      mu.reset(new TUntangleMu( sub ));
  }
    
    // define objective function
  IdealShapeTarget base_target;
  LambdaConstant target( lambda.average(), &base_target );
  TQualityMetric metric_0(&target, mu.get());
  ElementPMeanP metric( 1.0, &metric_0 );
  PMeanPTemplate objfunc( 1.0, &metric );
  
    // define termination criterion
  double eps = movementFactor * (edge_len.average() - edge_len.standard_deviation());
  TerminationCriterion inner("<type:untangle_inner>", TerminationCriterion::TYPE_INNER), 
    outer("<type:untangle_outer>", TerminationCriterion::TYPE_OUTER);
  outer.add_untangled_mesh();
  if (doCulling) 
    inner.cull_on_absolute_vertex_movement( eps );
  else
    outer.add_absolute_vertex_movement( eps );
  if (maxTime > 0.0) 
    outer.add_cpu_time( maxTime );
  inner.add_iteration_limit( NUM_INNER_ITERATIONS );
  if (maxIterations > 0)
    outer.add_iteration_limit(maxIterations);
  
    // construct solver
  SteepestDescent solver( &objfunc );
  solver.use_element_on_vertex_patch();
  solver.set_inner_termination_criterion( &inner );
  solver.set_outer_termination_criterion( &outer );
  if (doJacobi)
    solver.do_jacobi_optimization();
  else
    solver.do_gauss_optimization();
  
    // Run 
  qa->add_quality_assessment( &metric );
  InstructionQueue q;
  q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q.set_master_quality_improver( &solver, err ); MSQ_ERRRTN(err);
  q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q.run_common( mesh_and_domain, pmesh, settings, err ); MSQ_ERRRTN(err);
}


} // namespace MESQUITE_NS
