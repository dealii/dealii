/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file PaverMinEdgeLengthWrapper.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "PaverMinEdgeLengthWrapper.hpp"

#include "TerminationCriterion.hpp"
#include "InstructionQueue.hpp"
#include "QualityAssessor.hpp"
#include "MeshImpl.hpp"
#include "PlanarDomain.hpp"

#include "ElementPMeanP.hpp"
#include "PMeanPTemplate.hpp"
#include "TrustRegion.hpp"
#include "TQualityMetric.hpp"
#include "IdealShapeTarget.hpp"
#include "TShapeSizeB1.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "ReferenceMesh.hpp"

#include "EdgeLengthMetric.hpp"
#include "LambdaConstant.hpp"

#include "MsqFPE.hpp"
#include "MeshUtil.hpp"
#include "SimpleStats.hpp"

namespace MESQUITE_NS {

void PaverMinEdgeLengthWrapper::run_wrapper( MeshDomainAssoc* mesh_and_domain,
                                             ParallelMesh* pmesh,
                                             Settings* settings,
                                             QualityAssessor* qa,
                                             MsqError& err )
{
  InstructionQueue q;
  Mesh* mesh = mesh_and_domain->get_mesh();
 
    // calculate average lambda for mesh
  ReferenceMesh ref_mesh( mesh );
  RefMeshTargetCalculator W_0( &ref_mesh );
  SimpleStats lambda_stats;
  MeshUtil tool(mesh, settings);
  tool.lambda_distribution( lambda_stats, err ); MSQ_ERRRTN(err);
  double lambda = lambda_stats.average();
  
    // create objective function
  IdealShapeTarget W_i;
  LambdaConstant W( lambda, &W_i );
  TShapeSizeB1 tm;
  TQualityMetric mu_0( &W, &tm );
  ElementPMeanP mu( 1.0, &mu_0 );
  PMeanPTemplate of( 1.0, &mu );
  
    // create quality assessor
  EdgeLengthMetric len(0.0);
  qa->add_quality_assessment( &mu );
  qa->add_quality_assessment( &len );
  q.add_quality_assessor( qa, err );
  
    // create solver
  TrustRegion solver( &of );
  TerminationCriterion tc, ptc;
  tc.add_absolute_vertex_movement( maxVtxMovement );
  tc.add_iteration_limit( iterationLimit );
  ptc.add_iteration_limit( pmesh ? parallelIterations : 1 );
  solver.set_inner_termination_criterion( &tc );
  solver.set_outer_termination_criterion( &ptc );
  q.set_master_quality_improver( &solver, err ); MSQ_ERRRTN(err);
  q.add_quality_assessor( qa, err );

  // Optimize mesh
  q.run_common( mesh_and_domain, pmesh, settings, err ); MSQ_CHKERR(err);  
}

} // namespace MESQUITE_NS
