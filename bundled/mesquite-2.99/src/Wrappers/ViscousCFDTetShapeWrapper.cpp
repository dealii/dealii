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


/** \file ViscousCFDTetShapeWrapper.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ViscousCFDTetShapeWrapper.hpp"

#include "QualityAssessor.hpp"
#include "InstructionQueue.hpp"
#include "TagVertexMesh.hpp"
#include "TrustRegion.hpp"
#include "TerminationCriterion.hpp"

#include "PMeanPTemplate.hpp"
#include "TQualityMetric.hpp"
#include "AddQualityMetric.hpp"

#include "TShapeB1.hpp"
#include "TShapeNB1.hpp"
#include "TShapeSizeOrientB1.hpp"
#include "TShapeSizeOrientNB1.hpp"

#include "IdealShapeTarget.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "ReferenceMesh.hpp"
#include "TetDihedralWeight.hpp"
#include "RemainingWeight.hpp"

namespace MESQUITE_NS {

void ViscousCFDTetShapeWrapper::run_wrapper( MeshDomainAssoc* mesh_and_domain,
                                             ParallelMesh* pmesh,
                                             Settings* settings,
                                             QualityAssessor* qa,
                                             MsqError& err )
{
  InstructionQueue q;
  
  // Set up barrier metric to see if mesh contains inverted elements
  TShapeB1 mu_b;
  IdealShapeTarget w_ideal;
  TQualityMetric barrier( &w_ideal, &mu_b );
  
  // Check for inverted elements in the mesh
  QualityAssessor inv_check( &barrier );
  inv_check.disable_printing_results();
  q.add_quality_assessor( &inv_check, err );  MSQ_ERRRTN(err);
  q.run_common( mesh_and_domain, pmesh, settings, err ); MSQ_ERRRTN(err);
  q.remove_quality_assessor( 0, err ); MSQ_ERRRTN(err);
  const QualityAssessor::Assessor* inv_b = inv_check.get_results( &barrier );
  const bool use_barrier = (0 == inv_b->get_invalid_element_count());
  
  // Create remaining metric instances
  TShapeNB1 mu;
  TShapeSizeOrientNB1 mu_o;
  TShapeSizeOrientB1 mu_ob;
  
  // Select which target metrics to use
  TMetric *mu_p, *mu_op;
  if (use_barrier) {
    mu_p = &mu_b;
    mu_op = &mu_ob;
  }
  else {
    mu_p = &mu;
    mu_op = &mu_o;
  }
  

  // Set up target and weight calculators
  Mesh* mesh = mesh_and_domain->get_mesh();
  TagVertexMesh init_mesh( err, pmesh ? (Mesh*)pmesh : mesh );  MSQ_ERRRTN(err);
  ReferenceMesh ref_mesh( &init_mesh );
  RefMeshTargetCalculator w_init( &ref_mesh );
  TetDihedralWeight c_dihedral( &ref_mesh, dCutoff, aVal );
  RemainingWeight c_remaining( &c_dihedral );
  
  // Create objective function
  TQualityMetric metric1( &w_ideal, &c_dihedral,  mu_p  );
  TQualityMetric metric2( &w_init,  &c_remaining, mu_op );
  AddQualityMetric of_metric( &metric1, &metric2, err );  MSQ_ERRRTN(err);
  PMeanPTemplate obj_func( 1.0, &of_metric );
  
  // Create optimizer
  TrustRegion solver( &obj_func );
  TerminationCriterion term, ptc;
  term.add_iteration_limit( iterationLimit );
  term.add_absolute_vertex_movement( maxVtxMovement );
  ptc.add_iteration_limit( pmesh ? parallelIterations : 1 );
  solver.set_inner_termination_criterion( &term );
  solver.set_outer_termination_criterion( &ptc );
  
  // Create instruction queue
  qa->add_quality_assessment( &metric1 );
  qa->add_quality_assessment( &metric2 );
  qa->add_quality_assessment( &of_metric );
  q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  q.set_master_quality_improver( &solver, err ); MSQ_ERRRTN(err);
  q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);

  // Optimize mesh
  q.run_common( mesh_and_domain, pmesh, settings, err ); MSQ_CHKERR(err);  
}

} // namespace MESQUITE_NS
