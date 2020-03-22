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
/*!
  \file   VertexMover.cpp
  \brief  

  The VertexMover Class is the base class for all the smoothing algorythms 

  \author Thomas Leurent
  \date   2002-01-17
*/


#include "VertexMover.hpp"
#include "NonGradient.hpp"
#include "QualityMetric.hpp"
#include "TQualityMetric.hpp"
#include "TMetricBarrier.hpp"
#include "AWMetricBarrier.hpp"
#include "ElementMaxQM.hpp"
#include "ElementAvgQM.hpp"
#include "ElementPMeanP.hpp"
#include "PlanarDomain.hpp"
#include "ObjectiveFunctionTemplate.hpp"
#include "MaxTemplate.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"
#include "PatchSet.hpp"
#include "PatchData.hpp"
#include "ParallelHelperInterface.hpp"
#include <algorithm>

namespace MESQUITE_NS {

extern int get_parallel_rank();
extern int get_parallel_size();
  extern void parallel_barrier();

VertexMover::VertexMover( ObjectiveFunction* OF ) 
  : QualityImprover(),
    objFuncEval( OF ) ,
    jacobiOpt(false)
  {}


VertexMover::~VertexMover() {}

/*
  
    +-----------+
    |Reset Outer|
    |Criterion  |
    +-----------+
          |
          V
          +
        /   \
       /Outer\  YES
+--> <Criterion>-----> DONE
|      \Done?/
|       \   /
|         + 
|         |NO
|         V
|   +-----------+
1   |Reset Mesh |
|   | Iteration |
|   +-----------+
|         |
|         V
|         +
|       /  \
|   NO /Next\
+-----<Patch > <-----+
       \    /        |
        \  /         |
          +          |
       YES|          |
          V          |
    +-----------+    |
    |Reset Inner|    |
    |Criterion  |    2
    +-----------+    |
          |          |
          V          |  
          +          |               
        /   \        |
       /Inner\  YES  |
+--> <Criterion>-----+    --------------+
|      \Done?/                          |
|       \   /                           |
|         +                             |
|         |NO                           |
|         V                          Inside
3   +-----------+                    Smoother
|   |   Smooth  |                       |
|   |   Patch   |                       |
|   +-----------+                       |
|         |                             |
----------+               --------------+
                      
*/        

/*! \brief Improves the quality of the MeshSet, calling some
    methods specified in a class derived from VertexMover

    \param const MeshSet &: this MeshSet is looped over. Only the
    mutable data members are changed (such as currentVertexInd).
  */
double VertexMover::loop_over_mesh( MeshDomainAssoc* mesh_and_domain,
                                    const Settings* settings,
                                    MsqError& err )
{
  Mesh* mesh = mesh_and_domain->get_mesh();
  MeshDomain* domain = mesh_and_domain->get_domain();

  TagHandle coord_tag = 0; // store uncommitted coords for jacobi optimization 
  TagHandle* coord_tag_ptr = 0;

  TerminationCriterion* outer_crit = 0;
  TerminationCriterion* inner_crit = 0;

    // Clear culling flag, set hard fixed flag, etc on all vertices
  initialize_vertex_byte( mesh_and_domain, settings, err ); MSQ_ERRZERO(err);

    // Get the patch data to use for the first iteration
  OFEvaluator& obj_func = get_objective_function_evaluator();
 
    // Check for impromer use of MaxTemplate
  ObjectiveFunctionTemplate* maxt_ptr = dynamic_cast<MaxTemplate*>(obj_func.get_objective_function());
  if (maxt_ptr)
  {
  QualityImprover* ngqi_ptr = dynamic_cast<NonGradient*>(this);
  if (!ngqi_ptr)
      std::cout << "Warning: MaxTemplate results in non-differentiable objective function." << std::endl <<
                   "   Therefore, it is best to use the NonGradient solver. Other Mesquite" << std::endl <<
                   "   solvers require derivative information." << std::endl;
  }

  PatchData patch;
  patch.set_mesh( mesh );
  patch.set_domain( domain );
  if (settings)
    patch.attach_settings( settings );
  bool one_patch = false, inner_crit_terminated, all_culled;
  std::vector<Mesh::VertexHandle> patch_vertices;
  std::vector<Mesh::ElementHandle> patch_elements;
  bool valid;
  
  PatchSet* patch_set = get_patch_set();
  if (!patch_set) {
    MSQ_SETERR(err)("No PatchSet for QualityImprover!", MsqError::INVALID_STATE);
    return 0.0;
  }
  patch_set->set_mesh( mesh );
  
  std::vector<PatchSet::PatchHandle> patch_list;
  patch_set->get_patch_handles( patch_list, err ); MSQ_ERRZERO(err);

    // check for inverted elements when using a barrietr target metric
  TQualityMetric* tqm_ptr = NULL;
  TMetricBarrier* tm_ptr = NULL;
  AWMetricBarrier* awm_ptr = NULL;
  ElemSampleQM* sample_qm_ptr = NULL;
  QualityMetric* qm_ptr = NULL;
  TQualityMetric* pmeanp_ptr = NULL;
  ElementMaxQM* elem_max_ptr = NULL;
  ElementAvgQM* elem_avg_ptr = NULL;
  ElementPMeanP* elem_pmeanp_ptr = NULL;

  ObjectiveFunctionTemplate* of_ptr =  dynamic_cast<ObjectiveFunctionTemplate*>(obj_func.get_objective_function() );
  if (of_ptr)
    qm_ptr = of_ptr->get_quality_metric();
  if (qm_ptr)
  {
    pmeanp_ptr = dynamic_cast<TQualityMetric*>(qm_ptr);  // PMeanP case
    elem_max_ptr = dynamic_cast<ElementMaxQM*>(qm_ptr);
    elem_avg_ptr = dynamic_cast<ElementAvgQM*>(qm_ptr);
    elem_pmeanp_ptr = dynamic_cast<ElementPMeanP*>(qm_ptr);
  }
  if (elem_max_ptr)
    sample_qm_ptr = elem_max_ptr->get_quality_metric();
  else if (elem_pmeanp_ptr)
    sample_qm_ptr = elem_pmeanp_ptr->get_quality_metric();
  else if (elem_avg_ptr)
    sample_qm_ptr = elem_avg_ptr->get_quality_metric();
  else if (pmeanp_ptr)
  {
    tm_ptr =  dynamic_cast<TMetricBarrier*>(pmeanp_ptr->get_target_metric());
    awm_ptr =  dynamic_cast<AWMetricBarrier*>(pmeanp_ptr->get_target_metric());
  }
   
  if (sample_qm_ptr || pmeanp_ptr)
  { 
    if (!pmeanp_ptr)
    {
      tqm_ptr = dynamic_cast<TQualityMetric*>(sample_qm_ptr);
      if (tqm_ptr)
      {
        tm_ptr =  dynamic_cast<TMetricBarrier*>(tqm_ptr->get_target_metric());
        awm_ptr =  dynamic_cast<AWMetricBarrier*>(tqm_ptr->get_target_metric());
      }
    }
    else
    {
      tqm_ptr = dynamic_cast<TQualityMetric*>(pmeanp_ptr);
    }
    
    if (tqm_ptr && (tm_ptr || awm_ptr))
    {
        // check for inverted elements
      this->initialize(patch, err); 
      std::vector<size_t> handles;
 
        // set up patch data      
      std::vector<PatchSet::PatchHandle>::iterator patch_iter = patch_list.begin();
      while( patch_iter != patch_list.end() )
      {
        do 
          {
            patch_set->get_patch( *patch_iter, patch_elements, patch_vertices, err );
            if (MSQ_CHKERR(err)) goto ERROR;
            ++patch_iter;
          } while (patch_elements.empty() && patch_iter != patch_list.end()) ;

        patch.set_mesh_entities( patch_elements, patch_vertices, err );
        if (MSQ_CHKERR(err)) goto ERROR;
    
        qm_ptr->get_evaluations( patch, handles, true, err ); // MSQ_ERRFALSE(err);
  
          // do actual check for inverted elements
        std::vector<size_t>::const_iterator i;
        double tvalue;
        for (i = handles.begin(); i != handles.end(); ++i)
        {
          bool result = tqm_ptr->evaluate( patch, *i, tvalue, err );
          if (MSQ_CHKERR(err) || !result)
            return false;    // inverted element detected
        }
      }
    }
  }

    // Get termination criteria
  outer_crit=this->get_outer_termination_criterion();
  inner_crit=this->get_inner_termination_criterion();
  if(outer_crit == 0){
    MSQ_SETERR(err)("Termination Criterion pointer is Null", MsqError::INVALID_STATE);
    return 0.;
  }
  if(inner_crit == 0){
    MSQ_SETERR(err)("Termination Criterion pointer for inner loop is Null", MsqError::INVALID_STATE);
    return 0.;
  }

    // Set Termination Criterion defaults if no other Criterion is set
  if ( !outer_crit->criterion_is_set() )
    outer_crit->add_iteration_limit(1);
  if ( !inner_crit->criterion_is_set() )
    inner_crit->add_iteration_limit(10);

    // If using a local patch, suppress output of inner termination criterion
  if (patch_list.size() > 1) 
    inner_crit->set_debug_output_level(3);
  else
    one_patch = true;
  
  if (jacobiOpt) {
    coord_tag = get_jacobi_coord_tag(mesh, err);
    MSQ_ERRZERO(err);
    coord_tag_ptr = &coord_tag;
  }
  
    // Initialize outer loop
    
  this->initialize(patch, err);        
  if (MSQ_CHKERR(err)) goto ERROR;
  
  valid = obj_func.initialize( mesh_and_domain, settings, patch_set, err ); 
  if (MSQ_CHKERR(err)) goto ERROR;
  if (!valid) {
    MSQ_SETERR(err)("ObjectiveFunction initialization failed.  Mesh "
                    "invalid at one or more sample points.", 
                    MsqError::INVALID_MESH);
    goto ERROR;
  }
  
  outer_crit->reset_outer(mesh, domain, obj_func, settings, err); 
  if (MSQ_CHKERR(err)) goto ERROR;
  
 
    // if only one patch, get the patch now
  if (one_patch) {
    patch_set->get_patch( patch_list[0], patch_elements, patch_vertices, err );
    if (MSQ_CHKERR(err)) goto ERROR;
    patch.set_mesh_entities( patch_elements, patch_vertices, err );
    if (MSQ_CHKERR(err)) goto ERROR;
  }
  
   // Loop until outer termination criterion is met
  inner_crit_terminated = false;
  while (!outer_crit->terminate())
  {
    if (inner_crit_terminated) {
      MSQ_SETERR(err)("Inner termiation criterion satisfied for all patches "
                      "without meeting outer termination criterion.  This is "
                      "an infinite loop.  Aborting.", MsqError::INVALID_STATE);
      break;
    }
    inner_crit_terminated = true;
    all_culled = true;

        int num_patches=0;
      // Loop over each patch
    std::vector<PatchSet::PatchHandle>::iterator p_iter = patch_list.begin();
    while( p_iter != patch_list.end() )
    {
      if (!one_patch) { // if only one patch (global) re-use the previous one
          // loop until we get a non-empty patch.  patch will be empty
          // for culled vertices with element-on-vertex patches
        do {
          patch_set->get_patch( *p_iter, patch_elements, patch_vertices, err );
          if (MSQ_CHKERR(err)) goto ERROR;
          ++p_iter;
        } while (patch_elements.empty() && p_iter != patch_list.end()) ;
        
        if (patch_elements.empty()) { // no more non-culled vertices
                if (0)
                  std::cout << "P[" << get_parallel_rank() << "] tmp srk all vertices culled."  << std::endl;
          break;
        }
      
        all_culled = false;
        patch.set_mesh_entities( patch_elements, patch_vertices, err );
        if (MSQ_CHKERR(err)) goto ERROR;
      } else {
        ++p_iter;
        all_culled = false;
      }
        
            ++num_patches;

        // Initialize for inner iteration
        
      this->initialize_mesh_iteration(patch, err);
      if (MSQ_CHKERR(err)) goto ERROR;
      
      obj_func.reset();
      
      outer_crit->reset_patch( patch, err );
      if (MSQ_CHKERR(err)) goto ERROR;
      
      inner_crit->reset_inner( patch, obj_func, err );
      if (MSQ_CHKERR(err)) goto ERROR;
      
      inner_crit->reset_patch( patch, err );
      if (MSQ_CHKERR(err)) goto ERROR;
      
        // Don't even call optimizer if inner termination 
        // criterion has already been met.
      if (!inner_crit->terminate())
      {
        inner_crit_terminated = false;

          // Call optimizer - should loop on inner_crit->terminate()
        this->optimize_vertex_positions( patch, err );
        if (MSQ_CHKERR(err)) goto ERROR;
      
          // Update for changes during inner iteration 
          // (during optimizer loop)
        
        outer_crit->accumulate_patch( patch, err );
        if (MSQ_CHKERR(err)) goto ERROR;
        
        inner_crit->cull_vertices( patch, obj_func, err );
        if (MSQ_CHKERR(err)) goto ERROR;

                // FIXME
                if (0)
                  {
                    inner_crit->cull_vertices_global (patch, 
                                                      mesh, domain, settings,
                                                      obj_func, err);
                    if (MSQ_CHKERR(err)) goto ERROR;
                  }
        
        patch.update_mesh( err, coord_tag_ptr );
        if (MSQ_CHKERR(err)) goto ERROR;
      }
    } 

    if (jacobiOpt)
      commit_jacobi_coords( coord_tag, mesh, err );

    this->terminate_mesh_iteration(patch, err); 
    if (MSQ_CHKERR(err)) goto ERROR;
    
    outer_crit->accumulate_outer( mesh, domain, obj_func, settings, err );
    if (MSQ_CHKERR(err)) goto ERROR;
    
    if (all_culled)
      break;
  }


ERROR:  
  if (jacobiOpt)
    mesh->tag_delete( coord_tag, err );

    //call the criteria's cleanup funtions.
  if (outer_crit) 
    outer_crit->cleanup(mesh,domain,err);
  if (inner_crit)
  inner_crit->cleanup(mesh,domain,err);
    //call the optimization cleanup function.
  this->cleanup();

  return 0.;
}
  

static void checkpoint_bytes( Mesh* mesh, std::vector<unsigned char>& saved_bytes, MsqError& err)
{
  std::vector<Mesh::VertexHandle> vertexHandlesArray;
  mesh->get_all_vertices(vertexHandlesArray, err); MSQ_ERRRTN(err);
  saved_bytes.resize(vertexHandlesArray.size());
  mesh->vertices_get_byte( &vertexHandlesArray[0],
                           &saved_bytes[0],
                           vertexHandlesArray.size(),
                           err ); MSQ_ERRRTN(err);
}

static void restore_bytes( Mesh* mesh, std::vector<unsigned char>& saved_bytes, MsqError& err)
{
  std::vector<Mesh::VertexHandle> vertexHandlesArray;
  mesh->get_all_vertices(vertexHandlesArray, err); MSQ_ERRRTN(err);
  mesh->vertices_set_byte( &vertexHandlesArray[0],
                           &saved_bytes[0],
                           vertexHandlesArray.size(),
                           err ); MSQ_ERRRTN(err);
}

  static void save_or_restore_debug_state(bool save)
  {
    static bool debug[3] = {false,false,false};
    if (save) 
      {
        debug[0] = MsqDebug::get(1);
        debug[1] = MsqDebug::get(2);
        debug[2] = MsqDebug::get(3);
      }
    else
      {
        if (debug[0]) MsqDebug::enable(1);
        if (debug[1]) MsqDebug::enable(2);
        if (debug[2]) MsqDebug::enable(3);
      }
  }

/*! \brief Improves the quality of the MeshSet, calling some
    methods specified in a class derived from VertexMover

    \param const MeshSet &: this MeshSet is looped over. Only the
    mutable data members are changed (such as currentVertexInd).
  */
double VertexMover::loop_over_mesh( ParallelMesh* mesh,
                                    MeshDomain* domain,
                                    const Settings* settings,
                                    MsqError& err )
{
  std::vector<size_t> junk;
  Mesh::VertexHandle vertex_handle;
  TagHandle coord_tag = 0; // store uncommitted coords for jacobi optimization 
  TagHandle* coord_tag_ptr = 0;
  int outer_iter=0;
  int inner_iter=0;
    bool one_patch = false;

    // Clear culling flag, set hard fixed flag, etc on all vertices
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc((Mesh*)mesh, 0);
  initialize_vertex_byte( &mesh_and_domain, settings, err ); MSQ_ERRZERO(err);

    // Get the patch data to use for the first iteration
  OFEvaluator& obj_func = get_objective_function_evaluator();
  
  PatchData patch;
  patch.set_mesh( (Mesh*) mesh );
  patch.set_domain( domain );
  patch.attach_settings( settings );

  ParallelHelper* helper = mesh->get_parallel_helper();
  if (!helper) {
    MSQ_SETERR(err)("No ParallelHelper instance", MsqError::INVALID_STATE);
    return 0;
  }

  helper->smoothing_init(err);  MSQ_ERRZERO(err);

  bool inner_crit_terminated, all_culled;
  std::vector<Mesh::VertexHandle> patch_vertices;
  std::vector<Mesh::ElementHandle> patch_elements;
  std::vector<Mesh::VertexHandle> fixed_vertices;
  std::vector<Mesh::VertexHandle> free_vertices;
   
    // Get termination criteria
  TerminationCriterion* outer_crit=this->get_outer_termination_criterion();
  TerminationCriterion* inner_crit=this->get_inner_termination_criterion();
  if(outer_crit == 0){
    MSQ_SETERR(err)("Termination Criterion pointer is Null", MsqError::INVALID_STATE);
    return 0.;
  }
  if(inner_crit == 0){
    MSQ_SETERR(err)("Termination Criterion pointer for inner loop is Null", MsqError::INVALID_STATE);
    return 0.;
  }
  
  PatchSet* patch_set = get_patch_set();
  if (!patch_set) {
    MSQ_SETERR(err)("No PatchSet for QualityImprover!", MsqError::INVALID_STATE);
    return 0.0;
  }
  patch_set->set_mesh( (Mesh*)mesh );

  std::vector<PatchSet::PatchHandle> patch_list;
  patch_set->get_patch_handles( patch_list, err ); MSQ_ERRZERO(err);
  
    if (patch_list.size() > 1) 
      inner_crit->set_debug_output_level(3);
    else
      one_patch = true;

  if (jacobiOpt) {
    coord_tag = get_jacobi_coord_tag(mesh, err);
    MSQ_ERRZERO(err);
    coord_tag_ptr = &coord_tag;
  }
  
  // parallel error checking
  MsqError perr;
  bool pdone=false;

#define PERROR_COND continue

    // Initialize outer loop
    
  this->initialize(patch, err); 
  if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("initialize patch", MsqError::INVALID_STATE); } //goto ERROR;
  
  MeshDomainAssoc mesh_and_domain2 = MeshDomainAssoc((Mesh*)mesh, domain);
  obj_func.initialize( &mesh_and_domain2, settings, patch_set, err ); 
  if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("initialize obj_func", MsqError::INVALID_STATE);} //goto ERROR;
  
  outer_crit->reset_outer( (Mesh*)mesh, domain, obj_func, settings, err); 
  if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("reset_outer", MsqError::INVALID_STATE);} //goto ERROR;
   
   // Loop until outer termination criterion is met
  inner_crit_terminated = false;
  all_culled = false;
  for (;;)
  {
    if (0)
      std::cout << "P[" << get_parallel_rank() << "] tmp srk inner_iter= " << inner_iter << " outer_iter= " << outer_iter << std::endl;

    //PERROR:

    if (MSQ_CHKERR(perr))
      {
        std::cout << "P[" << get_parallel_rank() << "] VertexMover::loop_over_mesh found parallel error: " << perr << "\n quitting... pdone= " << pdone << std::endl;
        pdone = true;
      }

    helper->communicate_any_true (pdone, err);

    if (0)
      std::cout << "P[" << get_parallel_rank() << "] tmp srk inner_iter= " << inner_iter << " outer_iter= " << outer_iter << " pdone= " << pdone << std::endl;

    if (pdone)
      {
        std::cout << "P[" << get_parallel_rank() << "] VertexMover::loop_over_mesh found parallel error, quitting... pdone= " << pdone << std::endl;
        MSQ_SETERR(err)("PARALLEL ERROR", MsqError::PARALLEL_ERROR);
        break;
      }

    ++outer_iter;

    /// srkenno@sandia.gov 1/19/12: the logic here was changed so that all proc's must agree
    ///   on the values used for outer and inner termination before the iteration is stopped.
    ///   Previously, the ParallelHelper::communicate_all_true method returned true if any
    ///   proc sent it a true value, which seems to be a bug, at least in the name of the method.
    ///   The method has been changed to return true only if all proc's values are true.  
    ///   In the previous version, this meant that if any proc hit its inner or outer
    ///   termination criterion, the loop was exited, and thus some parts of the mesh
    ///   are potentially left unconverged.  Also, if the outer criterion was satisfied on
    ///   part of the mesh (say a uniform part), the iterations were not executed at all.
    /// Also, changed name of "did_some" to "inner_crit_terminated", and flipped its boolean 
    ///   value to be consistent with the name - for readability and for correctness since
    ///   we want to communicate a true value through the helper.

    bool outer_crit_terminated = outer_crit->terminate();
    bool outer_crit_terminated_local = outer_crit_terminated;
    helper->communicate_all_true( outer_crit_terminated, err ); 

    bool inner_crit_terminated_local = inner_crit_terminated;
    helper->communicate_all_true( inner_crit_terminated, err ); 

    bool all_culled_local = all_culled;
    helper->communicate_all_true( all_culled, err ); 

    bool done = all_culled || outer_crit_terminated;
    if (inner_crit_terminated) {
      MSQ_SETERR(err)("Inner termination criterion satisfied for all patches "
                      "without meeting outer termination criterion.  This is "
                      "an infinite loop.  Aborting.", MsqError::INVALID_STATE);
      done = true;
      helper->communicate_any_true( done, err ); 
    }

    bool local_done=done;

    helper->communicate_all_true( done, err ); 

    if (0)
      std::cout << "P[" << get_parallel_rank() << "] tmp srk done= " << done << " local_done= " << local_done 
                << " all_culled= " << all_culled 
                << " outer_crit->terminate()= " << outer_crit->terminate()
                << " outer_term= " << outer_crit_terminated
                << " outer_term_local= " << outer_crit_terminated_local
                << " inner_crit_terminated = " << inner_crit_terminated
                << " inner_crit_terminated_local = " << inner_crit_terminated_local
                << " all_culled = " << all_culled
                << " all_culled_local = " << all_culled_local
                << std::endl;

    if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("loop start", MsqError::INVALID_STATE);} //goto ERROR;


    if (done)
      break;
    
    
    inner_crit_terminated = true;
    all_culled = true;

    ///*** smooth the interior ***////

    // get the fixed vertices (i.e. the ones *not* part of the first independent set)
    helper->compute_first_independent_set(fixed_vertices); 

    // sort the fixed vertices
    std::sort(fixed_vertices.begin(), fixed_vertices.end());

      // Loop over each patch
        if (0 && MSQ_DBG(2))
          std::cout << "P[" << get_parallel_rank() << "] tmp srk number of patches = " << patch_list.size() 
                    << " inner_iter= " << inner_iter << " outer_iter= " << outer_iter 
                    << " inner.globalInvertedCount = " << inner_crit->globalInvertedCount
                    << " outer.globalInvertedCount = " << outer_crit->globalInvertedCount
                    << " inner.patchInvertedCount = " << inner_crit->patchInvertedCount
                    << " outer.patchInvertedCount = " << outer_crit->patchInvertedCount
                    << std::endl;

        save_or_restore_debug_state(true);
        //MsqDebug::disable_all();
      
        int num_patches=0;

    std::vector<PatchSet::PatchHandle>::iterator p_iter = patch_list.begin();
    while( p_iter != patch_list.end() )
    {


      // loop until we get a non-empty patch.  patch will be empty
      // for culled vertices with element-on-vertex patches
      do {
	patch_set->get_patch( *p_iter, patch_elements, patch_vertices, err );
	if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("get_patch", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
	++p_iter;
      } while (patch_elements.empty() && p_iter != patch_list.end()) ;
        
      if (patch_elements.empty()) { // no more non-culled vertices
              if (0)
                std::cout << "P[" << get_parallel_rank() << "] tmp srk all vertices culled."  << std::endl;
	break;
      }

      if (patch_vertices.empty()) // global patch hack (means all mesh vertices)
      {
	mesh->get_all_vertices(patch_vertices, err);
      }

      free_vertices.clear();

      for (size_t i = 0; i < patch_vertices.size(); ++i) 
	if (!std::binary_search(fixed_vertices.begin(), fixed_vertices.end(), patch_vertices[i]))
	  free_vertices.push_back(patch_vertices[i]);

      if (free_vertices.empty()) { // all vertices were fixed -> skip patch
	continue;
      }

            ++num_patches;
      all_culled = false;
      patch.set_mesh_entities( patch_elements, free_vertices, err );
      if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("set_mesh_entities", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
        
        // Initialize for inner iteration
        
      this->initialize_mesh_iteration(patch, err);
      if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("initialize_mesh_iteration", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
      
      obj_func.reset();
      
      outer_crit->reset_patch( patch, err );
      if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("reset_patch outer", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
      
      inner_crit->reset_inner( patch, obj_func, err );
      if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("reset_inner", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
      
      inner_crit->reset_patch( patch, err );
      if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("inner reset_patch", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
      
        // Don't even call optimizer if inner termination 
        // criterion has already been met.
      if (!inner_crit->terminate())
      {
        inner_crit_terminated = false;
        if (one_patch) ++inner_iter;

          // Call optimizer - should loop on inner_crit->terminate()
                size_t num_vert=patch.num_free_vertices();
                //std::cout << "P[" << get_parallel_rank() << "] tmp srk VertexMover num_vert= " << num_vert << std::endl;
              
        this->optimize_vertex_positions( patch, err );
        if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("optimize_vertex_positions", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
      
          // Update for changes during inner iteration 
          // (during optimizer loop)
        
        outer_crit->accumulate_patch( patch, err );
        if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("outer accumulate_patch", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
        
        inner_crit->cull_vertices( patch, obj_func, err );
        if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("inner cull_vertices", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;

                // experimental...
                if (0)
                  {
                    inner_crit->cull_vertices_global (patch, 
                                                      mesh, domain, settings,
                                                      obj_func, err);
                    if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("cull_vertices_global", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
                  }
        
        patch.update_mesh( err, coord_tag_ptr );
        if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("update_mesh", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
      }
    } // while(p_iter....

    save_or_restore_debug_state(false);

    if (1)
      {
        bool pdone_inner=false;

        if (MSQ_CHKERR(perr))
          {
            std::cout << "P[" << get_parallel_rank() << "] VertexMover::loop_over_mesh found parallel error: " << perr << "\n quitting... pdone_inner= " << pdone_inner << std::endl;
            pdone_inner = true;
          }

        helper->communicate_any_true (pdone_inner, err);

        if (0)
          std::cout << "P[" << get_parallel_rank() << "] tmp srk inner_iter= " << inner_iter << " outer_iter= " << outer_iter << " pdone_inner= " << pdone_inner << std::endl;

        if (pdone_inner)
          {
             if (0)
              std::cout << "P[" << get_parallel_rank() << "] tmp srk found parallel error, quitting... pdone_inner= " << pdone_inner << std::endl;
            MSQ_SETERR(err)("PARALLEL ERROR", MsqError::PARALLEL_ERROR);
            break;
          }
      }

    /// srkenno@sandia.gov save vertex bytes since boundary smoothing changes them
    std::vector<unsigned char> saved_bytes;
    checkpoint_bytes(mesh, saved_bytes, err); 
    if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("checkpoint_bytes ", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;

    helper->communicate_first_independent_set(err); 
    if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("communicate_first_independent_set ", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;

    ///*** smooth the boundary ***////
        save_or_restore_debug_state(true);
        MsqDebug::disable_all();

    while (helper->compute_next_independent_set())
    {
      // Loop over all boundary elements
      while(helper->get_next_partition_boundary_vertex(vertex_handle))
      {

	patch_vertices.clear();
	patch_vertices.push_back(vertex_handle);
	patch_elements.clear(); 
	mesh->vertices_get_attached_elements( &vertex_handle, 
                                              1,
                                              patch_elements, 
                                              junk, err );

	all_culled = false;
	patch.set_mesh_entities( patch_elements, patch_vertices, err );

	if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("set_mesh_entities 2 ", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
        
        // Initialize for inner iteration
        
	this->initialize_mesh_iteration(patch, err);
	if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)(" initialize_mesh_iteration 2", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
	
	obj_func.reset();
	
	outer_crit->reset_patch( patch, err );
	if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)("outer reset_patch 2 ", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
	
	inner_crit->reset_inner( patch, obj_func, err );
	if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)(" inner reset_inner 2", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
      
	inner_crit->reset_patch( patch, err );
	if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)(" inner_crit reset_patch 2", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
      
        // Don't even call optimizer if inner termination 
        // criterion has already been met.
	if (!inner_crit->terminate())
	{
	  inner_crit_terminated = false;
	  
          // Call optimizer - should loop on inner_crit->terminate()
	  this->optimize_vertex_positions( patch, err );
	  if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)(" optimize_vertex_positions 2", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
      
          // Update for changes during inner iteration 
          // (during optimizer loop)
	  
	  outer_crit->accumulate_patch( patch, err );
	  if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)(" outer accumulate_patch 2", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
	  
	  inner_crit->cull_vertices( patch, obj_func, err );
	  if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)(" inner cull_vertices", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;

                    // FIXME
                    if (0)
                      {
                        inner_crit->cull_vertices_global (patch, 
                                                          mesh, domain, settings,
                                                          obj_func, err);
                        if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)(" cull_vertices_global 2 ", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
                      }
        
          patch.update_mesh( err, coord_tag_ptr );
	  if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)(" update_mesh 2", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
	}
      }
      helper->communicate_next_independent_set(err);
      if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)(" communicate_next_independent_set 2", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
     } // while(helper->compute_next_independent_set())
 
    save_or_restore_debug_state(false);
 
    //if (!get_parallel_rank())
    //std::cout << "P[" << get_parallel_rank() << "] tmp srk num_patches= " << num_patches << std::endl;

    /// srkenno@sandia.gov restore vertex bytes since boundary smoothing changes them
    restore_bytes(mesh, saved_bytes, err);
    if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)(" restore_bytes", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;

    if (jacobiOpt)
      commit_jacobi_coords( coord_tag, mesh, err );

    this->terminate_mesh_iteration(patch, err); 
    if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)(" terminate_mesh_iteration", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
    
    outer_crit->accumulate_outer( mesh, domain, obj_func, settings, err );
    if (MSQ_CHKERR(err)) { MSQ_SETERR(perr)(" outer_crit accumulate_outer", MsqError::INVALID_STATE); PERROR_COND; } //goto ERROR;
  }

//ERROR: 

  if (MSQ_CHKERR(err)) {
    std::cout << "P[" << get_parallel_rank() << "] VertexMover::loop_over_mesh error = " << err.error_message() << std::endl;
  }

  if (jacobiOpt)
    mesh->tag_delete( coord_tag, err );

    //call the criteria's cleanup funtions.
  outer_crit->cleanup(mesh,domain,err); MSQ_CHKERR(err);
  inner_crit->cleanup(mesh,domain,err); MSQ_CHKERR(err);
    //call the optimization cleanup function.
  this->cleanup();
    // close the helper
  helper->smoothing_close(err); MSQ_CHKERR(err);

  return 0.;
}

    
void VertexMover::initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                    const Settings* settings,
                                    MsqError& err )
{
  QualityImprover::initialize_queue( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);
  objFuncEval.initialize_queue( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);
}

TagHandle VertexMover::get_jacobi_coord_tag( Mesh* mesh, MsqError& err )
{
    // Get tag handle
  const char tagname[] = "msq_jacobi_temp_coords";
  TagHandle tag = mesh->tag_create( tagname, Mesh::DOUBLE, 3, 0, err );  MSQ_ERRZERO(err);
    /* VertexMover will always delete the tag when it is done, so it is probably
       best not to accept an existing tag
  if (err.error_code() == TAG_ALREADY_EXISTS) {
    err.clear();
    tag = tag_get( tagname, err ); MSQ_ERRZERO(err);
    std::string name;
    Mesh::TagType type;
    unsigned length;
    mesh->tag_properties( tag, name, type, length, err ); MSQ_ERRZERO(err);
    if (type != Mesh::DOUBLE || length != 3) {
      MSQ_SETERR(err)(TAG_ALREADY_EXISTS,
                     "Tag \"%s\" already exists with invalid type", 
                     tagname);
      return 0;
    }
  }
    */
  
    // Initialize tag value with initial vertex coordinates so that
    // vertices that never get moved end up with the correct coords.
  std::vector<Mesh::VertexHandle> vertices;
  mesh->get_all_vertices( vertices, err ); MSQ_ERRZERO(err);
    // remove fixed vertices
    // to avoid really huge arrays (especially since we need to copy
    // coords out of annoying MsqVertex class to double array), work
    // in blocks
  const size_t blocksize = 4096;
  std::vector<bool> fixed(blocksize);
  MsqVertex coords[blocksize];
  double tag_data[3*blocksize];
  for (size_t i = 0; i < vertices.size(); i += blocksize) {
    size_t count = std::min( blocksize, vertices.size() - i );
      // remove fixed vertices
    mesh->vertices_get_fixed_flag( &vertices[i], fixed, count, err ); MSQ_ERRZERO(err);
    size_t w = 0;
    for (size_t j = 0; j < count; ++j)
      if (!fixed[j])
        vertices[i + w++] = vertices[i + j];
    count = w;
      // store tag data for free vertices
    mesh->vertices_get_coordinates( &vertices[i], coords, count, err ); MSQ_ERRZERO(err);
    for (size_t j = 0; j < count; ++j) {
      tag_data[3*j  ] = coords[j][0];
      tag_data[3*j+1] = coords[j][1];
      tag_data[3*j+2] = coords[j][2];
    }
    mesh->tag_set_vertex_data( tag, count, &vertices[i], tag_data, err ); MSQ_ERRZERO(err);
  }
  
  return tag;
}

void VertexMover::commit_jacobi_coords( TagHandle tag, Mesh* mesh, MsqError& err )
{
  Vector3D coords;
  std::vector<Mesh::VertexHandle> vertices;
  mesh->get_all_vertices( vertices, err ); MSQ_ERRRTN(err);
  std::vector<bool> fixed( vertices.size() );
  mesh->vertices_get_fixed_flag( &vertices[0], fixed, vertices.size(), err );
  for (size_t i = 0; i < vertices.size(); ++i) {
    if (!fixed[i]) {
      mesh->tag_get_vertex_data( tag, 1, &vertices[i], coords.to_array(), err ); MSQ_ERRRTN(err);
      mesh->vertex_set_coordinates( vertices[i], coords, err ); MSQ_ERRRTN(err);
    }
  }
}

} // namespace Mesquite
