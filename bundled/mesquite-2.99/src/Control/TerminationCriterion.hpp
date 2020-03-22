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

/*! \file TerminationCriterion.hpp

Header file for the TerminationCriterion classes.

  \author Michael Brewer
  \author Thomas Leurent
  \date   Feb. 14, 2003
 */


#ifndef TerminationCriterion_hpp
#define TerminationCriterion_hpp

#include "Mesquite.hpp"
#include "MsqTimer.hpp"
#include "Vector3D.hpp"

#include <string>
#include <vector>
#include <fstream>

namespace MESQUITE_NS
{
   class MsqError;
   class OFEvaluator;
   class PatchData;
   class PatchDataVerticesMemento;
   class Mesh;
   class MeshDomain;
   class MeshDomainAssoc;
   class Settings;
  class VertexMover;

  /*! \class TerminationCriterion

      \brief The TerminationCriterion class contains functionality to
      terminate the VertexMover's optimization.

      The Termination Criterion class has three roles.  It
      is used to terminate the optimization on a single patch; it
      is used to terminate the iterations over all patches in the
      mesh; and it is used to cull vertices from the optimization
      processes.  Thus, for each optimization, two TerminationCriterion
      objects are used.  The class contains five important member
      functions used in the VertexMover:  initialize(), reset(),
      terminate(), cull_vertices(), and cleanup().  These functions
      are each explained in detail below.  In general, the only one
      of these functions called directly from a concrete VertexMover
      is terminate() which allows the concrete VertexMover to determine
      when to stop producing new iterates on a given patch.  All other
      functionality is handled from the base VertexMover base class.

      There are several different types of termination criteria
      available. Multiple criteria types can be set on a given
      Termination Criterion object, and when this occurs, the
      optimization process will terminate whenever any of the
      criteria have been satisfied.
      
      The following is a brief description of how TerminationCriterion
      is used within Mesquite.  Functions called during QualityImprovement
      can be divided into three groups:
        reset_*      - Initialize data for an iteration
        accumulate_* - Update TC for changed data during iteration
        terminate    - Check if the termination criterion has been met.
      There are three different forms of the reset_* and accumulate_*
      functions which are called on the inner, outer, or both 
      TerminationCriterion classes:
        *_outer      - Called on outer termination criterion.
        *_inner      - Called on inner termination criterion.
        *_patch      - Called on outer termination criterion for
                       each patch and on inner termination criterion
                       for each inner iteration.
      
      If implementing a new TerminationCriterion, the following rules
      should be followed.  If the value must be calculated on a global
      patch for the outer TC, then:
        o The functionality should be added to *_inner (yes, INNER) 
        o The *_outer methods should be updated to call the *_inner 
            with a global patch when your TC is requested.
        o The internal data for any such TC should be initialized 
          in the reset_inner method.  
      If the value for the outer criterion can be calculated from each 
      local patch when iterating over the mesh with local patches, then:
        o The functionality should be added to *_patch
        o Any state values pertaining to the entire iteration must be 
           initialized in reset_inner(..) and cleared in terminate()
        o Any patch-specific data should be initialized in reset_patch
        o Care should be taken that terminate() does not check 
          uninitialized data if called before the first call to
          accumulate_patch()

  */
  class TerminationCriterion
  {
  public:
    
       //! checks the gradient \f$\nabla f \f$ of objective function 
       //! \f$f : I\!\!R^{3N} \rightarrow I\!\!R \f$ against a double \f$d\f$  
       //! and stops when \f$\sqrt{\sum_{i=1}^{3N}\nabla f_i^2}<d\f$  
    MESQUITE_EXPORT void add_absolute_gradient_L2_norm( double value );

       //! checks the gradient \f$\nabla f \f$ of objective function 
       //! \f$f : I\!\!R^{3N} \rightarrow I\!\!R \f$ against a double \f$d\f$  
       //! and stops when \f$ \max_{i=1}^{3N} \nabla f_i < d \f$  
    MESQUITE_EXPORT void add_absolute_gradient_inf_norm( double value );
    
         //!terminates on the j_th iteration when
         //! \f$\sqrt{\sum_{i=1}^{3N}\nabla f_{i,j}^2}<d\sqrt{\sum_{i=1}^{3N}\nabla f_{i,0}^2}\f$
         //!  That is, terminates when the norm of the gradient is smaller
         //!  than the specified fraction of the initial norm of the gradient. 
    MESQUITE_EXPORT void add_relative_gradient_L2_norm( double value );

       //!terminates on the j_th iteration when
         //! \f$\max_{i=1 \cdots 3N}\nabla f_{i,j}<d \max_{i=1 \cdots 3N}\nabla f_{i,0}\f$
         //!  That is, terminates when the norm of the gradient is small
         //! than some scaling factor times the norm of the original gradient.
         //! (Using the infinity norm.)
    MESQUITE_EXPORT void add_relative_gradient_inf_norm( double value );

         //!Terminates when the objective function value is smaller than
         //! the given scalar value.
    MESQUITE_EXPORT void add_absolute_quality_improvement( double value );

         //!Terminates when the objective function value is smaller than
         //! the given scalar value times the original objective function
         //! value.
    MESQUITE_EXPORT void add_relative_quality_improvement( double value );

         //!Terminates when a the maximum distance moved by any vertex
         //! during the previous iteration is below the given value.
    MESQUITE_EXPORT void add_absolute_vertex_movement( double value );

         //!Terminates when a the maximum distance moved by any vertex
         //! during the previous iteration is below the given value
         //! times the maximum distance moved by any vertex over the
         //! entire course of the optimization.
    MESQUITE_EXPORT void add_relative_vertex_movement( double value );

         //! Calculates a constant \f$ \delta = \beta \times (a - \sigma) \f$
         //! and terminates when the maximum vertex movement of an iteration
         //! is less than delta.
         //!
         //! \f$ \beta \f$ is the passed value, \f$ a \f$ 
         //! is the average edge length of the initial mesh and \f$ \sigma \f$
         //! is the standard deviation of the edge lengths of the initial
         //! mesh.  The initial mesh values are (re)calcualted for each 
         //! time the instruction queue is run. 
         //!
         //!\param value \f$ beta \f$.  Must be in range (0,1), exclusive.
    MESQUITE_EXPORT void add_absolute_vertex_movement_edge_length( double value );

         //!Terminates when the decrease in the objective function value since
         //! the previous iteration is below the given value.
    MESQUITE_EXPORT void add_absolute_successive_improvement( double value );

         //!Terminates when the decrease in the objective function value since
         //! the previous iteration is below the given value times the
         //! decrease in the objective function value since the beginning
         //! of this optimization process.
    MESQUITE_EXPORT void add_relative_successive_improvement( double value );
    
         //!Terminates when the algorithm exceeds an allotted time limit
         //! (given in seconds).
    MESQUITE_EXPORT void add_cpu_time( double seconds );
    
         //!Terminates when the number of iterations exceeds a given integer.
    MESQUITE_EXPORT void add_iteration_limit( unsigned int max_iterations );
    
         //!Terminates when any vertex leaves the bounding box, defined
         //! by the given value, d.  That is, when the absolute value of
         //! a single coordinate of vertex's position exceeds d.
    MESQUITE_EXPORT void add_bounded_vertex_movement( double value);
    
     //!Terminates when the mesh is detected to be untangled.
     //! Uses the same approach as QualityAssessor,
     //! checks the tau values at all the sample points.
    MESQUITE_EXPORT void add_untangled_mesh();
    
    MESQUITE_EXPORT void remove_all_criteria();
    
         //!Cull when the objective function value is smaller than
         //! the given scalar value.
    MESQUITE_EXPORT void cull_on_absolute_quality_improvement( double limit );
         //!Cull when the objective function value is smaller than
         //! the given scalar value times the original objective function
         //! value.
    MESQUITE_EXPORT void cull_on_relative_quality_improvement( double limit );
         //!Cull when a the maximum distance moved by any vertex
         //! during the previous iteration is below the given value.
    MESQUITE_EXPORT void cull_on_absolute_vertex_movement( double limit );
         //!Cull when a the maximum distance moved by any vertex
         //! during the previous iteration is below the given value
         //! times the maximum distance moved by any vertex over the
         //! entire course of the optimization.
    MESQUITE_EXPORT void cull_on_relative_vertex_movement( double limit );
         //!Cull when the decrease in the objective function value since
         //! the previous iteration is below the given value.
    MESQUITE_EXPORT void cull_on_absolute_successive_improvement( double limit );
         //! Calculates a constant \f$ \delta = \beta \times (a - \sigma) \f$
         //! and culls when the vertex movement  is less than delta.
         //!
         //! \f$ \beta \f$ is the passed value, \f$ a \f$ 
         //! is the average edge length of the initial mesh and \f$ \sigma \f$
         //! is the standard deviation of the edge lengths of the initial
         //! mesh.  The initial mesh values are (re)calcualted for each 
         //! time the instruction queue is run. 
         //!
         //!\param value \f$ beta \f$.  Must be in range (0,1), exclusive.
    MESQUITE_EXPORT void cull_on_absolute_vertex_movement_edge_length( double value );
         //!Cull when the decrease in the objective function value since
         //! the previous iteration is below the given value times the
         //! decrease in the objective function value since the beginning
         //! of this optimization process.
    MESQUITE_EXPORT void cull_on_relative_successive_improvement( double limit );

         //!Cull for a global patch - sets soft fixed flags for vertices
         //! that touch elements that are culled by above flags.
    MESQUITE_EXPORT void cull_for_global_patch( bool val=true );
    
     //!Cull when the mesh is detected to be untangled.
     //! Uses the same approach as QualityAssessor,
     //! checks the tau values at all the sample points.
    MESQUITE_EXPORT void cull_untangled_mesh();
    
    MESQUITE_EXPORT void remove_culling();
    
    enum InnerOuterType {
      TYPE_UNKNOWN,
      TYPE_INNER,
      TYPE_OUTER
    };
    
      //!Constructor which does not take any arguements
    MESQUITE_EXPORT TerminationCriterion(std::string name="", InnerOuterType innerOuterType=TYPE_UNKNOWN);
    
      //!Destructor
    MESQUITE_EXPORT ~TerminationCriterion(){};

      //!This function returns the current function value.
      /*! \todo Michael:  this function is not reliable.  It
        needs to be more robust.  How do we know whether
        currentOFValue got updated or not?  We may want to
        make sure that all the criteria get checked.*/
    MESQUITE_EXPORT double get_current_function_value()
       {return currentOFValue;}
       
    MESQUITE_EXPORT void set_debug_output_level( int i )
      { debugLevel = i; }
    
    enum TimeStepFileType { NOTYPE = 0, VTK, GNUPLOT };

    /**\brief Write mesh improvement animation 
     *
     * Write mesh at each iteration such that the sequence of mesh files can be used
     * to produce an animation of the mesh through the quality improvement process.
     * 
     * Files can be written either as VTK timesteps for viewing in a tool such as Paraview
     * or as GNU plot data files and a GNU plot script which, in combination with recent
     * versions of GNU plot, can be used to produce an animated GIF image.
     *
     * Writing of mesh steps can be disabled by calling this function with type == NOTYPE
     * and filename ==NULL.
     */
    MESQUITE_EXPORT void write_mesh_steps( const char* filename, TimeStepFileType type = VTK )
      { timeStepFileName = filename; timeStepFileType = type; }
    
    /*\brief Write data for generating plots of optimization process.
     *
     * Write data files in a format suitable for use with GNU plot and other applications.
     * The file will contain data corresponding to all termination criteria,
     * however, some values may be invalid if they are not calculated for
     * use in a termination criterion.
     */
    MESQUITE_EXPORT void write_iterations( const char* filename, MsqError& err );
    
    MESQUITE_EXPORT int get_iteration_count( ) const
      { return iterationCounter; }
      
    
      //! Clear any data accumulated during an outer iteration
    void reset_outer( Mesh* ms, MeshDomain* dm, OFEvaluator& of, 
                      const Settings* settings, MsqError& err );
    
      //! Clear any data accumulated during an inner iteration
    void reset_inner( PatchData& pd, OFEvaluator& of, MsqError& err );
    
      //! Shared inner and outer initialization during inner loop
    void reset_patch( PatchData& pd, MsqError& err );
    
      //! Accumulate data during inner iteration
    void accumulate_inner( PatchData& pd, OFEvaluator& eval, MsqError& err );
    
      //! Accumulate data during inner iteration
    void accumulate_inner( PatchData& pd, double of_value, Vector3D* of_grads, 
                           MsqError& err );
    
      //! Common code for both inner and outer termination 
      //! criteria during inner iteration.                       
    void accumulate_patch( PatchData& pd, MsqError& err );
    
    void accumulate_outer( Mesh* ms, MeshDomain* dm, OFEvaluator& eval, 
                           const Settings* settings, MsqError& err );
    
      //! Check if termination criterion has been met
    MESQUITE_EXPORT bool terminate();

      //! Check if at least one termination criterion is set
    MESQUITE_EXPORT bool criterion_is_set();
    
    
      //!Function which determines whether this patch should be 'culled'
    bool cull_vertices(PatchData &pd, OFEvaluator& obj_ptr, MsqError &err);

      //!experimental, first cut at culling for global patches - not finished
    bool cull_vertices_global(PatchData &global_patch,
                              Mesh *mesh, MeshDomain *domain, const Settings *settings,
                              OFEvaluator& of_eval,
                              MsqError &err);

      //!Cleans up after the TerminationCriterion is finished.
    void cleanup(Mesh* ms, MeshDomain* domain, MsqError &err);
    
    void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                           const Settings* settings,
                           MsqError& err );

    friend class VertexMover;

 protected:
    
    void write_timestep( PatchData& pd, const Vector3D* gradient, MsqError& err );
    
    static size_t count_inverted( PatchData& pd, MsqError& err );

    std::string par_string();  // for debug print of parallel rank
    
 private:
    //PRIVATE DATA MEMBERS
    long unsigned int terminationCriterionFlag;//!<Bit flag of termination crit
    long unsigned int cullingMethodFlag;/*!<Bit flag of criterion for culling*/
      //epsiloon used in culling methods.
    double cullingEps;

    bool cullingGlobalPatch;/*!<enable culling of pieces of a global patch*/

      //Data not specific to a single criterion
    double initialOFValue;
    double previousOFValue;
    double currentOFValue;
    double lowerOFBound;

      //Data specific to termination criterion 1 (gradient bounds)
    std::vector<Vector3D> mGrad;
    double initialGradL2NormSquared;
    double currentGradL2NormSquared;
    double gradL2NormAbsoluteEpsSquared;
    double gradL2NormRelativeEpsSquared;
    double initialGradInfNorm;
    double currentGradInfNorm;
    double gradInfNormAbsoluteEps;
    double gradInfNormRelativeEps;
      //Data specific to termination criterion 2 (KKT)
      //???????????????????????????????????????????
      //Data specific to termination criterion 3 (Quality Improvement)
    double qualityImprovementAbsoluteEps;
    double qualityImprovementRelativeEps;
      //Data specific to termination criterion 4 (inner iterations)
    int iterationBound;
    int iterationCounter;
      //Data specific to termination criterion 5 (cpu time)
    Timer mTimer;
    double timeBound;
      //Data specific to termination criterion 6 (vertex movement)
    PatchDataVerticesMemento* initialVerticesMemento;
    PatchDataVerticesMemento* previousVerticesMemento;//if we want relative
    double vertexMovementAbsoluteEps;
    double vertexMovementRelativeEps;
    double vertexMovementAvgBeta; //!< input beta value used to calculate \c vertexMovementAbsoluteAvg
    double vertexMovementAbsoluteAvgEdge; //!< calculated constant for \c VERTEX_MOVEMENT_ABS_EDGE_LENGTH
    double maxSquaredInitialMovement;
    double maxSquaredMovement;
    
      //Data specific to termination criterion 7 (successive improvement to F)
    double successiveImprovementsAbsoluteEps;
    double successiveImprovementsRelativeEps;
      //crit 8
    double boundedVertexMovementEps;
    int vertexMovementExceedsBound;
      
      // Data for untangled criterion
    size_t globalInvertedCount; //!< number of inverted elements in entire mesh
    size_t patchInvertedCount;  //!< number of inverted elements in previously tested patch
    
    int debugLevel;
    
    //! Plot data
    std::ofstream plotFile;
    
    //! Base name for timestep files
    std::string timeStepFileName;    
    TimeStepFileType timeStepFileType;
    std::string moniker;
    InnerOuterType innerOuterType;
  };

} //namespace


#endif // TerminationCriterion_hpp
