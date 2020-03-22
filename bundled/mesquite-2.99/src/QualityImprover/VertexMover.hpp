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
  \file   VertexMover.hpp
  \brief  

  The VertexMover Class is the base class for all the smoothing and 
  optimizing algorythms 

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_VertexMover_hpp 
#define Mesquite_VertexMover_hpp


#include "Mesquite.hpp"
#include "QualityImprover.hpp"
#include "OFEvaluator.hpp"
#include "MeshInterface.hpp"

namespace MESQUITE_NS
{
  class ObjectiveFunction;
  class PatchData;
  class ParallelMesh;

  /*! \class VertexMover
    Base class for all Vertex Movers.
   */  
  class MESQUITE_EXPORT VertexMover : public QualityImprover 
  {
  protected:
    VertexMover( ObjectiveFunction* OF = NULL );
    
  public:
    // virtual destructor ensures use of polymorphism during destruction
    virtual ~VertexMover();
    
    virtual void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                   const Settings* settings,
                                   MsqError& err );
    
    virtual double loop_over_mesh( MeshDomainAssoc* mesh_and_domain,
                                   const Settings* settings,
                                   MsqError &err);

    virtual double loop_over_mesh( ParallelMesh* mesh, 
                                   MeshDomain* domain,
                                   const Settings* settings,
                                   MsqError &err);

      /**\brief Do Nash-game type optimization
       *
       * Default, opposite of \c do_block_coordinate_descent().
       */
    inline
    void do_nash_game_optimization()
      { objFuncEval.do_nash_game(); }
      
      /**\brief Check if optmizer will do Nash-game type optimization
       *
       * Default, opposite of \c is_block_coordinate_descent_optimization().
       */
    inline
    bool is_nash_game_optimization() const
      { return objFuncEval.is_nash_game(); }

      /**\brief Do block coordinate descent optimization
       *
       * Opposite of \c do_nash_game().
       */
    inline
    void do_block_coordinate_descent_optimization()
      { objFuncEval.do_block_coordinate_descent(); }
      
      /**\brief Check if optmizer will do block coordinate descent type optimization
       *
       * Default, opposite of \c is_nash_game_optimization().
       */
    inline
    bool is_block_coordinate_descent_optimization() const
      { return objFuncEval.is_block_coordinate_descent(); }
    
      /**\brief Use Jacobi iteration for optimization
       *
       * Opposite of \c do_gauss_optimization()
       */
    inline
    void do_jacobi_optimization()
      { jacobiOpt = true; }

      /**\brief Check if optimization will use Jacobi iteration
       *
       * Opposite of \c is_gauss_optimization()
       */
    inline
    bool is_jacobi_optimization() const 
      { return jacobiOpt; }
    
      /**\brief Use Gauss-Seidel iteration for optimization
       *
       * Default, opposite of \c do_jacobi_optimization()
       */
    inline
    void do_gauss_optimization()
      { jacobiOpt = false; }

      /**\brief Check if optimization will use Gauss-Seidel iteration
       *
       * Default, opposite of \c is_jacobi_optimization()
       */
    inline
    bool is_gauss_optimization() const 
      { return !jacobiOpt; }

  protected:

    virtual void initialize(PatchData &pd, MsqError &err) = 0;
    virtual void cleanup() = 0;
    virtual void optimize_vertex_positions(PatchData &pd, 
                                           MsqError &err) = 0; // modifies the PatchData object

    virtual void initialize_mesh_iteration(PatchData &pd, 
                                         MsqError &err) = 0;
    virtual void terminate_mesh_iteration(PatchData &, 
                                         MsqError &err) = 0;

  
    OFEvaluator& get_objective_function_evaluator()
      { return objFuncEval; }
  
    static TagHandle get_jacobi_coord_tag( Mesh* mesh, MsqError& err );
    static void commit_jacobi_coords( TagHandle tag, Mesh* mesh, MsqError& err );
  
  private:
    OFEvaluator objFuncEval;
    bool jacobiOpt;
  };

  
} // namespace
#endif // Mesquite_VertexMover_hpp
