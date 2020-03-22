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

    (2006) kraftche@cae.wisc.edu    
  
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ObjectiveFunction.hpp

Header file for the Mesquite::ObjectiveFunction class

  \author Michael Brewer
  \author Thomas Leurent
  \date   2002-05-23
  \author Jason Kraftcheck
  \date   2006-04
 */


#ifndef OBJECTIVE_FUNCTION_HPP
#define OBJECTIVE_FUNCTION_HPP

#include "Mesquite.hpp"

#include <list>
#include <vector>
#include <cstddef>

namespace MESQUITE_NS
{
   class PatchData;
   class MsqHessian;
   class QualityMetric;
   class MsqVertex;
   class MsqError;
   class Vector3D;
   class SymMatrix3D;
   class Mesh;
   class MeshDomain;
   class PatchSet;
   class Settings;
   class MeshDomainAssoc;
   
  /*! \class ObjectiveFunction
       \brief Base class for concrete Objective Functions
       ObjectiveFunction contains a pointer to a QualityMetric.  If
       the ObjectiveFunction is associated with more than one
       QualityMetric (i.e., the Objective is a composite, and the
       composed ObjectiveFunctions are associated with different
       QualityMetrics), then the QualityMetric pointer is set
       to NULL..
      */
  class MESQUITE_EXPORT ObjectiveFunction
  {
   public:

    enum EvalType {
      /** Do not modify or use any accumulated value in the
       *  calculation of the objective function value.  
       *  Evaluate the objective function only over the passed 
       *  patch.  Used for NASH-type solutions.
       */
      CALCULATE,
      /** Incorporate the evaluation of the passed patch into
       *  the accumulated global objective function value.  This
       *  EvalType is used by the default initialize implemenation,
       *  and need not be considered by ObjectiveFunction 
       *  implementations that provide their own implementation of
       *  the initialize method.
       */
      ACCUMULATE,
      /** Save the evaluation over the passed patch such that it
       *  can later be removed from the accumulated global objective
       *  function value.  Assume that the current accumulated 
       *  value already includes the evaluation of this patch. Do
       *  *not* modify any accumulated, global value.
       */
      SAVE,
      /** Assume that the patch passed to this call is a modified
       *  version of the patch previously passed with either the
       *  SAVE or UPDATE EvalType specified.  Update the accumulated
       *  global value accordingly (remove previously saved local
       *  data from accumulated value, and incorporate evaluation
       *  over this patch into the accumulated value.)  Save the
       *  necessary data such that the results of incorporaring
       *  the evaluation of this patch into the accumulated global
       *  value can be removed upon a subsequent call with this
       *  EvalType.
       */
      UPDATE,
      /** For this EvalType, the passed back results should be the
       *  same as for the UPDATE EvalType, but any accumulated data
       *  or data corresponding to the previous local values should
       *  *not* be modified.
       */
      TEMPORARY
    };
    
      // virtual destructor ensures use of polymorphism during destruction
    virtual ~ObjectiveFunction();
     
      //!\brief Called at start of instruction queue processing
      //!
      //! Do any preliminary global initialization, consistency checking,
      //! etc.  This function is pure-virtual (abstract) in this class
      //! because every practical OF implementation at this time should
      //! have an implementation that at least recursively calls the same
      //! function on the underlying QualityMetric or ObjectiveFunction(s).
     virtual void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                    const Settings* settings,
                                    MsqError& err ) = 0;
    
      /**\brief Initial accumulated value for block coordinate descent algorithms
       *
       * Set accumulated value of objective function to the value for the
       * entire, unmodified mesh.  This is the initial state for a block
       * coordinate descent algorithm.  The ObjectiveFunction will asked
       * to add or remove values for a specific patch of the mesh during
       * the optimization.
       *\param mesh  The Mesh
       *\param domain The MeshDomain
       *\param user_set User-defined patch set - not relevant for most OF templates.
       */
    virtual bool initialize_block_coordinate_descent( MeshDomainAssoc* mesh_and_domain,
                                                      const Settings* settings,
                                                      PatchSet* user_set,
                                                      MsqError& err ) = 0;
    
      /**\brief Evaluate objective function for specified patch.
       *
       * Either evaluate the objective function over the passed patch
       * or update the accumulated, global objective function value for
       * changes in the passed patch, depending on the value of the
       * EvalType.
       *\param type  Evaluation type.
       *\param pd    The patch.
       *\param value_out  The passed-back value of the objective fuction. 
       *\param free  If true, incorporate the quality metric values only
       *             for those metric evaluations that depend on at least
       *             one free vertex
       *\return false if any QualityMetric evaluation returned false,
       *        true otherwise.
       */
    virtual bool evaluate( EvalType type, 
                           PatchData& pd,
                           double& value_out,
                           bool free,
                           MsqError& err ) = 0; 
    
    
      /**\brief Evaluate objective function and gradient for specified patch.
       *
       * Either evaluate the objective function over the passed patch
       * or update the accumulated, global objective function value for
       * changes in the passed patch, depending on the value of the
       * EvalType.  
       *
       * The default implementation of this function will
       * use the value-only variation of the evaluate method and numerical
       * approximation to calculate gradients.  Whenever possible, objective
       * function implementations should provide more efficient analyical
       * gradient calculations.
       * 
       *\param type  Evaluation type.
       *\param pd    The patch.
       *\param value_out  The passed-back value of the objective fuction. 
       *\param grad_out The gradient of the OF wrt the coordinates of
       *             each *free* vertex in the patch.
       *\return false if any QualityMetric evaluation returned false,
       *        true otherwise.
       */
    virtual bool evaluate_with_gradient( EvalType type, 
                                         PatchData& pd,
                                         double& value_out,
                                         std::vector<Vector3D>& grad_out,
                                         MsqError& err ); 
    
      /**\brief Evaluate objective function and diagonal blocks of Hessian for specified patch.
       *
       * Either evaluate the objective function over the passed patch
       * or update the accumulated, global objective function value for
       * changes in the passed patch, depending on the value of the
       * EvalType.  
       *
       * The default implementation of this function evaluate the
       * entire Hessian and discard non-diagonal portions.  Concrete
       * objective functions should provide a more efficient implementation
       * that evaluates and accumulates only the required terms.
       * 
       *\param type  Evaluation type.
       *\param pd    The patch.
       *\param value_out  The passed-back value of the objective fuction. 
       *\param grad_out The gradient of the OF wrt the coordinates of
       *             each *free* vertex in the patch.
       *\param hess_diag_out The diagonal blocks of a Hessian.  I.e. Decompose
       *             the Hessian into 3x3 submatrices and return only the 
       *             submatrices (blocks) along the diagonal.
       *\return false if any QualityMetric evaluation returned false,
       *        true otherwise.
       */
    virtual bool evaluate_with_Hessian_diagonal( EvalType type, 
                                        PatchData& pd,
                                        double& value_out,
                                        std::vector<Vector3D>& grad_out,
                                        std::vector<SymMatrix3D>& hess_diag_out,
                                        MsqError& err ); 
    
    
      /**\brief Evaluate objective function and Hessian for specified patch.
       *
       * Either evaluate the objective function over the passed patch
       * or update the accumulated, global objective function value for
       * changes in the passed patch, depending on the value of the
       * EvalType.  
       *
       * The default implementation of this function will fail.  
       * 
       *\param type  Evaluation type.
       *\param pd    The patch.
       *\param value_out  The passed-back value of the objective fuction. 
       *\param grad_out The gradient of the OF wrt the coordinates of
       *             each *free* vertex in the patch.
       *\param Hessian_out The Hessian of the OF wrt the coordinates of
       *             each *free* vertex in the patch.
       *\return false if any QualityMetric evaluation returned false,
       *        true otherwise.
       */
    virtual bool evaluate_with_Hessian( EvalType type, 
                                        PatchData& pd,
                                        double& value_out,
                                        std::vector<Vector3D>& grad_out,
                                        MsqHessian& Hessian_out,
                                        MsqError& err ); 

      /**\brief Create copy with same state
       *
       * Create a new instance of the objective function that
       * is a copy of the callee with the same accumulated
       * values, parameters, etc.
       */
    virtual ObjectiveFunction* clone() const = 0;
    
      /** Clear any values accumulated for BCD-related eval calls */       
    virtual void clear() = 0;
    
      /** Get the minimum number of layers of adjacent elements required
       *  in a patch to evaluate the objective function for a single 
       *  free vertex.  
       */
    virtual int min_patch_layers() const = 0;
    
   protected:
    
      /**\brief Returns eps used in the numerical gradient calculation.
       *
       * Returns an appropiate value (eps) to use as a delta step for
       * MsqVertex vertex in dimension k (i.e. k=0 -> x, k=1 -> y, k=2 -> z).
       * The objective function value at the perturbed vertex position is given
       * in local_val.
       */
    double get_eps(PatchData &pd, EvalType eval_type, double &local_val,
                   int k, size_t vertex_index, MsqError &err);
    
 private:
 
      /**\brief Compute numerical approx. of gradient for 1 vertex's coordinates.
       *
       * Compute the numerical approximation of the gradient of the objective
       * function for the coordinates of a single vertex given a patch for 
       * which that vertex is the only free vertex.
       *\param type   Evaluation type.
       *\param pd     A patch containing a single free vertex and the 
       *              adjacent elements necessary for complete evaluation
       *              of the dependence of the objective fuction on the
       *              coordinates of the free vertex.
       *\param flocal The objective function value for the unmodified
       *              subpatch.  (Output.)
       *\param grad   The gradient of the OF with respect to the coordinates
       *              of the free veretx.  (Output.)
       *\return       The result of calling the ObjectiveFunction::evaluate
       *              method.
       */
    bool compute_subpatch_numerical_gradient( EvalType type,
                                              EvalType get_eps_eval_type,
                                              PatchData& pd,
                                              double& flocal,
                                              Vector3D& grad,
                                              MsqError& err );
      
    bool compute_patch_numerical_gradient( EvalType type,
                                           EvalType get_eps_eval_type,
                                           PatchData& pd,
                                           double& flocal,
                                           std::vector<Vector3D>& grad,
                                           MsqError& err );
  };
  
} //namespace


#endif // ObjectiveFunction_hpp


