/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */

#ifndef MSQ_OFEVALUATOR_HPP
#define MSQ_OFEVALUATOR_HPP

/** \file OFEvaluator.hpp
 *  \brief 
 *  \author Jason Kraftcheck
 */

#include <vector>

#include "Vector3D.hpp"
#include "ObjectiveFunction.hpp"

namespace MESQUITE_NS {

class MsqError;
class Mesh;
class MeshDomain;
class PatchData;
class MsqHessian;
class PatchSet;

/**\brief Evaluate objective function
 *
 * This class handles the details of the differences between
 * Nash and block coordinate descent methods for evaluation
 * of the objective function, such that solvers need only
 * interact with this interface and need not be aware of the
 * Nash vs. BCD details.
 */
class MESQUITE_EXPORT OFEvaluator 
{
public:
  
  /**\brief Constructor
   *\param of The objective function (may be NULL for Laplacian-type solvers)
   *\param Nash True for Nash-type solutions, false for
   *       block coordinate descent.
   */
  OFEvaluator( ObjectiveFunction* of );
    
  void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                         const Settings* settings,
                         MsqError& err );
  
  /**\brief Initialize OFEvaluator 
   *
   * For a Nash-type algorithm, this method will initialize
   * some member variables.  For a block coordinate descent
   * solution, this method will calculate an initial value
   * of the objective function for the mesh.
   *
   *\param mesh The active mesh
   */
  bool initialize( MeshDomainAssoc* mesh_and_domain,
                   const Settings* settings,
                   PatchSet* patches,
                   MsqError& err );
  
  /**\brief Update accumulated values for changes to vertex positions
   *        in a patch.
   *
   * For a block coordinate descent solution, this method calculates
   * the updated global objective function value for any modifications
   * to the passed patch, as made by the solver.  The change to
   * the current patch state is considered relative to that of the
   * previous patch passed to any of the update methods, except when
   * the reset() method has been called by the solver to indicate
   * that a new inner iteration is starting.
   *
   * For a Nash-type solution, this method simply returns the evaluation
   * of the objective funtion for the local patch.  The behavior is identical
   * to calling the evaluate() method.
   *
   *\param pd  The mesh patch
   *\param value Output, the value of the objective function.
   */
  bool update( PatchData& pd, 
               double& value,
               MsqError& err );

  
  /**\brief Update accumulated values for changes to vertex positions
   *        in a patch.
   *
   * For a block coordinate descent solution, this method calculates
   * the updated global objective function value for any modifications
   * to the passed patch, as made by the solver.  The change to
   * the current patch state is considered relative to that of the
   * previous patch passed to any of the update methods, except when
   * the reset() method has been called by the solver to indicate
   * that a new inner iteration is starting.
   *
   * For a Nash-type solution, this method simply returns the evaluation
   * of the objective funtion for the local patch.  The behavior is identical
   * to calling the evaluate() method.
   *
   *\param pd  The mesh patch
   *\param value Output, the value of the objective function.
   *\param grad Output, the gradient of the objective function
   *             with respect to each FREE vertex in the patch.
   */
  bool update( PatchData& pd, double& value, 
               std::vector<Vector3D>& grad,
               MsqError& err );

  
  /**\brief Update accumulated values for changes to vertex positions
   *        in a patch.
   *
   * For a block coordinate descent solution, this method calculates
   * the updated global objective function value for any modifications
   * to the passed patch, as made by the solver.  The change to
   * the current patch state is considered relative to that of the
   * previous patch passed to any of the update methods, except when
   * the reset() method has been called by the solver to indicate
   * that a new inner iteration is starting.
   *
   * For a Nash-type solution, this method simply returns the evaluation
   * of the objective funtion for the local patch.  The behavior is identical
   * to calling the evaluate() method.
   *
   *\param pd  The mesh patch
   *\param value Output, the value of the objective function.
   *\param grad Output, the gradient of the objective function
   *             with respect to each FREE vertex in the patch.
   *\param Hessian_diag_blocks Output, 3x3 submatrices along diagonal of 
   *                           Hessian of objective function
   */
  bool update( PatchData& pd, double& value,
               std::vector<Vector3D>& grad, 
               std::vector<SymMatrix3D>& Hessian_diag_blocks,
               MsqError& err );
  
  /**\brief Update accumulated values for changes to vertex positions
   *        in a patch.
   *
   * For a block coordinate descent solution, this method calculates
   * the updated global objective function value for any modifications
   * to the passed patch, as made by the solver.  The change to
   * the current patch state is considered relative to that of the
   * previous patch passed to any of the update methods, except when
   * the reset() method has been called by the solver to indicate
   * that a new inner iteration is starting.
   *
   * For a Nash-type solution, this method simply returns the evaluation
   * of the objective funtion for the local patch.  The behavior is identical
   * to calling the evaluate() method.
   *
   *\param pd  The mesh patch
   *\param value Output, the value of the objective function.
   *\param grad Output, the gradient of the objective function
   *             with respect to each FREE vertex in the patch.
   *\param Hessian Output, the Hessian of the objective function.
   */
  bool update( PatchData& pd, double& value,
               std::vector<Vector3D>& grad, 
               MsqHessian& Hessian,
               MsqError& err );
               
  /**\brief Check if doing Nash game algorithm.*/
  bool is_nash_game() const { return !doBCD; }
  
  /**\brief Do Nash game algorithm.*/
  void do_nash_game() { doBCD = false; }
  
  /**\brief Check if doing block coordinate descent algorithm */
  bool is_block_coordinate_descent() const { return doBCD; }
  
  /**\brief Do block coordinate descent algorithm.*/
  void do_block_coordinate_descent() { doBCD = true; }
  
  /**\brief Evaluate the objective function without changing any 
   *        accumulated values.
   *
   * Evaluate the objective function for the specified patch
   * (or for the change to the specified patch for BCD).  This 
   * method does not change any internal state or accumulated
   * values.  It is provided for FeasibleNewton and other solvers
   * that need to obtain an OF value for some intermediate or 
   * temporary set of vertex positions.
   *
   *\param pd  The mesh patch
   *\param value Output, the value of the objective function.
   */
  bool evaluate( PatchData& pd, 
                 double& value,
                 MsqError& err ) const;

  /**\brief Evaluate the objective function without changing any 
   *        accumulated values.
   *
   * Evaluate the objective function for the specified patch
   * (or for the change to the specified patch for BCD).  This 
   * method does not change any internal state or accumulated
   * values.  It is provided for FeasibleNewton and other solvers
   * that need to obtain an OF value for some intermediate or 
   * temporary set of vertex positions.
   *
   *\param pd  The mesh patch
   *\param value Output, the value of the objective function.
   *\param grad Output, the gradient of the objective function
   *             with respect to each FREE vertex in the patch.
   */
  bool evaluate( PatchData& pd, 
                 double& value, 
                 std::vector<Vector3D>& grad,
                 MsqError& err ) const;

  /**\brief Evaluate the objective function without changing any 
   *        accumulated values.
   *
   * Evaluate the objective function for the specified patch
   * (or for the change to the specified patch for BCD).  This 
   * method does not change any internal state or accumulated
   * values.  It is provided for FeasibleNewton and other solvers
   * that need to obtain an OF value for some intermediate or 
   * temporary set of vertex positions.
   *
   *\param pd  The mesh patch
   *\param value Output, the value of the objective function.
   *\param grad Output, the gradient of the objective function
   *             with respect to each FREE vertex in the patch.
   *\param Hessian_diag_blocks Output, 3x3 submatrices along diagonal of 
   *                           Hessian of objective function
   */
  bool evaluate( PatchData& pd, 
                 double& value,
                 std::vector<Vector3D>& grad, 
                 std::vector<SymMatrix3D>& Hessian_diag_blocks,
                 MsqError& err ) const;

  /**\brief Evaluate the objective function without changing any 
   *        accumulated values.
   *
   * Evaluate the objective function for the specified patch
   * (or for the change to the specified patch for BCD).  This 
   * method does not change any internal state or accumulated
   * values.  It is provided for FeasibleNewton and other solvers
   * that need to obtain an OF value for some intermediate or 
   * temporary set of vertex positions.
   *
   *\param pd  The mesh patch
   *\param value Output, the value of the objective function.
   *\param grad Output, the gradient of the objective function
   *             with respect to each FREE vertex in the patch.
   *\param Hessian Output, the Hessian of the objective function.
   */
  bool evaluate( PatchData& pd, 
                 double& value,
                 std::vector<Vector3D>& grad, 
                 MsqHessian& Hessian,
                 MsqError& err ) const;
  
  /**\brief Reset for next inner iteration
   *
   * The control code for the vertex mover is responsible for
   * calling this method before the beginning of each inner
   * solver iteration so the necessary internal state can be
   * updated for the correct behavior of the first call to
   * the update() method for block coordinate descent algorithms.
   */
  bool reset();
  
  /**\brief Get ObjectiveFunction pointer */
  inline
  ObjectiveFunction* get_objective_function() const
    { return this->OF; }
    
  /**\brief Check if we have an objective function */
  inline
  bool have_objective_function() const
    { return 0 != get_objective_function(); }

private:

  /**\brief Disallow copying*/
  OFEvaluator( const OFEvaluator& );
  
  /**\brief Disallow assignment*/
  OFEvaluator& operator=( const OFEvaluator& );

  /** The ObjectiveFunction to evaluate */
  ObjectiveFunction *const OF;
  
  /** Nash or BCD */
  bool doBCD;
  
  /** Nash vs. BCD and state of BCD data */
  ObjectiveFunction::EvalType tempType, firstType, updateType, currUpdateType;
};

} // namespace Mesquite

#endif
