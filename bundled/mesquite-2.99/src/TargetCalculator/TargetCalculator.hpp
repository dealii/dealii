/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file TargetCalculator.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_CALCULATOR_HPP
#define MSQ_TARGET_CALCULATOR_HPP

#include "Mesquite.hpp"
#include "Sample.hpp"
#include "MsqMatrix.hpp"
#include "MeshInterface.hpp"
#include <stddef.h>

namespace MESQUITE_NS {

class PatchData;
class MsqError;
class ReferenceMeshInterface;
class Mesh;
class MeshDomain;
class Settings;

class MESQUITE_EXPORT TargetCalculator
{
public:
  virtual ~TargetCalculator();
     
   //!\brief Called at start of instruction queue processing
   //!
   //! Do any preliminary global initialization, consistency checking,
   //! etc.  Default implementation does nothing.
  virtual void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                 const Settings* settings,
                                 MsqError& err );

  /**\brief Get a target matrix
   *
   *\param pd      The current PatchData
   *\param element The index an element within the patch data.
   *\param sample  The sample point in the element.
   *\param W_out   The resulting target matrix.
   */
  virtual bool get_3D_target( PatchData& pd, 
                              size_t element,
                              Sample sample,
                              MsqMatrix<3,3>& W_out,
                              MsqError& err ) = 0;

  /**\brief Get a target matrix
   *
   *\param pd      The current PatchData
   *\param element The index an element within the patch data.
   *\param sample  The sample point in the element.
   *\param W_out   The resulting target matrix.
   */
  virtual bool get_surface_target( PatchData& pd, 
                                   size_t element,
                                   Sample sample,
                                   MsqMatrix<3,2>& W_out,
                                   MsqError& err ) = 0;

  /**\brief Get a target matrix
   *
   *\param pd      The current PatchData
   *\param element The index an element within the patch data.
   *\param sample  The sample point in the element.
   *\param W_out   The resulting target matrix.
   */
  virtual bool get_2D_target( PatchData& pd, 
                              size_t element,
                              Sample sample,
                              MsqMatrix<2,2>& W_out,
                              MsqError& err ) = 0;

    /**\brief Use 3x2 W for surface elements if true, 2x2 W if false
     *
     * If true, then the targets for surface elements attempt some
     * control of orientation and therefore get_surface_target must
     * be used to get the targets.  If false, then the target contains
     * no orientation data and is therefore the same as the corresponding
     * 2D target for surface elements.  In this case, get_2D_target should
     * be used.
     */
  virtual bool have_surface_orient() const = 0;

  /**\brief Factor some existing target or Jacobian matrix
   *
   * Utility method for subclasses to use in their implementation
   * of get_3D_target.
   *\param W       The matrix to factor
   *\param Lambda  Output: the size factor of W
   *\praam V       Output: the orientation factor of V
   *\param Q       Output: the skew factor of W
   *\param Delta   Output: the aspect ratio factor of W
   *\return bool   True if W can be factored, false otherwise.
   */
  static bool factor_3D( const MsqMatrix<3,3>& W,
                         double& Lambda,
                         MsqMatrix<3,3>& V,
                         MsqMatrix<3,3>& Q,
                         MsqMatrix<3,3>& Delta,
                         MsqError& err );

  /**\brief Factor some existing target or Jacobian matrix
   *
   * Utility method for subclasses to use in their implementation
   * of get_2D_target
   *\param W       The matrix to factor
   *\param Lambda  Output: the size factor of W
   *\praam V       Output: the orientation factor of V
   *\param Q       Output: the skew factor of W
   *\param Delta   Output: the aspect ratio factor of W
   *\return bool   True if W can be factored, false otherwise.
   */
  static bool factor_surface( const MsqMatrix<3,2>& W,
                         double& Lambda,
                         MsqMatrix<3,2>& V,
                         MsqMatrix<2,2>& Q,
                         MsqMatrix<2,2>& Delta,
                         MsqError& err );

  /**\brief Factor some existing target or Jacobian matrix
   *
   * Utility method for subclasses to use in their implementation
   * of get_2D_target
   *\param W       The matrix to factor
   *\param Lambda  Output: the size factor of W
   *\praam V       Output: the orientation factor of V
   *\param Q       Output: the skew factor of W
   *\param Delta   Output: the aspect ratio factor of W
   *\return bool   True if W can be factored, false otherwise.
   */
  static bool factor_2D( const MsqMatrix<2,2>& W,
                         double& Lambda,
                         MsqMatrix<2,2>& V,
                         MsqMatrix<2,2>& Q,
                         MsqMatrix<2,2>& Delta,
                         MsqError& err );

  /**\brief Get size component of W */
  static double size( const MsqMatrix<3,3>& W );
  /**\brief Get size component of W */
  static double size( const MsqMatrix<3,2>& W );
  /**\brief Get size component of W */
  static double size( const MsqMatrix<2,2>& W );

  /**\brief Get skew component of W */
  static MsqMatrix<3,3> skew( const MsqMatrix<3,3>& W );
  /**\brief Get skew component of W */
  static MsqMatrix<2,2> skew( const MsqMatrix<3,2>& W );
  /**\brief Get skew component of W */
  static MsqMatrix<2,2> skew( const MsqMatrix<2,2>& W );

  /**\brief Get skew component of W */
  static MsqMatrix<3,3> aspect( const MsqMatrix<3,3>& W );
  /**\brief Get skew component of W */
  static MsqMatrix<2,2> aspect( const MsqMatrix<3,2>& W );
  /**\brief Get skew component of W */
  static MsqMatrix<2,2> aspect( const MsqMatrix<2,2>& W );

  /**\brief Get shape (skew and aspect) component of W */
  static MsqMatrix<3,3> shape( const MsqMatrix<3,3>& W );
  /**\brief Get skew component of W */
  static MsqMatrix<2,2> shape( const MsqMatrix<3,2>& W );
  /**\brief Get skew component of W */
  static MsqMatrix<2,2> shape( const MsqMatrix<2,2>& W );
  

  /**\brief Create a new orientation matrix
   *
   * Create an orientation matrix such that
   * the first and second Jacobian columns of W are aligned to 
   * the passed vectors.
   */
  static MsqMatrix<3,3> new_orientation_3D( const MsqVector<3>& b1,
                                            const MsqVector<3>& b2 );
  /**\brief Create a new orientation matrix
   *
   * Create an orientation matrix such that
   * the first and second Jacobian columns of W are aligned to 
   * the passed vectors.
   */
  static MsqMatrix<3,2> new_orientation_2D( const MsqVector<3>& b1,
                                            const MsqVector<3>& b2 );

  /**\brief Get skew matrix for an ideally shaped element */
  static void ideal_skew_3D( EntityTopology element_type,
                             Sample s,
                             const PatchData& pd,
                             MsqMatrix<3,3>& W,
                             MsqError& err );

  /**\brief Get skew matrix for an ideally shaped element */
  static void ideal_skew_2D( EntityTopology element_type,
                             Sample s,
                             const PatchData& pd,
                             MsqMatrix<2,2>& W,
                             MsqError& err );

  /**\brief Get skew matrix for an ideally shaped element */
  static void ideal_shape_3D( EntityTopology element_type,
                              Sample s,
                              const PatchData& pd,
                              MsqMatrix<3,3>& W,
                              MsqError& err );

  /**\brief Get skew matrix for an ideally shaped element */
  static void ideal_shape_2D( EntityTopology element_type,
                              Sample s,
                              const PatchData& pd,
                              MsqMatrix<2,2>& W,
                              MsqError& err );

  /**\brief Create a new aspect ratio matrix
   *
   * Create an aspect ratio matrix such that the ratio of column
   * lengths is proportional to the ratio of the corresponding pair
   * of values in the passed vector.
   */
  static MsqMatrix<3,3> new_aspect_3D( const MsqVector<3>& r );

  /**\brief Create a new aspect ratio matrix
   *
   * Create an aspect ratio matrix such that the ratio of column
   * lengths is proportional to the ratio of the corresponding pair
   * of values in the passed vector.
   */
  static MsqMatrix<2,2> new_aspect_2D( const MsqVector<2>& r );

  /**\brief Create a new aspect ratio matrix
   *
   * Create an aspect ratio matrix such that the ratio of column
   * lengths is the passed value.
   */
  static MsqMatrix<2,2> new_aspect_2D( double rho );

  /**\brief Calculate the Jacobian given element vertex coordinates */
  static void jacobian_3D( PatchData& pd,  // for mapping function list
                           EntityTopology element_type,
                           int num_nodes,
                           Sample location,
                           const Vector3D* coords,
                           MsqMatrix<3,3>& W_out,
                           MsqError& err );

  /**\brief Calculate the Jacobian given element vertex coordinates */
  static void jacobian_2D( PatchData& pd,  // for mapping function list
                           EntityTopology element_type,
                           int num_nodes,
                           Sample location,
                           const Vector3D* coords,
                           MsqMatrix<3,2>& W_out,
                           MsqError& err );


  static void get_refmesh_Jacobian_3D( ReferenceMeshInterface* ref_mesh,
                                       PatchData& active_mesh,
                                       size_t element_no,
                                       Sample sample_no,
                                       MsqMatrix<3,3>& W_out,
                                       MsqError& err );

  static void get_refmesh_Jacobian_2D( ReferenceMeshInterface* ref_mesh,
                                       PatchData& active_mesh,
                                       size_t element_no,
                                       Sample sample_no,
                                       MsqMatrix<3,2>& W_out,
                                       MsqError& err );
};

} // namespace Mesquite

#endif
