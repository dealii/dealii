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


/** \file RefMeshTargetCalculator.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_REF_MESH_TARGET_CALCULATOR_HPP
#define MSQ_REF_MESH_TARGET_CALCULATOR_HPP

#include "Mesquite.hpp"
#include "Vector3D.hpp"
#include "JacobianCalculator.hpp"
#include "MeshInterface.hpp"
#include "TargetCalculator.hpp"
#include "NodeSet.hpp"

namespace MESQUITE_NS {

class ReferenceMeshInterface;

class RefMeshTargetCalculator : public TargetCalculator
{
public:
  MESQUITE_EXPORT RefMeshTargetCalculator( ReferenceMeshInterface* ref_mesh,
                                           bool orient_2d = false )
    : refMesh(ref_mesh), orient2D(orient_2d) {}


  MESQUITE_EXPORT virtual ~RefMeshTargetCalculator();

  /**\brief Get a target matrix
   *
   *\param pd      The current PatchData
   *\param element The index an element within the patch data.
   *\param sample  The sample point in the element.
   *\param W_out   The resulting target matrix.
   */
  MESQUITE_EXPORT virtual 
  bool get_3D_target( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<3,3>& W_out,
                      MsqError& err );

  /**\brief Get a target matrix
   *
   *\param pd      The current PatchData
   *\param element The index an element within the patch data.
   *\param sample  The sample point in the element.
   *\param W_out   The resulting target matrix.
   */
  MESQUITE_EXPORT virtual 
  bool get_2D_target( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<2,2>& W_out,
                      MsqError& err );

  /**\brief Get a target matrix
   *
   *\param pd      The current PatchData
   *\param element The index an element within the patch data.
   *\param sample  The sample point in the element.
   *\param W_out   The resulting target matrix.
   */
  MESQUITE_EXPORT virtual 
  bool get_surface_target( PatchData& pd, 
                           size_t element,
                           Sample sample,
                           MsqMatrix<3,2>& W_out,
                           MsqError& err );

  virtual bool have_surface_orient() const
    { return orient2D; }

private:
  ReferenceMeshInterface* refMesh;
  bool orient2D;
};

} // namespace Mesquite

#endif
