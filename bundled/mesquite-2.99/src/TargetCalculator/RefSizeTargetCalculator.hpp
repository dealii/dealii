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


/** \file RefSizeTargetCalculator.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_REF_SIZE_TARGET_CALCULATOR_HPP
#define MSQ_REF_SIZE_TARGET_CALCULATOR_HPP

#include "Mesquite.hpp"
#include "IdealShapeTarget.hpp"

namespace MESQUITE_NS {

class ReferenceMesh;

class RefSizeTargetCalculator : public TargetCalculator
{
public:
  RefSizeTargetCalculator( ReferenceMesh* reference_mesh,
                           TargetCalculator* tc );
  RefSizeTargetCalculator( ReferenceMesh* reference_mesh );
  
  virtual bool get_3D_target( PatchData& pd, 
                              size_t element,
                              Sample sample,
                              MsqMatrix<3,3>& W_out,
                              MsqError& err );

  virtual bool get_2D_target( PatchData& pd, 
                              size_t element,
                              Sample sample,
                              MsqMatrix<2,2>& W_out,
                              MsqError& err );

  virtual bool get_surface_target( PatchData& pd, 
                              size_t element,
                              Sample sample,
                              MsqMatrix<3,2>& W_out,
                              MsqError& err );
  
  virtual bool have_surface_orient() const 
    { return scaledTargets->have_surface_orient(); }
    
private:
  
  double average_edge_length( PatchData& pd, size_t element, MsqError& err );

  ReferenceMesh* refMesh;
  IdealShapeTarget defaultTargets;
  TargetCalculator* scaledTargets;
  
  /** Amount to scale average edge length by to achive correctly sized
   *  ideal area/volume target, indexed by element topology.
   */
  // double scaleFactor[MIXED];
};



} // namespace MESQUITE_NS

#endif
