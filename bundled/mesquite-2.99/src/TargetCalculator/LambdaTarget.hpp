/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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


/** \file LambdaTarget.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_LAMBDA_TARGET_HPP
#define MSQ_LAMBDA_TARGET_HPP

#include "Mesquite.hpp"
#include "TargetCalculator.hpp"

namespace MESQUITE_NS {

/**\brief Scale a target matrix by the size of another.
 *
 * Combine two target matrices by extracting a scalar representation of 
 * the size (lambda) of one and use that value to scale the other.  The
 * expected use of this target calculator is to preserve the size of
 * a reference mesh while targeting ideal elements for shape.
 */
class LambdaTarget : public TargetCalculator
{
public:
  /**
   *\param lambda_source Target calculator from which to extract
   *                     scaling factor (lambda).
   *\param composite_source Target calcualtor from which to obtain
   *                     a target that will be scaled by lambda.
   */
  LambdaTarget( TargetCalculator* lambda_source,
                TargetCalculator* composite_source );
  
  ~LambdaTarget();
  
  bool get_3D_target( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<3,3>& W_out,
                      MsqError& err );

  bool get_2D_target( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<2,2>& W_out,
                      MsqError& err );

  bool get_surface_target( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<3,2>& W_out,
                      MsqError& err );
  
  bool have_surface_orient() const
    { return compositeSource->have_surface_orient(); }
    
private:
  TargetCalculator* lambdaSource;
  TargetCalculator* compositeSource;
};


} // namespace MESQUITE_NS

#endif
