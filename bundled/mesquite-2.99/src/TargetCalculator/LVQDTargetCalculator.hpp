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


/** \file LVQDTargetCalculator.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_LVQD_TARGET_CALCULATOR_HPP
#define MSQ_LVQD_TARGET_CALCULATOR_HPP

#include "Mesquite.hpp"
#include "TargetCalculator.hpp"

namespace MESQUITE_NS {

/** \brief Construct target matrices from factors of guide matrices
 *
 * Given a set of guide matrices, extract certain properties
 * from them and combine those properties to produce a target
 * matrix.
 * \f$\mathbf{W}=\lambda(\mathbf{W_\lambda}) \times \mathbf{V}(\mathbf{W_V}) \times \mathbf{Q}(\mathbf{W_Q}) \times \mathbf{\Delta}(\mathbf{W_\Delta})\f$
 * - \f$\lambda\f$ - scalar - size component of \f$\mathbf{W_\lambda}\f$
 * - \f$\mathbf{V}\f$ - 3x3 or 3x2 - orientation component of \f$\mathbf{W_V}\f$
 * - \f$\mathbf{Q}\f$ - 3x3 or 2x2 - shape component of \f$\mathbf{W_Q}\f$
 * - \f$\mathbf{\Delta}\f$ - 3x3 or 2x2 - aspect ratio component of \f$\mathbf{W_\Delta}\f$
 */
class LVQDTargetCalculator : public TargetCalculator
{
public:
  LVQDTargetCalculator( TargetCalculator* lambda_source,
                        TargetCalculator* V_source,
                        TargetCalculator* Q_source,
                        TargetCalculator* delta_source );

  ~LVQDTargetCalculator();
  
  bool get_3D_target( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<3,3>& W_out,
                      MsqError& err );

  bool get_surface_target( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<3,2>& W_out,
                      MsqError& err );

  bool get_2D_target( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<2,2>& W_out,
                      MsqError& err );
  
  bool have_surface_orient() const;
  
private:
  int add_source( TargetCalculator* source );
  
  bool evaluate_guide_2D( PatchData& pd, 
                          size_t element,
                          Sample sample,
                          int idx,
                          double& lambda, 
                          MsqMatrix<3,2>& V,
                          MsqMatrix<2,2>& Q,
                          MsqMatrix<2,2>& delta,
                          MsqError& err );

  TargetCalculator* uniqueGuides[4]; //!< Up to 4 unique target sources
  int numUniqueGuides;               //!< Number of unique target sources
  int lambdaIdx, vIdx, qIdx, deltaIdx; //!< Index into uniqueGuides for respective source,
                                       //!< or -1 for no source.
};

} // namespace Mesquite

#endif
