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


/** \file StdDevTemplate.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_STD_DEV_TEMPLATE_HPP
#define MSQ_STD_DEV_TEMPLATE_HPP

#include "Mesquite.hpp"
#include "VarianceTemplate.hpp"

namespace MESQUITE_NS {

/**\brief standard deviation template
 *
 * This class implements an objective function that is the 
 * standard deviation of the quality metric evalutations.
 */
class StdDevTemplate : public VarianceTemplate
{
  public:
  
	MESQUITE_EXPORT
    StdDevTemplate( QualityMetric* qm = 0 ) : VarianceTemplate(qm)
      {}
    
	MESQUITE_EXPORT
    virtual ~StdDevTemplate() 
      {}
    
	MESQUITE_EXPORT
    virtual bool evaluate( EvalType type, 
                           PatchData& pd,
                           double& value_out,
                           bool free,
                           MsqError& err ); 

	MESQUITE_EXPORT
    virtual bool evaluate_with_gradient( EvalType type, 
                                         PatchData& pd,
                                         double& value_out,
                                         std::vector<Vector3D>& grad_out,
                                         MsqError& err ); 

	MESQUITE_EXPORT
    virtual bool evaluate_with_Hessian_diagonal( EvalType type, 
                                        PatchData& pd,
                                        double& value_out,
                                        std::vector<Vector3D>& grad_out,
                                        std::vector<SymMatrix3D>& hess_diag_out,
                                        MsqError& err ); 

	MESQUITE_EXPORT
    virtual ObjectiveFunction* clone() const;
  
};

} // namespace Mesquite

#endif
