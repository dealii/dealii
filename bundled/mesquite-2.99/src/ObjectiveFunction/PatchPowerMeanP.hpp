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


/** \file PatchPowerMeanP.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_PATCH_POWER_MEAN_P_HPP
#define MSQ_PATCH_POWER_MEAN_P_HPP

#include "Mesquite.hpp"
#include "PMeanPTemplate.hpp"
#include "Exponent.hpp"
#include "Matrix3D.hpp"

namespace MESQUITE_NS {

/**\brief Objective function: p-mean^p of p-mean^p of patch metric values
 */
class PatchPowerMeanP : public PMeanPTemplate
{
  public:
  
      /**
       *\param power   The exponent to use for the power-mean
       *\param qm      The quality metric.
       */
    MESQUITE_EXPORT 
    PatchPowerMeanP( double power, QualityMetric* qm = 0 ) 
      : PMeanPTemplate(power, qm) {}
    
      /**\brief copy constructor 
       *
       * Define a copy constructor because the compiler-provided 
       * default one would also copy the temporary arrays, which
       * would be a waste of time.
       */
    MESQUITE_EXPORT
    PatchPowerMeanP( const PatchPowerMeanP& copy )
      : PMeanPTemplate( copy ) {}
    
    MESQUITE_EXPORT
    virtual ~PatchPowerMeanP() 
      {}

    MESQUITE_EXPORT virtual 
    bool initialize_block_coordinate_descent( Mesh* mesh, 
                                              MeshDomain* domain, 
                                              const Settings* settings,
                                              PatchSet* user_set,
                                              MsqError& err );
    
    MESQUITE_EXPORT virtual 
    bool evaluate( EvalType type, 
                   PatchData& pd,
                   double& value_out,
                   bool free,
                   MsqError& err ); 

    MESQUITE_EXPORT virtual
    bool evaluate_with_gradient( EvalType type, 
                                 PatchData& pd,
                                 double& value_out,
                                 std::vector<Vector3D>& grad_out,
                                 MsqError& err ); 
    
    MESQUITE_EXPORT virtual 
    bool evaluate_with_Hessian( EvalType type, 
                                PatchData& pd,
                                double& value_out,
                                std::vector<Vector3D>& grad_out,
                                MsqHessian& Hessian_out,
                                MsqError& err ); 

    MESQUITE_EXPORT virtual ObjectiveFunction* clone() const;
};

} // namespace Mesquite

#endif
