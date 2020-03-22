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


/** \file PMeanPMetric.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_PMEAN_PMETRIC_HPP
#define MSQ_PMEAN_PMETRIC_HPP

#include "Mesquite.hpp"
#include "Exponent.hpp"

#include "Vector3D.hpp"
#include "Matrix3D.hpp"

#include <vector>

namespace MESQUITE_NS {

class PatchData;
class MsqError;
class QualityMetric;

class PMeanPMetric
{
  public:
    
    PMeanPMetric( double p )
      : P(p), P1(p-1.0), P2(p-2.0) {}

    const Exponent& get_power() const { return P; }

  protected:
      
    bool average( PatchData& pd, 
                  QualityMetric* qm,
                  const std::vector<size_t>& qm_handles, 
                  double& value, 
                  MsqError& );
    
    bool average_with_indices( PatchData& pd, 
                               QualityMetric* qm,
                               const std::vector<size_t>& qm_handles,
                               double& value, 
                               std::vector<size_t>& indices,
                               MsqError& err );
    
    bool average_with_gradient( PatchData& pd, 
                                QualityMetric* qm,
                                const std::vector<size_t>& qm_handles,
                                double& value, 
                                std::vector<size_t>& indices,
                                std::vector<Vector3D>& gradient,
                                MsqError& err );
   
    bool average_with_Hessian_diagonal( PatchData& pd,
                                        QualityMetric* metric,
                                        const std::vector<size_t>& qm_handles,
                                        double& value,
                                        std::vector<size_t>& indices,
                                        std::vector<Vector3D>& gradient,
                                        std::vector<SymMatrix3D>& Hessian_diagonal,
                                        MsqError& err );

    bool average_with_Hessian( PatchData& pd, 
                               QualityMetric* metric,
                               const std::vector<size_t>& qm_handles,
                               double& value, 
                               std::vector<size_t>& indices,
                               std::vector<Vector3D>& gradient,
                               std::vector<Matrix3D>& Hessian,
                               MsqError& err );
  private:
  
    Exponent P;
    Exponent P1;
    Exponent P2;
    std::vector<size_t> mIndices, mOffsets;
    std::vector<Vector3D> mGrad;
    std::vector<Matrix3D> mHess;
    std::vector<SymMatrix3D> mDiag;
    std::vector<double> mValues;
};

} // namespace Mesquite

#endif
