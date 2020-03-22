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


/** \file TQualityMetric.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_T_QUALITY_METRIC_HPP
#define MSQ_T_QUALITY_METRIC_HPP

#include "Mesquite.hpp"
#include "TMPQualityMetric.hpp"
#include "MsqMatrix.hpp"

namespace MESQUITE_NS {

class TMetric;

/**\brief Compare targets to mapping function Jacobian matrices
 *
 * A quality metric defined using 2D and 3D target metrics,
 * where the active (A) matrix compared to the target by
 * the underlying metrics is the Jacobian matrix of the
 * mapping function at a given sample point.  For surface
 * elements, A is rotated to align the normal with W, such that
 * both matrices can be reduced from 3x2 to 2x2.
 */
class TQualityMetric : public TMPQualityMetric
{
public:
  
  /** Used in tests and other templatized code */
  typedef TMetric MetricType;

  /**
   *\param tc   The target calculator 
   *\param wc   The weight calculator
   *\param target_metric The target metric
   */
  TQualityMetric( TargetCalculator* tc,
                  WeightCalculator* wc,
                  TMetric* target_metric ) 
    : TMPQualityMetric(tc,wc),
      targetMetric(target_metric)
   {}

  /**
   *\param tc   The target calculator 
   *\param target_metric The target metric
   */
  TQualityMetric( TargetCalculator* tc,
                  TMetric* target_metric ) 
    : TMPQualityMetric(tc,0),
      targetMetric(target_metric)
   {}
     
  MESQUITE_EXPORT virtual
  std::string get_name() const;
                 
  MESQUITE_EXPORT virtual
  bool evaluate_with_gradient( PatchData& pd,
                               size_t handle,
                               double& value,
                               std::vector<size_t>& indices,
                               std::vector<Vector3D>& gradient,
                               MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_Hessian_diagonal( PatchData& pd,
                                       size_t handle,
                                       double& value,
                                       std::vector<size_t>& indices,
                                       std::vector<Vector3D>& gradient,
                                       std::vector<SymMatrix3D>& Hessian_diagonal,
                                       MsqError& err );
                    
  MESQUITE_EXPORT virtual
  bool evaluate_with_Hessian( PatchData& pd,
                              size_t handle,
                              double& value,
                              std::vector<size_t>& indices,
                              std::vector<Vector3D>& gradient,
                              std::vector<Matrix3D>& Hessian,
                              MsqError& err );
  
  TMetric* get_target_metric() const { return targetMetric; }
  void set_target_metric( TMetric* m ) { targetMetric = m; }

protected:

  MESQUITE_EXPORT virtual
  bool evaluate_internal( PatchData& pd,
                 size_t handle,
                 double& value,
                 size_t* indices,
                 size_t& num_indices,
                 MsqError& err );

private:

  TMetric* targetMetric;
};

} // namespace Mesquite

#endif
