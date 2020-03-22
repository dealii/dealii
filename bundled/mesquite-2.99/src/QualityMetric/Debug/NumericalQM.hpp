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


/** \file NumericalQM.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_NUMERICAL_QM_HPP
#define MSQ_NUMERICAL_QM_HPP

#include "Mesquite.hpp"
#include "QualityMetric.hpp"

namespace MESQUITE_NS {

/**\brief Use finite difference rather than analytical derivative calculations.
 *
 * This class acts as a decorator for (wrapper around) an existing 
 * QualityMetric.  It bypasses any code provided by the underlying
 * metric for calculation of derivatives, forcing the use of finite
 * difference approximation of the derivatives.
 */
class NumericalQM : public QualityMetric
{
  public:
    
    /**
     *\param real_metric  The actual quality metric.
     *\param numerical_gradient Use finite difference to calculate first derivatives
     *\param numerical_hessian  Use finite difference to calculate second derivatives
     */
    NumericalQM( QualityMetric* real_metric,
                 bool numerical_gradient = true,
                 bool numerical_hessian = true );

    MetricType get_metric_type() const;

    std::string get_name() const;

    int get_negate_flag() const;

    void get_evaluations( PatchData& pd, 
                          std::vector<size_t>& handles, 
                          bool free_vertices_only,
                          MsqError& err );

    bool evaluate( PatchData& pd, 
                   size_t handle, 
                   double& value, 
                   MsqError& err );

    bool evaluate_with_indices( PatchData& pd,
                   size_t handle,
                   double& value,
                   std::vector<size_t>& indices,
                   MsqError& err );

    bool evaluate_with_gradient( PatchData& pd,
                   size_t handle,
                   double& value,
                   std::vector<size_t>& indices,
                   std::vector<Vector3D>& gradient,
                   MsqError& err );

    bool evaluate_with_Hessian_diagonal( PatchData& pd,
                   size_t handle,
                   double& value,
                   std::vector<size_t>& indices,
                   std::vector<Vector3D>& gradient,
                   std::vector<SymMatrix3D>& Hessian_diagonal,
                   MsqError& err );

    bool evaluate_with_Hessian( PatchData& pd,
                   size_t handle,
                   double& value,
                   std::vector<size_t>& indices,
                   std::vector<Vector3D>& gradient,
                   std::vector<Matrix3D>& Hessian,
                   MsqError& err );

  private:
    QualityMetric* realMetric;
    bool numericGrad, numericHess;

};


} // namespace MESQUITE_NS

#endif
