/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

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


/** \file AddQualityMetric.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_ADD_QUALITY_METRIC_HPP
#define MSQ_ADD_QUALITY_METRIC_HPP

#include "Mesquite.hpp"
#include "QualityMetric.hpp"

namespace MESQUITE_NS {

class AddQualityMetric : public QualityMetric
{
public:
  
  AddQualityMetric( QualityMetric* qm1, QualityMetric* qm2, MsqError& err );
  
  ~AddQualityMetric();
  
  MetricType get_metric_type() const;
  
  std::string get_name() const;
  
  int get_negate_flag() const;
  
  QualityMetric* get_first_metric() const { return &metric1; }
  QualityMetric* get_second_metric() const { return &metric2; }
  
  virtual
  void get_evaluations( PatchData& pd, 
                        std::vector<size_t>& handles, 
                        bool free_vertices_only,
                        MsqError& err );

  virtual
  bool evaluate( PatchData& pd, 
                 size_t handle, 
                 double& value, 
                 MsqError& err );


  virtual
  bool evaluate_with_indices( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 MsqError& err );


  virtual
  bool evaluate_with_gradient( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 MsqError& err );


  virtual
  bool evaluate_with_Hessian( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 std::vector<Matrix3D>& Hessian,
                 MsqError& err );
                 
private:
  QualityMetric &metric1, &metric2;
  mutable std::vector<size_t> mHandles;
  mutable std::vector<size_t> indices1, indices2;
  mutable std::vector<Vector3D> grad1, grad2;
  mutable std::vector<Matrix3D> Hess1, Hess2;
  
};




} // namespace Mesquite

#endif
