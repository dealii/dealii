/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ScalarMultiplyQualityMetric.hpp
\brief
Header file for the Mesquite::ScalarMultiplyQualityMetric class

  \author Todd Munson
  \date   2004-12-21
 */


#ifndef ScalarMultiplyQualityMetric_hpp
#define ScalarMultiplyQualityMetric_hpp

#include "Mesquite.hpp"
#include "QualityMetric.hpp"

namespace MESQUITE_NS
{
     /*! \class ScalarMultiplyQualityMetric
       \brief Multiplies quality metric value by a number (a double).
     */
class ScalarMultiplyQualityMetric : public QualityMetric
{
public:

  ScalarMultiplyQualityMetric( QualityMetric* metric, double scale )
    : mMetric(metric), 
      mScale( scale ) 
    {}
  
  ~ScalarMultiplyQualityMetric() {}
  
  MetricType get_metric_type() const
    { return mMetric->get_metric_type(); }
  
  std::string get_name() const;
  
  int get_negate_flag() const
    { return mMetric->get_negate_flag(); }
  
  void get_evaluations( PatchData& pd, 
                        std::vector<size_t>& handles, 
                        bool free_vertices_only,
                        MsqError& err );
  
  bool evaluate( PatchData& pd, size_t handle, double& value, MsqError& err );
  
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
  
  bool evaluate_with_Hessian( PatchData& pd,
                              size_t handle,
                              double& value,
                              std::vector<size_t>& indices,
                              std::vector<Vector3D>& gradient,
                              std::vector<Matrix3D>& Hessian,
                              MsqError& err );

private:
  
  QualityMetric* mMetric;
  double mScale;
};

} //namespace


#endif // ScalarMultiplyQualityMetric_hpp


