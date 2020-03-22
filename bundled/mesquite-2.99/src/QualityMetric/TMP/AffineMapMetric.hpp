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
 
    (2007) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file AffineMapMetric.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_AFFINE_MAP_METRIC_HPP
#define MSQ_AFFINE_MAP_METRIC_HPP

#include "Mesquite.hpp"
#include "ElemSampleQM.hpp"

namespace MESQUITE_NS {

class TargetCalculator;
class WeightCalculator;
class TMetric;
template <unsigned R, unsigned C> class MsqMatrix;

/**\brief Compare targets to affine map to ideal element.
 *
 * A quality metric defined using 2D and 3D target metrics,
 * where the active (A) matrix is an affine map to the ideal,
 * unit-side element.  
 */
class AffineMapMetric : public ElemSampleQM
{
public:

  MESQUITE_EXPORT
  AffineMapMetric( TargetCalculator* tc,
                   WeightCalculator* wc,
                   TMetric* target_metric );

  MESQUITE_EXPORT
  AffineMapMetric( TargetCalculator* tc,
                   TMetric* target_metric );
     
  MESQUITE_EXPORT virtual
  std::string get_name() const;
  
  MESQUITE_EXPORT virtual 
  int get_negate_flag() const;
  
  MESQUITE_EXPORT virtual
  void get_evaluations( PatchData& pd, 
                        std::vector<size_t>& handles, 
                        bool free_vertices_only,
                        MsqError& err );
  
  MESQUITE_EXPORT virtual 
  void get_element_evaluations( PatchData& pd, size_t elem_index,
                                std::vector<size_t>& handles,
                                MsqError& err );
  
  MESQUITE_EXPORT virtual
  bool evaluate( PatchData& pd, 
                 size_t handle, 
                 double& value, 
                 MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_indices( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 MsqError& err );
    
  void set_target_calculator( TargetCalculator* tc ) { targetCalc = tc; }
  void set_weight_calculator( WeightCalculator* wc ) { weightCalc = wc; }
  TargetCalculator* get_target_calculator() const { return targetCalc; }
  WeightCalculator* get_weight_calculator() const { return weightCalc; }
  
  TMetric* get_target_metric() const { return targetMetric; }
  void set_target_metric( TMetric* m ) { targetMetric = m; }
  
private:
  TargetCalculator* targetCalc;
  WeightCalculator* weightCalc;
  TMetric* targetMetric;
};

} // namespace Mesquite

#endif
