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
/*! \file VertexConditionNumberQualityMetric.hpp

Header file for the Mesquite::VertexConditionNumberQualityMetric class

  \author Michael Brewer
  \date   April 14, 2003
 */


#ifndef VertexConditionNumberQualityMetric_hpp
#define VertexConditionNumberQualityMetric_hpp

#include "Mesquite.hpp"
#include "VertexQM.hpp"
#include "AveragingQM.hpp"

namespace MESQUITE_NS
{
     /*! \class VertexConditionNumberQualityMetric
       \brief Computes the condition numbers of the corner's of elements
       connected to the given vertex and then averages those values.

       The metric does not use the sample point functionality or the
       compute_weighted_jacobian.  It uses the isotropic ideal
       element.  This metric does require a feasible region, and
       the metric needs to be minimized.
     */
   class VertexConditionNumberQualityMetric : public VertexQM, public AveragingQM
   {
  public:
     VertexConditionNumberQualityMetric();
     
       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~VertexConditionNumberQualityMetric()
        {}
     
     
     virtual std::string get_name() const;

      //! 1 if metric should be minimized, -1 if metric should be maximized.
     virtual int get_negate_flag() const;

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
  };
    
   

} //namespace

#endif // VertexConditionNumberQualityMetric_hpp

