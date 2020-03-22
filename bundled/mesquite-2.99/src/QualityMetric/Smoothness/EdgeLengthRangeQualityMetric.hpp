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

/*! \file EdgeLengthRangeQualityMetric.hpp

Header file for the Mesquite::EdgeLengthRangeQualityMetric class

  \author Michael Brewer
  \date   2002-06-13
 */


#ifndef EdgeLengthRangeQualityMetric_hpp
#define EdgeLengthRangeQualityMetric_hpp

#include "Mesquite.hpp"
#include "AveragingQM.hpp"
#include "VertexQM.hpp"

namespace MESQUITE_NS
{
     /*! \class EdgeLengthRangeQualityMetric
       \brief Computes the edge length range metric for a given vertex.
       
        EdgeLengthRangeQualityMetric is a vertex based metric which computes
        the lengths of the edges connected to a given vertex and then
        uses those values to form a metric.  The metric is created using
        two doubles A and B.  The value of the metric is zero (ideal) if
        the edge lengths are in the range [A,B].  Otherwise, the
        metric value is some positive number.  For a given vertex,
        v_i, with connected edges of lengths l_j for j=1...k, the metric
        value is the average (where the default average type is SUM) of
        u_j = ( | l_j - A | - (l_j - A) )^2 + ( | B - l_j | - (B - l_j) )^2.
     */
  class EdgeLengthRangeQualityMetric : public VertexQM, public AveragingQM
  {
   public:
    
    EdgeLengthRangeQualityMetric(double low_a, double high_a);

    virtual ~EdgeLengthRangeQualityMetric();
      
    virtual std::string get_name() const;
    
    virtual int get_negate_flag() const;

    virtual
    bool evaluate( PatchData& pd, 
                   size_t vertex, 
                   double& value, 
                   MsqError& err );

    virtual
    bool evaluate_with_indices( PatchData& pd,
                   size_t vertex,
                   double& value,
                   std::vector<size_t>& indices,
                   MsqError& err );
  
   private:
   
    bool evaluate_common( PatchData& pd, 
                          size_t vertex, 
                          double& value, 
                          std::vector<size_t>& vertices,
                          MsqError& err );

    double highVal;
    double lowVal;
    
  };


} //namespace


#endif // EdgeLengthRangeQualityMetric_hpp


