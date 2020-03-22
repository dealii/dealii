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


/** \file EdgeLengthMetric.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_EDGE_LENGTH_METRIC_HPP
#define MSQ_EDGE_LENGTH_METRIC_HPP

#include "EdgeQM.hpp"

namespace MESQUITE_NS {

class EdgeLengthMetric : public EdgeQM
{
public:
  MESQUITE_EXPORT EdgeLengthMetric( double gamma = 0.0 );

  MESQUITE_EXPORT virtual std::string get_name() const;

   //! 1 if metric should be minimized, -1 if metric should be maximized.
  MESQUITE_EXPORT virtual int get_negate_flag() const;

  MESQUITE_EXPORT virtual
  bool evaluate( PatchData& pd, 
                 size_t handle, 
                 double& value, 
                 MsqError& err );
  MESQUITE_EXPORT virtual
  bool evaluate_with_gradient( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 MsqError& err );

private:
  double mGamma;
};


} // namespace MESQUITE_NS

#endif
