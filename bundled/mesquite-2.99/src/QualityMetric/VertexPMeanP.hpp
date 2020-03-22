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


/** \file VertexPMeanP.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_VERTEX_PMEAN_P_HPP
#define MSQ_VERTEX_PMEAN_P_HPP

#include "Mesquite.hpp"
#include "PMeanPMetric.hpp"
#include "VertexQM.hpp"

namespace MESQUITE_NS {

class ElemSampleQM;

class VertexPMeanP : public VertexQM, public PMeanPMetric
{
public:

  MESQUITE_EXPORT VertexPMeanP( double p, ElemSampleQM* metric );
  
  MESQUITE_EXPORT virtual ~VertexPMeanP();
  
  MESQUITE_EXPORT ElemSampleQM* get_quality_metric() const 
    { return mMetric; }
  
  MESQUITE_EXPORT virtual std::string get_name() const;

  MESQUITE_EXPORT virtual int get_negate_flag() const;

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

  MESQUITE_EXPORT virtual
  bool evaluate_with_gradient( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_Hessian( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 std::vector<Matrix3D>& Hessian,
                 MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_Hessian_diagonal( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 std::vector<SymMatrix3D>& Hessian_diagonal,
                 MsqError& err );
private:

  ElemSampleQM* mMetric;
  mutable std::vector<size_t> mHandles;
};
  



} // namespace Mesquite

#endif
