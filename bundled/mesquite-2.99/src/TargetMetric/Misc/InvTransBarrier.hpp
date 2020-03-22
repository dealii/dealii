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


/** \file InvTransBarrier.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_INV_TRANS_BARRIER_HPP
#define MSQ_INV_TRANS_BARRIER_HPP

#include "Mesquite.hpp"
#include "TMetricBarrier.hpp"

namespace MESQUITE_NS {

/** Make a non-barrier metric into a barrier metric by passing it T^-t */
class InvTransBarrier : public TMetricBarrier
{
  public:
  InvTransBarrier( TMetric* metric ) : metricPtr(metric) {}

  MESQUITE_EXPORT virtual
  ~InvTransBarrier();

  MESQUITE_EXPORT virtual
  std::string get_name() const;
  
  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<2,2>& T, double& result, MsqError& err );
  
  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<3,3>& T, double& result, MsqError& err );

  private:
  TMetric* metricPtr;
};

} // namespace Mesquite

#endif
