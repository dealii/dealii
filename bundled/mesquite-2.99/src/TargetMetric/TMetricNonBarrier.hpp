/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TMetricNonBarrier.hpp
 *  \brief 
 *  \author Boyd Tidwell 
 */

#ifndef MSQ_T_METRIC_NON_BARRIER_HPP
#define MSQ_T_METRIC_NON_BARRIER_HPP

#include "Mesquite.hpp"
#include "TMetric.hpp"
#include <string>

namespace MESQUITE_NS {

class MsqError;
template <unsigned R, unsigned C> class MsqMatrix;
  
class TMetricNonBarrier : public TMetric 
{
public:

  MESQUITE_EXPORT virtual
  ~TMetricNonBarrier();

  MESQUITE_EXPORT virtual
  std::string get_name() const {return "TMetricNonBarrier";}

  static inline bool invalid_determinant( double d )
    { return d < 1e-12; }
};

class TMetricNonBarrier2D : public TMetric
{
public:

  MESQUITE_EXPORT virtual
  ~TMetricNonBarrier2D();

    /**\brief Evaluate \f$\mu(T)\f$
     *
     * This method always returns an error for 2D-only metrics
     */
  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<3,3>& T, 
                 double& result, 
                 MsqError& err );
};

class TMetricNonBarrier3D : public TMetric
{
public:

  MESQUITE_EXPORT virtual
  ~TMetricNonBarrier3D();

    /**\brief Evaluate \f$\mu(T)\f$
     *
     * This method always returns an error for 3D-only metrics
     */
  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<2,2>& T, 
                 double& result, 
                 MsqError& err );
};


} // namespace MESQUITE_NS

#endif
