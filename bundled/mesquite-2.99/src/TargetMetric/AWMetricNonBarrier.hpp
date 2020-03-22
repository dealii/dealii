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


/** \file AWMetricNonBarrier.hpp
 *  \brief 
 *  \author Boyd Tidwell 
 */

#ifndef MSQ_AW_METRIC_NON_BARRIER_HPP
#define MSQ_AW_METRIC_NON_BARRIER_HPP

#include "Mesquite.hpp"
#include "AWMetric.hpp"
#include <string>

namespace MESQUITE_NS {

class MsqError;
template <unsigned R, unsigned C> class MsqMatrix;

/**\brief The parent class for all AWMetricNonBarrier sub-classes 
 *
 */
class AWMetricNonBarrier : public AWMetric
{
public:
  MESQUITE_EXPORT virtual
  ~AWMetricNonBarrier();

  MESQUITE_EXPORT virtual
  std::string get_name() const {return "AWMetricNonBarrier";}

  static inline bool invalid_determinant( double d )
    { return d < 1e-12; }
};

class AWMetricNonBarrier2D : public AWMetricNonBarrier
{
public:

  MESQUITE_EXPORT virtual
  ~AWMetricNonBarrier2D();

    /**\brief Evaluate \f$\mu(A,W)\f$
     *
     * This method always returns an error for 2D-only metrics
     */
  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<3,3>& A, 
                 const MsqMatrix<3,3>& W, 
                 double& result, 
                 MsqError& err );
};

class AWMetricNonBarrier3D : public AWMetricNonBarrier
{
public:

  MESQUITE_EXPORT virtual
  ~AWMetricNonBarrier3D();

    /**\brief Evaluate \f$\mu(A,W)\f$
     *
     * This method always returns an error for 3D-only metrics
     */
  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<2,2>& A, 
                 const MsqMatrix<2,2>& W, 
                 double& result, 
                 MsqError& err );
};

} // namespace Mesquite

#endif
