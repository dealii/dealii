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


/** \file CurveDomain.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_CURVE_DOMAIN_HPP
#define MSQ_CURVE_DOMAIN_HPP

#include "Mesquite.hpp"

namespace MESQUITE_NS {




/**\brief Domain used for optional curve smoother 
 *
 * Provide necessary functionality for interogating a topologically
 * one-dimensional domain (a curve).
 */
class CurveDomain 
{
public:
  /**\brief Measure arc length along curve
   *
   * Get the length along a curve between two points on the curve.
   */
  virtual double arc_length( const double position1[3],
                             const double position2[3],
                             MsqError& err ) = 0;
  /**\brief Get a position on the curve given an arc length
   *
   * Measure a specified arc length along a curve.
   *\param from_here Point on curve from which to measure
   *\param length    Length along curve to move
   *\param result_point Output.  A point on the curve \c length units
   *                 from the point \c from_here in the forward
   *                 direction along the curve.
   */
  virtual void position_from_length( const double from_here[3],
                                     double length,
                                     double result_point[3],
                                     MsqError& err ) = 0;
};



} // namespace MESQUITE_NS

#endif
