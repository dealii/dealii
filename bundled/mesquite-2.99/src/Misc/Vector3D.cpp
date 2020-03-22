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
#include "Vector3D.hpp"
#include "MsqError.hpp"

#include <iostream>
#include <math.h>

#ifdef MSQ_HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

namespace MESQUITE_NS {

std::ostream& operator<<(std::ostream &s, const Mesquite::Vector3D &v)
{
    return s << v[0] << ' ' << v[1] << ' ' << v[2];
}

double Vector3D::interior_angle(const Vector3D &lhs,
                                const Vector3D &rhs,
                                MsqError& err)
  {
    double len1 = lhs.length();
    double len2 = rhs.length();
    double angle_cos = (lhs % rhs)/(len1 * len2);
    if (!finite( angle_cos ))
    {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
      return 0.0;
    }
    
    // Adjust the cosine if slightly out of range
    if ((angle_cos > 1.0) && (angle_cos < 1.0001))
      {
        angle_cos = 1.0;
      }
    else if (angle_cos < -1.0 && angle_cos > -1.0001)
      {
        angle_cos = -1.0;
      }
    
    return std::acos(angle_cos);
  }

}
