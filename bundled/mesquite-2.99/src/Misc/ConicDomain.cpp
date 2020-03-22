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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file ConicDomain.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ConicDomain.hpp"

namespace MESQUITE_NS {

ConicDomain::~ConicDomain() {}

void ConicDomain::evaluate( Mesh::VertexHandle,
                            const Vector3D& point,
                            Vector3D& closest,
                            Vector3D& normal ) const
{
  // translate such that cone point (mPoint) is at origin
  Vector3D pt = point - mPoint;

  // find the plane containing both the input point an the axis
  Vector3D n = mAxis * pt;
  double len = n.length();
  if (len < 1e-12) { 
      // point on axis
      // choose any plane that does't contain the axis
    if (fabs(mAxis[0]) <= fabs(mAxis[1]) && fabs(mAxis[0]) < fabs(mAxis[2]))
      n = mAxis * Vector3D(1,0,0);
    else if (fabs(mAxis[1]) <= fabs(mAxis[2]))
      n = mAxis * Vector3D(0,1,0);
    else
      n = mAxis * Vector3D(0,0,1);
  }
  else  {
    n /= len;
  }
  // Find two points that are the intersection of the plane with the
  // circular cross-section of the cone centered at mPoint
  Vector3D p1 = mRadius * (n * mAxis);
  Vector3D p2 = -p1;
  // Define two lines of intersect between the cone and the plane
  // as the two lines from the apex to each of p1 and p2.
  Vector3D apex = mHeight * mAxis;
  Vector3D v1 = p1 - apex;
  Vector3D v2 = p2 - apex;
  // Find closest point on each line to input position
  double t1 = v1 % (point - apex) / (v1 % v1);
  double t2 = v2 % (point - apex) / (v2 % v2);
  // Select the closest of the two
  p1 = apex + t1*v1;
  p2 = apex + t2*v2;
  double t;
  if ((p1 - point).length_squared() < (p2 - point).length_squared()) {
    normal = v1 * n;
    closest = p1;
    t = t1;
  }
  else {
    normal = v2 * n;
    closest = p2;
    t = t2;
  }
  // If we're below the apex (t > 0), then the normal
  // should be in the same direction as the axis.  Otherwise
  // it should be in the opposite direction.
  if (t * (normal % mAxis) < 0.0)
    normal = -normal;
  // normalize and translate
  normal *= outwardSign / normal.length();
  closest += mPoint;
}

void ConicDomain::snap_to( Mesh::VertexHandle h, Vector3D& v ) const
{
  Vector3D p(v), n;
  evaluate( h, p, v, n );
}

void ConicDomain::vertex_normal_at( Mesh::VertexHandle h, Vector3D& v ) const
{
  Vector3D p(v), l;
  evaluate( h, p, l, v );
}

void ConicDomain::element_normal_at( Mesh::ElementHandle h, Vector3D& v ) const
{
  Vector3D p(v), l;
    // NOTE: Explicitly invoke this class's evaluate method for elements.
    //       BoundedCylindarDomain overrides evaluate for vertices only.
  ConicDomain::evaluate( h, p, l, v );
}

void ConicDomain::vertex_normal_at( const Mesh::VertexHandle* h,
                                    Vector3D coords[],
                                    unsigned count, 
                                    MsqError& ) const
{
  for (unsigned i = 0; i < count; ++i)
    vertex_normal_at( h[i], coords[i] );
}
 
void ConicDomain::closest_point( Mesh::VertexHandle handle,
                                 const Vector3D& position,
                                 Vector3D& closest,
                                 Vector3D& normal,
                                 MsqError&  ) const
{
  evaluate( handle, position, closest, normal );
}

void ConicDomain::domain_DoF( const Mesh::VertexHandle* ,
                              unsigned short* dof_array,
                              size_t count,
                              MsqError&  ) const
{
  std::fill( dof_array, dof_array + count, (unsigned short)2 );
}



} // namespace Mesquite
