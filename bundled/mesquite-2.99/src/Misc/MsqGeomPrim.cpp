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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file MsqGeomPrim.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MsqGeomPrim.hpp"

namespace MESQUITE_NS {

bool MsqLine::intersect( const MsqLine& other, double& param, double epsilon ) const
{
  if (!closest( other, param ))
    return false;
  Vector3D p1 = point(param);
  Vector3D p2 = other.point(other.closest(p1));
  return (p1 - p2).length_squared() < epsilon*epsilon;
}
  

bool MsqLine::closest( const MsqLine& other, double& param ) const
{
  const Vector3D N = other.direction() * (direction() * other.direction());
  const double D = -(N % other.point());
  const double dot = N % direction();
  if (dot < DBL_EPSILON)
    return false; // parallel
  param = -(N % point() + D) / dot;
  return true;
}
  
bool MsqCircle::three_point( const Vector3D& p1, 
                             const Vector3D& p2, 
                             const Vector3D& p3, 
                             MsqCircle& result )
{
  Vector3D norm = (p1 - p2) * (p3 - p2);
  if (norm.length_squared() < DBL_EPSILON)
    return false;
  
  MsqLine line1( 0.5*(p1+p2), norm * (p1 - p2) );
  MsqLine line2( 0.5*(p2+p3), norm * (p3 - p2) );
  double t_xsect;
  if (!line1.closest( line2, t_xsect ))
    return false;
  
  Vector3D center = line1.point(t_xsect);
  double radius = ((center - p1).length() +
                   (center - p2).length() +
                   (center - p3).length()) / 3.0;
  result = MsqCircle( center, norm, radius );
  return true;
}
  
      
bool MsqCircle::two_point( const Vector3D& center, 
                           const Vector3D& p1, 
                           const Vector3D& p2, 
                           MsqCircle& result )
{
  Vector3D norm = (p1 - center) * (p2 - center);
  if (norm.length_squared() < DBL_EPSILON)
    return false;
  
  double radius = 0.5 * ((center - p1).length() + (center - p2).length());  
  result = MsqCircle( center, norm, radius );
  return true;
}

Vector3D MsqCircle::radial_vector() const
{
  int min_idx = 0;
  if (normal()[1] < normal()[min_idx])
    min_idx = 1;
  if (normal()[2] < normal()[min_idx])
    min_idx = 2;
  Vector3D vect(0,0,0);
  vect[min_idx] = 1;
  vect = normal() * vect;
  vect *= radius() / vect.length();
  return vect;
}


Vector3D MsqCircle::closest( const Vector3D& point ) const
{
  const Vector3D from_center = point - center();
  const Vector3D norm_proj = normal() * (normal() % from_center); // unit normal!
  const Vector3D in_plane = from_center - norm_proj;
  const double length = in_plane.length();
  if (length < DBL_EPSILON) 
    return center() + radial_vector();
  else
    return center() + in_plane * radius() / length;
}

bool MsqCircle::closest( const Vector3D& point,
                         Vector3D& result_pt,
                         Vector3D& result_tngt ) const
{
  const Vector3D from_center = point - center();
  Vector3D in_plane = from_center - (from_center % normal());
  if (in_plane.length_squared() < DBL_EPSILON)
    return false;
  
  result_pt = center() + in_plane * radius() / in_plane.length();
  result_tngt = in_plane * normal();
  return true;
}

MsqPlane::MsqPlane( const Vector3D& normal, double coeff )
{
  const double len = normal.length();
  mNormal = normal / len;
  mCoeff = coeff / len;
}

MsqPlane::MsqPlane( const Vector3D& normal, const Vector3D& point )
  : mNormal( normal / normal.length() ), mCoeff( -(mNormal % point) )
{
}

MsqPlane::MsqPlane( double a, double b, double c, double d )
  : mNormal(a,b,c), mCoeff(d)
{
  const double len = mNormal.length();
  mNormal /= len;
  mCoeff /= len;
}


bool MsqPlane::intersect( const MsqPlane& plane, MsqLine& result ) const
{
  const double dot = normal() % plane.normal();
  const double det = dot*dot - 1.0;
  if (fabs(det) < DBL_EPSILON) // parallel
    return false;
  
  const double s1 = (coefficient() - dot*plane.coefficient()) / det;
  const double s2 = (plane.coefficient() - dot*coefficient()) / det;
  result = MsqLine( s1*normal() + s2*plane.normal(), normal() * plane.normal() );
  return true;
}


bool MsqPlane::intersect( const MsqLine& line, double& result ) const
{
  const double dot = line.direction() % normal();
  if (fabs(dot) < DBL_EPSILON)
    return false;
  
  result = -(normal() % line.point() + coefficient()) / dot;
  return true;
}

Vector3D MsqSphere::closest( const Vector3D& point ) const
{
  Vector3D diff = point - center();
  double len = diff.length();
  if (len < DBL_EPSILON) {
    // pick any point
    diff = Vector3D(1,0,0);
    len = 1;
  }
  
  return center() + diff * radius() / len;
}

bool MsqSphere::closest( const Vector3D& point, 
                         Vector3D& point_on_sphere,
                         Vector3D& normal_at_point ) const
{
  normal_at_point = point - center();
  double len = normal_at_point.length();
  if (len < DBL_EPSILON) 
    return false;
  
  normal_at_point /= len;
  point_on_sphere = center() + radius() * normal_at_point;
  return true;
}
              
bool MsqSphere::intersect( const MsqPlane& plane, MsqCircle& result ) const
{
  const Vector3D plane_pt = plane.closest( center() );
  const Vector3D plane_dir = plane_pt - center();
  const double dir_sqr = plane_dir.length_squared();
  if (dir_sqr < DBL_EPSILON) { // plane passes through center of sphere
    result = MsqCircle( center(), plane.normal(), radius() );
    return true;
  }
  
  double rad_sqr = radius()*radius() - plane_dir.length_squared();
  if (rad_sqr < 0) // no intersection
    return false;
  
  result = MsqCircle( plane_pt, plane_dir, sqrt(rad_sqr) );
  return true;
}

bool MsqSphere::intersect( const MsqSphere& sphere, MsqCircle& result ) const
{
  const Vector3D d = sphere.center() - center();
  const double dist = d.length();
  if (dist > (radius() + sphere.radius()))
    return false;
  
  const double r1_sqr = radius()*radius();
  const double r2_sqr = sphere.radius() * sphere.radius();
  const double f = (d % d + r1_sqr - r2_sqr) / 2;
  //const double d1 = f / dist;
  
  const double rad = sqrt( r1_sqr - f*f/(d % d) );
  Vector3D c = center() + d * f / (d % d);
  result = MsqCircle( c, d, rad );
  return true;
}

} // namespace Mesquite
