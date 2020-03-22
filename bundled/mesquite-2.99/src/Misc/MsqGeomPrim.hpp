/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2008 Sandia National Laboratories.  Developed at the
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

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file MsqGeomPrim.hpp
 *  \brief Classes for primitive geometry
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_MSQ_GEOM_PRIM_HPP
#define MSQ_MSQ_GEOM_PRIM_HPP

#include "Mesquite.hpp"
#include "Vector3D.hpp"

namespace MESQUITE_NS {


/**\brief Line in R^3 */
class MESQUITE_EXPORT MsqLine 
{
  private:
    Vector3D mPoint;
    Vector3D mDirection; // unit direction
  
  public:
  
    MsqLine( const Vector3D& point, const Vector3D& dir )
      : mPoint(point), mDirection(dir/dir.length())
    {}
    
    MsqLine()
      : mPoint(0,0,0), mDirection(0,0,0)
    {}
    
    static
    MsqLine two_point( const Vector3D& p1, const Vector3D& p2 )
      { return MsqLine( p1, p2 - p1 ); }
    
    const Vector3D& point() const { return mPoint; }
    const Vector3D& direction() const { return mDirection; }
    
      //! Get point given parameter (mPoint + param * mDirection)
    Vector3D point( double param ) const 
      { return point() + param * direction(); }
    
      //! Get parameter value for location on line closest to input point.
    double closest( const Vector3D& from_point ) const
      { return mDirection % (from_point - mPoint); }
      
    double distance( const Vector3D& from_point ) const
      { return (point(closest(from_point)) - from_point).length(); }
      
      //! Find intersection between two lines
      //!\return false if no point intersection (skew, parallel, coincident, etc.)
    bool intersect( const MsqLine& other, double& param, double epsilon ) const;
    
      //! Find parameter of closest position on this line to another line
    bool closest( const MsqLine& other, double& param ) const;
};

/**\brief Circle in R^3 */
class MESQUITE_EXPORT MsqCircle
{
  private:
    Vector3D mCenter;
    Vector3D mNormal; //!< unit normal
    double mRadius;

  public:
  
    MsqCircle()
      : mCenter(0,0,0), mNormal(0,0,0), mRadius(0)
      {}
  
    MsqCircle( const Vector3D& center, const Vector3D& normal, double radius )
      : mCenter(center), mNormal(normal/normal.length()), mRadius(radius) 
      {}
      
      //! Find circle passing through three points
      //!\return false if points are colinear, true otherwise.
    static
    bool three_point( const Vector3D& p1, 
                      const Vector3D& p2, 
                      const Vector3D& p3, 
                      MsqCircle& result );
      
      //! Find circle with center and two points
      //!\return false if points are colinear or not equidistant from center, true otherwise.
    static
    bool two_point( const Vector3D& center, 
                    const Vector3D& p1, 
                    const Vector3D& p2, 
                    MsqCircle& result );
                                         
    const Vector3D& center() const { return mCenter; }
    const Vector3D& normal() const { return mNormal; }
    double radius() const          { return mRadius; }
    
      //! Get arbitrary radial vector (vector in plane with length equal to radius)
    Vector3D radial_vector() const;
    
      //! Find closest point on circle to input position.  Returns
      //! arbitrary location on circle if point is at center.
    Vector3D closest( const Vector3D& point ) const;
    
      //! Find closest point on circle to input position, and tangent
      //! at that point.  Fails and returns false if point is at center
      //! of circle.
    bool closest( const Vector3D& point,
                  Vector3D& point_on_circle,
                  Vector3D& tangent_at_point ) const;
};


/**\brief Plane */
class MESQUITE_EXPORT MsqPlane
{
  private:
    Vector3D mNormal; //!< unit normal
    double mCoeff;
  
  public:
  
    MsqPlane( const Vector3D& normal, double coeff );
    MsqPlane( const Vector3D& normal, const Vector3D& point );
      //! ax+by+cz+d=0
    MsqPlane( double a, double b, double c, double d );
    
      //! get unit normal vector for plane
    const Vector3D& normal() const { return mNormal; }
    
      //! get coefficient term for plane
    double coefficient() const { return mCoeff; }
    
      //! Get point on plane closest to origin.
    Vector3D point() const { return -mCoeff * mNormal; }
    
      //! Get distance from point to plane
    double distance( const Vector3D& point ) const
      { return fabs( normal() % point + coefficient() ); }
      
      //! Get closest location on plane to input position
    Vector3D closest( const Vector3D& point ) const
      { return point - normal() * (normal() % point + coefficient()); }
    
      //! Find intersection of this plane with another.
      //!\return false if parallel or invalid plane, true otherwise
    bool intersect( const MsqPlane& plane, MsqLine& result ) const;

      //! Find point of itersection with a line
      //!\return false if parallel, true otherwise
    bool intersect( const MsqLine& line, double& line_param ) const;
};

/**\brief Sphere */
class MESQUITE_EXPORT MsqSphere
{
  private:
    Vector3D mCenter;
    double mRadius;
  
  public:
  
    MsqSphere( const Vector3D& center, double radius )
      : mCenter(center), mRadius(radius) {}
    
    const Vector3D& center() const { return mCenter; }
    double radius() const { return mRadius; }
    
      //! Get location on sphere closest to input location.
      //! Returns arbitrary point if input location is at center.
    Vector3D closest( const Vector3D& point ) const;
    
      
      //! Get location on sphere closest to input location, and
      //! normal at closest point. 
      //!\return false if input point is at center, true otherwise.
    bool closest( const Vector3D& point, 
                  Vector3D& point_on_sphere,
                  Vector3D& normal_at_point ) const;
                  
    double distance( const Vector3D& point ) const
      { return fabs( (center() - point).length() - radius() ); }
    
      //! Intersect plane with sphere
      //!\return true if intersection exists, false otherwise.
    bool intersect( const MsqPlane& plane, MsqCircle& result ) const;
                  
      //! Intersect spheres
      //!\return true if intersection exists, false otherwise.
    bool intersect( const MsqSphere& sphere, MsqCircle& result ) const;
};


} // namespace Mesquite

#endif
