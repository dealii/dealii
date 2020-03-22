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


/** \file ConicDomain.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_CONIC_DOMAIN_HPP
#define MSQ_CONIC_DOMAIN_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"
#include "Vector3D.hpp"

namespace MESQUITE_NS {

/*! \class ConicDomain
     Define the geometry of an unbounded cone with circular cross-section.
  */
class MESQUITE_EXPORT ConicDomain : public Mesquite::MeshDomain
{
public:
    /**
     *\param radius_at_point    The radius of the cone at axis_point
     *\param height_to_apex     The distance in the direction of 
     *                          axis_direction from axis_point to the apex
     *\param axis_direction     Vector defining the direction of the axis
     *\param axis_point         A point through which the axis passes.
     *\NOTE Cone is not bounded at apex.  It extends infinitely in both
     *      directions.
     */
  inline ConicDomain( double radius_at_point,
                      double height_to_apex,
                      Vector3D axis_direction = Vector3D(0,0,1), 
                      Vector3D axis_point = Vector3D(0,0,0),
                      bool outward_normal = true)
    : mAxis( axis_direction / axis_direction.length() ),
      mPoint( axis_point ),
      mRadius( radius_at_point ),
      mHeight( height_to_apex ),
      outwardSign( outward_normal ? 1.0 : -1.0 )
    { }
    
  inline ConicDomain() {}

  virtual ~ConicDomain();

  virtual void snap_to(Mesh::VertexHandle entity_handle,
                       Vector3D &coordinate) const;

  virtual void vertex_normal_at(Mesh::VertexHandle entity_handle,
                                Vector3D &coordinate) const;


  virtual void element_normal_at(Mesh::ElementHandle entity_handle,
                                 Vector3D &coordinate) const;


  virtual void vertex_normal_at(const Mesh::VertexHandle* handle,
                                Vector3D coords[],
                                unsigned count,
                                MsqError& err) const;

  virtual void closest_point( Mesh::VertexHandle handle,
                              const Vector3D& position,
                              Vector3D& closest,
                              Vector3D& normal,
                              MsqError& err ) const;

  virtual void domain_DoF( const Mesh::VertexHandle* handle_array,
                           unsigned short* dof_array,
                           size_t count,
                           MsqError& err ) const;


  const Vector3D& axis() const { return mAxis; }
  const Vector3D& point() const { return mPoint; }
  double point_radius() const { return mRadius; }
  double height_from_point() const { return mHeight; }

protected:

  virtual void evaluate( Mesh::VertexHandle handle,
                         const Vector3D& point,
                         Vector3D& closest,
                         Vector3D& normal ) const;

private:

  Vector3D mAxis;   //!< Direction of central axis of cone. Unit vector.
  Vector3D mPoint;  //!< A point on the axis of the cone.
  double mRadius;   //!< Radius at mPoint
  double mHeight;   //!< Distance from mPoint to apex, in direction of mAxis
  double outwardSign; //!< 1.0 if normal points away from axis, -1.0 otherwise
};



} // namespace Mesquite

#endif
