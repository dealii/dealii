/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2005 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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

    kraftche@cae.wisc.edu    

  ***************************************************************** */

#ifndef MSQ_CYLINDER_DOMAIN_HPP
#define MSQ_CYLINDER_DOMAIN_HPP


#include "Mesquite.hpp"
#include "MeshInterface.hpp"
#include "Vector3D.hpp"

namespace MESQUITE_NS
{
  /*! \class CylinderDomain
       Define the geometry of an unbounded cylinder.
    */
  class MESQUITE_EXPORT CylinderDomain : public Mesquite::MeshDomain
  {
  public:
      /**
       *\param radius         - Radius of the cylinder
       *\param axis_direction - Vector defining the direction of the axis
       *\param axis_point     - A point through which the axis passes.
       */
    inline CylinderDomain( double radius,
                           Vector3D axis_direction = Vector3D(0,0,1), 
                           Vector3D axis_point = Vector3D(0,0,0),
                           bool outward_normal = true )
      : mAxis( axis_direction / axis_direction.length() ),
        mCenter( axis_point ),
        mRadius( radius ),
        outwardSign( outward_normal ? 1.0 : -1.0 )
      { }
    
    virtual ~CylinderDomain();

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
    const Vector3D& center() const { return mCenter; }
    double radius() const { return mRadius; }
    
  protected:
    
    virtual void evaluate( Mesh::VertexHandle handle,
                           const Vector3D& point,
                           Vector3D& closest,
                           Vector3D& normal ) const;
    
  private:
  
    Vector3D mAxis;
    Vector3D mCenter;
    double mRadius;
    double outwardSign;
  };


} // namespace Mesquite

#endif
