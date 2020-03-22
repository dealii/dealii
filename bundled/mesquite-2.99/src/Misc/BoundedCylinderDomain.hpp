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

#ifndef MSQ_BOUNDED_CYLINDER_DOMAIN_HPP
#define MSQ_BOUNDED_CYLINDER_DOMAIN_HPP

#include "CylinderDomain.hpp"
#include "Vector3D.hpp"

#include <list>

namespace MESQUITE_NS {

class BoundedCylinderDomain : public CylinderDomain
{
  public:
      /**
       *\param radius         - Radius of the cylinder
       *\param axis_direction - Vector defining the direction of the axis
       *\param axis_point     - A point through which the axis passes.
       */
    inline BoundedCylinderDomain( double radius,
                           Vector3D axis_direction = Vector3D(0,0,1), 
                           Vector3D axis_point = Vector3D(0,0,0) )
      : CylinderDomain( radius, axis_direction, axis_point ) {}

    
    virtual void domain_DoF( const Mesh::VertexHandle* handle_array,
                             unsigned short* dof_array,
                             size_t count,
                             MsqError& err ) const;
                             
      /**\brief define a circular curve bounding the cylinder
       *\param distance Location on cylinder at which to create
       *                  a circular curve, specified as the distance
       *                  along the cylinder axis from axis_point
       *                  specified to the constructor.
       *\param handles  A list of handles which are to be constidered
       *                  bound to the curve.
       */
    void create_curve( double distance, 
                       const std::vector<Mesh::VertexHandle>& handles );
                             
      /**\brief define a circular curve bounding the cylinder
       *\param distance  Location on cylinder at which to create
       *                  a circular curve, specified as the distance
       *                  along the cylinder axis from axis_point
       *                  specified to the constructor.
       *\param mesh      All vertices in this mesh within the specified
       *                  tolerance of the curve will be considered bound
       *                  to the curve.
       *\param tolerance The distance tolerance to use.
       */
    void create_curve( double distance, 
                       Mesh* mesh,
                       double tolerance,
                       MsqError& err );

  protected:
    
    void evaluate( double t,
                   const Vector3D& point,
                   Vector3D& closest,
                   Vector3D& normal ) const;
    
    virtual void evaluate( Mesh::VertexHandle handle,
                           const Vector3D& point,
                           Vector3D& closest,
                           Vector3D& normal ) const;
  
    bool find_curve( Mesh::VertexHandle handle, double& t ) const;
  
  private:
  
    struct Curve {
      double t;
      std::vector<Mesh::EntityHandle> handles;
    };
    
    std::list<Curve> curveList;
};


} // namespace Mesquite

#endif
