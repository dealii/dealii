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


/** \file XYRectangle.hpp
 *  \brief Define simple domain for 2D test problems
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_XY_RECTANGLE_HPP
#define MSQ_XY_RECTANGLE_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"

#include <map>

namespace MESQUITE_NS {

/**\brief Simple 2D Domain for free-smooth testing
 *
 * Define a simple bounded rectangular domain in the XY-plane.
 * Mesh vertices are classified as one of the following:
 *  - 0 DoF : On a corner of the rectangle
 *  - 1 DoF : On an edge of the rectangle
 *  - 2 DoF : In the interior of the rectangle
 */
class XYRectangle : public Mesquite::MeshDomain
{
  public:
    
    enum Plane { XY = 2, YZ = 0, ZX = 1 };
    
    
      /**\brief Define rectangular domain
       *
       *\param w Width of rectangle (X-range)
       *\param h Height of rectangle (Y-range)
       *\param x Minimum X coordinate of rectangle
       *\param y Minimum Y coordinate of rectangle
       *\param z Minimum Z coordinate of rectangle
       *\param plane Which plane (default is XY).
       *
       * Create w x h rectangle with, if plane is XY: X range of [x, x+w] and
       * Y range of [y, y+h].
       */
    MESQUITE_EXPORT
    XYRectangle( double w, double h, double x = 0, double y = 0, double z = 0, Plane plane = XY );
  
      /**\brief Classify mesh vertices against domain
       *
       * Figure out which input mesh vertices like on corners
       * or edges of the domain.  Will fail if any vertex is outside 
       * of the rectangle.
       */
    MESQUITE_EXPORT
    void setup( Mesquite::Mesh* mesh, Mesquite::MsqError& err );

    MESQUITE_EXPORT
    void snap_to( Mesquite::Mesh::VertexHandle entity_handle,
                  Mesquite::Vector3D &coordinate) const;
  
    MESQUITE_EXPORT
    void vertex_normal_at( Mesquite::Mesh::VertexHandle entity_handle,
                           Mesquite::Vector3D &coordinate) const;

    MESQUITE_EXPORT
    void element_normal_at( Mesquite::Mesh::ElementHandle entity_handle,
                            Mesquite::Vector3D &coordinate) const;

    MESQUITE_EXPORT
    void vertex_normal_at( const Mesquite::Mesh::VertexHandle* handles,
                           Mesquite::Vector3D coordinates[],
                           unsigned count,
                           Mesquite::MsqError& err ) const;
    
    MESQUITE_EXPORT
    void closest_point( Mesquite:: Mesh::VertexHandle handle,
                        const Mesquite::Vector3D& position,
                        Mesquite::Vector3D& closest,
                        Mesquite::Vector3D& normal,
                        Mesquite::MsqError& err ) const;
    
    MESQUITE_EXPORT
    void domain_DoF( const Mesquite::Mesh::VertexHandle* handle_array,
                     unsigned short* dof_array,
                     size_t num_handles,
                     Mesquite::MsqError& err ) const;
 
  private:
    double minCoords[3], maxCoords[3]; //!< corner coords
    const int normalDir, widthDir, heightDir;
  
    //! Single constraint on a vertex (reduces degrees of freedom by 1)
    struct VertexConstraint {
      enum Constraint { XC = 0, YC = 1, ZC = 2 };
      VertexConstraint(int a, double c)  : axis((Constraint)a), coord(c) {}
      Constraint axis;
      double coord;
    };
    
    //! Map vertex handles to constraints
    typedef std::multimap<Mesh::VertexHandle,VertexConstraint> constraint_t;
    constraint_t mConstraints;
};



} // namespace Mesquite

#endif
