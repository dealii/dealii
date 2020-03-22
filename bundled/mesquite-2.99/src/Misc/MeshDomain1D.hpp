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


/** \file MeshDomain1D.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_MESH_DOMAIN1D_HPP
#define MSQ_MESH_DOMAIN1D_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"
#include "MsqGeomPrim.hpp"
#include "CurveDomain.hpp"

namespace MESQUITE_NS {



class MESQUITE_EXPORT PointDomain : public MeshDomain
{
  private:
    Vector3D mGeom;
  
  public:
    
    PointDomain( const Vector3D& location ) : mGeom(location) {}

    const Vector3D& geom() const { return mGeom; }
    
     
    virtual void snap_to(Mesh::VertexHandle entity_handle,
                         Vector3D &coordinate) const;
    
    virtual void vertex_normal_at(Mesh::VertexHandle entity_handle,
                                  Vector3D &coordinate) const;
    virtual void element_normal_at(Mesh::ElementHandle entity_handle,
                                  Vector3D &coordinate) const;
                          
    virtual void vertex_normal_at( const Mesh::VertexHandle* handles,
                                   Vector3D coordinates[],
                                   unsigned count,
                                   MsqError& err ) const;
                            
    virtual void closest_point( Mesh::VertexHandle handle,
                                const Vector3D& position,
                                Vector3D& closest,
                                Vector3D& normal,
                                MsqError& err ) const;
                                
    virtual void domain_DoF( const Mesh::VertexHandle* handle_array,
                             unsigned short* dof_array,
                             size_t num_handles,
                             MsqError& err ) const;
};

class MESQUITE_EXPORT LineDomain : public MeshDomain, public CurveDomain
{
  private:
    MsqLine mGeom;
  
  public:
    
    LineDomain( const Vector3D& point, const Vector3D& dir )
      : mGeom( point, dir )
    {}
    
    LineDomain( const MsqLine& line )
      : mGeom( line )
    {}
    
    const MsqLine& geom() const { return mGeom; }
    
     
    virtual void snap_to(Mesh::VertexHandle entity_handle,
                         Vector3D &coordinate) const;
    
    virtual void vertex_normal_at(Mesh::VertexHandle entity_handle,
                                  Vector3D &coordinate) const;
    virtual void element_normal_at(Mesh::ElementHandle entity_handle,
                                  Vector3D &coordinate) const;
                          
    virtual void vertex_normal_at( const Mesh::VertexHandle* handles,
                                   Vector3D coordinates[],
                                   unsigned count,
                                   MsqError& err ) const;
                            
    virtual void closest_point( Mesh::VertexHandle handle,
                                const Vector3D& position,
                                Vector3D& closest,
                                Vector3D& normal,
                                MsqError& err ) const;
                                
    virtual void domain_DoF( const Mesh::VertexHandle* handle_array,
                             unsigned short* dof_array,
                             size_t num_handles,
                             MsqError& err ) const;
   
  virtual double arc_length( const double position1[3],
                             const double position2[3],
                             MsqError& err );

  virtual void position_from_length( const double from_here[3],
                                     double length,
                                     double result_point[3],
                                     MsqError& err );
};

class MESQUITE_EXPORT CircleDomain : public MeshDomain, CurveDomain
{
  private:
    MsqCircle mGeom;

  public:
  
    CircleDomain( const Vector3D& center, const Vector3D& normal, double radius )
      : mGeom( center, normal, radius )
      {}
  
    CircleDomain( const MsqCircle& circle )
      : mGeom( circle )
      {}
     
    const MsqCircle& geom() const { return mGeom; }
   
    
    virtual void snap_to(Mesh::VertexHandle entity_handle,
                         Vector3D &coordinate) const;
    
    virtual void vertex_normal_at(Mesh::VertexHandle entity_handle,
                                  Vector3D &coordinate) const;
    virtual void element_normal_at(Mesh::ElementHandle entity_handle,
                                  Vector3D &coordinate) const;
                          
    virtual void vertex_normal_at( const Mesh::VertexHandle* handles,
                                   Vector3D coordinates[],
                                   unsigned count,
                                   MsqError& err ) const;
                            
    virtual void closest_point( Mesh::VertexHandle handle,
                                const Vector3D& position,
                                Vector3D& closest,
                                Vector3D& normal,
                                MsqError& err ) const;
                                
    virtual void domain_DoF( const Mesh::VertexHandle* handle_array,
                             unsigned short* dof_array,
                             size_t num_handles,
                             MsqError& err ) const;
    
  virtual double arc_length( const double position1[3],
                             const double position2[3],
                             MsqError& err );

  virtual void position_from_length( const double from_here[3],
                                     double length,
                                     double result_point[3],
                                     MsqError& err );
};





} // namespace Mesquite

#endif
