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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov,
    kraftche@cae.wisc.edu         
   
  ***************************************************************** */
/*!
  \file   XYPlanarDomain.hpp
  \brief  


  \author Thomas Leurent
  \date   2002-01-17
  \author Boyd Tidwell
  \date   2012-03-07
*/


#ifndef MSQ_XYPLANAR_DOMAIN_HPP
#define MSQ_XYPLANAR_DOMAIN_HPP

#include "MeshInterface.hpp"
#include "Vector3D.hpp"

namespace MESQUITE_NS
{
  /*! \class XYPlanarDomain
       This is a template for a planar domain that exists in the xy-plane (z=0).
       It will provide the normal information necessary for surface mesh optimization.
    */
  class MESQUITE_EXPORT XYPlanarDomain : public Mesquite::MeshDomain
  {
  public:
   
    XYPlanarDomain() {};
    
    virtual ~XYPlanarDomain();
      
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
                             size_t num_vertices,
                             MsqError& err ) const;

  };
}

#endif
