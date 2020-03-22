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
#include "XYPlanarDomain.hpp"
#include "MsqError.hpp"
#include "MsqVertex.hpp"
#include "DomainUtil.hpp"

#include <algorithm>

Mesquite::XYPlanarDomain::~XYPlanarDomain() {}

void Mesquite::XYPlanarDomain::snap_to(Mesquite::Mesh::VertexHandle entity_handle,
                                       Vector3D &coordinate) const
{
  coordinate[2] = 0.0;
}


void Mesquite::XYPlanarDomain::vertex_normal_at(Mesquite::Mesh::VertexHandle /*entity_handle*/,
                                                Mesquite::Vector3D &coordinate) const
{
  coordinate = Vector3D(0.0, 0.0, 1.0);
}

void Mesquite::XYPlanarDomain::element_normal_at( Mesquite::Mesh::ElementHandle /*entity_handle*/,
                                                  Mesquite::Vector3D &coordinate) const
{
  coordinate = Vector3D(0.0, 0.0, 1.0);
}


void Mesquite::XYPlanarDomain::vertex_normal_at(const Mesquite::Mesh::VertexHandle* ,
                                                Vector3D coords[],
                                                unsigned count,
                                                Mesquite::MsqError& ) const
{
  for (unsigned i = 0; i < count; ++i)
    coords[i] = Vector3D(0.0, 0.0, 1.0);
}

void Mesquite::XYPlanarDomain::closest_point( Mesquite::Mesh::VertexHandle ,
                                              const Mesquite::Vector3D& position,
                                              Mesquite::Vector3D& closest,
                                              Mesquite::Vector3D& normal,
                                              Mesquite::MsqError& ) const
{
  closest = Vector3D(position[0], position[1], 0.0);
}

void Mesquite::XYPlanarDomain::domain_DoF( const Mesh::VertexHandle* ,
                                         unsigned short* dof_array,
                                         size_t num_vertices,
                                         MsqError&  ) const
{
  std::fill( dof_array, dof_array + num_vertices, 2 );
}


    
