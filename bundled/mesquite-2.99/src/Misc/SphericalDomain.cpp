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
#include "Mesquite.hpp"
#include "SphericalDomain.hpp"
#include "Vector3D.hpp"
#include "MsqError.hpp"
#include "MsqVertex.hpp"
#include "DomainUtil.hpp"
#include "MsqMatrix.hpp"

#ifdef MSQ_HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#include <algorithm>

Mesquite::SphericalDomain::~SphericalDomain() {}

void Mesquite::SphericalDomain::snap_to(Mesh::VertexHandle /*entity_handle*/,
                                        Vector3D &coordinate) const
{
    // Get vector center to coordinate, store in coordinate.
  coordinate -= mCenter;
    // Get distance from center of sphere
  double len = coordinate.length();
    // Scale vector to have length of radius
  coordinate *= mRadius / len;
    // If was at center, return arbitrary position on sphere
    // (all possitions are equally close)
  if (!finite(coordinate.x()))
    coordinate.set( mRadius, 0.0, 0.0 );
    // Get position from vector
  coordinate += mCenter;
}

void Mesquite::SphericalDomain::vertex_normal_at(Mesh::VertexHandle /*entity_handle*/,
                                                 Vector3D &coordinate) const
{
    // normal is vector from center to input position
  coordinate -= mCenter;
    // make it a unit vector
  double length = coordinate.length();
  coordinate /= length;
    // if input position was at center, choose same position
    // on sphere as snap_to.
  if (!finite(coordinate.x()))
    coordinate.set( 1.0, 0.0, 0.0 );
}
void Mesquite::SphericalDomain::element_normal_at(Mesh::ElementHandle h,
                                                 Vector3D &coordinate) const
{
  SphericalDomain::vertex_normal_at( h, coordinate );
}

void Mesquite::SphericalDomain::vertex_normal_at( 
                                const Mesquite::Mesh::VertexHandle* handle,
                                Mesquite::Vector3D coords[],
                                unsigned count,
                                Mesquite::MsqError& ) const
{
  for (unsigned i = 0; i < count; ++i)
    vertex_normal_at( handle[i], coords[i] );
}

void Mesquite::SphericalDomain::closest_point( Mesquite::Mesh::VertexHandle ,
                                               const Mesquite::Vector3D& position,
                                               Mesquite::Vector3D& closest,
                                               Mesquite::Vector3D& normal,
                                               Mesquite::MsqError& ) const
{
  normal = position - mCenter;
  normal.normalize();
  if (!finite(normal.x()))
    normal.set( 1.0, 0.0, 0.0 );
  closest = mCenter + mRadius * normal;
}


void Mesquite::SphericalDomain::domain_DoF( const Mesh::VertexHandle* ,
                                            unsigned short* dof_array,
                                            size_t num_vertices,
                                            MsqError&  ) const
{
  std::fill( dof_array, dof_array + num_vertices, 2 );
}

      
void Mesquite::SphericalDomain::fit_vertices( Mesh* mesh, MsqError& err, double epsilon )
{
  std::vector<Mesh::VertexHandle> verts;
  mesh->get_all_vertices( verts, err );
  if (!MSQ_CHKERR(err))
    fit_vertices( mesh, arrptr(verts), verts.size(), err, epsilon );
}

void Mesquite::SphericalDomain::fit_vertices( Mesh* mesh, 
                                              const Mesh::VertexHandle* verts,
                                              size_t num_verts,
                                              MsqError& err, 
                                              double epsilon )
{
  std::vector<MsqVertex> coords(num_verts);
  mesh->vertices_get_coordinates( verts, arrptr(coords), num_verts, err );
  MSQ_ERRRTN(err);
  
  if (epsilon <= 0.0)
    epsilon = DomainUtil::default_tolerance( arrptr(coords), num_verts );
  
  Vector3D pts[4];
  if (!DomainUtil::non_coplanar_vertices( arrptr(coords), num_verts, pts, epsilon )) {
    MSQ_SETERR(err)("All vertices are co-planar", MsqError::INVALID_MESH);
    return;
  }
  
    // solve deterinant form of four-point sphere
    
    // Define the bottom 4 rows of a 5x5 matrix.  The top
    // row contains the variables we are solving for, so just
    // fill it with ones.
  const double M_vals[25] = { 
    1,               1,         1,         1,         1,
    pts[0] % pts[0], pts[0][0], pts[0][1], pts[0][2], 1,
    pts[1] % pts[1], pts[1][0], pts[1][1], pts[1][2], 1,
    pts[2] % pts[2], pts[2][0], pts[2][1], pts[2][2], 1,
    pts[3] % pts[3], pts[3][0], pts[3][1], pts[3][2], 1 };
  MsqMatrix<5,5> M(M_vals);
  double M11 = det( MsqMatrix<4,4>(M,0,0) );
  double M12 = det( MsqMatrix<4,4>(M,0,1) );
  double M13 = det( MsqMatrix<4,4>(M,0,2) );
  double M14 = det( MsqMatrix<4,4>(M,0,3) );
  double M15 = det( MsqMatrix<4,4>(M,0,4) );

    // define the sphere
  Vector3D cent( 0.5*M12/M11, -0.5*M13/M11, 0.5*M14/M11 );
  this->set_sphere( cent, sqrt(cent%cent - M15/M11) );
}

