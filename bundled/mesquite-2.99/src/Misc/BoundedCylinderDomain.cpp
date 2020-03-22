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

#ifndef MSQ_BOUNDED_CYLINDER_DOMAIN_CPP
#define MSQ_BOUNDED_CYLINDER_DOMAIN_CPP

#include "BoundedCylinderDomain.hpp"
#include "MsqError.hpp"
#include "MsqVertex.hpp"

#include <limits>
#include <algorithm>

namespace MESQUITE_NS {

void BoundedCylinderDomain::domain_DoF( const Mesh::VertexHandle* handle_array,
                                        unsigned short* dof_array,
                                        size_t count,
                                        MsqError& err ) const
{
  double t;
  for (size_t i = 0; i < count; ++i)
    if (find_curve( handle_array[i], t ))
      dof_array[i] = 1;
    else
      dof_array[i] = 2;
}

void BoundedCylinderDomain::create_curve( double distance, 
                       const std::vector<Mesh::VertexHandle>& handles )
{
  Curve c;
  c.t = distance;
  c.handles = handles;
  std::sort( c.handles.begin(), c.handles.end() );
  curveList.push_back( c );
}

void BoundedCylinderDomain::create_curve( double distance,
                                          Mesh* mesh,
                                          double tolerance,
                                          MsqError& err )
{
  std::vector<Mesh::VertexHandle> handles;
  mesh->get_all_vertices( handles, err ); MSQ_ERRRTN(err);
  if (handles.empty()) {
    MSQ_SETERR(err)("No vertices in mesh.\n", MsqError::INVALID_ARG );
    return;
  }
  
  std::vector<MsqVertex> coords(handles.size());
  mesh->vertices_get_coordinates( arrptr(handles), arrptr(coords), handles.size(), err );
  MSQ_ERRRTN(err);
  
  std::vector<Mesh::EntityHandle> list;
  Vector3D close, normal;
  for (size_t i = 0; i < handles.size(); ++i)
  {
    evaluate( distance, coords[i], close, normal );
    if ((coords[i] - close).length() < tolerance)
      list.push_back( handles[i] );
  }
  
  if (list.empty())
  {
    MSQ_SETERR(err)("No vertices within specified tolerance.\n", MsqError::INVALID_ARG );
    return;
  }
  
  create_curve( distance, list );
}

void BoundedCylinderDomain::evaluate( double t, 
                                      const Vector3D& point,
                                      Vector3D& closest,
                                      Vector3D& normal ) const
{
  const double EPSILON = std::numeric_limits<double>::epsilon();
  double t2 = axis() % (point - center());
  const Vector3D circ_center = center() + t * axis();
  const Vector3D axis_point = center() + t2 * axis();
  
  normal = point - axis_point;
  const double len = normal.length();
  if (len < EPSILON)
  {
    this->CylinderDomain::evaluate( 0, axis_point, closest, normal );
  }
  else
  {
    normal /= len;
    closest = circ_center + radius() * normal;
  }
}

void BoundedCylinderDomain::evaluate( Mesh::VertexHandle handle,
                                      const Vector3D& point,
                                      Vector3D& closest,
                                      Vector3D& normal ) const
{
  double t;
  if (find_curve( handle, t ))
    evaluate( t, point, closest, normal );
  else
    this->CylinderDomain::evaluate( handle, point, closest, normal );
}

bool BoundedCylinderDomain::find_curve( Mesh::VertexHandle handle, double& t ) const
{
  for (std::list<Curve>::const_iterator i = curveList.begin(); i != curveList.end(); ++i)
    if (std::binary_search( i->handles.begin(), i->handles.end(), handle ))
    {
      t = i->t;
      return true;
    }
    
  return false;
}


} // namespace Mesquite

#endif
