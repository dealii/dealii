/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file IdealElements.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "IdealElements.hpp"
#include "Vector3D.hpp"

namespace MESQUITE_NS {

static Vector3D unit_quad[4] = { Vector3D( -0.5, -0.5, 0.0 ),
                                 Vector3D(  0.5, -0.5, 0.0 ),
                                 Vector3D(  0.5,  0.5, 0.0 ),
                                 Vector3D( -0.5,  0.5, 0.0 ) };

static Vector3D unit_hex[8] = { Vector3D(  0.5, -0.5, -0.5 ),
                                Vector3D(  0.5,  0.5, -0.5 ),
                                Vector3D( -0.5,  0.5, -0.5 ),
                                Vector3D( -0.5, -0.5, -0.5 ),
                                Vector3D(  0.5, -0.5,  0.5 ),
                                Vector3D(  0.5,  0.5,  0.5 ),
                                Vector3D( -0.5,  0.5,  0.5 ),
                                Vector3D( -0.5, -0.5,  0.5 ) };

static Vector3D unit_edge_tri[3];
static Vector3D unit_edge_tet[4];
static Vector3D unit_edge_pyr[5];
static Vector3D unit_edge_wdg[6];
static Vector3D unit_height_pyr[5];
                                 
static Vector3D unit_tri[3];
static Vector3D unit_tet[4];
static Vector3D unit_pyr[5];
static Vector3D unit_wdg[6];
static Vector3D unit_hex_pyr[5];

static void init_tri( Vector3D* coords, double side );
static void init_tet( Vector3D* coords, double side );
static void init_pyr( Vector3D* coords, double side );
static void init_wdg( Vector3D* coords, double side );
static void init_hex_pyr( Vector3D* coords, double height );

static const Vector3D* const* init_unit_edge( Vector3D** );
static const Vector3D* const* init_unit_elem( Vector3D** );

const Vector3D* unit_edge_element( EntityTopology type, bool unit_pyr )
{
  static Vector3D* values[MIXED+1];
  static const Vector3D* const* data = init_unit_edge( values );
  return (type == PYRAMID && unit_pyr) ? data[MIXED] : data[type];
}

const Vector3D* unit_element( EntityTopology type, bool unit_pyr )
{
  static Vector3D* values[MIXED+1];
  static const Vector3D* const* data = init_unit_elem( values );
  return (type == PYRAMID && unit_pyr) ? data[MIXED] : data[type];
}

static const Vector3D* const* init_unit_edge( Vector3D** ptr )
{
  for (unsigned i = 0; i < MIXED; ++i)
    ptr[i] = 0;
  
  init_tri( unit_edge_tri, 1.0 );
  init_tet( unit_edge_tet, 1.0 );
  init_pyr( unit_edge_pyr, 1.0 );
  init_wdg( unit_edge_wdg, 1.0 );
  init_hex_pyr( unit_height_pyr, 1.0 );
  
  ptr[TRIANGLE] = unit_edge_tri;
  ptr[QUADRILATERAL] = unit_quad;
  ptr[TETRAHEDRON] = unit_edge_tet;
  ptr[PYRAMID] = unit_edge_pyr;
  ptr[PRISM] = unit_edge_wdg;
  ptr[HEXAHEDRON] = unit_hex;
  ptr[MIXED] = unit_height_pyr;
  return ptr;
}

static const Vector3D* const* init_unit_elem( Vector3D** ptr )
{
  for (unsigned i = 0; i < MIXED; ++i)
    ptr[i] = 0;
  
  init_tri( unit_tri, 2.0 * pow( 3.0, -0.25 ) );
  init_tet( unit_tet, Mesquite::cbrt( 3.0 ) * sqrt(2.0) );
  init_pyr( unit_pyr, pow( 18.0, 1.0/6.0 ) );
  init_wdg( unit_wdg, Mesquite::cbrt( 4.0 ) * pow( 3.0, -1.0/6.0 ) );
  init_hex_pyr( unit_hex_pyr, Mesquite::cbrt( 3.0 ) );
  
  ptr[TRIANGLE] = unit_tri;
  ptr[QUADRILATERAL] = unit_quad;
  ptr[TETRAHEDRON] = unit_tet;
  ptr[PYRAMID] = unit_pyr;
  ptr[PRISM] = unit_wdg;
  ptr[HEXAHEDRON] = unit_hex;
  ptr[MIXED] = unit_hex_pyr;
  return ptr;
}

static void init_tri( Vector3D* coords, double side )
{
  const double third_height = side * sqrt(3.0) / 6.0;
  coords[1] = Vector3D( -0.5*side,  -third_height, 0.0 );
  coords[2] = Vector3D(  0.5*side,  -third_height, 0.0 );
  coords[0] = Vector3D(  0.0,      2*third_height, 0.0 );
}

static void init_tet( Vector3D* coords, double side )
{
  const double height = side * sqrt(2.0/3.0);
  const double third_base = side * sqrt(3.0) / 6.0;
  coords[0] = Vector3D( -0.5*side,  -third_base, -0.25*height );
  coords[1] = Vector3D(  0.5*side,  -third_base, -0.25*height );
  coords[2] = Vector3D(  0.0,      2*third_base, -0.25*height );
  coords[3] = Vector3D(  0.0,      0.0,           0.75*height );
}

static void init_pyr( Vector3D* coords, double side )
{
  const double height = side * sqrt(2.0) * 0.5;
  coords[0] = Vector3D(  0.5*side, -0.5*side, -0.25*height );
  coords[1] = Vector3D(  0.5*side,  0.5*side, -0.25*height );
  coords[2] = Vector3D( -0.5*side,  0.5*side, -0.25*height );
  coords[3] = Vector3D( -0.5*side, -0.5*side, -0.25*height );
  coords[4] = Vector3D(  0.0,       0.0,       0.75*height ); 
}

static void init_wdg( Vector3D* coords, double side )
{
  const double third_height = side * sqrt(3.0) / 6.0;
  coords[0] = Vector3D( -0.5*side,  -third_height, -0.5*side );
  coords[1] = Vector3D(  0.5*side,  -third_height, -0.5*side );
  coords[2] = Vector3D(  0.0,      2*third_height, -0.5*side );
  coords[3] = Vector3D( -0.5*side,  -third_height,  0.5*side );
  coords[4] = Vector3D(  0.5*side,  -third_height,  0.5*side );
  coords[5] = Vector3D(  0.0,      2*third_height,  0.5*side );
}

static void init_hex_pyr( Vector3D* coords, double side )
{
  coords[0] = Vector3D(  0.5*side, -0.5*side, -0.25*side );
  coords[1] = Vector3D(  0.5*side,  0.5*side, -0.25*side );
  coords[2] = Vector3D( -0.5*side,  0.5*side, -0.25*side );
  coords[3] = Vector3D( -0.5*side, -0.5*side, -0.25*side );
  coords[4] = Vector3D(  0.0,       0.0,       0.75*side ); 
}

} // namespace Mesquite
