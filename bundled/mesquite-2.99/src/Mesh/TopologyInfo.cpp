/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
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

#include "MsqError.hpp"
#include "TopologyInfo.hpp"

#include <string.h>
#include <assert.h>

namespace MESQUITE_NS {

TopologyInfo TopologyInfo::instance;

const char long_polygon_name[]       = "Polygon";
const char long_triangle_name[]      = "Triangle";
const char long_quadrilateral_name[] = "Quadrilateral";
const char long_polyhedron_name[]    = "Polyhedron";
const char long_tetrahedron_name[]   = "Tetrahedron";
const char long_hexahedron_name[]    = "Hexahedron";
const char long_prism_name[]         = "Prism";
const char long_pyramid_name[]       = "Pyramd";
const char long_septahedron_name[]   = "Septahedron";
const char short_polygon_name[]       = "Polygon";
const char short_triangle_name[]      = "Tri";
const char short_quadrilateral_name[] = "Quad";
const char short_polyhedron_name[]    = "Polyhedron";
const char short_tetrahedron_name[]   = "Tet";
const char short_hexahedron_name[]    = "Hex";
const char short_prism_name[]         = "Pri";
const char short_pyramid_name[]       = "Pyr";
const char short_septahedron_name[]   = "Sept";

TopologyInfo::TopologyInfo()
{
  memset( dimMap, 0, sizeof(dimMap) );
  memset( adjMap, 0, sizeof(adjMap) );
  memset( edgeMap, 0, sizeof(edgeMap) );
  memset( faceMap, 0, sizeof(faceMap) );
  memset( vertAdjMap, 0, sizeof(vertAdjMap) );
  memset( shortNames, 0, sizeof(shortNames) );
  memset( longNames, 0, sizeof(longNames) );
  
  longNames[POLYGON]       = long_polygon_name;
  longNames[TRIANGLE]      = long_triangle_name;
  longNames[QUADRILATERAL] = long_quadrilateral_name;
  longNames[POLYHEDRON]    = long_polyhedron_name;
  longNames[TETRAHEDRON]   = long_tetrahedron_name;
  longNames[HEXAHEDRON]    = long_hexahedron_name;
  longNames[PRISM]         = long_prism_name;
  longNames[PYRAMID]       = long_pyramid_name;
  longNames[SEPTAHEDRON]   = long_septahedron_name;
  
  shortNames[POLYGON]       = short_polygon_name;
  shortNames[TRIANGLE]      = short_triangle_name;
  shortNames[QUADRILATERAL] = short_quadrilateral_name;
  shortNames[POLYHEDRON]    = short_polyhedron_name;
  shortNames[TETRAHEDRON]   = short_tetrahedron_name;
  shortNames[HEXAHEDRON]    = short_hexahedron_name;
  shortNames[PRISM]         = short_prism_name;
  shortNames[PYRAMID]       = short_pyramid_name;
  shortNames[SEPTAHEDRON]   = short_septahedron_name;
  
  dimMap[POLYGON ]      = 2;
  dimMap[TRIANGLE]      = 2;
  dimMap[QUADRILATERAL] = 2;
  dimMap[POLYHEDRON]    = 3;
  dimMap[TETRAHEDRON]   = 3;
  dimMap[HEXAHEDRON]    = 3;
  dimMap[PRISM]         = 3;
  dimMap[PYRAMID]       = 3;
  dimMap[SEPTAHEDRON]   = 3;
  
  adjMap[TRIANGLE][0] = 3;
  adjMap[TRIANGLE][1] = 3;
  adjMap[TRIANGLE][2] = 1;
  adjMap[TRIANGLE][3] = 0;
  
  adjMap[QUADRILATERAL][0] = 4;
  adjMap[QUADRILATERAL][1] = 4;
  adjMap[QUADRILATERAL][2] = 1;
  adjMap[QUADRILATERAL][3] = 0;
  
  adjMap[TETRAHEDRON][0] = 4;
  adjMap[TETRAHEDRON][1] = 6;
  adjMap[TETRAHEDRON][2] = 4;
  adjMap[TETRAHEDRON][3] = 1;
  
  adjMap[HEXAHEDRON][0] = 8;
  adjMap[HEXAHEDRON][1] = 12;
  adjMap[HEXAHEDRON][2] = 6;
  adjMap[HEXAHEDRON][3] = 1;
  
  adjMap[PRISM][0] = 6;
  adjMap[PRISM][1] = 9;
  adjMap[PRISM][2] = 5;
  adjMap[PRISM][3] = 1;
  
  adjMap[PYRAMID][0] = 5;
  adjMap[PYRAMID][1] = 8;
  adjMap[PYRAMID][2] = 5;
  adjMap[PYRAMID][3] = 1;
  
  adjMap[SEPTAHEDRON][0] = 7;
  adjMap[SEPTAHEDRON][1] = 11;
  adjMap[SEPTAHEDRON][2] = 6;  /* See description in TSTT mesh interface doc */
  adjMap[SEPTAHEDRON][3] = 1;

  int side;
  for (side = 0; side < 3; ++side)
  {
    edgeMap[TRIANGLE-FIRST_FACE][side][0] = side;
    edgeMap[TRIANGLE-FIRST_FACE][side][1] = (side+1)%3;
  }
  for (side = 0; side < 4; ++side)
  {
    edgeMap[QUADRILATERAL-FIRST_FACE][side][0] = side;
    edgeMap[QUADRILATERAL-FIRST_FACE][side][1] = (side+1)%4;
  }
  for (side = 0; side < 3; ++side)
  {
    edgeMap[TETRAHEDRON-FIRST_FACE][side][0] = side;
    edgeMap[TETRAHEDRON-FIRST_FACE][side][1] = (side+1)%3;
  }
  for (side = 3; side < 6; ++side)
  {
    edgeMap[TETRAHEDRON-FIRST_FACE][side][0] = side -3 ;
    edgeMap[TETRAHEDRON-FIRST_FACE][side][1] = 3;
  }
  for (side = 0; side < 4; ++side)
  {
    edgeMap[HEXAHEDRON-FIRST_FACE][side][0] = side;
    edgeMap[HEXAHEDRON-FIRST_FACE][side][1] = (side+1)%4;
  }
  for (side = 4; side < 8; ++side)
  {
    edgeMap[HEXAHEDRON-FIRST_FACE][side][0] = side - 4;
    edgeMap[HEXAHEDRON-FIRST_FACE][side][1] = side;
  }
  for (side = 8; side < 12; ++side)
  {
    edgeMap[HEXAHEDRON-FIRST_FACE][side][0] = side - 4;
    edgeMap[HEXAHEDRON-FIRST_FACE][side][1] = 4+(side+1)%4;
  }
  for (side = 0; side < 3; ++side)
  {
    edgeMap[PRISM-FIRST_FACE][side][0] = side;
    edgeMap[PRISM-FIRST_FACE][side][1] = (side+1)%3;
  }
  for (side = 3; side < 6; ++side)
  {
    edgeMap[PRISM-FIRST_FACE][side][0] = side - 3;
    edgeMap[PRISM-FIRST_FACE][side][1] = side;
  }
  for (side = 6; side < 9; ++side)
  {
    edgeMap[PRISM-FIRST_FACE][side][0] = side-3;
    edgeMap[PRISM-FIRST_FACE][side][1] = 3+(side+1)%3;
  }
  for (side = 0; side < 4; ++side)
  {
    edgeMap[PYRAMID-FIRST_FACE][side][0] = side;
    edgeMap[PYRAMID-FIRST_FACE][side][1] = (side+1)%4;
  }
  for (side = 4; side < 8; ++side)
  {
    edgeMap[PYRAMID-FIRST_FACE][side][0] = side - 4;
    edgeMap[PYRAMID-FIRST_FACE][side][1] = 4;
  }
  
  for (side = 0; side < 3; ++side)
  {
    faceMap[TETRAHEDRON-FIRST_VOL][side][0] = 3;
    faceMap[TETRAHEDRON-FIRST_VOL][side][1] = side;
    faceMap[TETRAHEDRON-FIRST_VOL][side][2] = (side+1)%3;
    faceMap[TETRAHEDRON-FIRST_VOL][side][3] = 3;
  }
  faceMap[TETRAHEDRON-FIRST_VOL][3][0] = 3;
  faceMap[TETRAHEDRON-FIRST_VOL][3][1] = 2;
  faceMap[TETRAHEDRON-FIRST_VOL][3][2] = 1;
  faceMap[TETRAHEDRON-FIRST_VOL][3][3] = 0;

  for (side = 0; side < 4; ++side)
  {
    faceMap[HEXAHEDRON-FIRST_VOL][side][0] = 4;
    faceMap[HEXAHEDRON-FIRST_VOL][side][1] = side;
    faceMap[HEXAHEDRON-FIRST_VOL][side][2] = (side+1)%4;
    faceMap[HEXAHEDRON-FIRST_VOL][side][3] = 4+(side+1)%4;
    faceMap[HEXAHEDRON-FIRST_VOL][side][4] = side + 4;
  }
  faceMap[HEXAHEDRON-FIRST_VOL][4][0] = 4;
  faceMap[HEXAHEDRON-FIRST_VOL][4][1] = 3;
  faceMap[HEXAHEDRON-FIRST_VOL][4][2] = 2;
  faceMap[HEXAHEDRON-FIRST_VOL][4][3] = 1;
  faceMap[HEXAHEDRON-FIRST_VOL][4][4] = 0;
  faceMap[HEXAHEDRON-FIRST_VOL][5][0] = 4;
  faceMap[HEXAHEDRON-FIRST_VOL][5][1] = 4;
  faceMap[HEXAHEDRON-FIRST_VOL][5][2] = 5;
  faceMap[HEXAHEDRON-FIRST_VOL][5][3] = 6;
  faceMap[HEXAHEDRON-FIRST_VOL][5][4] = 7;
  
  for (side = 0; side < 4; ++side)
  {
    faceMap[PYRAMID-FIRST_VOL][side][0] = 3;
    faceMap[PYRAMID-FIRST_VOL][side][1] = side;
    faceMap[PYRAMID-FIRST_VOL][side][2] = (side+1)%4;
    faceMap[PYRAMID-FIRST_VOL][side][3] = 4;
  }
  faceMap[PYRAMID-FIRST_VOL][4][0] = 4;
  faceMap[PYRAMID-FIRST_VOL][4][1] = 3;
  faceMap[PYRAMID-FIRST_VOL][4][2] = 2;
  faceMap[PYRAMID-FIRST_VOL][4][3] = 1;
  faceMap[PYRAMID-FIRST_VOL][4][4] = 0;

  for (side = 0; side < 3; ++side)
  {
    faceMap[PRISM-FIRST_VOL][side][0] = 4;
    faceMap[PRISM-FIRST_VOL][side][1] = side;
    faceMap[PRISM-FIRST_VOL][side][2] = (side+1)%3;
    faceMap[PRISM-FIRST_VOL][side][3] = 3+(side+1)%3;
    faceMap[PRISM-FIRST_VOL][side][4] = side + 3;
  }
  faceMap[PRISM-FIRST_VOL][3][0] = 3;
  faceMap[PRISM-FIRST_VOL][3][1] = 2;
  faceMap[PRISM-FIRST_VOL][3][2] = 1;
  faceMap[PRISM-FIRST_VOL][3][3] = 0;
  faceMap[PRISM-FIRST_VOL][4][0] = 3;
  faceMap[PRISM-FIRST_VOL][4][1] = 3;
  faceMap[PRISM-FIRST_VOL][4][2] = 4;
  faceMap[PRISM-FIRST_VOL][4][3] = 5;
  
  int i;
  for (i = 0; i < 3; ++i)
  {
    vertAdjMap[TRIANGLE-FIRST_FACE][i][0] = 2;
    vertAdjMap[TRIANGLE-FIRST_FACE][i][1] = (i+1)%3;
    vertAdjMap[TRIANGLE-FIRST_FACE][i][2] = (i+2)%3;
  }
  
  for (i = 0; i < 4; ++i)
  {
    vertAdjMap[QUADRILATERAL-FIRST_FACE][i][0] = 2;
    vertAdjMap[QUADRILATERAL-FIRST_FACE][i][1] = (i+1)%4;
    vertAdjMap[QUADRILATERAL-FIRST_FACE][i][2] = (i+3)%4;
  }
  
  
  unsigned tet_corner_data[] = { 1, 2, 3, 
                                 0, 3, 2,
                                 3, 0, 1,
                                 2, 1, 0 };
  for (i = 0; i < 4; ++i)
  {
    vertAdjMap[TETRAHEDRON-FIRST_FACE][i][0] = 3;
    for (unsigned j = 0; j < 3; ++j)
      vertAdjMap[TETRAHEDRON-FIRST_FACE][i][j+1] = tet_corner_data[3*i+j];
  }
  
  for (i = 0; i < 4; ++i)
  {
    vertAdjMap[PYRAMID-FIRST_FACE][i][0] = 3;
    vertAdjMap[PYRAMID-FIRST_FACE][i][1] = (i+1)%4;
    vertAdjMap[PYRAMID-FIRST_FACE][i][2] = (i+3)%4;
    vertAdjMap[PYRAMID-FIRST_FACE][i][3] = 4;
  }
  vertAdjMap[PYRAMID-FIRST_FACE][4][0] = 4;
  for (i = 0; i < 4; i++)
    vertAdjMap[PYRAMID-FIRST_FACE][4][i+1] = 3 - i;
  
  for (i = 0; i < 4; ++i)
  {
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][0] = 3;
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][1] = (i+1)%4;
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][2] = (i+3)%4;
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][3] = i+4;
  }
  for (i = 4; i < 8; ++i)
  {
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][0] = 3;
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][1] = (i+3)%4+4;
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][2] = (i+1)%4+4;
    vertAdjMap[HEXAHEDRON-FIRST_FACE][i][3] = i-4;
  }
  
  for (i = 0; i < 3; ++i)
  {
    vertAdjMap[PRISM-FIRST_FACE][i][0] = 3;
    vertAdjMap[PRISM-FIRST_FACE][i][1] = (i+1)%3;
    vertAdjMap[PRISM-FIRST_FACE][i][2] = (i+2)%3;
    vertAdjMap[PRISM-FIRST_FACE][i][3] = i+3;
  }
  for (i = 3; i < 6; ++i)
  {
    vertAdjMap[PRISM-FIRST_FACE][i][0] = 3;
    vertAdjMap[PRISM-FIRST_FACE][i][1] = (i+2)%3+3;
    vertAdjMap[PRISM-FIRST_FACE][i][2] = (i+1)%3+3;
    vertAdjMap[PRISM-FIRST_FACE][i][3] = i-3;
  }
  
    // Build reverse vertex-vertex adjacency index map
  const EntityTopology types[] = { TRIANGLE, 
                                   QUADRILATERAL, 
                                   TETRAHEDRON,
                                   PYRAMID,
                                   PRISM, 
                                   HEXAHEDRON };
  const int num_types = sizeof(types)/sizeof(types[0]);
  for (i = 0; i < num_types; ++i)
  {
    const unsigned num_vert = corners( types[i] );
    for (unsigned v = 0; v < num_vert; ++v)
    {
      unsigned num_v_adj;
      const unsigned* v_adj = adjacent_vertices( types[i], v, num_v_adj );
      unsigned* reverse = revVertAdjIdx[types[i]-FIRST_FACE][v];
      reverse[0] = num_v_adj;
      
      for (unsigned j = 0; j < num_v_adj; ++j)
      {
        unsigned num_j_adj, k;
        const unsigned* j_adj = adjacent_vertices( types[i], v_adj[j], num_j_adj );
        for (k = 0; k < num_j_adj && j_adj[k] != v; ++k);
        assert( k < num_j_adj ); // If this fails, vertAdjMap is corrupt!
        reverse[j+1] = k;
      }
    }
  }
}

void TopologyInfo::higher_order( EntityTopology topo,
                                     unsigned num_nodes,
                                     bool& midedge,
                                     bool& midface,
                                     bool& midvol,
                                     MsqError& err )
{
  int ho = higher_order( topo, num_nodes, err );
  midedge = (bool)( (ho & (1<<1)) >> 1);
  midface = (bool)( (ho & (1<<2)) >> 2);
  midvol  = (bool)( (ho & (1<<3)) >> 3);
}


int TopologyInfo::higher_order( EntityTopology topo, 
                                unsigned num_nodes, 
                                MsqError& err )
{
  int result = 0;
  if (topo == POLYGON)  // polygons currently do not have higher order elements
    return 0;

  if (topo >= MIXED || num_nodes < instance.adjMap[topo][0])
  {
    MSQ_SETERR(err)("Invalid element topology", MsqError::INVALID_ARG);
    return 0;
  }
  
  unsigned dim = instance.dimMap[topo];
  assert( num_nodes >= instance.adjMap[topo][0] );
  unsigned nodes = num_nodes - instance.adjMap[topo][0];
  unsigned edges = instance.adjMap[topo][1];
  unsigned faces = instance.adjMap[topo][2];
  if (edges && nodes >= edges)
  {
    nodes -= edges;
    result |= 1<<1;
  }
  if (faces && nodes >= faces)
  {
    nodes -= faces;
    result |= 1<<2;
  }
  if (1 == nodes)
  {
    if (2 == dim)
    {
      nodes -= 1;
      result |= 1<<2;
    }
    else if(3 == dim)
    {
      nodes -= 1;
      result |= 1<<3;
    }
  }
  
  if (nodes)
  {
    MSQ_SETERR(err)("Invalid element topology", MsqError::INVALID_STATE);
  }
  
  return result;
}
      
int TopologyInfo::higher_order_from_side( EntityTopology topo,
                                          unsigned num_nodes,
                                          unsigned side_dimension,
                                          unsigned side_number,
                                          MsqError& err )
{
  bool mids[4] = { true };
  higher_order( topo, num_nodes, mids[1], mids[2], mids[3], err );
  MSQ_ERRZERO(err);
  
  if (side_dimension > dimension(topo) || 
      side_number > adjacent(topo, side_dimension)) {
    MSQ_SETERR(err)(MsqError::INVALID_ARG,"Invalid side number: %u\n", side_number );
    return 0;
  }
  
  if (!mids[side_dimension])
    return -1;
  
  int result = side_number;
  switch (side_dimension) {
    case 3: if (mids[2]) result += faces(topo);
    case 2: if (mids[1]) result += edges(topo);
    case 1: result += corners(topo);
    case 0: break;
    default: 
      MSQ_SETERR(err)(MsqError::INVALID_ARG,"Invalid dimension: %u\n", side_dimension );
      return 0;
  }
  return result;
}

void TopologyInfo::side_from_higher_order( EntityTopology topo,
                                           unsigned num_nodes,
                                           unsigned node_number,
                                           unsigned& side_dim_out,
                                           unsigned& side_num_out,
                                           MsqError& err )
{
  bool midedge, midface, midvol;
  higher_order( topo, num_nodes, midedge, midface, midvol, err );
  MSQ_ERRRTN(err);
  side_num_out = node_number;
  
  if (side_num_out < corners(topo)) {
    side_dim_out = 0;
    return;
  }
  side_num_out -= corners(topo);
  
  if (midedge) {
    if (side_num_out < edges(topo)) {
      side_dim_out = 1;
      return;
    }
    side_num_out -= edges(topo);
  }
  
  if (midface) {
    if (side_num_out < faces(topo)) {
      side_dim_out = 2;
      return;
    }
    side_num_out -= faces(topo);
  }
  
  if (midvol && side_num_out == 0) {
    side_dim_out = 3;
    return;
  }
  
  MSQ_SETERR(err)(MsqError::INVALID_ARG,"Invalid node index\n");
}

const unsigned*  TopologyInfo::edge_vertices( EntityTopology topo,
                                              unsigned edge, 
                                              MsqError& err)
{
  if (topo < (EntityTopology)FIRST_FACE || 
      topo > (EntityTopology)LAST_VOL || 
      edge >= edges( topo ) )
  {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    topo = (EntityTopology)FIRST_FACE;
    edge = 0;
  }
  
  return instance.edgeMap[topo-FIRST_FACE][edge];
}

const unsigned*  TopologyInfo::edge_vertices( EntityTopology topo,
                                              unsigned edge )
{
  if (topo < (EntityTopology)FIRST_FACE || 
      topo > (EntityTopology)LAST_VOL || 
      edge >= edges( topo ) )
  {
    return 0;
  }
  return instance.edgeMap[topo-FIRST_FACE][edge];
}

const unsigned* TopologyInfo::face_vertices( EntityTopology topo,
                                             unsigned face,
                                             unsigned& length,
                                             MsqError& err )
{
  if (topo < (EntityTopology)FIRST_VOL || 
      topo > (EntityTopology)LAST_VOL || 
      face >= faces( topo ) )
  {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    topo = (EntityTopology)FIRST_VOL;
    face = 0;
  }
  
  length = instance.faceMap[topo-FIRST_VOL][face][0];
  return instance.faceMap[topo-FIRST_VOL][face] + 1;
}
const unsigned* TopologyInfo::face_vertices( EntityTopology topo,
                                             unsigned face,
                                             unsigned& length )
{
  if (topo < (EntityTopology)FIRST_VOL || 
      topo > (EntityTopology)LAST_VOL || 
      face >= faces( topo ) )
  {
    return 0;
  }
  
  length = instance.faceMap[topo-FIRST_VOL][face][0];
  return instance.faceMap[topo-FIRST_VOL][face] + 1;
}




const unsigned* TopologyInfo::side_vertices( EntityTopology topo,
                                             unsigned dim,
                                             unsigned side,
                                             unsigned& count_out,
                                             MsqError& err )
{
  static const unsigned all[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  const unsigned* result;
  
  if (dim != 0 && dim == dimension(topo))
  {
    count_out = corners( topo );
    result = all;
  }
  else if (dim == 1)
  {
    count_out = 2;
    result = edge_vertices( topo, side, err );
  }
  else if( dim == 2)
  {
    result = face_vertices( topo, side, count_out, err );
  } 
  else
  {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    count_out = 0;
    result = 0;
  }    
  return result;
}
const unsigned* TopologyInfo::side_vertices( EntityTopology topo,
                                             unsigned dim,
                                             unsigned side,
                                             unsigned& count_out  )
{
  static const unsigned all[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  const unsigned* result;
  
  if (dim != 0 && dim == dimension(topo))
  {
    count_out = corners( topo );
    result = all;
  }
  else if (dim == 1)
  {
    count_out = 2;
    result = edge_vertices( topo, side );
  }
  else if( dim == 2)
  {
    result = face_vertices( topo, side, count_out );
  } 
  else
  {
    result = 0;
  }    
  return result;
}
      


void TopologyInfo::side_number( EntityTopology topo,
                                    unsigned num_nodes,
                                    unsigned node_index,
                                    unsigned& side_dim_out,
                                    unsigned& side_num_out,
                                    MsqError& err )
{
  if (topo >= (EntityTopology)MIXED || num_nodes < instance.adjMap[topo][0])
  {
    MSQ_SETERR(err)("Invalid element topology", MsqError::INVALID_ARG);
    return;
  }
  
  unsigned nodes = instance.adjMap[topo][0];
  unsigned edges = instance.adjMap[topo][1];
  unsigned faces = instance.adjMap[topo][2];
  side_num_out = node_index;

  if (side_num_out < nodes)
  {
    side_dim_out = 0;
    return;
  }
  num_nodes -= nodes;
  side_num_out -= nodes;
  
  if (edges && num_nodes >= edges)
  {
    if (side_num_out < edges)
    {
      side_dim_out = 1;
      return;
    }
    num_nodes -= edges;
    side_num_out -= edges;
  }
  if (faces && num_nodes >= faces)
  {
    if (side_num_out < faces)
    {
      side_dim_out = 2;
      return;
    }
    num_nodes -= faces;
    side_num_out -= faces;
  }
  if (side_num_out == 0)
  {
    side_dim_out = instance.dimMap[topo];
    side_num_out = 0;
    return;
  }
  
  MSQ_SETERR(err)(MsqError::INVALID_ARG);
}
  
  

const unsigned* TopologyInfo::adjacent_vertices( EntityTopology topo,
                                              unsigned index,
                                              unsigned& num_adj_out )
{
  const unsigned count = corners( topo );
  if (!count || index >= count)
  {
    num_adj_out = 0;
    return 0;
  }
  
  const unsigned* vect = instance.vertAdjMap[topo-FIRST_FACE][index];
  num_adj_out = vect[0];
  return vect + 1;
}

const unsigned* TopologyInfo::reverse_vertex_adjacency_offsets(
                                              EntityTopology topo,
                                              unsigned index,
                                              unsigned& num_adj_out )
{
  const unsigned count = corners( topo );
  if (!count || index >= count)
  {
    num_adj_out = 0;
    return 0;
  }
  
  const unsigned* vect = instance.revVertAdjIdx[topo-FIRST_FACE][index];
  num_adj_out = vect[0];
  return vect + 1;
}

bool TopologyInfo::compare_sides( const size_t* verts1,
                               EntityTopology type1,
                               unsigned side1,
                               const size_t* verts2,
                               EntityTopology type2,
                               unsigned side2,
                               unsigned side_dim,
                               MsqError& err )
{
  const unsigned *conn1, *conn2;
  unsigned len1, len2;
  
  conn1 = side_vertices( type1, side_dim, side1, len1, err );
  MSQ_ERRZERO(err);
  conn2 = side_vertices( type2, side_dim, side2, len2, err );
  MSQ_ERRZERO(err);
  
    // obviously not the same if different number of vertices
    // (triangular face cannot match quadrilateral face)
  if (len1 != len2)
    return false;
  
    // Find location (i) in vertices of side of second element
    // that matches the first vertex in the side of the first 
    // element.
  unsigned i, j;
  for (i = 0; i < len2; ++i)
    if (verts1[conn1[0]] == verts2[conn2[i]])
      break;
    // If not found, then no match
  if (i == len2)
    return false;
  
    // Try comparing side connectivity in forward order
  for (j = 1; j < len1; ++j)
    if (verts1[conn1[j]] != verts2[conn2[(i+j)%len2]])
      break;
    // If they match, we're done
  if (j == len1)
    return true;
  
    // Try comparing in reverse order
  for (j = 1; j < len1; ++j)
    if (verts1[conn1[j]] != verts2[conn2[(i+len2-j)%len2]])
      return false;
    // If here, matched in reverse order
  return true;
}

unsigned TopologyInfo::find_edge( EntityTopology topo,
                                  const unsigned* side_vertices,
                                  bool& reversed_out,
                                  MsqError& err )
{
  if (dimension(topo) <= 1) {
    MSQ_SETERR(err)(MsqError::INVALID_ARG,"Invalid element dimension");
    return (unsigned)-1;
  }

  for (unsigned i = 0; i < edges(topo); ++i) {
    const unsigned* edge = edge_vertices( topo, i, err );
    MSQ_ERRZERO(err);

    if (edge[0] == side_vertices[0] &&
        edge[1] == side_vertices[1]) {
      reversed_out = false;
      return i;
    }

    if (edge[0] == side_vertices[1] &&
        edge[1] == side_vertices[0]) {
      reversed_out = true;
      return i;
    }
  }
  
  MSQ_SETERR(err)(MsqError::INVALID_ARG,"No such edge");
  return (unsigned)-1;
}

unsigned TopologyInfo::find_face( EntityTopology topo,
                                  const unsigned* side_vertices,
                                  unsigned num_vertices,
                                  bool& reversed_out,
                                  MsqError& err )
{
  if (dimension(topo) <= 2) {
    MSQ_SETERR(err)(MsqError::INVALID_ARG,"Invalid element dimension");
    return (unsigned)-1;
  }

  for (unsigned i = 0; i < faces(topo); ++i) {
    unsigned j, n, offset;
    const unsigned* face = face_vertices( topo, i, n, err );
    MSQ_ERRZERO(err);
    if (n != num_vertices)
      continue;

    for (offset = 0; offset < num_vertices; ++offset)
      if (face[offset] == side_vertices[0])
        break;
    if (offset == num_vertices)
      continue;

    for (j = 1; j < num_vertices; ++j)
      if (side_vertices[j] != face[(offset + j)%num_vertices])
        break;
    if (j == num_vertices) {
      reversed_out = false;
      return i;
    }

    for (j = 1; j < num_vertices; ++j)
      if (side_vertices[j] != face[(offset + num_vertices - j)%num_vertices])
        break;
    if (j == num_vertices) {
      reversed_out = true;
      return i;
    }
  }
  
  MSQ_SETERR(err)(MsqError::INVALID_ARG,"No such face");
  return (unsigned)-1;
}


void TopologyInfo::find_side( EntityTopology topo, 
                              const unsigned* side_vertices,
                              unsigned num_vertices,
                              unsigned& dimension_out,
                              unsigned& number_out,
                              bool& reversed_out,
                              MsqError& err )
{
  switch (num_vertices) {
  case 1:
    dimension_out = 0;
    number_out = *side_vertices;
    reversed_out = false;
    if (*side_vertices >= corners(topo)) 
      MSQ_SETERR(err)(MsqError::INVALID_ARG,"Invalid corner number: %u\n", *side_vertices);
    break;
  case 2:
    dimension_out = 1;
    number_out = find_edge( topo, side_vertices, reversed_out, err );
    MSQ_CHKERR(err);
    break;
  case 3:
  case 4:
    dimension_out = 2;
    number_out = find_face( topo, side_vertices, num_vertices, reversed_out, err );
    MSQ_CHKERR(err);
    break;
  default:
    MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT, "Invalid number of side vertices: %u\n", num_vertices );
    break;
  }
}
  
      


} //namepsace Mesquite
