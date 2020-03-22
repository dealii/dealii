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

#include "VtkTypeInfo.hpp"
#include "MsqError.hpp"
#include <assert.h>

namespace MESQUITE_NS {

const unsigned vtk_pixel_order[] = { 0, 1, 3, 2 };

const unsigned vtk_voxel_order[] = { 0, 1, 3, 2, 4, 5, 7, 6 };

const unsigned vtk_wedge_order[] = { 0, 2, 1,   // bottom corners
                                     3, 5, 4,   // top corners
                                     8, 7, 6,   // bottom edges
                                    12,14,13,   // lateral edges
                                    11,10, 9,   // top edges
                                    17,16,15 }; // quadrilateral faces

const unsigned vtk_hex_order[] = { 0,  1,  2,  3, // corners (same)
                                   4,  5,  6,  7,
                                   8,  9, 10, 11, // mid-edge (top before lateral)
                                  16, 17, 18, 19,
                                  12, 13, 14, 15,
                                  22, 21, 23, 20, // mid-face (mixed up) & mid-region (same)
                                  24, 25, 26 };

static const VtkTypeInfo typeInfoList[] = {
      { 0,                         0, MIXED,         0, 0 },
      { "vertex",                  1, MIXED,         1, 0 },
      { "polyvertex",              2, MIXED,         0, 0 },
      { "line",                    3, MIXED,         2, 0 },
      { "polyline",                4, MIXED,         0, 0 },
      { "triangle",                5, TRIANGLE,      3, 0 },
      { "triangle strip",          6, MIXED,         0, 0 },
      { "polygon",                 7, POLYGON,       0, 0 },
      { "pixel",                   8, QUADRILATERAL, 4, vtk_pixel_order },
      { "quadrilateral",           9, QUADRILATERAL, 4, 0 }, 
      { "tetrahedron",            10, TETRAHEDRON,   4, 0 }, 
      { "voxel",                  11, HEXAHEDRON,    8, vtk_voxel_order }, 
      { "hexahedron",             12, HEXAHEDRON,    8, 0 }, 
      { "wedge",                  13, PRISM,         6, vtk_wedge_order }, 
      { "pyramid",                14, PYRAMID,       5, 0 },
      { "pentagonal prism",       15, MIXED,        10, 0 }, // not supported
      { "hexagonal prism",        16, MIXED,        12, 0 }, // not supported
      { "invalid (17)",           17, MIXED,         0, 0 },
      { "invalid (18)",           18, MIXED,         0, 0 },
      { "invalid (19)",           19, MIXED,         0, 0 },
      { "invalid (20)",           20, MIXED,         0, 0 },
      { "quadratic edge",         21, MIXED,         3, 0 },
      { "quadratic tri",          22, TRIANGLE,      6, 0 },
      { "quadratic quad",         23, QUADRILATERAL, 8, 0 },
      { "quadratic tet",          24, TETRAHEDRON,  10, 0 },
      { "quadratic hex",          25, HEXAHEDRON,   20, vtk_hex_order },
      { "quadratic wedge",        26, PRISM,        15, vtk_wedge_order },
      { "quadratic pyramid",      27, PYRAMID,      13, 0 },
      { "bi-quadratic quad",      28, QUADRILATERAL, 9, 0 },
      { "tri-quadratic hex",      29, HEXAHEDRON,   27, vtk_hex_order },
      { "quadratic-linear quad",  30, MIXED,         6, 0 },     // not supported
      { "quadratic-linear wedge", 31, MIXED,        12, vtk_wedge_order }, // not supported
      { "bi-quadratic wedge",     32, MIXED,        18, vtk_wedge_order }, // not supported
      { "bi-quadratic hex",       33, MIXED,        24, vtk_hex_order },  // not supported
      { 0,                        34, MIXED,         0, 0 }
    };

const unsigned typeInfoListLen = sizeof(typeInfoList)/sizeof(typeInfoList[0]);

const unsigned reverseIndexList[][3] = {
  {  0,  0,  0 }, // 0
  {  0,  0,  0 }, // 1
  {  0,  0,  0 }, // 2
  {  0,  0,  0 }, // 3
  {  0,  0,  0 }, // 4
  {  0,  0,  0 }, // 5
  {  0,  0,  0 }, // 6
  {  7,  0,  0 }, // POLYGON
  {  5, 22,  0 }, // TRIANGLE
  {  9, 23, 28 }, // QUADRILATERAL
  {  0,  0,  0 }, // POLYHEDRON
  { 10, 24,  0 }, // TETRAHEDRON
  { 12, 25, 29 }, // HEXAHEDRON
  { 13, 26,  0 }, // PRISM
  { 14, 27,  0 }, // PYRAMID
  {  0,  0,  0 }, // SEPTAHEDRON
  {  0,  0,  0 } };

const VtkTypeInfo* VtkTypeInfo::find_type( unsigned vtk_type, MsqError& err )
{
  if (vtk_type >= typeInfoListLen)
  {
    MSQ_SETERR(err)("Type out of bounds", MsqError::INVALID_ARG);
    return 0;
  }
  
  return &typeInfoList[vtk_type];
}

const VtkTypeInfo* VtkTypeInfo::find_type( EntityTopology msq_type,
                                           unsigned num_nodes,
                                           MsqError& err )
{
  if      ( typeInfoList[ reverseIndexList[msq_type][0] ].numNodes == num_nodes )
    return &typeInfoList[ reverseIndexList[msq_type][0] ];
  else if ( typeInfoList[ reverseIndexList[msq_type][1] ].numNodes == num_nodes )
    return &typeInfoList[ reverseIndexList[msq_type][1] ];
  else if ( typeInfoList[ reverseIndexList[msq_type][2] ].numNodes == num_nodes )
    return &typeInfoList[ reverseIndexList[msq_type][2] ];
 
  if (msq_type == POLYGON && num_nodes >= 3 && num_nodes <= 12)
    return &typeInfoList[ reverseIndexList[msq_type][0] ];


  MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT, "VTK file does not support element type %d with %u nodes", (int)msq_type, num_nodes  );
  return 0;
}

void VtkTypeInfo::mesquiteToVtkOrder( std::vector<size_t>& conn_list ) const
{
  assert(conn_list.size() == numNodes);
  if (vtkConnOrder)
  {
    std::vector<size_t> temp_list(numNodes);
    std::swap( temp_list, conn_list );
    for (size_t i = 0; i < numNodes; ++i)
      conn_list[vtkConnOrder[i]] = temp_list[i];
  }
}

} // namespace Mesquite
