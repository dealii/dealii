//----------------------------  geometry_info.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  geometry_info.cc  ---------------------------


#include <grid/geometry_info.h>


const unsigned int GeometryInfo<1>::dim;
const unsigned int GeometryInfo<1>::children_per_cell;
const unsigned int GeometryInfo<1>::faces_per_cell;
const unsigned int GeometryInfo<1>::subfaces_per_face;
const unsigned int GeometryInfo<1>::vertices_per_cell;
const unsigned int GeometryInfo<1>::vertices_per_face;
const unsigned int GeometryInfo<1>::lines_per_face;
const unsigned int GeometryInfo<1>::quads_per_face;
const unsigned int GeometryInfo<1>::lines_per_cell;
const unsigned int GeometryInfo<1>::quads_per_cell;
const unsigned int GeometryInfo<1>::hexes_per_cell;

const unsigned int GeometryInfo<2>::dim;
const unsigned int GeometryInfo<2>::children_per_cell;
const unsigned int GeometryInfo<2>::faces_per_cell;
const unsigned int GeometryInfo<2>::subfaces_per_face;
const unsigned int GeometryInfo<2>::vertices_per_cell;
const unsigned int GeometryInfo<2>::vertices_per_face;
const unsigned int GeometryInfo<2>::lines_per_face;
const unsigned int GeometryInfo<2>::quads_per_face;
const unsigned int GeometryInfo<2>::lines_per_cell;
const unsigned int GeometryInfo<2>::quads_per_cell;
const unsigned int GeometryInfo<2>::hexes_per_cell;

const unsigned int GeometryInfo<3>::dim;
const unsigned int GeometryInfo<3>::children_per_cell;
const unsigned int GeometryInfo<3>::faces_per_cell;
const unsigned int GeometryInfo<3>::subfaces_per_face;
const unsigned int GeometryInfo<3>::vertices_per_cell;
const unsigned int GeometryInfo<3>::vertices_per_face;
const unsigned int GeometryInfo<3>::lines_per_face;
const unsigned int GeometryInfo<3>::quads_per_face;
const unsigned int GeometryInfo<3>::lines_per_cell;
const unsigned int GeometryInfo<3>::quads_per_cell;
const unsigned int GeometryInfo<3>::hexes_per_cell;



const unsigned int GeometryInfo<1>::opposite_face[GeometryInfo<1>::faces_per_cell]
= { 0, 1 };


const unsigned int GeometryInfo<2>::opposite_face[GeometryInfo<2>::faces_per_cell]
= { 2, 3, 0, 1 };


const unsigned int GeometryInfo<3>::opposite_face[GeometryInfo<3>::faces_per_cell]
= { 1, 0, 4, 5, 2, 3 };





unsigned int
GeometryInfo<2>::child_cell_on_face (const unsigned int face,
                                     const unsigned int subface)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface<subfaces_per_face, ExcIndexRange(subface, 0, subfaces_per_face));
  
  const unsigned subcells[faces_per_cell][subfaces_per_face] = {{0,1},
								{1,2},
								{3,2},
								{0,3}};
  return subcells[face][subface];
};




unsigned int
GeometryInfo<3>::child_cell_on_face (const unsigned int face,
				     const unsigned int subface)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface<subfaces_per_face, ExcIndexRange(subface, 0, subfaces_per_face));
  
  const unsigned subcells[faces_per_cell][subfaces_per_face] = {{0, 1, 2, 3},
								{4, 5, 6, 7},
								{0, 1, 5, 4},
								{1, 5, 6, 2},
								{3, 2, 6, 7},
								{0, 4, 7, 3}};
  return subcells[face][subface];
};
