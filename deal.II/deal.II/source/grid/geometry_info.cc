//----------------------------  geometry_info.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998-2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  geometry_info.cc  ---------------------------


#include <grid/geometry_info.h>



template <int dim> const unsigned int GeometryInfo<dim>::children_per_cell;
template <int dim> const unsigned int GeometryInfo<dim>::faces_per_cell;
template <int dim> const unsigned int GeometryInfo<dim>::subfaces_per_face;
template <int dim> const unsigned int GeometryInfo<dim>::vertices_per_cell;
template <int dim> const unsigned int GeometryInfo<dim>::vertices_per_face;
template <int dim> const unsigned int GeometryInfo<dim>::lines_per_face;
template <int dim> const unsigned int GeometryInfo<dim>::quads_per_face;
template <int dim> const unsigned int GeometryInfo<dim>::lines_per_cell;
template <int dim> const unsigned int GeometryInfo<dim>::quads_per_cell;
template <int dim> const unsigned int GeometryInfo<dim>::hexes_per_cell;


template <>
const unsigned int GeometryInfo<1>::opposite_face[GeometryInfo<1>::faces_per_cell]
= { 0, 1 };


template <>
const unsigned int GeometryInfo<2>::opposite_face[GeometryInfo<2>::faces_per_cell]
= { 2, 3, 0, 1 };


template <>
const unsigned int GeometryInfo<3>::opposite_face[GeometryInfo<3>::faces_per_cell]
= { 1, 0, 4, 5, 2, 3 };



template <>
const unsigned int GeometryInfo<1>::dx_to_deal[GeometryInfo<1>::vertices_per_cell]
= { 0, 1};


template <>
const unsigned int GeometryInfo<2>::dx_to_deal[GeometryInfo<2>::vertices_per_cell]
= { 0, 3, 1, 2};


template <>
const unsigned int GeometryInfo<3>::dx_to_deal[GeometryInfo<3>::vertices_per_cell]
= { 0, 3, 4, 7, 1, 2, 5, 6};

//TODO: Use these values in mappings.
template <>
const unsigned int
GeometryInfo<1>::unit_normal_direction[GeometryInfo<1>::faces_per_cell]
= { 0, 0 };
template <>
const int
GeometryInfo<1>::unit_normal_orientation[GeometryInfo<1>::faces_per_cell]
= { -1, 1 };

template <>
const unsigned int
GeometryInfo<2>::unit_normal_direction[GeometryInfo<2>::faces_per_cell]
= { 1, 0, 1, 0 };
template <>
const int
GeometryInfo<2>::unit_normal_orientation[GeometryInfo<2>::faces_per_cell]
= { -1, 1, 1, -1 };

template <>
const unsigned int
GeometryInfo<3>::unit_normal_direction[GeometryInfo<3>::faces_per_cell]
= { 1, 1, 2, 0, 2, 0 };
template <>
const int
GeometryInfo<3>::unit_normal_orientation[GeometryInfo<3>::faces_per_cell]
= { -1, 1, -1, 1, 1, -1 };



template <>
unsigned int
GeometryInfo<1>::child_cell_on_face (const unsigned int face,
                                     const unsigned int subface)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface==subfaces_per_face, ExcIndexRange(subface, 0, 1));

  return face;
}



template <>
unsigned int
GeometryInfo<2>::child_cell_on_face (const unsigned int face,
                                     const unsigned int subface)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface<subfaces_per_face, ExcIndexRange(subface, 0, subfaces_per_face));
  
  static const unsigned
    subcells[faces_per_cell][subfaces_per_face] = {{0,1},
                                                   {1,2},
                                                   {3,2},
                                                   {0,3}};
  return subcells[face][subface];
}



template <>
unsigned int
GeometryInfo<3>::child_cell_on_face (const unsigned int face,
				     const unsigned int subface)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface<subfaces_per_face, ExcIndexRange(subface, 0, subfaces_per_face));
  
  static const unsigned
    subcells[faces_per_cell][subfaces_per_face] = {{0, 1, 2, 3},
                                                   {4, 5, 6, 7},
                                                   {0, 1, 5, 4},
                                                   {1, 5, 6, 2},
                                                   {3, 2, 6, 7},
                                                   {0, 4, 7, 3}};
  return subcells[face][subface];
}



template <>
unsigned int
GeometryInfo<1>::line_to_cell_vertices (const unsigned int line,
					const unsigned int vertex)
{
  Assert (line<lines_per_cell, ExcIndexRange(line, 0, lines_per_cell));
  Assert (vertex<2, ExcIndexRange(vertex, 0, 2));

  return vertex;
}


template <>
unsigned int
GeometryInfo<2>::line_to_cell_vertices (const unsigned int line,
					const unsigned int vertex)
{
  return child_cell_on_face(line, vertex);
}



template <>
unsigned int
GeometryInfo<3>::line_to_cell_vertices (const unsigned int line,
					const unsigned int vertex)
{
  Assert (line<lines_per_cell, ExcIndexRange(line, 0, lines_per_cell));
  Assert (vertex<2, ExcIndexRange(vertex, 0, 2));
  
  static const unsigned
    vertices[lines_per_cell][2] = {{0, 1},  // front face
				   {1, 2},
				   {3, 2},
				   {0, 3},
				   {4, 5},  // back face
				   {5, 6},
				   {7, 6},
				   {4, 7},
				   {0, 4},  // connects of front and back face
				   {1, 5},
				   {2, 6},
				   {3, 7}};
  return vertices[line][vertex];
}



template <>
unsigned int
GeometryInfo<3>::face_to_cell_lines (const unsigned int face,
				     const unsigned int line)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (line<lines_per_face, ExcIndexRange(line, 0, lines_per_face));

  static const unsigned
    lines[faces_per_cell][lines_per_face] = {{0, 1, 2, 3},
					     {4, 5, 6, 7},
					     {0, 9, 4, 8},
					     {9, 5,10, 1},
					     {2,10, 6,11},
					     {8, 7,11, 3}};

  return lines[face][line];
}



template <int dim>
inline
unsigned int
GeometryInfo<dim>::face_to_cell_vertices (const unsigned int face,
					  const unsigned int vertex)
{
  return child_cell_on_face(face, vertex);
}




template class GeometryInfo<0>;
template class GeometryInfo<1>;
template class GeometryInfo<2>;
template class GeometryInfo<3>;
template class GeometryInfo<4>;
