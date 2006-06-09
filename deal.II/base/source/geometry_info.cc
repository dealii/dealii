//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/geometry_info.h>



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


using namespace deal_II_numbers;

// make sure that also the icc compiler defines (and not only declares) 
// these variables
namespace internal
{
  void foo (const unsigned int *) {}

  template <int dim>
  void define_variables () 
  { 
    foo(&GeometryInfo<dim>::vertices_per_cell);
  }

  template void define_variables<2> ();
  template void define_variables<3> ();
  template void define_variables<4> ();
}



template <>
const unsigned int
GeometryInfo<1>::unit_normal_direction[faces_per_cell]
= { 0, 0 };

template <>
const unsigned int
GeometryInfo<2>::unit_normal_direction[faces_per_cell]
= { 0, 0, 1, 1 };

template <>
const unsigned int
GeometryInfo<3>::unit_normal_direction[faces_per_cell]
= { 0, 0, 1, 1, 2, 2 };

template <>
const unsigned int
GeometryInfo<4>::unit_normal_direction[faces_per_cell]
= { 0, 0, 1, 1, 2, 2, 3, 3 };



template <>
const int
GeometryInfo<1>::unit_normal_orientation[faces_per_cell]
= { -1, 1 };

template <>
const int
GeometryInfo<2>::unit_normal_orientation[faces_per_cell]
= { -1, 1, -1, 1 };

template <>
const int
GeometryInfo<3>::unit_normal_orientation[faces_per_cell]
= { -1, 1, -1, 1, -1, 1 };

template <>
const int
GeometryInfo<4>::unit_normal_orientation[faces_per_cell]
= { -1, 1, -1, 1, -1, 1, -1, 1 };



template <>
const unsigned int
GeometryInfo<1>::opposite_face[faces_per_cell]
= { 1, 0 };

template <>
const unsigned int
GeometryInfo<2>::opposite_face[faces_per_cell]
= { 1, 0, 3, 2 };

template <>
const unsigned int
GeometryInfo<3>::opposite_face[faces_per_cell]
= { 1, 0, 3, 2, 5, 4 };

template <>
const unsigned int
GeometryInfo<4>::opposite_face[faces_per_cell]
= { 1, 0, 3, 2, 5, 4, 7, 6 };



template <>
const unsigned int GeometryInfo<1>::ucd_to_deal[GeometryInfo<1>::vertices_per_cell]
= { 0, 1};

template <>
const unsigned int GeometryInfo<2>::ucd_to_deal[GeometryInfo<2>::vertices_per_cell]
= { 0, 1, 3, 2};

template <>
const unsigned int GeometryInfo<3>::ucd_to_deal[GeometryInfo<3>::vertices_per_cell]
= { 0, 1, 5, 4, 2, 3, 7, 6};

template <>
const unsigned int GeometryInfo<4>::ucd_to_deal[GeometryInfo<4>::vertices_per_cell]
= {  invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int};


template <>
const unsigned int GeometryInfo<1>::dx_to_deal[GeometryInfo<1>::vertices_per_cell]
= { 0, 1};

template <>
const unsigned int GeometryInfo<2>::dx_to_deal[GeometryInfo<2>::vertices_per_cell]
= { 0, 2, 1, 3};

template <>
const unsigned int GeometryInfo<3>::dx_to_deal[GeometryInfo<3>::vertices_per_cell]
= { 0, 4, 2, 6, 1, 5, 3, 7};

template <>
const unsigned int GeometryInfo<4>::dx_to_deal[GeometryInfo<4>::vertices_per_cell]
= {  invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int};

template <>
const unsigned int GeometryInfo<1>::vertex_to_face
          [GeometryInfo<1>::vertices_per_cell][1]
= { { 0 },
    { 1 } };

template <>
const unsigned int GeometryInfo<2>::vertex_to_face
          [GeometryInfo<2>::vertices_per_cell][2]
= { { 0, 2 },
    { 1, 2 },
    { 0, 3 },
    { 1, 3 } };

template <>
const unsigned int GeometryInfo<3>::vertex_to_face
          [GeometryInfo<3>::vertices_per_cell][3]
= { { 0, 2, 4 },
    { 1, 2, 4 },
    { 0, 3, 4 },
    { 1, 3, 4 },
    { 0, 2, 5 },
    { 1, 2, 5 },
    { 0, 3, 5 },
    { 1, 3, 5 } };

template <>
const unsigned int GeometryInfo<4>::vertex_to_face
          [GeometryInfo<4>::vertices_per_cell][4]
= { { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
    { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int }};

template <>
unsigned int
GeometryInfo<1>::child_cell_on_face (const unsigned int face,
                                     const unsigned int subface,
				     const bool)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface<subfaces_per_face,
	  ExcIndexRange(subface, 0, subfaces_per_face));

  return face;
}



template <>
unsigned int
GeometryInfo<2>::child_cell_on_face (const unsigned int face,
                                     const unsigned int subface,
				     const bool)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface<subfaces_per_face, ExcIndexRange(subface, 0, subfaces_per_face));
  
  static const unsigned
    subcells[faces_per_cell][subfaces_per_face] = {{0,2},
                                                   {1,3},
                                                   {0,1},
                                                   {2,3}};

  return subcells[face][subface];
}



template <>
unsigned int
GeometryInfo<3>::child_cell_on_face (const unsigned int face,
				     const unsigned int subface,
				     const bool face_orientation)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface<subfaces_per_face, ExcIndexRange(subface, 0, subfaces_per_face));

  static const unsigned int flip[subfaces_per_face] = { 0, 2, 1, 3 };

  static const unsigned
    subcells[faces_per_cell][subfaces_per_face] = {{0, 2, 4, 6},
                                                   {1, 3, 5, 7},
                                                   {0, 4, 1, 5},
                                                   {2, 6, 3, 7},
                                                   {0, 1, 2, 3},
                                                   {4, 5, 6, 7}};
  return face_orientation ? subcells[face][subface] : subcells[face][flip[subface]];
}



template <>
unsigned int
GeometryInfo<4>::child_cell_on_face (const unsigned int,
                                   const unsigned int,
                                   const bool)
{
  Assert(false, ExcNotImplemented());
  return invalid_unsigned_int;
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
    vertices[lines_per_cell][2] = {{0, 2},  // bottom face
				   {1, 3},
				   {0, 1},
				   {2, 3},
				   {4, 6},  // top face
				   {5, 7},
				   {4, 5},
				   {6, 7},
				   {0, 4},  // connects of bottom
				   {1, 5},  //   top face
				   {2, 6},
				   {3, 7}};
  return vertices[line][vertex];
}


template <>
unsigned int
GeometryInfo<4>::line_to_cell_vertices (const unsigned int,
                                      const unsigned int)
{
  Assert(false, ExcNotImplemented());
  return invalid_unsigned_int;
}


template <>
unsigned int
GeometryInfo<1>::face_to_cell_lines (const unsigned int face,
				     const unsigned int line,
				     const bool)
{
  Assert (face+1<faces_per_cell+1, ExcIndexRange(face, 0, faces_per_cell));
  Assert (line+1<lines_per_face+1, ExcIndexRange(line, 0, lines_per_face));
  
				   // There is only a single line, so
				   // it must be this.
  return 0;
}



template <>
unsigned int
GeometryInfo<2>::face_to_cell_lines (const unsigned int face,
				     const unsigned int line,
				     const bool)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (line<lines_per_face, ExcIndexRange(line, 0, lines_per_face));

				   // The face is a line itself.
  return face;
}



template <>
unsigned int
GeometryInfo<3>::face_to_cell_lines (const unsigned int face,
				     const unsigned int line,
				     const bool face_orientation)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (line<lines_per_face, ExcIndexRange(line, 0, lines_per_face));

				   // also represented by (line+2)%4
  static const unsigned int flip[lines_per_face] = { 2, 3, 0, 1 };

  static const unsigned
    lines[faces_per_cell][lines_per_face] = {{8,10, 0, 4}, // left face
					     {9,11, 1, 5}, // right face
					     {2, 6, 8, 9}, // front face
					     {3, 7,10,11}, // back face
					     {0, 1, 2, 3}, // bottom face
					     {4, 5, 6, 7}};// top face
  return face_orientation ? lines[face][line] : lines[face][flip[line]];
}



template<int dim>
unsigned int
GeometryInfo<dim>::face_to_cell_lines (const unsigned int,
				       const unsigned int,
				       const bool)
{
  Assert(false, ExcNotImplemented());
  return invalid_unsigned_int;
}



template <int dim>
unsigned int
GeometryInfo<dim>::face_to_cell_vertices (const unsigned int face,
					  const unsigned int vertex,
					  const bool face_orientation)
{
  return child_cell_on_face(face, vertex, face_orientation);
}



template <int dim>
Point<dim>
GeometryInfo<dim>::project_to_unit_cell (const Point<dim> &q)
{
  Point<dim> p = q;
  for(unsigned i=0; i<dim; i++)
    if      (p[i] < 0.)  p[i] = 0.;
    else if (p[i] > 1.)  p[i] = 1.;

  return p;
}



template <int dim>
double
GeometryInfo<dim>::distance_to_unit_cell (const Point<dim> &p)
{
   double result = 0.0;
   
   for(unsigned i=0; i<dim; i++)
      if ((-p[i]) > result)
         result = -p[i];
      else if ((p[i]-1.) > result)
         result = (p[i] - 1.);

   return result;
}

template class GeometryInfo<1>;
template class GeometryInfo<2>;
template class GeometryInfo<3>;
template class GeometryInfo<4>;
