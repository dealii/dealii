/* $Id$ */

#include <grid/geometry_info.h>


// enable these lines for gcc2.95
//
//const unsigned int GeometryInfo<deal_II_dimension>::vertices_per_cell;
//const unsigned int GeometryInfo<deal_II_dimension>::lines_per_cell;
//const unsigned int GeometryInfo<deal_II_dimension>::quads_per_cell;
//const unsigned int GeometryInfo<deal_II_dimension>::hexes_per_cell;
//const unsigned int GeometryInfo<deal_II_dimension>::children_per_cell;


template <>
const unsigned int GeometryInfo<1>::opposite_face[GeometryInfo<1>::faces_per_cell]
= { 0, 1 };


template <>
const unsigned int GeometryInfo<2>::opposite_face[GeometryInfo<2>::faces_per_cell]
= { 2, 3, 0, 1 };


template <>
const unsigned int GeometryInfo<3>::opposite_face[GeometryInfo<3>::faces_per_cell]
= { 1, 0, 4, 5, 2, 3 };

