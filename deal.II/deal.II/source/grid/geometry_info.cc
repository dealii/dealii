//----------------------------  geometry_info.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  geometry_info.cc  ---------------------------


#include <grid/geometry_info.h>


// egcs 1.1 does not need these definitions of static member
// variables, as later compilers should need. on the other hand, if we
// define them, then egcs1.1 wants initialization, which we would have
// to mirror from the .h file. rather, we omit it here

#if ! ((__GNUC__==2) && (__GNUC_MINOR__ < 95))
const unsigned int GeometryInfo<deal_II_dimension>::vertices_per_cell;
const unsigned int GeometryInfo<deal_II_dimension>::lines_per_cell;
const unsigned int GeometryInfo<deal_II_dimension>::quads_per_cell;
const unsigned int GeometryInfo<deal_II_dimension>::hexes_per_cell;
const unsigned int GeometryInfo<deal_II_dimension>::children_per_cell;
#endif


template <>
const unsigned int GeometryInfo<1>::opposite_face[GeometryInfo<1>::faces_per_cell]
= { 0, 1 };


template <>
const unsigned int GeometryInfo<2>::opposite_face[GeometryInfo<2>::faces_per_cell]
= { 2, 3, 0, 1 };


template <>
const unsigned int GeometryInfo<3>::opposite_face[GeometryInfo<3>::faces_per_cell]
= { 1, 0, 4, 5, 2, 3 };

