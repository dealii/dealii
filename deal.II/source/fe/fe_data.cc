// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2001 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/geometry_info.h>
#include <deal.II/fe/fe.h>

DEAL_II_NAMESPACE_OPEN

template<int dim>
FiniteElementData<dim>::FiniteElementData ()
  :
  dofs_per_vertex(0),
  dofs_per_line(0),
  dofs_per_quad(0),
  dofs_per_hex(0),
  first_line_index(0),
  first_quad_index(0),
  first_hex_index(0),
  first_face_line_index(0),
  first_face_quad_index(0),
  dofs_per_face(0),
  dofs_per_cell (0),
  components(0),
  degree(0),
  conforming_space(unknown),
  cached_primitivity(false)
{}



template <int dim>
FiniteElementData<dim>::
FiniteElementData (const std::vector<unsigned int> &dofs_per_object,
                   const unsigned int n_components,
                   const unsigned int degree,
                   const Conformity conformity,
                   const unsigned int)
  :
  dofs_per_vertex(dofs_per_object[0]),
  dofs_per_line(dofs_per_object[1]),
  dofs_per_quad(dim>1? dofs_per_object[2]:0),
  dofs_per_hex(dim>2? dofs_per_object[3]:0),
  first_line_index(GeometryInfo<dim>::vertices_per_cell
                   * dofs_per_vertex),
  first_quad_index(first_line_index+
                   GeometryInfo<dim>::lines_per_cell
                   * dofs_per_line),
  first_hex_index(first_quad_index+
                  GeometryInfo<dim>::quads_per_cell
                  * dofs_per_quad),
  first_face_line_index(GeometryInfo<dim-1>::vertices_per_cell
                        * dofs_per_vertex),
  first_face_quad_index((dim==3 ?
                         GeometryInfo<dim-1>::vertices_per_cell
                         * dofs_per_vertex :
                         GeometryInfo<dim>::vertices_per_cell
                         * dofs_per_vertex) +
                        GeometryInfo<dim-1>::lines_per_cell
                        * dofs_per_line),
  dofs_per_face(GeometryInfo<dim>::vertices_per_face * dofs_per_vertex +
                GeometryInfo<dim>::lines_per_face * dofs_per_line +
                GeometryInfo<dim>::quads_per_face *dofs_per_quad),
  dofs_per_cell (GeometryInfo<dim>::vertices_per_cell * dofs_per_vertex +
                 GeometryInfo<dim>::lines_per_cell * dofs_per_line +
                 GeometryInfo<dim>::quads_per_cell * dofs_per_quad +
                 GeometryInfo<dim>::hexes_per_cell *dofs_per_hex),
  components(n_components),
  degree(degree),
  conforming_space(conformity),
  block_indices_data(1, dofs_per_cell)
{
  Assert(dofs_per_object.size()==dim+1, ExcDimensionMismatch(dofs_per_object.size()-1,dim));
}



template<int dim>
bool FiniteElementData<dim>::operator== (const FiniteElementData<dim> &f) const
{
  return ((dofs_per_vertex == f.dofs_per_vertex) &&
          (dofs_per_line == f.dofs_per_line) &&
          (dofs_per_quad == f.dofs_per_quad) &&
          (dofs_per_hex == f.dofs_per_hex) &&
          (components == f.components) &&
          (degree == f.degree) &&
          (conforming_space == f.conforming_space));
}

template<int dim>
unsigned int
FiniteElementData<dim>::
face_to_cell_index (const unsigned int face_index,
                    const unsigned int face,
                    const bool face_orientation,
                    const bool face_flip,
                    const bool face_rotation) const
{
  Assert (face_index < this->dofs_per_face,
          ExcIndexRange(face_index, 0, this->dofs_per_face));
  Assert (face < GeometryInfo<dim>::faces_per_cell,
          ExcIndexRange(face, 0, GeometryInfo<dim>::faces_per_cell));

  // see the function's documentation for an explanation of this
  // assertion -- in essence, derived classes have to implement
  // an overloaded version of this function if we are to use any
  // other than standard orientation
  if ((face_orientation != true) || (face_flip != false) || (face_rotation != false))
    Assert ((this->dofs_per_line <= 1) && (this->dofs_per_quad <= 1),
            ExcMessage ("The function in this base class can not handle this case. "
                        "Rather, the derived class you are using must provide "
                        "an overloaded version but apparently hasn't done so. See "
                        "the documentation of this function for more information."));

  // DoF on a vertex
  if (face_index < this->first_face_line_index)
    {
      // Vertex number on the face
      const unsigned int face_vertex = face_index / this->dofs_per_vertex;
      return face_index % this->dofs_per_vertex
             + GeometryInfo<dim>::face_to_cell_vertices(face, face_vertex,
                                                        face_orientation,
                                                        face_flip,
                                                        face_rotation)
             * this->dofs_per_vertex;
    }
  // Else, DoF on a line?
  if (face_index < this->first_face_quad_index)
    {
      // Ignore vertex dofs
      const unsigned int index = face_index - this->first_face_line_index;
      // Line number on the face
      const unsigned int face_line = index / this->dofs_per_line;
      return this->first_line_index + index % this->dofs_per_line
             + GeometryInfo<dim>::face_to_cell_lines(face, face_line,
                                                     face_orientation,
                                                     face_flip,
                                                     face_rotation)
             * this->dofs_per_line;
    }
  // Else, DoF is on a quad

  // Ignore vertex and line dofs
  const unsigned int index = face_index - this->first_face_quad_index;
  return this->first_quad_index + index
         + face * this->dofs_per_quad;
}


template class FiniteElementData<1>;
template class FiniteElementData<2>;
template class FiniteElementData<3>;

DEAL_II_NAMESPACE_CLOSE
