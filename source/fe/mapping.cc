// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
Mapping<dim, spacedim>::get_vertices(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell> vertices;
  for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
    {
      vertices[i] = cell->vertex(i);
    }
  return vertices;
}


template <int dim, int spacedim>
Point<dim - 1>
Mapping<dim, spacedim>::project_real_point_to_unit_point_on_face(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int &                                        face_no,
  const Point<spacedim> &                                     p) const
{
  // The function doesn't make physical sense for dim=1
  Assert(dim > 1, ExcNotImplemented());
  // Not implemented for higher dimensions
  Assert(dim <= 3, ExcNotImplemented());

  Point<dim> unit_cell_pt = transform_real_to_unit_cell(cell, p);

  Point<dim - 1> unit_face_pt;

  if (dim == 2)
    {
      if (GeometryInfo<dim>::unit_normal_direction[face_no] == 0)
        unit_face_pt = Point<dim - 1>(unit_cell_pt(1));
      else if (GeometryInfo<dim>::unit_normal_direction[face_no] == 1)
        unit_face_pt = Point<dim - 1>(unit_cell_pt(0));
    }
  else if (dim == 3)
    {
      if (GeometryInfo<dim>::unit_normal_direction[face_no] == 0)
        unit_face_pt = Point<dim - 1>(unit_cell_pt(1), unit_cell_pt(2));
      else if (GeometryInfo<dim>::unit_normal_direction[face_no] == 1)
        unit_face_pt = Point<dim - 1>(unit_cell_pt(0), unit_cell_pt(2));
      else if (GeometryInfo<dim>::unit_normal_direction[face_no] == 2)
        unit_face_pt = Point<dim - 1>(unit_cell_pt(0), unit_cell_pt(1));
    }

  return unit_face_pt;
}

/* ---------------------------- InternalDataBase --------------------------- */


template <int dim, int spacedim>
Mapping<dim, spacedim>::InternalDataBase::InternalDataBase() :
  update_each(update_default)
{}



template <int dim, int spacedim>
std::size_t
Mapping<dim, spacedim>::InternalDataBase::memory_consumption() const
{
  return sizeof(*this);
}


/*------------------------------ InternalData ------------------------------*/



// explicit instantiations
#include "mapping.inst"


DEAL_II_NAMESPACE_CLOSE
