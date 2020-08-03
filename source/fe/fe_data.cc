// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
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

#include <deal.II/base/geometry_info.h>

#include <deal.II/fe/fe.h>

DEAL_II_NAMESPACE_OPEN

template <int dim>
FiniteElementData<dim>::FiniteElementData(
  const std::vector<unsigned int> &dofs_per_object,
  const unsigned int               n_components,
  const unsigned int               degree,
  const Conformity                 conformity,
  const BlockIndices &             block_indices)
  : FiniteElementData(dofs_per_object,
                      dim == 0 ?
                        ReferenceCell::Type::Vertex :
                        (dim == 1 ? ReferenceCell::Type::Line :
                                    (dim == 2 ? ReferenceCell::Type::Quad :
                                                ReferenceCell::Type::Hex)),
                      n_components,
                      degree,
                      conformity,
                      block_indices)
{}

template <int dim>
FiniteElementData<dim>::FiniteElementData(
  const std::vector<unsigned int> &dofs_per_object,
  const ReferenceCell::Type        cell_type,
  const unsigned int               n_components,
  const unsigned int               degree,
  const Conformity                 conformity,
  const BlockIndices &             block_indices)
  : cell_type(cell_type)
  , dofs_per_vertex(dofs_per_object[0])
  , dofs_per_line(dofs_per_object[1])
  , dofs_per_quad(dim > 1 ? dofs_per_object[2] : 0)
  , dofs_per_quad_max(dofs_per_quad)
  , dofs_per_hex(dim > 2 ? dofs_per_object[3] : 0)
  , non_local_dofs_per_cell(
      dofs_per_object.size() == dim + 2 ? dofs_per_object[dim + 1] : 0)
  , first_line_index(
      ReferenceCell::internal::Info::get_cell(cell_type).n_vertices() *
      dofs_per_vertex)
  , first_quad_index(
      first_line_index +
      ReferenceCell::internal::Info::get_cell(cell_type).n_lines() *
        dofs_per_line)
  , first_hex_index(
      first_quad_index +
      (dim == 2 ?
         1 :
         (dim == 3 ?
            ReferenceCell::internal::Info::get_cell(cell_type).n_faces() :
            0)) *
        dofs_per_quad)
  , first_face_line_index(
      ReferenceCell::internal::Info::get_face(cell_type, 0).n_vertices() *
      dofs_per_vertex)
  , first_face_quad_index(
      (dim == 3 ?
         ReferenceCell::internal::Info::get_face(cell_type, 0).n_vertices() *
           dofs_per_vertex :
         ReferenceCell::internal::Info::get_cell(cell_type).n_vertices() *
           dofs_per_vertex) +
      ReferenceCell::internal::Info::get_face(cell_type, 0).n_lines() *
        dofs_per_line)
  , dofs_per_face(
      ReferenceCell::internal::Info::get_face(cell_type, 0).n_vertices() *
        dofs_per_vertex +
      ReferenceCell::internal::Info::get_face(cell_type, 0).n_lines() *
        dofs_per_line +
      (dim == 3 ? 1 : 0) * dofs_per_quad)
  , dofs_per_face_max(dofs_per_face)
  , dofs_per_cell(
      ReferenceCell::internal::Info::get_cell(cell_type).n_vertices() *
        dofs_per_vertex +
      ReferenceCell::internal::Info::get_cell(cell_type).n_lines() *
        dofs_per_line +
      (dim == 2 ?
         1 :
         (dim == 3 ?
            ReferenceCell::internal::Info::get_cell(cell_type).n_faces() :
            0)) *
        dofs_per_quad +
      (dim == 3 ? 1 : 0) * dofs_per_hex + non_local_dofs_per_cell)
  , components(n_components)
  , degree(degree)
  , conforming_space(conformity)
  , block_indices_data(block_indices.size() == 0 ?
                         BlockIndices(1, dofs_per_cell) :
                         block_indices)
{
  Assert(dofs_per_object.size() == dim + 1 || dofs_per_object.size() == dim + 2,
         ExcMessage("dofs_per_object should have size of either " +
                    std::to_string(dim + 1) + " or " +
                    std::to_string(dim + 2)));
}



template <int dim>
bool
FiniteElementData<dim>::operator==(const FiniteElementData<dim> &f) const
{
  return ((dofs_per_vertex == f.dofs_per_vertex) &&
          (dofs_per_line == f.dofs_per_line) &&
          (dofs_per_quad == f.dofs_per_quad) &&
          (dofs_per_hex == f.dofs_per_hex) && (components == f.components) &&
          (degree == f.degree) && (conforming_space == f.conforming_space) &&
          (non_local_dofs_per_cell == f.non_local_dofs_per_cell));
}


template class FiniteElementData<1>;
template class FiniteElementData<2>;
template class FiniteElementData<3>;

DEAL_II_NAMESPACE_CLOSE
