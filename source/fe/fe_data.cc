// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/geometry_info.h>

#include <deal.II/fe/fe.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  internal::GenericDoFsPerObject
  expand(const unsigned int               dim,
         const std::vector<unsigned int> &dofs_per_object,
         const ReferenceCell              cell_type)
  {
    internal::GenericDoFsPerObject result;

    const unsigned int face_no = 0;

    result.dofs_per_object_exclusive.resize(4, std::vector<unsigned int>(1));
    result.dofs_per_object_inclusive.resize(4, std::vector<unsigned int>(1));
    result.object_index.resize(4, std::vector<unsigned int>(1));
    result.first_object_index_on_face.resize(3, std::vector<unsigned int>(1));

    // dofs_per_vertex
    const unsigned int dofs_per_vertex     = dofs_per_object[0];
    result.dofs_per_object_exclusive[0][0] = dofs_per_vertex;

    // dofs_per_line
    const unsigned int dofs_per_line       = dofs_per_object[1];
    result.dofs_per_object_exclusive[1][0] = dofs_per_line;

    // dofs_per_quad
    const unsigned int dofs_per_quad       = dim > 1 ? dofs_per_object[2] : 0;
    result.dofs_per_object_exclusive[2][0] = dofs_per_quad;

    // dofs_per_hex
    const unsigned int dofs_per_hex        = dim > 2 ? dofs_per_object[3] : 0;
    result.dofs_per_object_exclusive[3][0] = dofs_per_hex;


    // first_line_index
    const unsigned int first_line_index =
      (cell_type.n_vertices() * dofs_per_vertex);
    result.object_index[1][0] = first_line_index;

    // first_quad_index
    const unsigned int first_quad_index =
      (first_line_index + cell_type.n_lines() * dofs_per_line);
    result.object_index[2][0] = first_quad_index;

    // first_hex_index
    result.object_index[3][0] =
      (first_quad_index +
       (dim == 2 ? 1 : (dim == 3 ? cell_type.n_faces() : 0)) * dofs_per_quad);

    // first_face_line_index
    result.first_object_index_on_face[1][0] =
      (cell_type.face_reference_cell(face_no).n_vertices() * dofs_per_vertex);

    // first_face_quad_index
    result.first_object_index_on_face[2][0] =
      ((dim == 3 ? cell_type.face_reference_cell(face_no).n_vertices() *
                     dofs_per_vertex :
                   cell_type.n_vertices() * dofs_per_vertex) +
       cell_type.face_reference_cell(face_no).n_lines() * dofs_per_line);

    // dofs_per_face
    result.dofs_per_object_inclusive[dim - 1][0] =
      (cell_type.face_reference_cell(face_no).n_vertices() * dofs_per_vertex +
       cell_type.face_reference_cell(face_no).n_lines() * dofs_per_line +
       (dim == 3 ? 1 : 0) * dofs_per_quad);


    // dofs_per_cell
    result.dofs_per_object_inclusive[dim][0] =
      (cell_type.n_vertices() * dofs_per_vertex +
       cell_type.n_lines() * dofs_per_line +
       (dim == 2 ? 1 : (dim == 3 ? cell_type.n_faces() : 0)) * dofs_per_quad +
       (dim == 3 ? 1 : 0) * dofs_per_hex);

    return result;
  }

  unsigned int
  number_unique_entries(const std::vector<unsigned int> &vector)
  {
    if (std::all_of(vector.begin(), vector.end(), [&](const auto &e) {
          return e == vector.front();
        }))
      {
        return 1;
      }
    else
      return vector.size();
  }
} // namespace internal



template <int dim>
FiniteElementData<dim>::FiniteElementData(
  const std::vector<unsigned int> &dofs_per_object,
  const unsigned int               n_components,
  const unsigned int               degree,
  const Conformity                 conformity,
  const BlockIndices              &block_indices)
  : FiniteElementData(dofs_per_object,
                      dim == 0 ?
                        ReferenceCells::Vertex :
                        (dim == 1 ? ReferenceCells::Line :
                                    (dim == 2 ? ReferenceCells::Quadrilateral :
                                                ReferenceCells::Hexahedron)),
                      n_components,
                      degree,
                      conformity,
                      block_indices)
{}

template <int dim>
FiniteElementData<dim>::FiniteElementData(
  const std::vector<unsigned int> &dofs_per_object,
  const ReferenceCell              cell_type,
  const unsigned int               n_components,
  const unsigned int               degree,
  const Conformity                 conformity,
  const BlockIndices              &block_indices)
  : FiniteElementData(internal::expand(dim, dofs_per_object, cell_type),
                      cell_type,
                      n_components,
                      degree,
                      conformity,
                      block_indices)
{}



template <int dim>
FiniteElementData<dim>::FiniteElementData(
  const internal::GenericDoFsPerObject &data,
  const ReferenceCell                   reference_cell,
  const unsigned int                    n_components,
  const unsigned int                    degree,
  const Conformity                      conformity,
  const BlockIndices                   &block_indices)
  : reference_cell_kind(reference_cell)
  , number_of_unique_2d_subobjects(
      internal::number_unique_entries(data.dofs_per_object_inclusive[2]))
  , number_unique_faces(
      internal::number_unique_entries(data.dofs_per_object_inclusive[dim - 1]))
  , dofs_per_vertex(data.dofs_per_object_exclusive[0][0])
  , dofs_per_line(data.dofs_per_object_exclusive[1][0])
  , n_dofs_on_quad(data.dofs_per_object_exclusive[2])
  , dofs_per_quad(n_dofs_on_quad[0])
  , dofs_per_quad_max(
      *max_element(n_dofs_on_quad.begin(), n_dofs_on_quad.end()))
  , dofs_per_hex(data.dofs_per_object_exclusive[3][0])
  , first_line_index(data.object_index[1][0])
  , first_index_of_quads(data.object_index[2])
  , first_quad_index(first_index_of_quads[0])
  , first_hex_index(data.object_index[3][0])
  , first_line_index_of_faces(data.first_object_index_on_face[1])
  , first_face_line_index(first_line_index_of_faces[0])
  , first_quad_index_of_faces(data.first_object_index_on_face[2])
  , first_face_quad_index(first_quad_index_of_faces[0])
  , n_dofs_on_face(data.dofs_per_object_inclusive[dim - 1])
  , dofs_per_face(n_dofs_on_face[0])
  , dofs_per_face_max(
      *max_element(n_dofs_on_face.begin(), n_dofs_on_face.end()))
  , dofs_per_cell(data.dofs_per_object_inclusive[dim][0])
  , components(n_components)
  , degree(degree)
  , conforming_space(conformity)
  , block_indices_data(block_indices.size() == 0 ?
                         BlockIndices(1, dofs_per_cell) :
                         block_indices)
{}



template <int dim>
bool
FiniteElementData<dim>::operator==(const FiniteElementData<dim> &f) const
{
  return ((dofs_per_vertex == f.dofs_per_vertex) &&
          (dofs_per_line == f.dofs_per_line) &&
          (dofs_per_quad == f.dofs_per_quad) &&
          (dofs_per_hex == f.dofs_per_hex) && (components == f.components) &&
          (degree == f.degree) && (conforming_space == f.conforming_space));
}


template class FiniteElementData<1>;
template class FiniteElementData<2>;
template class FiniteElementData<3>;

DEAL_II_NAMESPACE_CLOSE
