// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_p1.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include <algorithm>
#include <iostream>
#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace
{
  namespace VTKCellType
  {
    // Define VTK constants for linear, quadratic and
    // high-order Lagrange geometrices
    enum : unsigned int
    {
      VTK_VERTEX = 1,
      // Linear cells
      VTK_LINE       = 3,
      VTK_TRIANGLE   = 5,
      VTK_QUAD       = 9,
      VTK_TETRA      = 10,
      VTK_HEXAHEDRON = 12,
      VTK_WEDGE      = 13,
      VTK_PYRAMID    = 14,
      // Quadratic cells
      VTK_QUADRATIC_EDGE       = 21,
      VTK_QUADRATIC_TRIANGLE   = 22,
      VTK_QUADRATIC_QUAD       = 23,
      VTK_QUADRATIC_TETRA      = 24,
      VTK_QUADRATIC_HEXAHEDRON = 25,
      VTK_QUADRATIC_WEDGE      = 26,
      VTK_QUADRATIC_PYRAMID    = 27,
      // Lagrange cells
      VTK_LAGRANGE_CURVE         = 68,
      VTK_LAGRANGE_TRIANGLE      = 69,
      VTK_LAGRANGE_QUADRILATERAL = 70,
      VTK_LAGRANGE_TETRAHEDRON   = 71,
      VTK_LAGRANGE_HEXAHEDRON    = 72,
      VTK_LAGRANGE_WEDGE         = 73,
      VTK_LAGRANGE_PYRAMID       = 74,
      // Invalid code
      VTK_INVALID = numbers::invalid_unsigned_int
    };

  } // namespace VTKCellType

} // namespace

constexpr ndarray<unsigned int, 2, 2> ReferenceCell::line_vertex_permutations;

constexpr ndarray<unsigned int, 6, 3>
  ReferenceCell::triangle_vertex_permutations;

constexpr ndarray<unsigned int, 8, 4>
  ReferenceCell::quadrilateral_vertex_permutations;

std::string
ReferenceCell::to_string() const
{
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return "Vertex";
      case ReferenceCells::Line:
        return "Line";
      case ReferenceCells::Triangle:
        return "Tri";
      case ReferenceCells::Quadrilateral:
        return "Quad";
      case ReferenceCells::Tetrahedron:
        return "Tet";
      case ReferenceCells::Pyramid:
        return "Pyramid";
      case ReferenceCells::Wedge:
        return "Wedge";
      case ReferenceCells::Hexahedron:
        return "Hex";
      case ReferenceCells::Invalid:
        return "Invalid";
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }

  return "Invalid";
}



template <int dim>
std::pair<unsigned int, RefinementCase<dim - 1>>
ReferenceCell::equivalent_refinement_case(
  const types::geometric_orientation combined_face_orientation,
  const internal::SubfaceCase<dim>   subface_case,
  const unsigned int                 subface_no) const
{
  if constexpr (dim == 3)
    {
      static const RefinementCase<dim - 1>
        equivalent_refine_case[internal::SubfaceCase<dim>::case_isotropic + 1]
                              [GeometryInfo<3>::max_children_per_face] = {
                                // case_none. there should be only
                                // invalid values here. However, as
                                // this function is also called (in
                                // tests) for cells which have no
                                // refined faces, use isotropic
                                // refinement instead
                                {RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy},
                                // case_x
                                {RefinementCase<dim - 1>::cut_x,
                                 RefinementCase<dim - 1>::cut_x,
                                 RefinementCase<dim - 1>::no_refinement,
                                 RefinementCase<dim - 1>::no_refinement},
                                // case_x1y
                                {RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_x,
                                 RefinementCase<dim - 1>::no_refinement},
                                // case_x2y
                                {RefinementCase<dim - 1>::cut_x,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::no_refinement},
                                // case_x1y2y
                                {RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy},
                                // case_y
                                {RefinementCase<dim - 1>::cut_y,
                                 RefinementCase<dim - 1>::cut_y,
                                 RefinementCase<dim - 1>::no_refinement,
                                 RefinementCase<dim - 1>::no_refinement},
                                // case_y1x
                                {RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_y,
                                 RefinementCase<dim - 1>::no_refinement},
                                // case_y2x
                                {RefinementCase<dim - 1>::cut_y,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::no_refinement},
                                // case_y1x2x
                                {RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy},
                                // case_xy (case_isotropic)
                                {RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy,
                                 RefinementCase<dim - 1>::cut_xy}};

      constexpr unsigned int X = numbers::invalid_unsigned_int;
      static const unsigned int
        equivalent_subface_number[internal::SubfaceCase<dim>::case_isotropic +
                                  1][GeometryInfo<3>::max_children_per_face] = {
          // case_none, see above
          {0, 1, 2, 3},
          // case_x
          {0, 1, X, X},
          // case_x1y
          {0, 2, 1, X},
          // case_x2y
          {0, 1, 3, X},
          // case_x1y2y
          {0, 2, 1, 3},
          // case_y
          {0, 1, X, X},
          // case_y1x
          {0, 1, 1, X},
          // case_y2x
          {0, 2, 3, X},
          // case_y1x2x
          {0, 1, 2, 3},
          // case_xy (case_isotropic)
          {0, 1, 2, 3}};

      static const RefinementCase<dim - 1> rotated_refinement_case[4] = {
        RefinementCase<dim - 1>::no_refinement,
        RefinementCase<dim - 1>::cut_y,
        RefinementCase<dim - 1>::cut_x,
        RefinementCase<dim - 1>::cut_xy};
      const auto [face_orientation, face_rotation, face_flip] =
        internal::split_face_orientation(combined_face_orientation);

      const auto equivalent_refinement_case =
        equivalent_refine_case[subface_case][subface_no];
      const unsigned int equivalent_subface_no =
        equivalent_subface_number[subface_case][subface_no];
      // make sure, that we got a valid subface and RefineCase
      Assert(equivalent_refinement_case != RefinementCase<dim>::no_refinement,
             ExcInternalError());
      Assert(equivalent_subface_no != X, ExcInternalError());
      // now, finally respect non-standard faces
      const RefinementCase<dim - 1> final_refinement_case =
        (face_orientation == face_rotation ?
           rotated_refinement_case[equivalent_refinement_case] :
           equivalent_refinement_case);

      const unsigned int final_subface_no =
        GeometryInfo<dim>::child_cell_on_face(RefinementCase<dim>(
                                                final_refinement_case),
                                              /*face_no = */ 4,
                                              equivalent_subface_no,
                                              face_orientation,
                                              face_flip,
                                              face_rotation,
                                              equivalent_refinement_case);

      return std::make_pair(final_subface_no, final_refinement_case);
    }
  else
    {
      (void)combined_face_orientation;
      (void)subface_case;
      (void)subface_no;

      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
ReferenceCell::get_default_mapping(const unsigned int degree) const
{
  AssertDimension(dim, get_dimension());

  if (is_hyper_cube())
    return std::make_unique<MappingQ<dim, spacedim>>(degree);
  else if (is_simplex())
    return std::make_unique<MappingFE<dim, spacedim>>(
      FE_SimplexP<dim, spacedim>(degree));
  else if (*this == ReferenceCells::Pyramid)
    return std::make_unique<MappingFE<dim, spacedim>>(
      FE_PyramidP<dim, spacedim>(degree));
  else if (*this == ReferenceCells::Wedge)
    return std::make_unique<MappingFE<dim, spacedim>>(
      FE_WedgeP<dim, spacedim>(degree));
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
    }

  return std::make_unique<MappingQ<dim, spacedim>>(degree);
}



template <int dim, int spacedim>
const Mapping<dim, spacedim> &
ReferenceCell::get_default_linear_mapping() const
{
  AssertDimension(dim, get_dimension());

  if (is_hyper_cube())
    {
      return StaticMappingQ1<dim, spacedim>::mapping;
    }
  else if (is_simplex())
    {
      static const MappingP1<dim, spacedim> mapping;
      return mapping;
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      static const MappingFE<dim, spacedim> mapping(
        FE_PyramidP<dim, spacedim>(1));
      return mapping;
    }
  else if (*this == ReferenceCells::Wedge)
    {
      static const MappingFE<dim, spacedim> mapping(
        FE_WedgeP<dim, spacedim>(1));
      return mapping;
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
    }

  return StaticMappingQ1<dim, spacedim>::mapping; // never reached
}



template <int dim>
Quadrature<dim>
ReferenceCell::get_gauss_type_quadrature(const unsigned n_points_1d) const
{
  AssertDimension(dim, get_dimension());

  if (is_hyper_cube())
    return QGauss<dim>(n_points_1d);
  else if (is_simplex())
    return QGaussSimplex<dim>(n_points_1d);
  else if (*this == ReferenceCells::Pyramid)
    return QGaussPyramid<dim>(n_points_1d);
  else if (*this == ReferenceCells::Wedge)
    return QGaussWedge<dim>(n_points_1d);
  else
    DEAL_II_NOT_IMPLEMENTED();

  return Quadrature<dim>(); // never reached
}



template <int dim>
const Quadrature<dim> &
ReferenceCell::get_nodal_type_quadrature() const
{
  AssertDimension(dim, get_dimension());

  // A function that is used to fill a quadrature object of the
  // desired type the first time we encounter a particular
  // reference cell
  const auto create_quadrature = [](const ReferenceCell &reference_cell) {
    std::vector<Point<dim>> vertices(reference_cell.n_vertices());
    for (const unsigned int v : reference_cell.vertex_indices())
      vertices[v] = reference_cell.vertex<dim>(v);

    return Quadrature<dim>(vertices);
  };

  if (is_hyper_cube())
    {
      static const Quadrature<dim> quadrature = create_quadrature(*this);
      return quadrature;
    }
  else if (is_simplex())
    {
      static const Quadrature<dim> quadrature = create_quadrature(*this);
      return quadrature;
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      static const Quadrature<dim> quadrature = create_quadrature(*this);
      return quadrature;
    }
  else if (*this == ReferenceCells::Wedge)
    {
      static const Quadrature<dim> quadrature = create_quadrature(*this);
      return quadrature;
    }
  else
    DEAL_II_NOT_IMPLEMENTED();

  static const Quadrature<dim> dummy;
  return dummy; // never reached
}



unsigned int
ReferenceCell::exodusii_vertex_to_deal_vertex(const unsigned int vertex_n) const
{
  AssertIndexRange(vertex_n, n_vertices());

  switch (this->kind)
    {
      case ReferenceCells::Line:
      case ReferenceCells::Triangle:
        return vertex_n;
      case ReferenceCells::Quadrilateral:
        {
          constexpr std::array<unsigned int, 4> exodus_to_deal{{0, 1, 3, 2}};
          return exodus_to_deal[vertex_n];
        }
      case ReferenceCells::Tetrahedron:
        return vertex_n;
      case ReferenceCells::Hexahedron:
        {
          constexpr std::array<unsigned int, 8> exodus_to_deal{
            {0, 1, 3, 2, 4, 5, 7, 6}};
          return exodus_to_deal[vertex_n];
        }
      case ReferenceCells::Wedge:
        {
          constexpr std::array<unsigned int, 6> exodus_to_deal{
            {2, 1, 0, 5, 4, 3}};
          return exodus_to_deal[vertex_n];
        }
      case ReferenceCells::Pyramid:
        {
          constexpr std::array<unsigned int, 5> exodus_to_deal{{0, 1, 3, 2, 4}};
          return exodus_to_deal[vertex_n];
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }

  return numbers::invalid_unsigned_int;
}



unsigned int
ReferenceCell::exodusii_face_to_deal_face(const unsigned int face_n) const
{
  AssertIndexRange(face_n, n_faces());

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return 0;
      case ReferenceCells::Line:
      case ReferenceCells::Triangle:
        return face_n;
      case ReferenceCells::Quadrilateral:
        {
          constexpr std::array<unsigned int, 4> exodus_to_deal{{2, 1, 3, 0}};
          return exodus_to_deal[face_n];
        }
      case ReferenceCells::Tetrahedron:
        {
          constexpr std::array<unsigned int, 4> exodus_to_deal{{1, 3, 2, 0}};
          return exodus_to_deal[face_n];
        }
      case ReferenceCells::Hexahedron:
        {
          constexpr std::array<unsigned int, 6> exodus_to_deal{
            {2, 1, 3, 0, 4, 5}};
          return exodus_to_deal[face_n];
        }
      case ReferenceCells::Wedge:
        {
          constexpr std::array<unsigned int, 6> exodus_to_deal{{3, 4, 2, 0, 1}};
          return exodus_to_deal[face_n];
        }
      case ReferenceCells::Pyramid:
        {
          constexpr std::array<unsigned int, 5> exodus_to_deal{{3, 2, 4, 1, 0}};
          return exodus_to_deal[face_n];
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }

  return numbers::invalid_unsigned_int;
}



unsigned int
ReferenceCell::unv_vertex_to_deal_vertex(const unsigned int vertex_n) const
{
  AssertIndexRange(vertex_n, n_vertices());
  // Information on this file format isn't easy to find - the documents here
  //
  // https://www.ceas3.uc.edu/sdrluff/
  //
  // Don't actually explain anything about the sections we care about (2412) in
  // any detail. For node numbering I worked backwards from what is actually in
  // our test files (since that's supposed to work), which all use some
  // non-standard clockwise numbering scheme which starts at the bottom right
  // vertex.
  if (*this == ReferenceCells::Line)
    {
      return vertex_n;
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      constexpr std::array<unsigned int, 4> unv_to_deal{{1, 0, 2, 3}};
      return unv_to_deal[vertex_n];
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      constexpr std::array<unsigned int, 8> unv_to_deal{
        {6, 7, 5, 4, 2, 3, 1, 0}};
      return unv_to_deal[vertex_n];
    }

  DEAL_II_NOT_IMPLEMENTED();

  return numbers::invalid_unsigned_int;
}



unsigned int
ReferenceCell::vtk_linear_type() const
{
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return VTKCellType::VTK_VERTEX;
      case ReferenceCells::Line:
        return VTKCellType::VTK_LINE;
      case ReferenceCells::Triangle:
        return VTKCellType::VTK_TRIANGLE;
      case ReferenceCells::Quadrilateral:
        return VTKCellType::VTK_QUAD;
      case ReferenceCells::Tetrahedron:
        return VTKCellType::VTK_TETRA;
      case ReferenceCells::Pyramid:
        return VTKCellType::VTK_PYRAMID;
      case ReferenceCells::Wedge:
        return VTKCellType::VTK_WEDGE;
      case ReferenceCells::Hexahedron:
        return VTKCellType::VTK_HEXAHEDRON;
      case ReferenceCells::Invalid:
        return VTKCellType::VTK_INVALID;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }

  return VTKCellType::VTK_INVALID;
}



unsigned int
ReferenceCell::vtk_quadratic_type() const
{
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return VTKCellType::VTK_VERTEX;
      case ReferenceCells::Line:
        return VTKCellType::VTK_QUADRATIC_EDGE;
      case ReferenceCells::Triangle:
        return VTKCellType::VTK_QUADRATIC_TRIANGLE;
      case ReferenceCells::Quadrilateral:
        return VTKCellType::VTK_QUADRATIC_QUAD;
      case ReferenceCells::Tetrahedron:
        return VTKCellType::VTK_QUADRATIC_TETRA;
      case ReferenceCells::Pyramid:
        return VTKCellType::VTK_QUADRATIC_PYRAMID;
      case ReferenceCells::Wedge:
        return VTKCellType::VTK_QUADRATIC_WEDGE;
      case ReferenceCells::Hexahedron:
        return VTKCellType::VTK_QUADRATIC_HEXAHEDRON;
      case ReferenceCells::Invalid:
        return VTKCellType::VTK_INVALID;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }

  return VTKCellType::VTK_INVALID;
}



unsigned int
ReferenceCell::vtk_lagrange_type() const
{
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return VTKCellType::VTK_VERTEX;
      case ReferenceCells::Line:
        return VTKCellType::VTK_LAGRANGE_CURVE;
      case ReferenceCells::Triangle:
        return VTKCellType::VTK_LAGRANGE_TRIANGLE;
      case ReferenceCells::Quadrilateral:
        return VTKCellType::VTK_LAGRANGE_QUADRILATERAL;
      case ReferenceCells::Tetrahedron:
        return VTKCellType::VTK_LAGRANGE_TETRAHEDRON;
      case ReferenceCells::Pyramid:
        return VTKCellType::VTK_LAGRANGE_PYRAMID;
      case ReferenceCells::Wedge:
        return VTKCellType::VTK_LAGRANGE_WEDGE;
      case ReferenceCells::Hexahedron:
        return VTKCellType::VTK_LAGRANGE_HEXAHEDRON;
      case ReferenceCells::Invalid:
        return VTKCellType::VTK_INVALID;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }

  return VTKCellType::VTK_INVALID;
}



template <>
unsigned int
ReferenceCell::vtk_lexicographic_to_node_index<0>(
  const std::array<unsigned, 0> &,
  const std::array<unsigned, 0> &,
  const bool) const
{
  DEAL_II_NOT_IMPLEMENTED();
  return 0;
}



/**
 * Modified from
 * https://github.com/Kitware/VTK/blob/265ca48a79a36538c95622c237da11133608bbe5/Common/DataModel/vtkLagrangeCurve.cxx#L478
 */
template <>
unsigned int
ReferenceCell::vtk_lexicographic_to_node_index<1>(
  const std::array<unsigned, 1> &node_indices,
  const std::array<unsigned, 1> &nodes_per_direction,
  const bool) const
{
  const unsigned int i = node_indices[0];

  const bool ibdy = (i == 0 || i == nodes_per_direction[0]);
  // How many boundaries do we lie on at once?
  const int nbdy = (ibdy ? 1 : 0);

  if (nbdy == 1) // Vertex DOF
    { // ijk is a corner node. Return the proper index (somewhere in [0,7]):
      return i ? 1 : 0;
    }

  const int offset = 2;
  return (i - 1) + offset;
}



/**
 * Modified from
 * https://github.com/Kitware/VTK/blob/265ca48a/Common/DataModel/vtkLagrangeQuadrilateral.cxx#L558
 */
template <>
unsigned int
ReferenceCell::vtk_lexicographic_to_node_index<2>(
  const std::array<unsigned, 2> &node_indices,
  const std::array<unsigned, 2> &nodes_per_direction,
  const bool) const
{
  Assert(*this == ReferenceCells::Quadrilateral, ExcNotImplemented());

  const unsigned int i = node_indices[0];
  const unsigned int j = node_indices[1];

  const bool ibdy = (i == 0 || i == nodes_per_direction[0]);
  const bool jbdy = (j == 0 || j == nodes_per_direction[1]);
  // How many boundaries do we lie on at once?
  const int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0);

  if (nbdy == 2) // Vertex DOF
    { // ijk is a corner node. Return the proper index (somewhere in [0,3]):
      return (i != 0u ? (j != 0u ? 2 : 1) : (j != 0u ? 3 : 0));
    }

  int offset = 4;
  if (nbdy == 1) // Edge DOF
    {
      if (!ibdy)
        { // On i axis
          return (i - 1) +
                 (j != 0u ?
                    nodes_per_direction[0] - 1 + nodes_per_direction[1] - 1 :
                    0) +
                 offset;
        }

      if (!jbdy)
        { // On j axis
          return (j - 1) +
                 (i != 0u ? nodes_per_direction[0] - 1 :
                            2 * (nodes_per_direction[0] - 1) +
                              nodes_per_direction[1] - 1) +
                 offset;
        }
    }

  offset += 2 * (nodes_per_direction[0] - 1 + nodes_per_direction[1] - 1);
  // nbdy == 0: Face DOF
  return offset + (i - 1) + (nodes_per_direction[0] - 1) * ((j - 1));
}



/**
 * Modified from
 * https://github.com/Kitware/VTK/blob/265ca48a/Common/DataModel/vtkLagrangeHexahedron.cxx#L734
 * (legacy_format=true) and from
 * https://github.com/Kitware/VTK/blob/256fe70de00e3441f126276ca4a8c5477d0bcb86/Common/DataModel/vtkHigherOrderHexahedron.cxx#L593
 * (legacy_format=false). The two versions differ regarding the ordering of
 * lines 10 and 11 (clockwise vs. anti-clockwise). See also:
 * https://github.com/Kitware/VTK/blob/7a0b92864c96680b1f42ee84920df556fc6ebaa3/Documentation/release/dev/node-numbering-change-for-VTK_LAGRANGE_HEXAHEDRON.md
 *
 */
template <>
unsigned int
ReferenceCell::vtk_lexicographic_to_node_index<3>(
  const std::array<unsigned, 3> &node_indices,
  const std::array<unsigned, 3> &nodes_per_direction,
  const bool                     legacy_format) const
{
  Assert(*this == ReferenceCells::Hexahedron, ExcNotImplemented());

  const unsigned int i = node_indices[0];
  const unsigned int j = node_indices[1];
  const unsigned int k = node_indices[2];

  const bool ibdy = (i == 0 || i == nodes_per_direction[0]);
  const bool jbdy = (j == 0 || j == nodes_per_direction[1]);
  const bool kbdy = (k == 0 || k == nodes_per_direction[2]);
  // How many boundaries do we lie on at once?
  const int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);

  if (nbdy == 3) // Vertex DOF
    { // ijk is a corner node. Return the proper index (somewhere in [0,7]):
      return (i != 0u ? (j != 0u ? 2 : 1) : (j != 0u ? 3 : 0)) +
             (k != 0u ? 4 : 0);
    }

  int offset = 8;
  if (nbdy == 2) // Edge DOF
    {
      if (!ibdy)
        { // On i axis
          return (i - 1) +
                 (j != 0u ?
                    nodes_per_direction[0] - 1 + nodes_per_direction[1] - 1 :
                    0) +
                 (k != 0u ? 2 * (nodes_per_direction[0] - 1 +
                                 nodes_per_direction[1] - 1) :
                            0) +
                 offset;
        }
      if (!jbdy)
        { // On j axis
          return (j - 1) +
                 (i != 0u ? nodes_per_direction[0] - 1 :
                            2 * (nodes_per_direction[0] - 1) +
                              nodes_per_direction[1] - 1) +
                 (k != 0u ? 2 * (nodes_per_direction[0] - 1 +
                                 nodes_per_direction[1] - 1) :
                            0) +
                 offset;
        }
      // !kbdy, On k axis
      offset +=
        4 * (nodes_per_direction[0] - 1) + 4 * (nodes_per_direction[1] - 1);
      if (legacy_format)
        return (k - 1) +
               (nodes_per_direction[2] - 1) *
                 (i != 0u ? (j != 0u ? 3 : 1) : (j != 0u ? 2 : 0)) +
               offset;
      else
        return (k - 1) +
               (nodes_per_direction[2] - 1) *
                 (i != 0u ? (j != 0u ? 2 : 1) : (j != 0u ? 3 : 0)) +
               offset;
    }

  offset += 4 * (nodes_per_direction[0] - 1 + nodes_per_direction[1] - 1 +
                 nodes_per_direction[2] - 1);
  if (nbdy == 1) // Face DOF
    {
      if (ibdy) // On i-normal face
        {
          return (j - 1) + ((nodes_per_direction[1] - 1) * (k - 1)) +
                 (i != 0u ? (nodes_per_direction[1] - 1) *
                              (nodes_per_direction[2] - 1) :
                            0) +
                 offset;
        }
      offset += 2 * (nodes_per_direction[1] - 1) * (nodes_per_direction[2] - 1);
      if (jbdy) // On j-normal face
        {
          return (i - 1) + ((nodes_per_direction[0] - 1) * (k - 1)) +
                 (j != 0u ? (nodes_per_direction[2] - 1) *
                              (nodes_per_direction[0] - 1) :
                            0) +
                 offset;
        }
      offset += 2 * (nodes_per_direction[2] - 1) * (nodes_per_direction[0] - 1);
      // kbdy, On k-normal face
      return (i - 1) + ((nodes_per_direction[0] - 1) * (j - 1)) +
             (k != 0u ?
                (nodes_per_direction[0] - 1) * (nodes_per_direction[1] - 1) :
                0) +
             offset;
    }

  // nbdy == 0: Body DOF
  offset += 2 * ((nodes_per_direction[1] - 1) * (nodes_per_direction[2] - 1) +
                 (nodes_per_direction[2] - 1) * (nodes_per_direction[0] - 1) +
                 (nodes_per_direction[0] - 1) * (nodes_per_direction[1] - 1));
  return offset + (i - 1) +
         (nodes_per_direction[0] - 1) *
           ((j - 1) + (nodes_per_direction[1] - 1) * ((k - 1)));
}



unsigned int
ReferenceCell::vtk_vertex_to_deal_vertex(const unsigned int vertex_index) const
{
  AssertIndexRange(vertex_index, n_vertices());

  // For some of the following, deal.II uses the same ordering as VTK
  // and in that case, we only need to return 'vertex_index' (i.e.,
  // use the identity mapping). For some others, we need to translate.
  //
  // For the ordering, see the VTK manual (for example at
  // http://www.princeton.edu/~efeibush/viscourse/vtk.pdf, page 9).
  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return vertex_index;
      case ReferenceCells::Line:
        return vertex_index;
      case ReferenceCells::Triangle:
        return vertex_index;
      case ReferenceCells::Quadrilateral:
        {
          static constexpr std::array<unsigned int, 4> index_translation_table =
            {{0, 1, 3, 2}};
          return index_translation_table[vertex_index];
        }
      case ReferenceCells::Tetrahedron:
        return vertex_index;
      case ReferenceCells::Pyramid:
        {
          static constexpr std::array<unsigned int, 5> index_translation_table =
            {{0, 1, 3, 2, 4}};
          return index_translation_table[vertex_index];
        }
      case ReferenceCells::Wedge:
        return vertex_index;
      case ReferenceCells::Hexahedron:
        {
          static constexpr std::array<unsigned int, 8> index_translation_table =
            {{0, 1, 3, 2, 4, 5, 7, 6}};
          return index_translation_table[vertex_index];
        }
      case ReferenceCells::Invalid:
        {
          DEAL_II_NOT_IMPLEMENTED();
          return numbers::invalid_unsigned_int;
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }

  return numbers::invalid_unsigned_int;
}



unsigned int
ReferenceCell::gmsh_element_type() const
{
  /*
    From the GMSH documentation:

    elm-type
    defines the geometrical type of the n-th element:

    1
    Line (2 nodes).

    2
    Triangle (3 nodes).

    3
    Quadrangle (4 nodes).

    4
    Tetrahedron (4 nodes).

    5
    Hexahedron (8 nodes).

    6
    Prism (6 nodes).

    7
    Pyramid (5 nodes).

    8
    Second order line (3 nodes: 2 associated with the vertices and 1 with the
    edge).

    9
    Second order triangle (6 nodes: 3 associated with the vertices and 3 with
    the edges).

    10 Second order quadrangle (9 nodes: 4 associated with the
    vertices, 4 with the edges and 1 with the face).

    11 Second order tetrahedron (10 nodes: 4 associated with the vertices and 6
    with the edges).

    12 Second order hexahedron (27 nodes: 8 associated with the vertices, 12
    with the edges, 6 with the faces and 1 with the volume).

    13 Second order prism (18 nodes: 6 associated with the vertices, 9 with the
    edges and 3 with the quadrangular faces).

    14 Second order pyramid (14 nodes: 5 associated with the vertices, 8 with
    the edges and 1 with the quadrangular face).

    15 Point (1 node).
  */

  switch (this->kind)
    {
      case ReferenceCells::Vertex:
        return 15;
      case ReferenceCells::Line:
        return 1;
      case ReferenceCells::Triangle:
        return 2;
      case ReferenceCells::Quadrilateral:
        return 3;
      case ReferenceCells::Tetrahedron:
        return 4;
      case ReferenceCells::Pyramid:
        return 7;
      case ReferenceCells::Wedge:
        return 6;
      case ReferenceCells::Hexahedron:
        return 5;
      case ReferenceCells::Invalid:
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }

  return numbers::invalid_unsigned_int;
}



namespace
{
  // Compute the nearest point to @p on the line segment and the square of its
  // distance to @p.
  template <int dim>
  std::pair<Point<dim>, double>
  project_to_line(const Point<dim> &x0,
                  const Point<dim> &x1,
                  const Point<dim> &p)
  {
    Assert(x0 != x1, ExcInternalError());
    // t is the convex combination coefficient (x = (1 - t) * x0 + t * x1)
    // defining the position of the closest point on the line (not line segment)
    // to p passing through x0 and x1. This formula is equivalent to the
    // standard 'project a vector onto another vector', where each vector is
    // shifted to start at x0.
    const double t = ((x1 - x0) * (p - x0)) / ((x1 - x0).norm_square());

    // Only consider points between x0 and x1
    if (t <= 0)
      return std::make_pair(x0, x0.distance_square(p));
    else if (t <= 1)
      {
        const auto p2 = x0 + t * (x1 - x0);
        return std::make_pair(p2, p2.distance_square(p));
      }
    else
      return std::make_pair(x1, x1.distance_square(p));
  }

  // template base-case
  template <int dim>
  std::pair<Point<dim>, double>
  project_to_quad(const std::array<Point<dim>, 3> & /*vertices*/,
                  const Point<dim> & /*p*/,
                  const ReferenceCell /*reference_cell*/)
  {
    DEAL_II_ASSERT_UNREACHABLE();
    return std::make_pair(Point<dim>(),
                          std::numeric_limits<double>::signaling_NaN());
  }

  /**
   * Compute the nearest point on a 2d subobject (something with structdim = 2
   * and spacedim = 3) and the square of that point's distance to @p p. Here,
   * the subobject is described with three vertices: either the three vertices
   * of a Triangle or the first three of a Quadrilateral (as the fourth one
   * can be computed, in that case, from the first three).
   *
   * If the given point cannot be projected via a normal vector (i.e., if the
   * line parallel to the normal vector intersecting @p p does not intersect the
   * subobject) then this function returns the origin and the largest double
   * precision number. distance_to_line_square() is in charge of computing the
   * distance to lines.
   *
   * @note This function is for Quadrilaterals and Triangles because it is
   * intended to find the shortest distance to the face of a reference cell
   * (which must be a Triangle or Quadrilateral) in 3d.
   */
  template <>
  std::pair<Point<3>, double>
  project_to_quad(const std::array<Point<3>, 3> &vertices,
                  const Point<3>                &p,
                  const ReferenceCell            face_reference_cell)
  {
    Assert(face_reference_cell == ReferenceCells::Triangle ||
             face_reference_cell == ReferenceCells::Quadrilateral,
           ExcNotImplemented());

    // Make the problem slightly easier by shifting everything to avoid a point
    // at the origin (this way we can invert the matrix of vertices). Use 2.0 so
    // that the bottom left vertex of a Pyramid is now at x = 1.
    std::array<Point<3>, 3> shifted_vertices = vertices;
    const Tensor<1, 3>      shift{{2.0, 2.0, 2.0}};
    for (Point<3> &shifted_vertex : shifted_vertices)
      shifted_vertex += shift;
    const Point<3> shifted_p = p + shift;

    // As we are projecting onto a face of a reference cell, the vectors
    // describing its local coordinate system should be orthogonal. We don't
    // know which of the three vectors computed from p are mutually orthogonal
    // for triangles so that case requires an extra check.
    Tensor<1, 3>   e0;
    Tensor<1, 3>   e1;
    const Point<3> vertex = shifted_vertices[0];
    // Triangles are difficult because of two cases:
    // 1. the top face of a Tetrahedron, which does not have a right angle
    // 2. wedges and pyramids, whose faces do not lie on the reference simplex
    //
    // Deal with both by creating a locally orthogonal (but not necessarily
    // orthonormal) coordinate system and testing if the projected point is in
    // the triangle by expressing it as a convex combination of the vertices.
    if (face_reference_cell == ReferenceCells::Triangle)
      {
        e0 = shifted_vertices[1] - shifted_vertices[0];
        e1 = shifted_vertices[2] - shifted_vertices[0];
        e1 -= (e0 * e1) * e0 / (e0.norm_square());
      }
    else
      {
        e0 = shifted_vertices[1] - shifted_vertices[0];
        e1 = shifted_vertices[2] - shifted_vertices[0];
      }
    Assert(std::abs(e0 * e1) <= 1e-14, ExcInternalError());
    // the quadrilaterals on pyramids and wedges don't necessarily have edge
    // lengths of 1 so we cannot skip the denominator
    const double   c0 = e0 * (shifted_p - vertex) / e0.norm_square();
    const double   c1 = e1 * (shifted_p - vertex) / e1.norm_square();
    const Point<3> projected_shifted_p = vertex + c0 * e0 + c1 * e1;

    bool in_quad = false;
    if (face_reference_cell == ReferenceCells::Triangle)
      {
        Tensor<2, 3> shifted_vertex_matrix;
        for (unsigned int i = 0; i < 3; ++i)
          shifted_vertex_matrix[i] = shifted_vertices[i];
        const Tensor<1, 3> combination_coordinates =
          invert(transpose(shifted_vertex_matrix)) * projected_shifted_p;
        bool is_convex_combination = true;
        for (unsigned int i = 0; i < 3; ++i)
          is_convex_combination = is_convex_combination &&
                                  (0.0 <= combination_coordinates[i]) &&
                                  (combination_coordinates[i] <= 1.0);
        in_quad = is_convex_combination;
      }
    else
      in_quad = (0.0 <= c0 && c0 <= 1.0 && 0.0 <= c1 && c1 <= 1.0);

    if (in_quad)
      return std::make_pair(projected_shifted_p - shift,
                            shifted_p.distance_square(projected_shifted_p));
    else
      return std::make_pair(Point<3>(), std::numeric_limits<double>::max());
  }
} // namespace



template <int dim>
Point<dim>
ReferenceCell::closest_point(const Point<dim> &p) const
{
  AssertDimension(dim, get_dimension());

  // Handle simple cases first:
  if (dim == 0)
    return p;
  if (contains_point(p, 0.0))
    return p;
  if (dim == 1)
    return project_to_line(vertex<dim>(0), vertex<dim>(1), p).first;

  // Find the closest vertex so that we only need to check adjacent faces and
  // lines.
  Point<dim>   result;
  unsigned int closest_vertex_no        = 0;
  double closest_vertex_distance_square = vertex<dim>(0).distance_square(p);
  for (unsigned int i = 1; i < n_vertices(); ++i)
    {
      const double new_vertex_distance_square =
        vertex<dim>(i).distance_square(p);
      if (new_vertex_distance_square < closest_vertex_distance_square)
        {
          closest_vertex_no              = i;
          closest_vertex_distance_square = new_vertex_distance_square;
        }
    }

  double min_distance_square = std::numeric_limits<double>::max();
  if (dim == 2)
    {
      for (const unsigned int face_no :
           faces_for_given_vertex(closest_vertex_no))
        {
          const Point<dim> v0 = vertex<dim>(line_to_cell_vertices(face_no, 0));
          const Point<dim> v1 = vertex<dim>(line_to_cell_vertices(face_no, 1));

          auto pair = project_to_line(v0, v1, p);
          if (pair.second < min_distance_square)
            {
              result              = pair.first;
              min_distance_square = pair.second;
            }
        }
    }
  else
    {
      // Check faces and then lines.
      //
      // For reference cells with sloped faces (i.e., all 3D shapes except
      // Hexahedra), we might be able to do a valid normal projection to a face
      // with a different slope which is on the 'other side' of the reference
      // cell. To catch that case we have to unconditionally check lines after
      // checking faces.
      //
      // For pyramids the closest vertex might not be on the closest face: for
      // example, the origin is closest to vertex 4 which is not on the bottom
      // plane. Get around that by just checking all faces for pyramids.
      const std::array<unsigned int, 5> all_pyramid_faces{{0, 1, 2, 3, 4}};
      const auto &faces = *this == ReferenceCells::Pyramid ?
                            ArrayView<const unsigned int>(all_pyramid_faces) :
                            faces_for_given_vertex(closest_vertex_no);
      for (const unsigned int face_no : faces)
        {
          auto face_cell = face_reference_cell(face_no);
          // We only need the first three points since for quads the last point
          // is redundant
          std::array<Point<dim>, 3> vertices;
          for (unsigned int vertex_no = 0; vertex_no < 3; ++vertex_no)
            vertices[vertex_no] = vertex<dim>(face_to_cell_vertices(
              face_no, vertex_no, numbers::default_geometric_orientation));

          auto pair = project_to_quad(vertices, p, face_cell);
          if (pair.second < min_distance_square)
            {
              result              = pair.first;
              min_distance_square = pair.second;
            }
        }

      for (const unsigned int face_no :
           faces_for_given_vertex(closest_vertex_no))
        {
          auto face_cell = face_reference_cell(face_no);
          for (const unsigned int face_line_no : face_cell.line_indices())
            {
              const auto cell_line_no =
                face_to_cell_lines(face_no,
                                   face_line_no,
                                   numbers::default_geometric_orientation);
              const auto v0 =
                vertex<dim>(line_to_cell_vertices(cell_line_no, 0));
              const auto v1 =
                vertex<dim>(line_to_cell_vertices(cell_line_no, 1));
              auto pair = project_to_line(v0, v1, p);
              if (pair.second < min_distance_square)
                {
                  result              = pair.first;
                  min_distance_square = pair.second;
                }
            }
        }
    }

  Assert(min_distance_square < std::numeric_limits<double>::max(),
         ExcInternalError());

  // If necessary, slightly adjust the computed point so that it is closer to
  // being on the surface of the reference cell. Due to roundoff it is difficult
  // to place points on sloped surfaces (e.g., for Pyramids) so this check isn't
  // perfect but does improve the accuracy of the projected point.
  if (!contains_point(result, 0.0))
    {
      constexpr unsigned int x_index = 0;
      constexpr unsigned int y_index = (dim >= 2 ? 1 : 0);
      constexpr unsigned int z_index = (dim >= 3 ? 2 : 0);
      switch (this->kind)
        {
          case ReferenceCells::Vertex:
            DEAL_II_ASSERT_UNREACHABLE();
            break;
            // the bounds for each dimension of a hypercube are mutually
            // independent:
          case ReferenceCells::Line:
          case ReferenceCells::Quadrilateral:
          case ReferenceCells::Hexahedron:
            for (unsigned int d = 0; d < dim; ++d)
              result[d] = std::clamp(result[d], 0.0, 1.0);
            // simplices can use the standard definition of a simplex:
            break;
          case ReferenceCells::Triangle:
            result[x_index] = std::clamp(result[x_index], 0.0, 1.0);
            result[y_index] =
              std::clamp(result[y_index], 0.0, 1.0 - result[x_index]);
            break;
          case ReferenceCells::Tetrahedron:
            result[x_index] = std::clamp(result[x_index], 0.0, 1.0);
            result[y_index] =
              std::clamp(result[y_index], 0.0, 1.0 - result[x_index]);
            result[z_index] =
              std::clamp(result[z_index],
                         0.0,
                         1.0 - result[x_index] - result[y_index]);
            break;
          // wedges and pyramids are more ad-hoc:
          case ReferenceCells::Wedge:
            result[x_index] = std::clamp(result[x_index], 0.0, 1.0);
            result[y_index] =
              std::clamp(result[y_index], 0.0, 1.0 - result[x_index]);
            result[z_index] = std::clamp(result[z_index], 0.0, 1.0);
            break;
          case ReferenceCells::Pyramid:
            {
              result[x_index] = std::clamp(result[x_index], -1.0, 1.0);
              result[y_index] = std::clamp(result[y_index], -1.0, 1.0);
              // It suffices to transform everything to the first quadrant to
              // adjust z:
              const auto x_abs = std::abs(result[x_index]);
              const auto y_abs = std::abs(result[y_index]);

              if (y_abs <= x_abs)
                result[z_index] = std::clamp(result[z_index], 0.0, 1.0 - x_abs);
              else
                result[z_index] = std::clamp(result[z_index], 0.0, 1.0 - y_abs);
            }
            break;
          default:
            DEAL_II_NOT_IMPLEMENTED();
        }
    }
  // We should be within 4 * eps of the cell by this point. The roundoff error
  // comes from, e.g., computing (1 - x) + x when moving points onto the top of
  // a Pyramid.
  Assert(contains_point(result, 4.0 * std::numeric_limits<double>::epsilon()),
         ExcInternalError());

  return result;
}



std::ostream &
operator<<(std::ostream &out, const ReferenceCell &reference_cell)
{
  AssertThrow(out.fail() == false, ExcIO());

  // Output as an integer to avoid outputting it as a character with
  // potentially non-printing value:
  out << static_cast<unsigned int>(reference_cell.kind);
  return out;
}



std::istream &
operator>>(std::istream &in, ReferenceCell &reference_cell)
{
  AssertThrow(in.fail() == false, ExcIO());

  // Read the information as an integer and convert it to the correct type
  unsigned int value;
  in >> value;
  reference_cell.kind = static_cast<decltype(reference_cell.kind)>(value);

  // Ensure that the object we read is valid
  Assert(
    (reference_cell == ReferenceCells::Vertex) ||
      (reference_cell == ReferenceCells::Line) ||
      (reference_cell == ReferenceCells::Triangle) ||
      (reference_cell == ReferenceCells::Quadrilateral) ||
      (reference_cell == ReferenceCells::Tetrahedron) ||
      (reference_cell == ReferenceCells::Hexahedron) ||
      (reference_cell == ReferenceCells::Wedge) ||
      (reference_cell == ReferenceCells::Pyramid) ||
      (reference_cell == ReferenceCells::Invalid),
    ExcMessage(
      "The reference cell kind just read does not correspond to one of the "
      "valid choices. There must be an error."));

  return in;
}

// explicitly instantiate dimension 0 quadrature in addition to the standard
// dimensions
template Quadrature<0>
ReferenceCell::get_gauss_type_quadrature(const unsigned n_points_1D) const;
template const Quadrature<0> &
ReferenceCell::get_nodal_type_quadrature() const;

#include "grid/reference_cell.inst"

DEAL_II_NAMESPACE_CLOSE
