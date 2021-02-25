// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN


std::string
ReferenceCell::to_string() const
{
  if (*this == ReferenceCells::Vertex)
    return "Vertex";
  else if (*this == ReferenceCells::Line)
    return "Line";
  else if (*this == ReferenceCells::Triangle)
    return "Tri";
  else if (*this == ReferenceCells::Quadrilateral)
    return "Quad";
  else if (*this == ReferenceCells::Tetrahedron)
    return "Tet";
  else if (*this == ReferenceCells::Pyramid)
    return "Pyramid";
  else if (*this == ReferenceCells::Wedge)
    return "Wedge";
  else if (*this == ReferenceCells::Hexahedron)
    return "Hex";
  else if (*this == ReferenceCells::Invalid)
    return "Invalid";

  Assert(false, ExcNotImplemented());

  return "Invalid";
}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
ReferenceCell::get_default_mapping(const unsigned int degree) const
{
  AssertDimension(dim, get_dimension());

  if (is_hyper_cube())
    return std::make_unique<MappingQGeneric<dim, spacedim>>(degree);
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
      Assert(false, ExcNotImplemented());
    }

  return std::make_unique<MappingQGeneric<dim, spacedim>>(degree);
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
      static const MappingFE<dim, spacedim> mapping(
        FE_SimplexP<dim, spacedim>(1));
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
      Assert(false, ExcNotImplemented());
    }

  return StaticMappingQ1<dim, spacedim>::mapping; // never reached
}



template <int dim>
Quadrature<dim>
ReferenceCell::get_gauss_type_quadrature(const unsigned n_points_1D) const
{
  AssertDimension(dim, get_dimension());

  if (is_hyper_cube())
    return QGauss<dim>(n_points_1D);
  else if (is_simplex())
    return QGaussSimplex<dim>(n_points_1D);
  else if (*this == ReferenceCells::Pyramid)
    return QGaussPyramid<dim>(n_points_1D);
  else if (*this == ReferenceCells::Wedge)
    return QGaussWedge<dim>(n_points_1D);
  else
    Assert(false, ExcNotImplemented());

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
    Triangulation<dim> tria;
    GridGenerator::reference_cell(reference_cell, tria);

    return Quadrature<dim>(tria.get_vertices());
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
    Assert(false, ExcNotImplemented());

  static const Quadrature<dim> dummy;
  return dummy; // never reached
}



unsigned int
ReferenceCell::exodusii_vertex_to_deal_vertex(const unsigned int vertex_n) const
{
  AssertIndexRange(vertex_n, n_vertices());

  if (*this == ReferenceCells::Line)
    {
      return vertex_n;
    }
  else if (*this == ReferenceCells::Triangle)
    {
      return vertex_n;
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      constexpr std::array<unsigned int, 4> exodus_to_deal{{0, 1, 3, 2}};
      return exodus_to_deal[vertex_n];
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      return vertex_n;
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      constexpr std::array<unsigned int, 8> exodus_to_deal{
        {0, 1, 3, 2, 4, 5, 7, 6}};
      return exodus_to_deal[vertex_n];
    }
  else if (*this == ReferenceCells::Wedge)
    {
      constexpr std::array<unsigned int, 6> exodus_to_deal{{2, 1, 0, 5, 4, 3}};
      return exodus_to_deal[vertex_n];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      constexpr std::array<unsigned int, 5> exodus_to_deal{{0, 1, 3, 2, 4}};
      return exodus_to_deal[vertex_n];
    }

  Assert(false, ExcNotImplemented());

  return numbers::invalid_unsigned_int;
}



unsigned int
ReferenceCell::exodusii_face_to_deal_face(const unsigned int face_n) const
{
  AssertIndexRange(face_n, n_faces());

  if (*this == ReferenceCells::Vertex)
    {
      return 0;
    }
  if (*this == ReferenceCells::Line)
    {
      return face_n;
    }
  else if (*this == ReferenceCells::Triangle)
    {
      return face_n;
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      constexpr std::array<unsigned int, 4> exodus_to_deal{{2, 1, 3, 0}};
      return exodus_to_deal[face_n];
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      constexpr std::array<unsigned int, 4> exodus_to_deal{{1, 3, 2, 0}};
      return exodus_to_deal[face_n];
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      constexpr std::array<unsigned int, 6> exodus_to_deal{{2, 1, 3, 0, 4, 5}};
      return exodus_to_deal[face_n];
    }
  else if (*this == ReferenceCells::Wedge)
    {
      constexpr std::array<unsigned int, 6> exodus_to_deal{{3, 4, 2, 0, 1}};
      return exodus_to_deal[face_n];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      constexpr std::array<unsigned int, 5> exodus_to_deal{{3, 2, 4, 1, 0}};
      return exodus_to_deal[face_n];
    }

  Assert(false, ExcNotImplemented());

  return numbers::invalid_unsigned_int;
}



#include "reference_cell.inst"

DEAL_II_NAMESPACE_CLOSE
