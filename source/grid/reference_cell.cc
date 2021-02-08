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
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include <deal.II/simplex/fe_lib.h>
#include <deal.II/simplex/polynomials.h>
#include <deal.II/simplex/quadrature_lib.h>

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
      Simplex::FE_P<dim, spacedim>(degree));
  else if (*this == ReferenceCells::Pyramid)
    return std::make_unique<MappingFE<dim, spacedim>>(
      Simplex::FE_PyramidP<dim, spacedim>(degree));
  else if (*this == ReferenceCells::Wedge)
    return std::make_unique<MappingFE<dim, spacedim>>(
      Simplex::FE_WedgeP<dim, spacedim>(degree));
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
        Simplex::FE_P<dim, spacedim>(1));
      return mapping;
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      static const MappingFE<dim, spacedim> mapping(
        Simplex::FE_PyramidP<dim, spacedim>(1));
      return mapping;
    }
  else if (*this == ReferenceCells::Wedge)
    {
      static const MappingFE<dim, spacedim> mapping(
        Simplex::FE_WedgeP<dim, spacedim>(1));
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
    return Simplex::QGauss<dim>(n_points_1D);
  else if (*this == ReferenceCells::Pyramid)
    return Simplex::QGaussPyramid<dim>(n_points_1D);
  else if (*this == ReferenceCells::Wedge)
    return Simplex::QGaussWedge<dim>(n_points_1D);
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

#include "reference_cell.inst"

DEAL_II_NAMESPACE_CLOSE
