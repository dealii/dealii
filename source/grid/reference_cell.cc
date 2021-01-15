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

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include <deal.II/simplex/fe_lib.h>
#include <deal.II/simplex/polynomials.h>
#include <deal.II/simplex/quadrature_lib.h>

DEAL_II_NAMESPACE_OPEN

namespace ReferenceCell
{
  template <int dim, int spacedim>
  void
  make_triangulation(const Type &                  reference_cell,
                     Triangulation<dim, spacedim> &tria)
  {
    AssertDimension(dim, get_dimension(reference_cell));

    if (reference_cell == get_hypercube(dim))
      {
        GridGenerator::hyper_cube(tria, 0, 1);
      }
    else if (reference_cell == Type::Tri)
      {
        static const std::vector<Point<spacedim>> vertices = {
          {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}};

        std::vector<CellData<dim>> cells(1);
        cells[0].vertices = {0, 1, 2};

        tria.create_triangulation(vertices, cells, {});
      }
    else if (reference_cell == Type::Tet)
      {
        AssertDimension(spacedim, 3);

        static const std::vector<Point<spacedim>> vertices = {
          {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};

        std::vector<CellData<dim>> cells(1);
        cells[0].vertices = {0, 1, 2, 3};

        tria.create_triangulation(vertices, cells, {});
      }
    else if (reference_cell == Type::Pyramid)
      {
        AssertDimension(spacedim, 3);

        static const std::vector<Point<spacedim>> vertices = {
          {{-1.0, -1.0, 0.0},
           {+1.0, -1.0, 0.0},
           {-1.0, +1.0, 0.0},
           {+1.0, +1.0, 0.0},
           {+0.0, +0.0, 1.0}}};

        std::vector<CellData<dim>> cells(1);
        cells[0].vertices = {0, 1, 2, 3, 4};

        tria.create_triangulation(vertices, cells, {});
      }
    else if (reference_cell == Type::Wedge)
      {
        AssertDimension(spacedim, 3);

        static const std::vector<Point<spacedim>> vertices = {
          {{1.0, 0.0, 0.0},
           {0.0, 1.0, 0.0},
           {0.0, 0.0, 0.0},
           {1.0, 0.0, 1.0},
           {0.0, 1.0, 1.0},
           {0.0, 0.0, 1.0}}};

        std::vector<CellData<dim>> cells(1);
        cells[0].vertices = {0, 1, 2, 3, 4, 5};

        tria.create_triangulation(vertices, cells, {});
      }
    else
      {
        Assert(false, ExcNotImplemented());
      }
  }



  template <int dim, int spacedim>
  const Mapping<dim, spacedim> &
  get_default_linear_mapping(const Type &reference_cell)
  {
    AssertDimension(dim, get_dimension(reference_cell));

    if (reference_cell == get_hypercube(dim))
      {
        return StaticMappingQ1<dim, spacedim>::mapping;
      }
    else if (reference_cell == Type::Tri || reference_cell == Type::Tet)
      {
        static const MappingFE<dim, spacedim> mapping(
          Simplex::FE_P<dim, spacedim>(1));
        return mapping;
      }
    else if (reference_cell == Type::Pyramid)
      {
        static const MappingFE<dim, spacedim> mapping(
          Simplex::FE_PyramidP<dim, spacedim>(1));
        return mapping;
      }
    else if (reference_cell == Type::Wedge)
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



  template <int dim, int spacedim>
  const Mapping<dim, spacedim> &
  get_default_linear_mapping(const Triangulation<dim, spacedim> &triangulation)
  {
    const auto &reference_cell_types = triangulation.get_reference_cell_types();
    Assert(reference_cell_types.size() == 1,
           ExcMessage(
             "This function can only work for triangulations that "
             "use only a single cell type -- for example, only triangles "
             "or only quadrilaterals. For mixed meshes, there is no "
             "single linear mapping object that can be used for all "
             "cells of the triangulation. The triangulation you are "
             "passing to this function uses multiple cell types."));

    return get_default_linear_mapping<dim, spacedim>(
      reference_cell_types.front());
  }



  template <int dim>
  Quadrature<dim>
  get_gauss_type_quadrature(const Type &   reference_cell,
                            const unsigned n_points_1D)
  {
    AssertDimension(dim, get_dimension(reference_cell));

    if (reference_cell == get_hypercube(dim))
      return QGauss<dim>(n_points_1D);
    else if (reference_cell == Type::Tri || reference_cell == Type::Tet)
      return Simplex::QGauss<dim>(n_points_1D);
    else if (reference_cell == Type::Pyramid)
      return Simplex::QGaussPyramid<dim>(n_points_1D);
    else if (reference_cell == Type::Wedge)
      return Simplex::QGaussWedge<dim>(n_points_1D);
    else
      Assert(false, ExcNotImplemented());

    return Quadrature<dim>(); // never reached
  }

  template <int dim>
  Quadrature<dim> &
  get_nodal_type_quadrature(const Type &reference_cell)
  {
    AssertDimension(dim, get_dimension(reference_cell));

    const auto create_quadrature = [](const Type &reference_cell) {
      Triangulation<dim> tria;
      make_triangulation(reference_cell, tria);

      return Quadrature<dim>(tria.get_vertices());
    };

    if (reference_cell == get_hypercube(dim))
      {
        static Quadrature<dim> quadrature = create_quadrature(reference_cell);
        return quadrature;
      }
    else if (reference_cell == Type::Tri || reference_cell == Type::Tet)
      {
        static Quadrature<dim> quadrature = create_quadrature(reference_cell);
        return quadrature;
      }
    else if (reference_cell == Type::Pyramid)
      {
        static Quadrature<dim> quadrature = create_quadrature(reference_cell);
        return quadrature;
      }
    else if (reference_cell == Type::Wedge)
      {
        static Quadrature<dim> quadrature = create_quadrature(reference_cell);
        return quadrature;
      }
    else
      Assert(false, ExcNotImplemented());

    static Quadrature<dim> dummy;

    return dummy; // never reached
  }

#include "reference_cell.inst"

} // namespace ReferenceCell

DEAL_II_NAMESPACE_CLOSE
