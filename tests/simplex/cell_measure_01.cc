// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test TriaAccessor::measure(), TriaAccessor::diameter(), and
// TriaAccessor::barycenter().

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim, int spacedim>
void
process(const std::vector<Point<spacedim>> &vertices,
        const std::vector<CellData<dim>>   &cells)
{
  Triangulation<dim, spacedim> tria;
  tria.create_triangulation(vertices, cells, SubCellData());

  std::shared_ptr<MappingFE<dim>> mapping;
  const auto                      reference_cells = tria.get_reference_cells();
  AssertDimension(reference_cells.size(), 1);

  if (reference_cells[0] == ReferenceCells::get_simplex<dim>())
    mapping = std::make_shared<MappingFE<dim>>(FE_SimplexP<dim>(1));
  else if (reference_cells[0] == ReferenceCells::Wedge)
    mapping = std::make_shared<MappingFE<dim>>(FE_WedgeP<dim>(1));
  else if (reference_cells[0] == ReferenceCells::Pyramid)
    mapping = std::make_shared<MappingFE<dim>>(FE_PyramidP<dim>(1));
  else
    AssertThrow(false, ExcNotImplemented());

  const Quadrature<dim> quad =
    reference_cells.front().template get_gauss_type_quadrature<dim>(2);
  FE_Nothing<dim, spacedim> fe(reference_cells.front());
  FEValues<dim, spacedim>   fe_values(*mapping, fe, quad, update_JxW_values);

  deallog << "dim=" << dim << " spacedim=" << spacedim << ':' << std::endl;
  for (const auto &cell : tria.active_cell_iterators())
    {
      fe_values.reinit(cell);
      const double measure = std::accumulate(fe_values.get_JxW_values().begin(),
                                             fe_values.get_JxW_values().end(),
                                             0.0);

      deallog << "measure: " << cell->measure() << std::endl;
      deallog << "measure (via quadrature): " << measure << std::endl;
      deallog << "diameter " << cell->diameter() << std::endl;
      // Only implemented for tensor-product and simplices at the moment
      if (cell->reference_cell() == ReferenceCells::get_simplex<dim>())
        deallog << "barycenter " << cell->barycenter() << std::endl;
    }

  deallog << "diameter_min " << GridTools::minimal_cell_diameter(tria, *mapping)
          << std::endl;
  deallog << "diameter_max " << GridTools::maximal_cell_diameter(tria, *mapping)
          << std::endl;
  deallog << std::endl;
}

template <int dim>
void
test()
{
  DEAL_II_NOT_IMPLEMENTED();
}

template <>
void
test<2>()
{
  const int dim      = 2;
  const int spacedim = 2;

  std::vector<Point<spacedim>> vertices;
  vertices.emplace_back(0, 0);
  vertices.emplace_back(1, 0);
  vertices.emplace_back(0, 1);

  std::vector<CellData<dim>> cells;
  CellData<dim>              cell;
  cell.vertices = {0, 1, 2};
  cells.push_back(cell);

  process(vertices, cells);
}

template <>
void
test<3>()
{
  const int dim      = 3;
  const int spacedim = 3;

  deallog.push("basic tetrahedron");
  {
    std::vector<Point<spacedim>> vertices;
    vertices.emplace_back(0, 0, 0);
    vertices.emplace_back(1, 0, 0);
    vertices.emplace_back(0, 1, 0);
    vertices.emplace_back(0, 0, 1);

    std::vector<CellData<dim>> cells;
    CellData<dim>              cell;
    cell.vertices = {0, 1, 2, 3};
    cells.push_back(cell);

    process(vertices, cells);
  }
  deallog.pop();

  deallog.push("basic wedge");
  {
    std::vector<Point<spacedim>> vertices;
    vertices.emplace_back(0, 0, 0);
    vertices.emplace_back(1, 0, 0);
    vertices.emplace_back(0, 1, 0);
    vertices.emplace_back(0, 0, 1);
    vertices.emplace_back(1, 0, 1);
    vertices.emplace_back(0, 1, 1);

    std::vector<CellData<dim>> cells;
    CellData<dim>              cell;
    cell.vertices = {0, 1, 2, 3, 4, 5};
    cells.push_back(cell);

    process(vertices, cells);
  }
  deallog.pop();

  deallog.push("basic pyramid");
  {
    std::vector<Point<spacedim>> vertices;
    vertices.emplace_back(0, 0, 0);
    vertices.emplace_back(1, 0, 0);
    vertices.emplace_back(0, 1, 0);
    vertices.emplace_back(1, 1, 0);
    vertices.emplace_back(0.5, 0.5, 1.0);

    std::vector<CellData<dim>> cells;
    CellData<dim>              cell;
    cell.vertices = {0, 1, 2, 3, 4};
    cells.push_back(cell);

    process(vertices, cells);
  }
  deallog.pop();

  deallog.push("general tetrahedron");
  {
    std::vector<Point<spacedim>> vertices;
    vertices.emplace_back(0, 0, 0);
    vertices.emplace_back(0.5, 0.25, 0.1);
    vertices.emplace_back(0, 1, 0.2);
    vertices.emplace_back(0.25, 0.5, 0.5);

    std::vector<CellData<dim>> cells;
    CellData<dim>              cell;
    cell.vertices = {0, 1, 2, 3};
    cells.push_back(cell);

    process(vertices, cells);
  }
  deallog.pop();

  deallog.push("general wedge");
  {
    std::vector<Point<spacedim>> vertices;
    vertices.emplace_back(0, 0.2, 0);
    vertices.emplace_back(1, 0, 0.1);
    vertices.emplace_back(0.1, 1.3, 0.2);
    vertices.emplace_back(-1, 0, 1);
    vertices.emplace_back(1, 1, 1.3);
    vertices.emplace_back(0, 2, 1.4);

    std::vector<CellData<dim>> cells;
    CellData<dim>              cell;
    cell.vertices = {0, 1, 2, 3, 4, 5};
    cells.push_back(cell);

    process(vertices, cells);
  }
  deallog.pop();

  deallog.push("general pyramid");
  {
    std::vector<Point<spacedim>> vertices;
    vertices.emplace_back(-2, -1, -1.5);
    vertices.emplace_back(1.2, -0.5, -1);
    vertices.emplace_back(-0.3, 1.2, -0.5);
    vertices.emplace_back(1, 1, 0.0);
    vertices.emplace_back(0.5, 0.5, 1.0);

    std::vector<CellData<dim>> cells;
    CellData<dim>              cell;
    cell.vertices = {0, 1, 2, 3, 4};
    cells.push_back(cell);

    process(vertices, cells);
  }
  deallog.pop();
}

int
main()
{
  initlog();

  test<2>();
  test<3>();
}
