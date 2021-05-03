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



// Test TriaAccessor::measure(), TriaAccessor::diameter(), and
// TriaAccessor::barycenter().

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim, int spacedim>
void
process(const std::vector<Point<spacedim>> &vertices,
        const std::vector<CellData<dim>> &  cells)
{
  Triangulation<dim, spacedim> tria;
  tria.create_triangulation(vertices, cells, SubCellData());

  deallog << "dim=" << dim << " spacedim=" << spacedim << ":" << std::endl;
  for (const auto &cell : tria.active_cell_iterators())
    {
      deallog << "measure:  " << cell->measure() << std::endl;
      deallog << "diameter: " << cell->diameter() << std::endl;
      // Only implemented for tensor-product and simplices at the moment
      if (cell->reference_cell() == ReferenceCells::get_simplex<dim>())
        deallog << "barycenter: " << cell->barycenter() << std::endl;
    }

  std::shared_ptr<MappingFE<dim>> mapping;

  const auto reference_cells = tria.get_reference_cells();

  AssertDimension(reference_cells.size(), 1);

  if (reference_cells[0] == ReferenceCells::get_simplex<dim>())
    mapping = std::make_shared<MappingFE<dim>>(FE_SimplexP<dim>(1));
  else if (reference_cells[0] == ReferenceCells::Wedge)
    mapping = std::make_shared<MappingFE<dim>>(FE_WedgeP<dim>(1));
  else
    AssertThrow(false, ExcNotImplemented());

  deallog << "diameter_min: "
          << GridTools::minimal_cell_diameter(tria, *mapping) << std::endl;
  deallog << "diameter_max: "
          << GridTools::maximal_cell_diameter(tria, *mapping) << std::endl;
  deallog << std::endl;
}

template <int dim>
void
test()
{
  Assert(false, ExcNotImplemented());
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
}

int
main()
{
  initlog();

  test<2>();
  test<3>();
}
