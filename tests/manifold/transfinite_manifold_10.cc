// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// This test verifies that the transfinite interpolation works on a torus

#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <memory>

#include "../tests.h"


class InnerTorusManifold : public ChartManifold<3, 3, 3>
{
public:
  static const unsigned int dim      = 3;
  static const unsigned int spacedim = 3;
  static const unsigned int chartdim = 3;

  InnerTorusManifold()
    : ChartManifold<dim, spacedim>(Tensor<1, 3>({0., 0., 2 * numbers::PI}))
  {}

  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override
  {
    return std::make_unique<InnerTorusManifold>();
  }

  virtual Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const override
  {
    const double phi    = std::atan2(space_point[2], space_point[0]);
    const double radius = std::sqrt(space_point[0] * space_point[0] +
                                    space_point[2] * space_point[2]);
    return Point<3>(radius, space_point[1], phi);
  }

  virtual Point<spacedim>
  push_forward(const Point<chartdim> &chart_point) const override
  {
    return Point<3>(chart_point[0] * std::cos(chart_point[2]),
                    chart_point[1],
                    chart_point[0] * std::sin(chart_point[2]));
  }
};



int
main()
{
  initlog();
  deallog << std::setprecision(9);

  const unsigned int dim = 3;
  Triangulation<dim> tria;

  const double       R       = 2;
  const double       r       = 0.5;
  const unsigned int n_cells = 3;

  // the first 8 vertices are in the x-y-plane
  Point<3> const        p = Point<3>(R, 0.0, 0.0);
  double const          a = 1. / (1 + std::sqrt(2.0));
  std::vector<Point<3>> vertices(8 * n_cells);
  vertices[0] = p + Point<3>(-1, -1, 0) * (r / std::sqrt(2.0)),
  vertices[1] = p + Point<3>(+1, -1, 0) * (r / std::sqrt(2.0)),
  vertices[2] = p + Point<3>(-1, -1, 0) * (r / std::sqrt(2.0) * a),
  vertices[3] = p + Point<3>(+1, -1, 0) * (r / std::sqrt(2.0) * a),
  vertices[4] = p + Point<3>(-1, +1, 0) * (r / std::sqrt(2.0) * a),
  vertices[5] = p + Point<3>(+1, +1, 0) * (r / std::sqrt(2.0) * a),
  vertices[6] = p + Point<3>(-1, +1, 0) * (r / std::sqrt(2.0)),
  vertices[7] = p + Point<3>(+1, +1, 0) * (r / std::sqrt(2.0));

  // create remaining vertices by rotating around positive y-axis
  double phi_cell = 2.0 * numbers::PI / n_cells;
  for (unsigned int c = 1; c < n_cells; ++c)
    {
      for (unsigned int v = 0; v < 8; ++v)
        {
          double const r_2d      = vertices[v][0];
          vertices[8 * c + v][0] = r_2d * std::cos(phi_cell * c);
          vertices[8 * c + v][1] = vertices[v][1];
          vertices[8 * c + v][2] = r_2d * std::sin(phi_cell * c);
        }
    }

  // cell connectivity
  std::vector<CellData<3>> cells(5 * n_cells);
  for (unsigned int c = 0; c < n_cells; ++c)
    {
      for (unsigned int j = 0; j < 2; ++j)
        {
          unsigned int offset = (8 * (c + j)) % (8 * n_cells);
          // cell 0 in x-y-plane
          cells[5 * c].vertices[0 + j * 4] = offset + 0;
          cells[5 * c].vertices[1 + j * 4] = offset + 1;
          cells[5 * c].vertices[2 + j * 4] = offset + 2;
          cells[5 * c].vertices[3 + j * 4] = offset + 3;
          // cell 1 in x-y-plane
          cells[5 * c + 1].vertices[0 + j * 4] = offset + 2;
          cells[5 * c + 1].vertices[1 + j * 4] = offset + 3;
          cells[5 * c + 1].vertices[2 + j * 4] = offset + 4;
          cells[5 * c + 1].vertices[3 + j * 4] = offset + 5;
          // cell 2 in x-y-plane
          cells[5 * c + 2].vertices[0 + j * 4] = offset + 4;
          cells[5 * c + 2].vertices[1 + j * 4] = offset + 5;
          cells[5 * c + 2].vertices[2 + j * 4] = offset + 6;
          cells[5 * c + 2].vertices[3 + j * 4] = offset + 7;
          // cell 3 in x-y-plane
          cells[5 * c + 3].vertices[0 + j * 4] = offset + 0;
          cells[5 * c + 3].vertices[1 + j * 4] = offset + 2;
          cells[5 * c + 3].vertices[2 + j * 4] = offset + 6;
          cells[5 * c + 3].vertices[3 + j * 4] = offset + 4;
          // cell 4 in x-y-plane
          cells[5 * c + 4].vertices[0 + j * 4] = offset + 3;
          cells[5 * c + 4].vertices[1 + j * 4] = offset + 1;
          cells[5 * c + 4].vertices[2 + j * 4] = offset + 5;
          cells[5 * c + 4].vertices[3 + j * 4] = offset + 7;
        }

      cells[5 * c].material_id     = 0;
      cells[5 * c + 1].material_id = 0;
      cells[5 * c + 2].material_id = 0;
      cells[5 * c + 3].material_id = 0;
      cells[5 * c + 4].material_id = 0;
    }

  tria.create_triangulation(vertices, cells, SubCellData());

  tria.set_all_manifold_ids(0);

  for (auto &cell : tria.cell_iterators())
    {
      bool cell_at_boundary = false;
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        if (cell->at_boundary(f))
          cell_at_boundary = true;
      if (cell_at_boundary == false)
        cell->set_all_manifold_ids(2);
    }
  tria.set_all_manifold_ids_on_boundary(1);
  tria.set_manifold(1, TorusManifold<dim>(2, 0.5));
  tria.set_manifold(2, InnerTorusManifold());
  TransfiniteInterpolationManifold<dim> transfinite;
  transfinite.initialize(tria);
  tria.set_manifold(0, transfinite);

  for (unsigned int cycle = 0; cycle < 2; ++cycle)
    {
      tria.refine_global(1);
      deallog << "Cycle " << cycle << std::endl;
      deallog << "Number of active cells: " << tria.n_active_cells()
              << std::endl;
      deallog << "Vertices: " << std::endl;
      for (auto &v : tria.get_vertices())
        deallog << v << std::endl;
      deallog << std::endl;
    }
  return 0;
}
