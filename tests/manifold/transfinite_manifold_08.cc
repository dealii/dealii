// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// This test verifies that the transfinite interpolation works on a cylinder
// with oblique faces

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(9);

  const double       radius = 1.;
  const double       length = 3.;
  constexpr int      dim    = 3;
  Triangulation<dim> tria;

  // position of auxiliary point to achieve an angle of 120 degrees in corner
  // of inner cell
  const double ycord =
    0.55 * radius * std::cos(numbers::PI / 12) /
    (std::sin(numbers::PI / 12) + std::cos(numbers::PI / 12));
  // vertices for quarter of circle
  std::vector<Point<3>> vertices{{0, 0, 0},
                                 {0.55 * radius, 0, 0},
                                 {ycord, ycord, 0},
                                 {radius, 0, 0},
                                 {radius * std::sqrt(0.5),
                                  radius * std::sqrt(0.5),
                                  0}};

  // create additional vertices for other three quarters of circle -> gives 17
  // vertices in total
  for (unsigned int a = 1; a < 4; ++a)
    {
      Tensor<2, 3> transform;
      transform[0][0] = a == 2 ? -1. : 0;
      transform[1][0] = a == 2 ? 0 : (a == 1 ? 1 : -1);
      transform[0][1] = -transform[1][0];
      transform[1][1] = transform[0][0];
      transform[2][2] = 1;
      for (unsigned int i = 1; i < 5; ++i)
        vertices.push_back(Point<3>(transform * vertices[i]));
    }
  const unsigned int n_base_vertices = vertices.size();
  for (unsigned int i = 0; i < n_base_vertices; ++i)
    {
      Point<3> new_point;
      for (unsigned int d = 0; d < 2; ++d)
        new_point[d] = vertices[i][d];
      new_point[2] =
        length - std::tan(numbers::PI / 8) * std::abs(vertices[i][0]) / radius;
      vertices.push_back(new_point);
    }

  // create cells for mesh on z=0; the first four elements are at the
  // center of the circle
  std::vector<CellData<3>> cell_data(12);
  cell_data[0].vertices[0] = 0;
  cell_data[0].vertices[1] = 1;
  cell_data[0].vertices[2] = 5;
  cell_data[0].vertices[3] = 2;
  cell_data[1].vertices[0] = 9;
  cell_data[1].vertices[1] = 0;
  cell_data[1].vertices[2] = 6;
  cell_data[1].vertices[3] = 5;
  cell_data[2].vertices[0] = 10;
  cell_data[2].vertices[1] = 13;
  cell_data[2].vertices[2] = 9;
  cell_data[2].vertices[3] = 0;
  cell_data[3].vertices[0] = 13;
  cell_data[3].vertices[1] = 14;
  cell_data[3].vertices[2] = 0;
  cell_data[3].vertices[3] = 1;

  // the next 8 elements describe the rim; we take one quarter of the circle
  // in each loop iteration
  for (unsigned int a = 0; a < 4; ++a)
    {
      cell_data[4 + a * 2].vertices[0] = 1 + a * 4;
      cell_data[4 + a * 2].vertices[1] = 3 + a * 4;
      cell_data[4 + a * 2].vertices[2] = 2 + a * 4;
      cell_data[4 + a * 2].vertices[3] = 4 + a * 4;
      cell_data[5 + a * 2].vertices[0] = 2 + a * 4;
      cell_data[5 + a * 2].vertices[1] = 4 + a * 4;
      AssertIndexRange(4 + a * 4, vertices.size());
      cell_data[5 + a * 2].vertices[2] = a == 3 ? 1 : 5 + a * 4;
      cell_data[5 + a * 2].vertices[3] = a == 3 ? 3 : 7 + a * 4;
    }

  // create extrusion of the base mesh
  for (unsigned int i = 0; i < cell_data.size(); ++i)
    for (unsigned int v = 0; v < 4; ++v)
      cell_data[i].vertices[4 + v] = n_base_vertices + cell_data[i].vertices[v];

  SubCellData subcell_data;
  tria.create_triangulation(vertices, cell_data, subcell_data);

  // set manifolds on circumferential direction (face ids between 0 and 4)
  tria.set_all_manifold_ids(1);
  for (auto &cell : tria.active_cell_iterators())
    for (unsigned int f = 0; f < 4; ++f)
      if (cell->at_boundary(f))
        cell->face(f)->set_all_manifold_ids(2);

  // attach 3 cylindrical manifolds to mesh
  tria.set_manifold(2,
                    CylindricalManifold<dim>(Point<dim>{0, 0, 1},
                                             Point<dim>{0, 0, 0}));
  TransfiniteInterpolationManifold<dim> transfinite;
  transfinite.initialize(tria);
  tria.set_manifold(1, transfinite);

  for (unsigned int l = 0; l < 2; ++l)
    {
      for (auto &cell : tria.active_cell_iterators())
        {
          deallog << "Face centers:" << std::endl;
          for (const unsigned int f : GeometryInfo<dim>::face_indices())
            deallog << cell->face(f)->center(true) << std::endl;
          deallog << std::endl;
        }
      tria.refine_global();
    }
}
