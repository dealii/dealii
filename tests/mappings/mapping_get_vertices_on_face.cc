// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Checks if vertices on face are given in correct order

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
std::vector<Point<dim>>
get_vertices(const typename Triangulation<dim>::cell_iterator &cell,
             const typename Triangulation<dim>::face_iterator &face,
             const Mapping<dim>                               &mapping)
{
  std::vector<Point<dim>> vertices(face->n_vertices());

  for (unsigned int i = 0; i < face->n_vertices(); ++i)
    {
      vertices[i] = mapping.transform_unit_to_real_cell(
        cell, mapping.transform_real_to_unit_cell(cell, face->vertex(i)));
    }
  return vertices;
}

template <int dim>
void
test_vertex_order()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1.0, 1.0);

  MappingQ<dim> mapping(1);

  for (const auto &cell : tria.active_cell_iterators())
    for (unsigned int f = 0; f < cell->n_faces(); ++f)
      {
        const auto &vertices_detour =
          get_vertices(cell, cell->face(f), mapping);
        const auto &vertices = mapping.get_vertices(cell, f);

        AssertDimension(vertices_detour.size(), vertices.size());

        for (unsigned int i = 0; i < vertices.size(); ++i)
          {
            const auto diff = vertices[i] - vertices_detour[i];
            AssertThrow(diff.norm() < 1e-12, ExcMessage("Vertices differ!"));
          }
      }

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(12);

  test_vertex_order<1>();
  test_vertex_order<2>();
  test_vertex_order<3>();

  return 0;
}
