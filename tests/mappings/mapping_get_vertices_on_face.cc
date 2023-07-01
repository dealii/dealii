// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2023 by the deal.II authors
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


// Checks if vertices on face are given in correct oder

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
std::vector<Point<dim>>
get_vertices(const typename Triangulation<dim>::cell_iterator &cell,
             const typename Triangulation<dim>::face_iterator &face,
             const Mapping<dim> &                              mapping)
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
        auto const &vertices_detour =
          get_vertices(cell, cell->face(f), mapping);
        auto const &vertices = mapping.get_vertices(cell, f);

        AssertDimension(vertices_detour.size(), vertices.size());

        for (unsigned int i = 0; i < vertices.size(); ++i)
          {
            auto const diff = vertices[i] - vertices_detour[i];
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
