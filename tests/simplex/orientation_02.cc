// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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



// Test a mesh with two tetrahedra for all possible orientations.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

void
test(const unsigned int orientation)
{
  const unsigned int face_no = 0;


  Triangulation<3> dummy, tria;
  GridGenerator::reference_cell(dummy, ReferenceCells::Tetrahedron);

  auto vertices = dummy.get_vertices();

  std::vector<CellData<3>> cells;

  {
    CellData<3> cell;
    cell.vertices    = {0, 1, 2, 3};
    cell.material_id = 0;
    cells.push_back(cell);
  }

  {
    const auto &face = dummy.begin()->face(face_no);
    const auto  permuted =
      ReferenceCell(ReferenceCells::Triangle)
        .permute_according_orientation(
          std::array<unsigned int, 3>{{face->vertex_index(0),
                                       face->vertex_index(1),
                                       face->vertex_index(2)}},
          orientation);

    auto direction =
      cross_product_3d(vertices[permuted[1]] - vertices[permuted[0]],
                       vertices[permuted[2]] - vertices[permuted[0]]);
    direction = direction / direction.norm();

    vertices.push_back(face->center() + direction);

    CellData<3> cell;
    cell.vertices = {permuted[0], permuted[1], permuted[2], 4u};

    cell.material_id = 1;
    cells.push_back(cell);
  }

  tria.create_triangulation(vertices, cells, {});

  auto cell = tria.begin();
  cell++;

  // check orientation
  deallog << "face orientation: " << orientation << " " << std::endl;
  AssertDimension(orientation,
                  (cell->face_orientation(0) * 1 + cell->face_rotation(0) * 2 +
                   cell->face_flip(0) * 4));

  // check vertices
  {
    for (unsigned int v = 0; v < cell->n_vertices(); ++v)
      deallog << cell->vertex_index(v) << " ";
    deallog << " vs. ";
    for (unsigned int v = 0; v < cell->n_vertices(); ++v)
      deallog << cells[1].vertices[v] << " ";
    deallog << std::endl;

    for (const auto v : cell->vertex_indices())
      AssertDimension(cell->vertex_index(v), cells[1].vertices[v]);
  }

  const auto face = cell->face(0);

  // check lines
  for (const auto l : face->line_indices())
    {
      const unsigned int l_ =
        ReferenceCells::Tetrahedron.standard_to_real_face_line(l,
                                                               face_no,
                                                               orientation);

      std::array<unsigned int, 2> a = {
        {face->line(l_)->vertex_index(0), face->line(l_)->vertex_index(1)}};
      std::sort(a.begin(), a.end());

      std::array<unsigned int, 2> b = {
        {cells[1].vertices[l],
         cells[1].vertices[(l + 1) % face->n_vertices()]}};
      std::sort(b.begin(), b.end());

      deallog << a[0] << "-" << a[1] << " vs. " << b[0] << "-" << b[1]
              << std::endl;

      AssertDimension(a[0], b[0]);
      AssertDimension(a[1], b[1]);
    }

  deallog << std::endl;
}

int
main()
{
  initlog();

  for (unsigned int o = 0; o < 6; ++o)
    test(o);
}
