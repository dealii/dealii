// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test a mesh with two tetrahedra for all possible orientations.

#include <deal.II/base/array_view.h>

#include <deal.II/grid/cell_data.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

void
test_tri(const types::geometric_orientation orientation)
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
      ReferenceCells::Triangle.permute_by_combined_orientation(
        make_array_view(std::array<unsigned int, 3>{{0, 1, 2}}), orientation);

    // add new vertex we need for our new triangle
    vertices.push_back(face->center() + Point<3>(0, 0, 1));

    CellData<3> cell;
    cell.vertices    = {face->vertex_index(permuted[0]),
                        face->vertex_index(permuted[1]),
                        face->vertex_index(permuted[2]),
                        4};
    cell.material_id = 1;
    cells.push_back(cell);
  }

  tria.create_triangulation(vertices, cells, {});

  auto cell = tria.begin();
  cell++;

  // check orientation
  deallog << "face orientation: "
          << int(cell->combined_face_orientation(face_no)) << std::endl;

  // check vertices
  {
    // we do know which vertices should be part of triangle 0 of our cell
    // (0, 1, 2). Lets' compare their indices.
    deallog << "expected indices:   ";
    for (unsigned int v = 0; v < 3; ++v)
      deallog << cell->vertex_index(v) << ' ';

    deallog << std::endl << "calculated indices: ";

    for (unsigned int v = 0; v < 3; ++v)
      {
        const auto orientation = cell->combined_face_orientation(face_no);

        unsigned int corrected_vertex_index =
          ReferenceCell::standard_to_real_quad_vertex(v,
                                                      orientation,
                                                      ReferenceCells::Triangle);

        deallog << cell->face(face_no)->vertex_index(corrected_vertex_index)
                << ' ';
      }
  }

  deallog << std::endl;

  // check lines
  {
    // we do know which lines should be part of triangle 0 of our cell
    // (0, 1, 2). Lets' compare their indices.
    deallog << "expected lines:   ";
    for (unsigned int l = 0; l < 3; ++l)
      deallog << cell->line_index(l) << ' ';

    deallog << std::endl << "calculated lines: ";

    for (unsigned int v = 0; v < 3; ++v)
      {
        const auto orientation = cell->combined_face_orientation(face_no);

        unsigned int corrected_line_index =
          ReferenceCell::standard_to_real_quad_line(v,
                                                    orientation,
                                                    ReferenceCells::Triangle);

        deallog << cell->face(face_no)->line_index(corrected_line_index) << ' ';
      }
  }

  deallog << std::endl << std::endl;
}

void
test_quad(const types::geometric_orientation orientation)
{
  const unsigned int face_no = 0;

  Triangulation<3> dummy, tria;
  GridGenerator::reference_cell(dummy, ReferenceCells::Pyramid);

  // for (auto f : dummy.begin()->face_indices()) {
  //   deallog << f << ": ";
  //   for (auto v : dummy.begin()->face(f)->vertex_indices())
  //     deallog << dummy.begin()->face(f)->vertex_index(v) << " ";
  //   deallog << std::endl;
  // }

  auto vertices = dummy.get_vertices();

  std::vector<CellData<3>> cells;

  {
    CellData<3> cell;
    cell.vertices    = {0, 1, 2, 3, 4};
    cell.material_id = 0;
    cells.push_back(cell);
  }

  {
    const auto &face = dummy.begin()->face(face_no);
    const auto  permuted =
      ReferenceCells::Quadrilateral.permute_by_combined_orientation(
        make_array_view(std::array<unsigned int, 4>{{0, 1, 2, 3}}),
        orientation);

    CellData<3> cell;
    cell.vertices    = {face->vertex_index(permuted[0]),
                        face->vertex_index(permuted[1]),
                        face->vertex_index(permuted[2]),
                        face->vertex_index(permuted[3]),
                        4};
    cell.material_id = 1;
    cells.push_back(cell);
  }

  tria.create_triangulation(vertices, cells, {});

  auto cell = tria.begin();
  cell++;

  // check orientation
  deallog << "face orientation: "
          << int(cell->combined_face_orientation(face_no)) << std::endl;

  // check vertices
  {
    // we do know which vertices should be part of triangle 0 of our cell
    // (0, 1, 2). Lets' compare their indices.
    deallog << "expected indices:   ";
    for (unsigned int v = 0; v < 4; ++v)
      deallog << cell->vertex_index(v) << ' ';

    deallog << std::endl << "calculated indices: ";

    for (unsigned int v = 0; v < 4; ++v)
      {
        const auto orientation = cell->combined_face_orientation(face_no);

        unsigned int corrected_vertex_index =
          ReferenceCell::standard_to_real_quad_vertex(
            v, orientation, ReferenceCells::Quadrilateral);

        deallog << cell->face(face_no)->vertex_index(corrected_vertex_index)
                << ' ';
      }
  }

  deallog << std::endl;

  // check lines
  {
    // we do know which lines should be part of triangle 0 of our cell
    // (0, 1, 2). Lets' compare their indices.
    deallog << "expected lines:   ";
    for (unsigned int l = 0; l < 4; ++l)
      deallog << cell->line_index(l) << ' ';

    deallog << std::endl << "calculated lines: ";

    for (unsigned int v = 0; v < 4; ++v)
      {
        const auto orientation = cell->combined_face_orientation(face_no);

        unsigned int corrected_line_index =
          ReferenceCell::standard_to_real_quad_line(
            v, orientation, ReferenceCells::Quadrilateral);

        deallog << cell->face(face_no)->line_index(corrected_line_index) << ' ';
      }
  }

  deallog << std::endl << std::endl;
}

int
main()
{
  initlog();

  deallog << "--- Triangle ---" << std::endl;
  for (unsigned int o = 0; o < 6; ++o)
    test_tri(o);

  deallog << "--- Quad ---" << std::endl;
  for (unsigned int o = 0; o < 8; ++o)
    test_quad(o);
}
