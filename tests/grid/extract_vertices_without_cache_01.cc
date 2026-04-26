// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/grid/cell_data.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_topology.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

// Verify that vertices set by the new extract_vertices_without_cache() function
// are identical to existing vertex indices.

// Also verify against the former algorithm we used to extract vertices.
template <int dim, int spacedim>
void
extract_vertices(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  ArrayView<unsigned int>                                     vertex_indices)
{
  const auto reference_cell = cell->reference_cell();
  for (const unsigned int i : cell->vertex_indices())
    {
      const auto [face_index, vertex_index] =
        reference_cell.standard_vertex_to_face_and_vertex_index(i);
      // Since the last lookup is on the face (not the cell!) we
      // need to extract this index relative to the orientation (not
      // reverse orientation)
      const auto vertex_within_face_index =
        reference_cell.standard_to_real_face_vertex(
          vertex_index,
          face_index,
          cell->combined_face_orientation(face_index));
      vertex_indices[i] =
        cell->face(face_index)->vertex_index(vertex_within_face_index);
    }
}

template <int dim, int spacedim>
void
test(Triangulation<dim, spacedim> &tria)
{
  for (unsigned int i = 0; i < 5 - dim; ++i)
    {
      tria.begin_active()->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  std::vector<unsigned int> cell_vertices, cell_vertices_2;
  for (const auto &cell : tria.cell_iterators())
    {
      cell_vertices.resize(cell->n_vertices());
      cell_vertices_2.resize(cell->n_vertices());
      GridTools::internal::extract_vertices_without_cache<dim, spacedim>(
        cell, make_array_view(cell_vertices));
      extract_vertices<dim, spacedim>(cell, make_array_view(cell_vertices_2));

      deallog << cell << ": " << cell->vertex_index(0);
      for (unsigned int vertex_no = 1; vertex_no < cell_vertices.size();
           ++vertex_no)
        {
          deallog << ", " << cell->vertex_index(vertex_no);
          AssertThrow(cell->vertex_index(vertex_no) == cell_vertices[vertex_no],
                      ExcInternalError());
          AssertThrow(cell_vertices_2[vertex_no] == cell_vertices[vertex_no],
                      ExcInternalError());
        }
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  {
    Triangulation<1> tria;
    GridGenerator::reference_cell(tria, ReferenceCells::Line);
    deallog << "Line:" << std::endl;
    test(tria);
  }

  {
    Triangulation<2> tria;
    GridGenerator::reference_cell(tria, ReferenceCells::Triangle);
    deallog << "Triangle:" << std::endl;
    test(tria);
  }

  {
    Triangulation<2> tria;
    GridGenerator::reference_cell(tria, ReferenceCells::Quadrilateral);
    deallog << "Quadrilateral:" << std::endl;
    test(tria);
  }

  // Test tets with full permutations:
  {
    Triangulation<3>      tria;
    std::vector<Point<3>> vertices{{0.0, 0.0, 0.0},
                                   {1.0, 0.0, 0.0},
                                   {0.0, 1.0, 0.0},
                                   {1.0, 1.0, 0.0},
                                   {0.0, 0.0, 1.0},
                                   {1.0, 1.0, 1.0}};

    std::vector<CellData<3>> cells(2);
    cells[0].vertices = {0u, 1u, 2u, 3u};
    cells[1].vertices = {1u, 2u, 3u, 4u};

    unsigned int permutation_no = 0;
    do
      {
        tria.clear();
        tria.create_triangulation(vertices, cells, SubCellData());
        deallog << "Tetrahedron 2 permutation " << permutation_no << ":"
                << std::endl;
        test(tria);
        ++permutation_no;
      }
    while (std::next_permutation(cells.back().vertices.begin(),
                                 cells.back().vertices.end()));
  }

  {
    Triangulation<3> tria;
    GridGenerator::reference_cell(tria, ReferenceCells::Wedge);
    deallog << "Wedge:" << std::endl;
    test(tria);
  }

  // Test hexes with rotations:
  {
    Triangulation<3>      tria;
    std::vector<Point<3>> vertices{{0.0, 0.0, 0.0},
                                   {1.0, 0.0, 0.0},
                                   {0.0, 1.0, 0.0},
                                   {1.0, 1.0, 0.0},
                                   {0.0, 0.0, 1.0},
                                   {1.0, 0.0, 1.0},
                                   {0.0, 1.0, 1.0},
                                   {1.0, 1.0, 1.0},
                                   {0.0, 0.0, 2.0},
                                   {1.0, 0.0, 2.0},
                                   {0.0, 1.0, 2.0},
                                   {1.0, 1.0, 2.0}};

    // Some of these hexes don't make sense but we are just checking
    // permutations
    for (types::geometric_orientation rotation_no = 0; rotation_no < 8;
         ++rotation_no)
      {
        std::vector<CellData<3>> cells(2);
        std::iota(cells[0].vertices.begin(), cells[0].vertices.end(), 0u);
        std::iota(cells[1].vertices.begin(), cells[1].vertices.end(), 4u);
        const auto new_lower_face =
          ReferenceCells::Quadrilateral.permute_by_combined_orientation(
            ArrayView<const unsigned int>(cells[1].vertices.data(), 4),
            rotation_no);
        std::copy(new_lower_face.begin(),
                  new_lower_face.end(),
                  cells[1].vertices.begin());
        const auto new_upper_face =
          ReferenceCells::Quadrilateral.permute_by_combined_orientation(
            ArrayView<const unsigned int>(cells[1].vertices.data() + 4, 4),
            rotation_no);
        std::copy(new_upper_face.begin(),
                  new_upper_face.end(),
                  cells[1].vertices.begin() + 4);

        tria.clear();
        tria.create_triangulation(vertices, cells, SubCellData());
        deallog << "Hexahedron 2 rotation " << int(rotation_no) << ":"
                << std::endl;
        test(tria);
      }
  }
}
