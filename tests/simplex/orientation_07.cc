/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2023 - 2026 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Test a mesh with two pyramids for all possible orientations. Similar to
// orientation_04

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

void
create_triangulation(const std::vector<Point<3>>    &vertices_,
                     const std::vector<CellData<3>> &cell_data_,
                     const unsigned int              face_n,
                     const unsigned int              final_orientation,
                     Triangulation<3>               &tria)
{
  const ReferenceCell ref_cell  = ReferenceCells::Pyramid;
  auto                cell_data = cell_data_;
  cell_data.emplace_back();

  if (face_n == 0)
    {
      auto vertices = vertices_;
      vertices.push_back(Point<3>(0., 0., -1.0));
      cell_data.back().vertices.resize(0);

      std::vector<unsigned int>     vertices_index = {{0, 1, 2, 3}};
      ArrayView<const unsigned int> const_view(&vertices_index[0], 4);

      const auto permuted_indices =
        ReferenceCells::Quadrilateral.permute_by_combined_orientation(
          const_view, final_orientation);

      for (unsigned int i = 0; i < vertices_index.size(); ++i)
        cell_data.back().vertices.push_back(permuted_indices[i]);

      cell_data.back().vertices.push_back(ref_cell.n_vertices());

      tria.clear();
      tria.create_triangulation(vertices, cell_data, SubCellData());
    }
  else
    {
      bool found_config = false;

      for (unsigned int face_index = 1; face_index < 5 && found_config == false;
           ++face_index)
        for (unsigned int o = 0; o < 6 && found_config == false; ++o)
          {
            auto vertices = vertices_;

            vertices.push_back(
              ref_cell.vertex<3>(
                ref_cell.face_to_cell_vertices(face_index, 0, o)) +
              ref_cell.template face_normal_vector<3>(face_index));
            vertices.push_back(
              ref_cell.vertex<3>(
                ref_cell.face_to_cell_vertices(face_index, 1, o)) +
              ref_cell.template face_normal_vector<3>(face_index));

            std::vector<unsigned int> vertices_index;
            vertices_index.emplace_back(
              ref_cell.face_to_cell_vertices(face_index, 0, o));
            vertices_index.emplace_back(
              ref_cell.face_to_cell_vertices(face_index, 1, o));
            vertices_index.emplace_back(ref_cell.n_vertices());
            vertices_index.emplace_back(ref_cell.n_vertices() + 1);

            ArrayView<const unsigned int> const_vertices_index(
              &vertices_index[0], 4);


            for (unsigned int counter = 0; counter < 8 && found_config == false;
                 ++counter)
              {
                cell_data.back().vertices.resize(0);

                const auto permuted_indices =
                  ReferenceCells::Quadrilateral.permute_by_combined_orientation(
                    const_vertices_index, counter);

                for (unsigned int i = 0; i < vertices_index.size(); ++i)
                  cell_data.back().vertices.push_back(permuted_indices[i]);

                cell_data.back().vertices.push_back(
                  ref_cell.face_to_cell_vertices(face_index, 2, o));

                tria.clear();
                tria.create_triangulation(vertices, cell_data, SubCellData());

                const auto cell = tria.begin();
                const auto face = cell->face(face_index);

                auto ncell = tria.begin();
                ncell++;

                unsigned int nf;
                for (unsigned int i = 0; i < ref_cell.n_faces(); ++i)
                  if (ncell->face(i) == face)
                    nf = i;

                const auto current_orientation =
                  ncell->combined_face_orientation(nf);

                if (current_orientation == final_orientation && nf == face_n)
                  {
                    found_config = true;
                  }
              }
            ++o;
          }
    }
}



void
test()
{
  const unsigned int dim = 3;

  for (unsigned int f = 0; f < 5; ++f)
    {
      const unsigned int n_orientations = f == 0 ? 8 : 6;
      for (unsigned int r = 0; r < n_orientations; ++r)
        {
          deallog << "Face " << f << " in orientation " << r;
          const unsigned int orientation = r;
          const unsigned int face_no     = f;

          Triangulation<dim> tria;

          Triangulation<3> dummy;
          GridGenerator::reference_cell(dummy, ReferenceCells::Pyramid);

          auto vertices = dummy.get_vertices();

          std::vector<CellData<3>> cells;

          {
            CellData<3> cell;
            cell.vertices    = {0, 1, 2, 3, 4};
            cell.material_id = 0;
            cells.push_back(cell);
          }
          create_triangulation(vertices, cells, face_no, orientation, tria);

          bool success = true;

          auto cell = tria.begin();
          cell++;

          std::vector<unsigned int> verticess;

          for (const auto v : cell->vertex_indices())
            verticess.emplace_back(cell->vertex_index(v));

          const unsigned int lines_per_face = face_no == 0 ? 4 : 3;
          for (unsigned int ll = 0; ll < lines_per_face; ++ll)
            {
              const unsigned int l =
                cell->reference_cell().face_to_cell_lines(face_no,
                                                          ll,
                                                          orientation);

              // calls face_to_cell_line_orientation()
              const auto orientation_exp = cell->line_orientation(l);

              std::pair<unsigned int, unsigned int> p0;
              p0.first =
                verticess[cell->reference_cell().line_to_cell_vertices(l, 0)];
              p0.second =
                verticess[cell->reference_cell().line_to_cell_vertices(l, 1)];

              std::pair<unsigned int, unsigned int> p1;
              p1.first  = cell->line(l)->vertex_index(0);
              p1.second = cell->line(l)->vertex_index(1);

              if (orientation_exp == numbers::reverse_line_orientation)
                std::swap(p1.first, p1.second);

              success &= (p0 == p1);
            }

          if (success)
            deallog << " success" << std::endl;
          else
            deallog << " failure" << std::endl;
        }
    }

  deallog << std::endl;
}

int
main()
{
  initlog();

  test();
}
