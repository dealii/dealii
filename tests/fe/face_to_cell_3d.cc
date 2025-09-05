/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2020 - 2025 by the deal.II authors
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

// Test that the new function (FETools::face_to_cell_index()) computes correct
// values for twisted hexahedral grids with periodic boundaries.

#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/tria.h>

#include <numeric>
#include <vector>

#include "../tests.h"

int
main()
{
  initlog();

  const auto reference_cell = ReferenceCells::Hexahedron;

  const Point<3>                  corner_p1(-0.5, 1.0, 0.0);
  const Point<3>                  corner_p2(0.5, 2.0, 1.2);
  const std::vector<unsigned int> sub_divisions{2, 2, 2};

  const std::vector<Point<3>> base_vertices{{0.0, 0.0, 0.0},
                                            {1.0, 0.0, 0.0},
                                            {0.0, 1.0, 0.0},
                                            {1.0, 1.0, 0.0}};
  for (types::geometric_orientation combined_orientation = 0;
       combined_orientation < reference_cell.n_face_orientations(0);
       ++combined_orientation)
    {
      deallog << "orientation = " << int(combined_orientation) << std::endl;
      // set up vertices in lexical order:
      std::vector<Point<3>> vertices(base_vertices.begin(),
                                     base_vertices.end());
      for (const auto &base_vertex : base_vertices)
        vertices.emplace_back(base_vertex[0], base_vertex[1], 1.0);
      for (const auto &base_vertex : base_vertices)
        vertices.emplace_back(base_vertex[0], base_vertex[1], 2.0);

      // set up the first cell:
      std::vector<CellData<3>> cells(2);
      std::iota(cells[0].vertices.begin(), cells[0].vertices.end(), 0u);

      // set up the second cell:
      const std::vector<unsigned int> middle_vertices{4u, 5u, 6u, 7u};
      const auto                      permuted_vertices =
        reference_cell.face_reference_cell(0).permute_by_combined_orientation(
          make_array_view(middle_vertices), combined_orientation);
      std::copy(permuted_vertices.begin(),
                permuted_vertices.end(),
                cells[1].vertices.begin());
      for (unsigned int vertex_no = 0; vertex_no < 4; ++vertex_no)
        cells[1].vertices[vertex_no + 4] = cells[1].vertices[vertex_no] + 4;

      // set up the FE things for the rest of the test:
      Triangulation<3> tria;
      tria.create_triangulation(vertices, cells, SubCellData());

      const std::array<unsigned int, 3> fe_degrees{{1, 2, 3}};

      for (const auto &fe_degree : fe_degrees)
        {
          FE_Q<3> fe(fe_degree);
          deallog << "FE = " << fe.get_name() << std::endl;
          DoFHandler<3> dof_handler(tria);
          dof_handler.distribute_dofs(fe);
          const auto support_points = DoFTools::map_dofs_to_support_points(
            reference_cell.get_default_linear_mapping<3>(), dof_handler);
          const auto cell_0    = dof_handler.begin_active();
          const auto cell_1    = ++dof_handler.begin_active();
          const auto face_0_no = 5;
          const auto face_1_no = 4;
          // these assertions verify the test is set up correctly
          Assert(cell_0->face(face_0_no) == cell_1->face(face_1_no),
                 ExcInternalError());
          Assert(cell_0->combined_face_orientation(face_0_no) ==
                   numbers::default_geometric_orientation,
                 ExcInternalError());
          Assert(cell_1->combined_face_orientation(face_1_no) ==
                   combined_orientation,
                 ExcInternalError());

          // now for the actual test:
          std::vector<types::global_dof_index> dofs_0(fe.n_dofs_per_cell()),
            dofs_1(fe.n_dofs_per_cell());
          cell_0->get_dof_indices(dofs_0);
          cell_1->get_dof_indices(dofs_1);

          // We don't yet implement rotation of dofs on quads in this function
          // so, in that case, ignore them
          const auto max_face_dof =
            fe.n_dofs_per_quad() == 1 ?
              fe.n_dofs_per_face() :
              fe.n_dofs_per_face() - fe.n_dofs_per_quad();
          for (unsigned int i = 0; i < max_face_dof; ++i)
            {
              const auto i0 = dofs_0[fe.face_to_cell_index(
                i, face_0_no, cell_0->combined_face_orientation(face_0_no))];
              const auto i1 = dofs_1[fe.face_to_cell_index(
                i, face_1_no, cell_1->combined_face_orientation(face_1_no))];
              Assert(i0 == i1, ExcInternalError());
              deallog << std::setw(2) << i0 << ":" << (i0 == i1) << ":"
                      << support_points.at(i0) << std::endl;
            }
        }

      deallog << std::endl;
    }
  deallog << "OK" << std::endl;
}
