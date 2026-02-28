// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test ReferenceCell::face_to_cell_line_orientation()
#include <deal.II/grid/reference_cell.h>

#include "../tests.h"


void
test(const ReferenceCell reference_cell)
{
  for (const auto &face_no : reference_cell.face_indices())
    for (const auto &face_line_no :
         reference_cell.face_reference_cell(face_no).line_indices())
      for (unsigned int combined_face_orientation = 0;
           combined_face_orientation <
           reference_cell.n_face_orientations(face_no);
           ++combined_face_orientation)
        for (const auto &line_orientation :
             {numbers::default_geometric_orientation,
              numbers::reverse_line_orientation})
          {
            types::geometric_orientation face_to_cell_line_orientation = 0;

            if ((reference_cell == ReferenceCells::Tetrahedron &&
                 (face_no == 0 ||
                  (face_no == 1 && (face_line_no == 1 || face_line_no == 2)) ||
                  (face_no == 2 && face_line_no == 1))) ||
                (reference_cell == ReferenceCells::Pyramid &&
                 (face_no == 0 ||
                  ((face_no == 1 || face_no == 2) &&
                   (face_line_no == 1 || face_line_no == 2)))) ||
                reference_cell.is_hyper_cube())

              {
                face_to_cell_line_orientation =
                  reference_cell.face_to_cell_line_orientation(
                    face_line_no,
                    face_no,
                    combined_face_orientation,
                    line_orientation);
                deallog << "Face number: " << face_no
                        << " face line number: " << face_line_no
                        << " combined_face_orientation "
                        << combined_face_orientation << " line orientation "
                        << std::to_string(line_orientation)
                        << " face_to_cell_line_orientation: "
                        << std::to_string(face_to_cell_line_orientation)
                        << std::endl;
              }
          }
  deallog << std::endl;
}

int
main()
{
  initlog();

  deallog.push("Hexahedron");
  test(ReferenceCells::Hexahedron);
  deallog.pop();

  deallog.push("Tetrahedron");
  test(ReferenceCells::Tetrahedron);
  deallog.pop();

  deallog.push("Pyramid");
  test(ReferenceCells::Pyramid);
  deallog.pop();
}
