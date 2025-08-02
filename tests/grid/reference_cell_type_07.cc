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


// Test ReferenceCell::subface_vertex_location()

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test(const ReferenceCell &reference_cell)
{
  const bool check_vs_tria =
    reference_cell.is_simplex() || reference_cell.is_hyper_cube();
  Triangulation<dim> tria;
  if (check_vs_tria)
    {
      GridGenerator::reference_cell(tria, reference_cell);
      tria.refine_global(1);
    }

  deallog << reference_cell.to_string() << std::endl;
  for (const unsigned int face_no : reference_cell.face_indices())
    {
      deallog << "face_no = " << face_no << std::endl;
      for (unsigned int subface_no = 0;
           subface_no <
           (reference_cell.face_reference_cell(face_no).n_isotropic_children());
           ++subface_no)
        {
          deallog << "  subface_no = " << subface_no << std::endl;
          std::vector<Point<dim>> subface_vertices;
          for (const unsigned int subface_vertex_no :
               reference_cell.face_reference_cell(face_no).vertex_indices())
            subface_vertices.push_back(
              reference_cell.subface_vertex_location<dim>(
                face_no,
                subface_no,
                subface_vertex_no,
                RefinementCase<dim- 1>::isotropic_refinement));

          for (const auto &vertex : subface_vertices)
            deallog << "    " << vertex << std::endl;

          if (check_vs_tria)
            {
              const auto subface =
                tria.begin(0)->face(face_no)->child(subface_no);

              for (const unsigned int vertex_no : subface->vertex_indices())
                Assert(subface->vertex(vertex_no) ==
                         subface_vertices[vertex_no],
                       ExcMessage("should match"));
            }
        }
    }
}

int
main()
{
  initlog();

  {
    deallog.push("2D");
    test<2>(ReferenceCells::Quadrilateral);
    test<2>(ReferenceCells::Triangle);
    deallog.pop();
  }

  {
    deallog.push("3D");
    test<3>(ReferenceCells::Tetrahedron);
    test<3>(ReferenceCells::Pyramid);
    test<3>(ReferenceCells::Wedge);
    test<3>(ReferenceCells::Hexahedron);
    deallog.pop();
  }
}
