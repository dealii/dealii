// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check CylindricalManifold<3>::normal_vector()

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"


void
check()
{
  constexpr int dim = 3;

  Triangulation<dim> triangulation;
  GridGenerator::truncated_cone(triangulation);
  GridTools::transform(
    [](const Point<3> &p) { return Point<3>(-p[1], p[0], p[2]); },
    triangulation);
  static const CylindricalManifold<dim> boundary(1);
  triangulation.set_manifold(0, boundary);

  triangulation.refine_global(2);

  for (const auto &cell : triangulation.active_cell_iterators())
    for (const unsigned int face_no : GeometryInfo<dim>::face_indices())
      {
        const auto face = cell->face(face_no);
        if (face->boundary_id() == 0)
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face;
               ++v)
            {
              const Point<dim>     vertex = face->vertex(v);
              const Tensor<1, dim> tangent_1({-vertex[2], 0., vertex[0]});
              const Tensor<1, dim> tangent_2 = vertex - Point<dim>(0, 3, 0);

              // get the normal vector and test it
              const Tensor<1, 3> normal_vector =
                boundary.normal_vector(face, vertex);

              Assert(std::fabs(normal_vector.norm() - 1) < 1e-12,
                     ExcInternalError());
              Assert(std::fabs(normal_vector * tangent_1) < 1e-12,
                     ExcInternalError());
              Assert(std::fabs(normal_vector * tangent_2) < 1e-12,
                     ExcInternalError());
            }
      }
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  check();
}
