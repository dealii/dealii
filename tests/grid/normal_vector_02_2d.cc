// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



int
main()
{
  initlog();
  deallog << std::setprecision(3);
  deallog << std::fixed;

  SphericalManifold<2> boundary(Point<2>(1, 0));

  Triangulation<2>               tria;
  Manifold<2>::FaceVertexNormals normals;

  GridGenerator::hyper_ball(tria, Point<2>(1, 0), 3);

  Triangulation<2>::active_cell_iterator cell = tria.begin_active();
  for (; cell != tria.end(); ++cell)
    for (const unsigned int face_no : GeometryInfo<2>::face_indices())
      if (cell->at_boundary(face_no))
        {
          Triangulation<2>::face_iterator face = cell->face(face_no);
          boundary.get_normals_at_vertices(face, normals);
          for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_face; ++v)
            AssertThrow((boundary.normal_vector(face, face->vertex(v)) -
                         normals[v] / normals[v].norm())
                            .norm() < 1e-12,
                        ExcInternalError());
        }

  deallog << "OK" << std::endl;
}
