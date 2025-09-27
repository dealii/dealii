// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test that at the vertices, Manifold::normal_vector returns the same as
// Manifold::get_normals_at_vertices once the latter vectors are normalized



#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



void
create_triangulation(const unsigned int case_no, Triangulation<2> &tria)
{
  switch (case_no)
    {
      case 0:
        GridGenerator::hyper_cube(tria, 1., 3.);
        break;
      case 1:
        {
          GridGenerator::hyper_cube(tria, 1., 3.);
          Point<2> &v0 = tria.begin_active()->vertex(0);
          v0           = Point<2>(-0.5, -1);
          Point<2> &v1 = tria.begin_active()->vertex(1);
          v1           = Point<2>(0.25, 0.25);
          break;
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
    };
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);
  deallog << std::fixed;

  Triangulation<2>               tria;
  FlatManifold<2>                boundary;
  Manifold<2>::FaceVertexNormals normals;
  for (unsigned int case_no = 0; case_no < 2; ++case_no)
    {
      deallog << "Case" << case_no << std::endl;
      create_triangulation(case_no, tria);
      const Triangulation<2>::active_cell_iterator cell = tria.begin_active();
      Triangulation<2>::face_iterator              face;
      for (const unsigned int face_no : GeometryInfo<2>::face_indices())
        {
          face = cell->face(face_no);
          boundary.get_normals_at_vertices(face, normals);
          for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_face; ++v)
            AssertThrow((boundary.normal_vector(face, face->vertex(v)) -
                         normals[v] / normals[v].norm())
                            .norm() < 1e-12,
                        ExcInternalError());
        }
      tria.clear();
    }

  deallog << "OK" << std::endl;
}
