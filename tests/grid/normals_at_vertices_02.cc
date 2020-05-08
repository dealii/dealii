// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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

  SphericalManifold<3> boundary(Point<3>(1, 0, 0));

  Triangulation<3>               tria;
  Manifold<3>::FaceVertexNormals normals;

  GridGenerator::hyper_ball(tria, Point<3>(1, 0, 0), 3);

  Triangulation<3>::active_cell_iterator cell = tria.begin_active();
  for (; cell != tria.end(); ++cell)
    for (const unsigned int face_no : GeometryInfo<3>::face_indices())
      if (cell->at_boundary(face_no))
        {
          deallog << " Face" << face_no << std::endl;
          Triangulation<3>::face_iterator face = cell->face(face_no);
          boundary.get_normals_at_vertices(face, normals);
          for (unsigned int v = 0; v < GeometryInfo<3>::vertices_per_face; ++v)
            {
              deallog << "  vertex=" << face->vertex(v)
                      << ",  normal=" << normals[v] << std::endl;

              // note that we can't check
              // here that the normal vector
              // is, in fact, normalized,
              // since the function does not
              // actually guarantee that
            }
        }
}
