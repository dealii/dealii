// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// check CylindricalManifold<3>::normal_vector()

#include "../tests.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

void
check()
{
  constexpr int dim = 3;

  Triangulation<dim> triangulation;
  GridGenerator::truncated_cone(triangulation);
  GridTools::transform(
    [](const Point<3>& p) { return Point<3>(-p[1], p[0], p[2]); },
    triangulation);
  static const CylindricalManifold<dim> boundary(1);
  triangulation.set_manifold(0, boundary);

  triangulation.refine_global(2);

  for(const auto cell : triangulation.active_cell_iterators())
    for(unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
        ++face_no)
      {
        const auto face = cell->face(face_no);
        if(face->boundary_id() == 0)
          for(unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
            {
              const Point<dim>     vertex = face->vertex(v);
              const Tensor<1, dim> tangent_1({-vertex(2), 0., vertex(0)});
              const Tensor<1, dim> tangent_2 = vertex - Point<dim>(0, 3, 0);

              // get the normal vector and test it
              const Tensor<1, 3> normal_vector
                = boundary.normal_vector(face, vertex);

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
