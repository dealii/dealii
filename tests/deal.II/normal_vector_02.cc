// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>

#include <fstream>
#include <iomanip>




int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  HyperBallBoundary<3> boundary (Point<3>(1,0,0));

  Triangulation<3> tria;
  Boundary<3>::FaceVertexNormals normals;

  GridGenerator::hyper_ball (tria, Point<3>(1,0,0), 3);

  Triangulation<3>::active_cell_iterator cell=tria.begin_active();
  for (; cell!=tria.end(); ++cell)
    for (unsigned int face_no=0;
         face_no<GeometryInfo<3>::faces_per_cell; ++face_no)
      if (cell->at_boundary(face_no))
        {
          Triangulation<3>::face_iterator face = cell->face(face_no);
          boundary.get_normals_at_vertices(face, normals);
          for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_face; ++v)
            Assert ((boundary.normal_vector (face,
                                             face->vertex(v))
                     -
                     normals[v] / normals[v].norm()).norm()
                    <
                    1e-12,
                    ExcInternalError());
        }

  deallog << "OK" << std::endl;
}
