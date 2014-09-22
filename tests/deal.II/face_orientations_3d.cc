// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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



// just like grid_in_3d, but count the number of misoriented faces in
// the meshes that can be oriented

#include "../tests.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iomanip>


std::ofstream logfile("output");


void test (const char *filename)
{
  Triangulation<3> tria;
  GridIn<3> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (filename);

  gi.read_xda (in);
  deallog << "  " << tria.n_active_hexs() << " active cells" << std::endl;
  deallog << "  " << tria.n_active_quads() << " active faces" << std::endl;

  unsigned int misoriented_faces = 0;
  for (Triangulation<3>::active_cell_iterator cell=tria.begin_active();
       cell!=tria.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
      if (cell->face_orientation (f) == false)
        {
          ++misoriented_faces;

          // check that the face is
          // correctly oriented from
          // the other side at
          // least. note that if this
          // face is misoriented,
          // then there must be a
          // neighbor over there
          Assert (cell->neighbor(f)
                  ->face_orientation(cell->neighbor_of_neighbor(f))
                  == true,
                  ExcInternalError());
        }
  deallog << "  " << misoriented_faces << " misoriented faces" << std::endl;
}


int main ()
{
  deallog << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test (SOURCE_DIR "/grid_in_3d/1.in");
  test (SOURCE_DIR "/grid_in_3d/2.in");
  test (SOURCE_DIR "/grid_in_3d/3.in");
  test (SOURCE_DIR "/grid_in_3d/4.in");

  test (SOURCE_DIR "/grid_in_3d/evil_0.in");
  test (SOURCE_DIR "/grid_in_3d/evil_4.in");
}

