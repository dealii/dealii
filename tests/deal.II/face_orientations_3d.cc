//----------------------------  face_orientations_3d.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  face_orientations_3d.cc  ---------------------------


// just like grid_in_3d, but count the number of misoriented faces in
// the meshes that can be oriented

#include "../tests.h"
#include <dofs/dof_handler.h>
#include <grid/tria.h>
#include <grid/tria_boundary.h>
#include <grid/tria_boundary_lib.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_out.h>
#include <grid/grid_in.h>
#include <grid/grid_generator.h>
#include <base/logstream.h>

#include <fstream>


std::ofstream logfile("face_orientations_3d.output");


void test (const char *filename)
{
  deallog << "Reading " << filename << std::endl;
  
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
  logfile.precision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ("grid_in_3d_1.in");
  test ("grid_in_3d_2.in");
  test ("grid_in_3d_3.in");
  test ("grid_in_3d_4.in");

  test ("grid_in_3d_evil_0.in");
  test ("grid_in_3d_evil_4.in");
}

