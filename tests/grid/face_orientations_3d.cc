// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// just like grid_in_3d, but count the number of misoriented faces in
// the meshes that can be oriented

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



void
test(const char *filename)
{
  Triangulation<3> tria;
  GridIn<3>        gi;
  gi.attach_triangulation(tria);
  std::ifstream in(filename);

  gi.read_xda(in);
  deallog << "  " << tria.n_active_hexs() << " active cells" << std::endl;
  deallog << "  " << tria.n_active_quads() << " active faces" << std::endl;

  unsigned int misoriented_faces = 0;
  for (Triangulation<3>::active_cell_iterator cell = tria.begin_active();
       cell != tria.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<3>::face_indices())
      if (cell->face_orientation(f) == false)
        {
          ++misoriented_faces;

          // check that the face is
          // correctly oriented from
          // the other side at
          // least. note that if this
          // face is misoriented,
          // then there must be a
          // neighbor over there
          AssertThrow(cell->neighbor(f)->face_orientation(
                        cell->neighbor_of_neighbor(f)) == true,
                      ExcInternalError());
        }
  deallog << "  " << misoriented_faces << " misoriented faces" << std::endl;
}


int
main()
{
  deallog << std::setprecision(2);
  initlog();

  test(SOURCE_DIR "/grid_in_3d/1.in");
  test(SOURCE_DIR "/grid_in_3d/2.in");
  test(SOURCE_DIR "/grid_in_3d/3.in");
  test(SOURCE_DIR "/grid_in_3d/4.in");

  test(SOURCE_DIR "/grid_in_3d/evil_0.in");
  test(SOURCE_DIR "/grid_in_3d/evil_4.in");
}
