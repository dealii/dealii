// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// similar to the hyper_ball_3d test, but for the cylinder grid. here,
// the cause for the failure was different, though: the description of
// the cells was wrong, and they were not sent through the
// GridTools::consistently_order_cells() function which would have cured the
// problem.


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



int
main()
{
  initlog();
  deallog << std::setprecision(2);

  // generate a cylinder
  Triangulation<3> tria;
  GridGenerator::cylinder(tria, std::sqrt(2.));

  // output all faces. here, we
  // should have 18 (two layers of
  // cells, each with 4 outer faces,
  // plus 5 faces each for the top
  // and bottom of the cylinder)
  unsigned int external_faces = 0;
  for (Triangulation<3>::face_iterator face = tria.begin_face();
       face != tria.end_face();
       ++face)
    {
      deallog << face << "   " << face->boundary_id() << "  " << '<'
              << face->vertex(0) << '>' << std::endl
              << "           <" << face->vertex(1) << '>' << std::endl
              << "           <" << face->vertex(2) << '>' << std::endl
              << "           <" << face->vertex(3) << '>' << std::endl;
      if (face->at_boundary())
        ++external_faces;
    }

  deallog << "External faces: " << external_faces << std::endl;

  Assert(external_faces == 18, ExcInternalError());

  return 0;
}
