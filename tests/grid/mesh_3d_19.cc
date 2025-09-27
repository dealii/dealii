// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// an attempt to understand the failure of mesh_3d_18

char logname[] = "output";


#include <deal.II/base/function.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <vector>

#include "../tests.h"



void
test_with_wrong_face_orientation()
{
  Triangulation<3> triangulation;
  GridGenerator::hyper_ball(triangulation);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  Triangulation<3>::active_cell_iterator cell = triangulation.begin_active();
  ++cell;
  ++cell;

  deallog << "cell=" << cell << std::endl;
  deallog << "cell->neighbor(3)=" << cell->neighbor(3) << std::endl;
  deallog << "cell->face(3)=" << cell->face(3) << std::endl;

  for (unsigned int i = 0; i < 6; ++i)
    deallog << "cell->neighbor(3)->face(" << i
            << ")=" << cell->neighbor(3)->face(i) << std::endl;
}



int
main()
{
  std::ofstream logfile(logname);
  deallog << std::setprecision(3);

  deallog.attach(logfile);

  test_with_wrong_face_orientation();

  deallog << "OK" << std::endl;
}
