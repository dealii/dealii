// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// adapted from hp/crash_06, which itself is from
// make_hanging_node_constraints for hp-elements. used to crash. triggers the
// crash that at the time of writing the test afflicts all
// hp/hp_constraints_*_03 tests

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
  deallog << "cell->face_orientation(3)="
          << (cell->face_orientation(3) ? "true" : "false") << std::endl;

  const Triangulation<3>::active_cell_iterator neighbor_child =
    cell->neighbor_child_on_subface(3, 1);

  deallog << "cell->neighbor_child_on_subface(3,1)=" << neighbor_child
          << std::endl;
  deallog << "cell->neighbor_child_on_subface(3,1)->neighbor(5)="
          << neighbor_child->neighbor(5) << std::endl;

  deallog << "cell->neighbor_child_on_subface(3,1)->face_orientation(5)="
          << (neighbor_child->face_orientation(5) ? "true" : "false")
          << std::endl;

  deallog << "cell->face(3)=" << cell->face(3) << std::endl;
  for (unsigned int i = 0; i < 4; ++i)
    deallog << "cell->face(3)->child(" << i << ")=" << cell->face(3)->child(i)
            << std::endl;

  for (unsigned int i = 0; i < 6; ++i)
    deallog << "cell->neighbor(3)->face(" << i
            << ")=" << cell->neighbor(3)->face(i) << std::endl;

  for (unsigned int i = 0; i < 6; ++i)
    deallog << "cell->neighbor_child_on_subface(3,1)->face(" << i
            << ")=" << cell->neighbor_child_on_subface(3, 1)->face(i)
            << std::endl;

  // The following assertion was originally in
  // make_hanging_node_constraints and
  // triggered for the mesh and cell here at
  // hand. however, if one carefully reads the
  // new comment for
  // CellAccessor::neighbor_child_on_subface,
  // one realizes that it should be true in any
  // case, no matter whether
  // cell->face_orientation()==false or not for
  // the face we are presently considering. the
  // original assertion failed,
  // whereas with the new functionality of
  // CellAccessor::neighbor_child_on_subface it
  // should work. let's make sure we get the
  // status we expect.
  AssertThrow(cell->face(3)->child(1) ==
                neighbor_child->face(cell->neighbor_of_neighbor(3)),
              ExcInternalError());
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
