// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2020 by the deal.II authors
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


// generate and refine a hyper shell in 3d with 96 cells.
//
// this mesh ran into an interesting problem (somewhat hard to find
// :-( ) in that it violates the condition
//   cell->face(f)->at_boundary() == cell->at_boundary(f)
// when looping over all cells. The underlying cause was that the
// colorize option for the mesh was causing trouble. Verify that this
// problem has been fixed

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
check(const unsigned int n)
{
  deallog << "n=" << n << std::endl;

  Point<dim>         center;
  Triangulation<dim> tria(Triangulation<dim>::none);
  GridGenerator::hyper_shell(tria, center, 0.5, 1, n, true);

  // this is the test that failed
  // before
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      AssertThrow(cell->face(f)->at_boundary() == cell->at_boundary(f),
                  ExcInternalError());

  // also output something slightly
  // more useful
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        deallog << cell->face(f) << ' ' << (int)cell->face(f)->boundary_id()
                << ' ' << cell->face(f)->center().norm() << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<3>(6);
  check<3>(12);
  check<3>(96);
}
