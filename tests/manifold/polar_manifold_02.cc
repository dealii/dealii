// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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


// Test that the flat manifold does what it should on a sphere.

#include "../tests.h"



// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

// Helper function
template <int dim, int spacedim>
void test(unsigned int ref=1)
{
  PolarManifold<dim,spacedim> manifold;

  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_ball (tria);
  tria.reset_manifold(0);

  typename Triangulation<dim,spacedim>::active_cell_iterator cell;

  for (cell = tria.begin_active(); cell != tria.end(); ++cell)
    cell->set_all_manifold_ids(1);

  for (cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      if (cell->center().distance(Point<spacedim>()) < 1e-10)
        cell->set_all_manifold_ids(0);
    }

  tria.set_manifold(1, manifold);
  tria.refine_global(2);

  GridOut gridout;
  gridout.write_msh(tria, deallog.get_file_stream());
}

int main ()
{
  initlog();

  test<2,2>();
  test<3,3>();

  return 0;
}

