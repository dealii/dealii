// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2017 by the deal.II authors
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



// like find_cell_2, but with the strange meshes from the mesh_3d_* tests

#include "../tests.h"
#include "../grid/mesh_3d.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/fe/mapping_q.h>



void check (Triangulation<3> &tria)
{
  MappingQ<3> map(3);

  Point<3> p(1./3.,1./2.,-1./5.);

  std::pair<Triangulation<3>::active_cell_iterator, Point<3> > cell
    = GridTools::find_active_cell_around_point (map, tria, p);

  deallog << cell.first << std::endl;
  for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_cell; ++v)
    deallog << "<" << cell.first->vertex(v) << "> ";
  deallog << "[ " << cell.second << "] ";
  deallog << std::endl;

  Assert (p.distance (cell.first->center()) < cell.first->diameter()/2,
          ExcInternalError());
}


int main ()
{
  initlog();

  {
    Triangulation<3> coarse_grid;
    create_two_cubes (coarse_grid);
    coarse_grid.refine_global (1);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_L_shape (coarse_grid);
    coarse_grid.refine_global (1);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball (coarse_grid);
    coarse_grid.reset_all_manifolds();
    coarse_grid.refine_global (1);
    check (coarse_grid);
  }

}



