// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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



// Test extracted orginally from coarsening_02 on the
// branch_distributed_grids. we have a triangulation with
// Triangulation<dim>::limit_level_difference_at_vertices set, but produce a
// mesh in which this isn't honored
//
// here, test only the 2d/3d cases and do the 1d case in a separate test since
// it uses different code paths


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/base/utilities.h>

#include <fstream>
#include <cstdlib>


template <int dim>
bool
satisfies_level1_at_vertex_rule (const Triangulation<dim> &tr)
{
  std::vector<unsigned int> min_adjacent_cell_level (tr.n_vertices(),
                                                     tr.n_levels());
  std::vector<unsigned int> max_adjacent_cell_level (tr.n_vertices(),
                                                     0);

  for (typename Triangulation<dim>::active_cell_iterator
       cell = tr.begin_active();
       cell != tr.end(); ++cell)
    for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
      {
        min_adjacent_cell_level[cell->vertex_index(v)]
          = std::min<unsigned int>
            (min_adjacent_cell_level[cell->vertex_index(v)],
             cell->level());
        max_adjacent_cell_level[cell->vertex_index(v)]
          = std::max<unsigned int> (min_adjacent_cell_level[cell->vertex_index(v)],
                                    cell->level());
      }

  for (unsigned int k=0; k<tr.n_vertices(); ++k)
    if (tr.vertex_used(k))
      if (max_adjacent_cell_level[k] -
          min_adjacent_cell_level[k] > 1)
        return false;
  return true;
}


template<int dim>
void test()
{
  Triangulation<dim> triangulation (Triangulation<dim>::limit_level_difference_at_vertices);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global (2);

  const unsigned int n_refinements[] = { 0, 13, 11, 6 };
  for (unsigned int i=0; i<n_refinements[dim]; ++i)
    {
      // refine one-fifth of cells randomly
      std::vector<bool> flags (triangulation.n_active_cells(), false);
      for (unsigned int k=0; k<flags.size()/5 + 1; ++k)
        flags[Testing::rand() % flags.size()] = true;
      // make sure there's at least one that
      // will be refined
      flags[0] = true;

      // refine triangulation
      unsigned int index=0;
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active();
           cell != triangulation.end(); ++cell, ++index)
        if (flags[index])
          cell->set_refine_flag();
      Assert (index == triangulation.n_active_cells(), ExcInternalError());

      // flag all other cells for coarsening
      // (this should ensure that at least
      // some of them will actually be
      // coarsened)
      index=0;
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active();
           cell != triangulation.end(); ++cell, ++index)
        if (!flags[index])
          cell->set_coarsen_flag();

      triangulation.execute_coarsening_and_refinement ();

      // verify that none of the cells
      // violates the level-1-at-vertex rule
      Assert (satisfies_level1_at_vertex_rule (triangulation),
              ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int main()
{
  initlog();
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
