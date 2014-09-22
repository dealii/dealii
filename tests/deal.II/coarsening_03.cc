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



// a simplified version of coarsening_03. what's happening is that we have a
// uniformly 4x4 refined mesh, mark the bottom left 4 cells for refinement,
// mark the top right 4 for coarsening, and then end up with a mesh that
// violates the Triangulation<dim>::limit_level_difference_at_vertices
// constraint we have requested

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

  for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
    triangulation.begin(1)->child(c)->set_refine_flag();
  for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
    (--triangulation.begin(2))->child(c)->set_coarsen_flag();

  triangulation.prepare_coarsening_and_refinement ();

  for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    deallog << cell << ' '
            << (cell->refine_flag_set() ?
                "to be refined" :
                (cell->coarsen_flag_set() ?
                 "to be coarsened" :
                 ""))
            << std::endl;

  triangulation.execute_coarsening_and_refinement ();

  // verify that none of the cells
  // violates the level-1-at-vertex rule
  Assert (satisfies_level1_at_vertex_rule (triangulation),
          ExcInternalError());
}


int main()
{
  initlog();
  deallog.threshold_double(1.e-10);

  deallog.push("1d");
  test<1>();
  deallog.pop();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
