// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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


// test copy_marked_cells

#include "../tests.h"
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iomanip>
#include <cstdio>

std::ofstream logfile("output");



template <int dim, int spacedim>
void test ()
{
  Triangulation<dim, spacedim> tria;
  std::vector<unsigned int> reps(dim, 1);
  reps[0]=3;
  Point<dim> p1;
  Point<dim> p2;
  for (unsigned int d=0; d<dim; ++d)
    p2[d] = 1.0;
  p2[0] = 3.0;

  GridGenerator::subdivided_hyper_rectangle (tria, reps,
                                             p1, p2, true);
  tria.begin_active()->set_refine_flag();
  (++tria.begin_active())->set_refine_flag();
  tria.execute_coarsening_and_refinement ();

  typename Triangulation<dim>::active_cell_iterator
  cell = tria.begin_active();

  for (unsigned int i=0; i < tria.n_active_cells(); ++i)
    {
      if (i < tria.n_active_cells()/2-1)
        {
          cell->set_material_id(1+i);
          cell->set_user_flag();
        }
      else
        cell->set_material_id(0);

      if (i==0)
        cell->set_all_manifold_ids(17);

      cell->set_subdomain_id(i%2);
      cell->set_level_subdomain_id(i%2);
      ++cell;
    }
  GridOut gout;
  GridOutFlags::Svg svg_flags;

  svg_flags.coloring = GridOutFlags::Svg::level_number;
  svg_flags.label_cell_index = false;
  svg_flags.label_level_number = false;
  svg_flags.label_material_id = true;
  svg_flags.label_subdomain_id = true;
  svg_flags.label_level_subdomain_id = true;

  gout.set_flags(svg_flags);

  {
    std::ofstream f("in.svg");
    gout.write_svg(tria, f);
    gout.write_svg(tria, deallog.get_file_stream());
  }

  Triangulation<dim, spacedim> tria2;
  GridGenerator::create_mesh_from_marked_cells(tria2, tria);

  {
    std::ofstream f("out.svg");
    gout.write_svg(tria2, f);
    gout.write_svg(tria2, deallog.get_file_stream());
  }


  deallog << dim
          << tria2.n_levels()
          << tria2.n_active_cells()
          << std::endl;
}


int main ()
{
  initlog();

  //test<1>();
  test<2,2>();
  //test<3>();

  return 0;
}
