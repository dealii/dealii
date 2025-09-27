// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Plot some GNUPLOT output to make sure that all of the flags work as
// intended. This verifies that the new collinear point removal algorithm works
// and also that we can do output in dim = 2, spacedim = 3.


#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

// overloads to get multiple grids for multiple dim and spacedim combinations
void
make_grid(Triangulation<2, 2> &triangulation)
{
  GridGenerator::hyper_shell(triangulation, Point<2>(), 2.0, 6.0, 12);
}

void
make_grid(Triangulation<2, 3> &triangulation)
{
  GridGenerator::hyper_sphere(triangulation, Point<3>(), 6.0);
  triangulation.refine_global(1); // need more cells
}

void
make_grid(Triangulation<3, 3> &triangulation)
{
  GridGenerator::hyper_shell(triangulation, Point<3>(), 2.0, 6.0, 12);
  triangulation.refine_global(0);
}



template <int dim, int spacedim = dim>
void
gnuplot_output(const GridOutFlags::Gnuplot &flags)
{
  deallog << "Triangulation<" << dim << ", " << spacedim << '>' << std::endl;

  Triangulation<dim, spacedim> triangulation;
  make_grid(triangulation);

  MappingQ<dim, spacedim> mapping(3);

  auto cell = triangulation.begin_active();
  cell->set_refine_flag(); // 0
  ++cell;
  ++cell;
  cell->set_refine_flag(); // 2
  ++cell;
  ++cell;
  cell->set_refine_flag(); // 4
  triangulation.execute_coarsening_and_refinement();

  GridOut grid_out;
  grid_out.set_flags(flags);
  grid_out.write_gnuplot(triangulation, deallog.get_file_stream(), &mapping);
  deallog << std::endl << std::endl;
}



int
main()
{
  initlog();

  {
    deallog << "don't curve anything" << std::endl;
    GridOutFlags::Gnuplot flags;
    flags.n_extra_curved_line_points = 0;

    flags.curved_inner_cells              = false;
    flags.write_additional_boundary_lines = false;
    gnuplot_output<2>(flags);
    gnuplot_output<2, 3>(flags);
    gnuplot_output<3>(flags);
  }

  {
    deallog << "don't curve inner lines and don't output extra boundary lines"
            << std::endl;
    GridOutFlags::Gnuplot flags;
    flags.n_extra_curved_line_points      = 1;
    flags.curved_inner_cells              = false;
    flags.write_additional_boundary_lines = false;
    gnuplot_output<2>(flags);
    gnuplot_output<2, 3>(flags);
    gnuplot_output<3>(flags);
  }

  {
    GridOutFlags::Gnuplot flags;
    deallog << "don't curve inner lines but output extra boundary lines"
            << std::endl;
    flags.n_extra_curved_line_points      = 1;
    flags.curved_inner_cells              = false;
    flags.write_additional_boundary_lines = true;
    gnuplot_output<2>(flags);
    gnuplot_output<2, 3>(flags); // not implemented, but produces output
    gnuplot_output<3>(flags);
  }

  {
    GridOutFlags::Gnuplot flags;
    deallog << "curve lines and output extra boundary lines" << std::endl;
    flags.n_extra_curved_line_points      = 1;
    flags.curved_inner_cells              = true;
    flags.write_additional_boundary_lines = true;
    gnuplot_output<2>(flags);
    gnuplot_output<2, 3>(flags);
    gnuplot_output<3>(flags);
  }

  deallog << "OK" << std::endl;
}
