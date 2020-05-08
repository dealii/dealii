// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2019 by the deal.II authors
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


// test the draq_*-flags for GridOut::write_vtk

#include <deal.II/base/geometry_info.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim, int spacedim>
void
do_test(std::ostream &logfile,
        const bool    output_cells,
        const bool    output_faces,
        const bool    output_edges,
        const bool    output_only_relevant)
{
  Triangulation<dim, spacedim> tria_1;
  GridGenerator::hyper_cube(tria_1, 0, 1, true);
  tria_1.refine_global(1);

  tria_1.begin_active()->set_refine_flag();
  tria_1.execute_coarsening_and_refinement();

  GridOutFlags::Vtk flags;
  flags.output_cells         = output_cells;
  flags.output_faces         = output_faces;
  flags.output_edges         = output_edges;
  flags.output_only_relevant = output_only_relevant;

  GridOut grid_out_1;
  grid_out_1.set_flags(flags);
  grid_out_1.write_vtk(tria_1, logfile);

  // write to buffer
  std::ostringstream buf;
  grid_out_1.write_vtk(tria_1, buf);


  std::istringstream buf_in;
  buf_in.str(buf.str());


  // read from buffer and create a new triangulation
  Triangulation<dim, spacedim> tria_2;
  GridIn<dim, spacedim>        grid_in;
  grid_in.attach_triangulation(tria_2);
  grid_in.read_vtk(buf_in);

  // output the new triangulation
  GridOut grid_out_2;
  grid_out_2.write_vtk(tria_2, logfile);
}



template <int dim, int spacedim>
void
test(std::ostream &logfile)
{
  do_test<dim, spacedim>(logfile, true, false, false, true);
  do_test<dim, spacedim>(logfile, true, false, true, true);
  do_test<dim, spacedim>(logfile, true, true, false, true);
  do_test<dim, spacedim>(logfile, true, true, true, true);

  do_test<dim, spacedim>(logfile, true, false, false, false);
  do_test<dim, spacedim>(logfile, true, false, true, false);
  do_test<dim, spacedim>(logfile, true, true, false, false);
  do_test<dim, spacedim>(logfile, true, true, true, false);
}



int
main()
{
  initlog("output");
  test<1, 1>(deallog.get_file_stream());
  test<1, 2>(deallog.get_file_stream());
  test<2, 2>(deallog.get_file_stream());
  test<2, 3>(deallog.get_file_stream());
  test<3, 3>(deallog.get_file_stream());
}
