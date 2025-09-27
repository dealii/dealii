// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// read a 3d file in the VTK format

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <string>

#include "../tests.h"

template <int dim>
void
check_file(const std::string name, typename GridIn<dim>::Format format)
{
  Triangulation<dim> tria;
  GridIn<dim>        gi;
  gi.attach_triangulation(tria);
  gi.read(name, format);
  deallog << '\t' << tria.n_vertices() << '\t' << tria.n_cells() << std::endl;

  GridOut grid_out;
  grid_out.write_gnuplot(tria, deallog.get_file_stream());

  // now also output boundary ids for all boundary faces. this is
  // where the bug was for which this test was written
  deallog << "Face boundary_ids:" << std::endl;
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        deallog << cell << ' ' << f << ": " << (int)cell->face(f)->boundary_id()
                << std::endl;
}

void
filename_resolution()
{
  check_file<3>(std::string(SOURCE_DIR "/grid_in_vtk_3d_04.vtk"),
                GridIn<3>::vtk);
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  filename_resolution();
}
