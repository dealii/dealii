// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// in implementing GMSH version 2 input, we forgot reading vertices in
// version 2 format as well. A fix by Victor Prosolin helped this
// problem

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
  std::string source_dir(SOURCE_DIR "/");
  std::string relative_name(name.begin() + source_dir.size(), name.end());

  deallog << relative_name << '\t' << tria.n_vertices() << '\t'
          << tria.n_cells() << std::endl;

  GridOut grid_out;
  grid_out.write_gnuplot(tria, deallog.get_file_stream());
}

void
filename_resolution()
{
  check_file<2>(std::string(SOURCE_DIR "/grid_in_msh_version_1/input_v1.msh"),
                GridIn<2>::msh);
  check_file<2>(std::string(SOURCE_DIR "/grid_in_msh_version_1/input_v2.msh"),
                GridIn<2>::msh);
}


int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(5);

  filename_resolution();
}
