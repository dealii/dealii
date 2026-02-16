// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <string>

#include "../tests.h"

template <int dim>
void
test(const std::string &infilename)
{
  Triangulation<dim> tria;
  GridIn<dim>        gi;
  gi.attach_triangulation(tria);
  gi.read(infilename);

  deallog.get_file_stream()
    << "------------------------------------------" << std::endl;

  GridOut grid_out;
  grid_out.set_flags(GridOutFlags::Ucd(true));
  grid_out.write_ucd(tria, deallog.get_file_stream());
}

int
main()
{
  initlog();
  test<2>(std::string(SOURCE_DIR "/grid_in_tecplot/1.dat"));
  test<2>(std::string(SOURCE_DIR "/grid_in_tecplot/2.dat"));
  test<2>(std::string(SOURCE_DIR "/grid_in_tecplot/3.dat"));
  test<2>(std::string(SOURCE_DIR "/grid_in_tecplot/4.dat"));
}
