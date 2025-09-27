// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// make sure that we write boundary lines marked with a non-zero boundary
// indicator correctly in MSH format

#include <deal.II/base/geometry_info.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.begin_active()->line(0)->set_boundary_id(1);
  tria.begin_active()->face(2 * dim - 1)->set_boundary_id(2);

  GridOut           grid_out;
  GridOutFlags::Msh flags;
  flags.write_lines = flags.write_faces = true;
  grid_out.set_flags(flags);
  grid_out.write_msh(tria, deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test<2>();
  test<3>();
}
