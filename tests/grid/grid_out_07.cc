// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test GridOut::write_vtu

#include <deal.II/base/geometry_info.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim, int spacedim>
void
test(std::ostream &logfile)
{
  Triangulation<dim, spacedim> tria;
  std::vector<unsigned int>    legs(2 * dim, 1);
  if (dim > 1)
    GridGenerator::hyper_cross(tria, legs, true);
  else
    GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  GridOut           grid_out;
  GridOutFlags::Vtu vtu_flags;
  vtu_flags.compression_level = DataOutBase::CompressionLevel::best_compression;
  grid_out.set_flags(vtu_flags);
  grid_out.write_vtu(tria, logfile);
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
