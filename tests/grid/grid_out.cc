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



#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim>                  tria;
  static const SphericalManifold<dim> x;
  if (dim == 2)
    {
      tria.set_manifold(0, x);
      GridGenerator::hyper_ball(tria);
    }
  else
    GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  GridOut              grid_out;
  GridOutFlags::Eps<2> eps2(
    GridOutFlags::EpsFlagsBase::width, 300, .5, false, 5, true);
  grid_out.set_flags(eps2);

  if (dim != 1)
    grid_out.write_eps(tria, deallog.get_file_stream());
  grid_out.write_gnuplot(tria, deallog.get_file_stream());
  grid_out.set_flags(GridOutFlags::Ucd(true));
  grid_out.write_ucd(tria, deallog.get_file_stream());
  if (dim != 1)
    grid_out.write_dx(tria, deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(2);

  test<1>();
  test<2>();
  test<3>();
}
