// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check that refinement of the flat surfaces of a
// GridGenerator::half_hyper_ball works correctly using SphericalManifold as
// boundary description.

#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test(std::ostream &out)
{
  Point<dim> p1;
  p1[0] = 2.;
  if (dim > 1)
    p1[1] = -1.;
  if (dim > 2)
    p1[2] = 0.;
  Point<dim> p2;
  p2[0] = 3.;
  if (dim > 1)
    p2[1] = 2.;
  if (dim > 2)
    p2[2] = 4.;
  Point<dim> p3;
  p3[0] = 2.;
  if (dim > 1)
    p3[1] = 1.;
  if (dim > 2)
    p3[2] = 4.;

  SphericalManifold<dim> boundary_description(p1);
  GridOut                go;
  GridOut::OutputFormat  format = GridOut::gnuplot;

  {
    deallog << "half_hyper_ball" << std::endl;
    Triangulation<dim> tr;
    GridGenerator::half_hyper_ball(tr, p1, 3.);
    tr.set_manifold(0, boundary_description);

    tr.refine_global(2);
    go.write(tr, out, format);
  }
}


int
main()
{
  initlog();

  deallog.push("3d");
  test<3>(deallog.get_file_stream());
  deallog.pop();
}
