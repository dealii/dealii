// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Create a Triangulation<2,3> from an implicit function

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


class ImplicitFunction : public Function<3>
{
public:
  virtual double
  value(const Point<3> &p, const unsigned int component = 0) const override
  {
    return std::pow(1 - std::sqrt(p[0] * p[0] + p[1] * p[1]), 2) + p[2] * p[2] -
           .25;
  }
};

int
main()
{
  initlog();
  // Build a deal.II triangulation
  Triangulation<2, 3>             tria;
  ImplicitFunction                implicit_function;
  CGALWrappers::AdditionalData<2> data;
  data.angular_bound  = 30.;
  data.radius_bound   = .1;
  data.distance_bound = .1;
  GridGenerator::implicit_function(
    tria, implicit_function, data, Point<3>(1, 0, 0), 10.);
  {
    GridOut       go;
    std::ofstream of("tria.vtk");
    go.write_vtk(tria, of);
  }
  remove("tria.vtk");
  //  If we got here, everything was ok, including writing the grid.
  deallog << "OK" << std::endl;
}
