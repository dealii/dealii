// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test GridGenerator::moebius

#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
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


  // loop without rotation
  if (true)
    {
      deallog << "moebius, no rotation" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::moebius(tr, 20, 0, 10.0, 2.0);
      GridOut go;
      go.set_flags(GridOutFlags::Ucd(true));
      go.write_ucd(tr, out);
    }

  // loop with quarter rotation (1 * pi/2)
  if (true)
    {
      deallog << "---------------------------" << std::endl
              << "moebius, quarter rotation (1* Pi/2)" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::moebius(tr, 20, 1, 10.0, 2.0);
      GridOut go;
      go.set_flags(GridOutFlags::Ucd(true));
      go.write_ucd(tr, out);
    }

  // loop with half rotation (2 * pi/2)
  if (true)
    {
      deallog << "---------------------------" << std::endl
              << "moebius, half rotation (2* Pi/2)" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::moebius(tr, 20, 2, 10.0, 2.0);
      GridOut go;
      go.set_flags(GridOutFlags::Ucd(true));
      go.write_ucd(tr, out);
    }

  // loop with three quarter rotation (3 * pi/2)
  if (true)
    {
      deallog << "---------------------------" << std::endl
              << "moebius, three quarter rotation (3* Pi/2)" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::moebius(tr, 20, 3, 10.0, 2.0);
      GridOut go;
      go.set_flags(GridOutFlags::Ucd(true));
      go.write_ucd(tr, out);
    }

  // loop with full rotation (1 * pi/2)
  if (true)
    {
      deallog << "---------------------------" << std::endl
              << "moebius, full rotation (2* Pi)" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::moebius(tr, 20, 4, 10.0, 2.0);
      GridOut go;
      go.set_flags(GridOutFlags::Ucd(true));
      go.write_ucd(tr, out);
    }
}


int
main()
{
  initlog();

  test<3>(deallog.get_file_stream());
}
