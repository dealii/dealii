// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <array>

#include "../tests.h"

/*
 * Verify that the correct exception (ExcInvalidInputOrientation) is thrown
 * when this function is called for invalid input.
 */

template <int dim>
void
check_parallelepiped(std::ostream &logfile)
{
  // Data structure defining dim coordinates that make up a
  // parallelepiped.
  std::array<Tensor<1, dim>, dim> edges;

  switch (dim)
    {
      case 1:
        edges[0][0] = -0.5;
        break;

      case 2:
        edges[0][1] = 0.5;
        edges[1][0] = 0.5;
        break;

      case 3:
        edges[0][0] = 1.0;
        edges[1][1] = 1.0;
        edges[2][2] = -1.0;
        break;

      default:
        Assert(false, ExcInternalError());
    }

  Point<dim> origin;

  Triangulation<dim> triangulation(Triangulation<dim>::MeshSmoothing::none,
                                   /*check_for_distorted_cells*/ true);
  std::vector<unsigned int> subdivisions;

  try
    {
      GridGenerator::subdivided_parallelepiped<dim>(
        triangulation, origin, edges, subdivisions, false);
    }
  catch (ExceptionBase &exc)
    {
      logfile << exc.get_exc_name() << std::endl;
    }
}

int
main()
{
  deal_II_exceptions::disable_abort_on_exception();
  initlog();
  std::ostream &logfile = deallog.get_file_stream();

  check_parallelepiped<1>(logfile);
  check_parallelepiped<2>(logfile);
  check_parallelepiped<3>(logfile);
  logfile << "OK" << std::endl;
}
