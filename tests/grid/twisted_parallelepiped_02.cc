// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/tensor.h>

#include <fstream>

/*
 * Verify that the correct exception (ExcGridHasInvalidCell) is thrown when
 * this function is called for invalid input.
 */

template<int dim>
void check_parallelepiped (std::ofstream &logfile)
{
  // Data structure defining dim coordinates that make up a
  // parallelepiped.
  std_cxx11::array<Tensor<1, dim>, dim> edges;

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
      Assert (false, ExcInternalError ());
    }

  Point<dim> origin;

  Triangulation<dim> triangulation;
  std::vector<unsigned int> subdivisions;

  try
    {
      GridGenerator::subdivided_parallelepiped<dim>(triangulation, origin, edges,
                                                    subdivisions, false);
    }
  catch (ExceptionBase &exc)
    {
      logfile << exc.get_exc_name() << std::endl;
    }
}

int main ()
{
  deal_II_exceptions::disable_abort_on_exception();
  std::ofstream logfile ("output");

  check_parallelepiped<1>(logfile);
  check_parallelepiped<2>(logfile);
  check_parallelepiped<3>(logfile);
  logfile << "OK" << std::endl;
}
