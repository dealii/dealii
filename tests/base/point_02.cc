// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2017 by the deal.II authors
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


// test for Point::distance_square()

#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
void
check()
{
  Point<dim> p1, p2;
  for (unsigned int i = 0; i < dim; ++i)
    {
      p1[i] = 10.0 + 0.12345 * i;
      p1[i] = 0.5 + 0.6789 * i;
    }

  const double d  = p1.distance(p2);
  const double d2 = p1.distance_square(p2);

  AssertThrow(std::abs(d - std::sqrt(d2)) < 1e-10, ExcInternalError())

      deallog
    << "Ok" << std::endl;
}

int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);

  check<1>();
  check<2>();
  check<3>();
}
