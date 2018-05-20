// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II authors
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

// check Table<2,double> with elements numbers::signaling_nan<double>
//
// the test only checks that the function can be called. It would have
// been nicer to output the tensors (which would also verify the
// correct output type, as well as that indeed *each* element is
// correctly filled), but outputting a sNaN triggers a floating point
// exception as well

#include "../tests.h"
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/table.h>
#include <limits>

template <typename T>
void
check()
{
  Table<2, T> t;

  t.reinit(2, 3);
  t.fill(numbers::signaling_nan<T>());

  deallog << "OK" << std::endl;
}

int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);

  check<float>();
  check<double>();

  return 0;
}
