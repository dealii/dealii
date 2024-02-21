// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check numbers::signaling_nan<Point>
//
// the test only checks that the function can be called. It would have
// been nicer to output the tensors (which would also verify the
// correct output type, as well as that indeed *each* element is
// correctly filled), but outputting a sNaN triggers a floating point
// exception as well

#include <deal.II/base/signaling_nan.h>

#include <limits>

#include "../tests.h"


template <typename T>
void
check()
{
  numbers::signaling_nan<Point<1, T>>();
  numbers::signaling_nan<Point<2, T>>();
  numbers::signaling_nan<Point<3, T>>();

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<float>();
  check<double>();

  return 0;
}
