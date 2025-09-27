// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check Table<2,double> with elements numbers::signaling_nan<double>
//
// the test only checks that the function can be called. It would have
// been nicer to output the tensors (which would also verify the
// correct output type, as well as that indeed *each* element is
// correctly filled), but outputting a sNaN triggers a floating point
// exception as well

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/table.h>

#include <limits>

#include "../tests.h"


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
  initlog();
  deallog << std::setprecision(3);

  check<float>();
  check<double>();

  return 0;
}
