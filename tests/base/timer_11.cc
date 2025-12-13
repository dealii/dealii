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

// test that a timer cannot be stopped while not running

#include <deal.II/base/timer.h>

#include "../tests.h"

int
main()
{
  initlog();

  Timer t;
  t.stop();

  // this will throw an exception
  t.stop();

  deallog << "OK" << std::endl;
}
