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

// the same as timer_02.cc but test the function n_laps().

#include <deal.II/base/timer.h>

#include "../tests.h"


int
main()
{
  initlog();

  Timer t;

  AssertThrow(t.n_laps() == 1, ExcInternalError());
  t.stop();
  AssertThrow(t.n_laps() == 1, ExcInternalError());
  t.start();
  AssertThrow(t.n_laps() == 2, ExcInternalError());
  t.stop();
  AssertThrow(t.n_laps() == 2, ExcInternalError());
  t.reset();
  AssertThrow(t.n_laps() == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}
